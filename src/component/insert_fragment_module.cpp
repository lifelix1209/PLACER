#include "pipeline.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/sam.h>

namespace placer {
namespace {

struct ClipInfo {
    int32_t leading = 0;
    int32_t trailing = 0;
    int32_t right_anchor_after_leading = 0;
    int32_t left_anchor_before_trailing = 0;
    int32_t ref_end = 0;
};

struct InsOp {
    int32_t start = 0;
    int32_t len = 0;
    int32_t ref_pos = -1;
    int32_t left_anchor = 0;
    int32_t right_anchor = 0;
};

struct SAEntry {
    std::string rname;
    int32_t pos = 0;  // 1-based in SA tag
    char strand = '+';
    std::string cigar;
    int32_t mapq = 0;
    int32_t nm = 0;
};

struct CigarStringOp {
    int32_t len = 0;
    char op = '\0';
};

struct QueryInterval {
    int32_t qstart = 0;
    int32_t qend = 0;
    int32_t q_aligned = 0;
    int32_t leading_s = 0;
    int32_t trailing_s = 0;
};

struct NormalizedAln {
    std::string chrom;
    int32_t ref_start = -1;
    int32_t ref_end = -1;
    int32_t qstart = -1;
    int32_t qend = -1;
    bool is_reverse = false;
    int32_t mapq = 0;
    int32_t nm = -1;
    bool is_primary = false;
    bool reliable = true;
    int32_t anchor_len = 0;
    bool has_large_indel_near_bp = false;
};

bool is_match_like(int op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}

bool is_match_like_char(char op) {
    return op == 'M' || op == '=' || op == 'X';
}

bool consumes_query_char(char op) {
    return op == 'M' || op == '=' || op == 'X' || op == 'I' || op == 'S';
}

bool consumes_ref_char(char op) {
    return op == 'M' || op == '=' || op == 'X' || op == 'D' || op == 'N';
}

int find_first_non_hard_clip(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

int find_last_non_hard_clip(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = n_cigar - 1; i >= 0; --i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

int32_t compute_ref_end(const ReadView& read) {
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return read.pos();
    }
    int32_t ref_pos = read.pos();
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        if ((bam_cigar_type(op) & 2) != 0) {
            ref_pos += len;
        }
    }
    return ref_pos;
}

ClipInfo analyze_clip_info(const ReadView& read) {
    ClipInfo info;

    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        info.ref_end = read.pos();
        return info;
    }

    const int first = find_first_non_hard_clip(cigar, n_cigar);
    const int last = find_last_non_hard_clip(cigar, n_cigar);
    if (first < 0 || last < 0 || first > last) {
        info.ref_end = compute_ref_end(read);
        return info;
    }

    if (bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
        info.leading = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
        int32_t flank = 0;
        for (int i = first + 1; i <= last; ++i) {
            const int op = bam_cigar_op(cigar[i]);
            if (!is_match_like(op)) {
                break;
            }
            flank += static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        }
        info.right_anchor_after_leading = flank;
    }

    if (bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
        info.trailing = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
        int32_t flank = 0;
        for (int i = last - 1; i >= first; --i) {
            const int op = bam_cigar_op(cigar[i]);
            if (!is_match_like(op)) {
                break;
            }
            flank += static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        }
        info.left_anchor_before_trailing = flank;
    }

    info.ref_end = compute_ref_end(read);
    return info;
}

int32_t contiguous_match_run_after(const uint32_t* cigar, int32_t n_cigar, int32_t idx) {
    int32_t run = 0;
    for (int32_t j = idx + 1; j < n_cigar; ++j) {
        const int op = bam_cigar_op(cigar[j]);
        if (!is_match_like(op)) {
            break;
        }
        run += static_cast<int32_t>(bam_cigar_oplen(cigar[j]));
    }
    return run;
}

std::vector<InsOp> find_long_insertions(const ReadView& read, int32_t min_long_ins) {
    std::vector<InsOp> ops;

    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return ops;
    }

    int32_t qpos = 0;
    int32_t rpos = read.pos();
    int32_t prev_match_run = 0;
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));

        if (is_match_like(op)) {
            prev_match_run += len;
            qpos += len;
            rpos += len;
            continue;
        }

        if (op == BAM_CINS) {
            if (len >= min_long_ins) {
                InsOp row;
                row.start = qpos;
                row.len = len;
                row.ref_pos = rpos;
                row.left_anchor = prev_match_run;
                row.right_anchor = contiguous_match_run_after(cigar, n_cigar, i);
                ops.push_back(row);
            }
            qpos += len;
            prev_match_run = 0;
            continue;
        }

        if ((bam_cigar_type(op) & 1) != 0) {
            qpos += len;
        }
        if ((bam_cigar_type(op) & 2) != 0) {
            rpos += len;
        }
        prev_match_run = 0;
    }

    return ops;
}

bool parse_int32(const char* text, int32_t& out) {
    if (!text || !*text) {
        return false;
    }
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || (end && *end != '\0')) {
        return false;
    }
    out = static_cast<int32_t>(value);
    return true;
}

std::vector<SAEntry> parse_sa_tag_z(const std::string& sa_z) {
    std::vector<SAEntry> out;
    if (sa_z.empty()) {
        return out;
    }

    size_t start = 0;
    while (start < sa_z.size()) {
        const size_t rec_end = sa_z.find(';', start);
        const size_t rec_len =
            (rec_end == std::string::npos) ? (sa_z.size() - start) : (rec_end - start);
        if (rec_len == 0) {
            break;
        }

        const std::string rec = sa_z.substr(start, rec_len);
        std::vector<std::string> fields;
        fields.reserve(6);
        size_t field_start = 0;
        while (field_start <= rec.size()) {
            size_t comma = rec.find(',', field_start);
            if (comma == std::string::npos) {
                comma = rec.size();
            }
            fields.push_back(rec.substr(field_start, comma - field_start));
            if (comma == rec.size()) {
                break;
            }
            field_start = comma + 1;
        }

        if (fields.size() >= 6) {
            SAEntry entry;
            entry.rname = fields[0];
            parse_int32(fields[1].c_str(), entry.pos);
            entry.strand = fields[2].empty() ? '+' : fields[2][0];
            entry.cigar = fields[3];
            parse_int32(fields[4].c_str(), entry.mapq);
            parse_int32(fields[5].c_str(), entry.nm);
            if (!entry.cigar.empty()) {
                out.push_back(std::move(entry));
            }
        }

        if (rec_end == std::string::npos) {
            break;
        }
        start = rec_end + 1;
    }

    return out;
}

bool parse_cigar_ops(const std::string& cigar, std::vector<CigarStringOp>& ops) {
    ops.clear();
    if (cigar.empty()) {
        return false;
    }
    int32_t num = 0;
    bool has_num = false;
    for (char c : cigar) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            has_num = true;
            num = num * 10 + (c - '0');
            continue;
        }
        if (!has_num || num <= 0) {
            return false;
        }
        ops.push_back({num, c});
        has_num = false;
        num = 0;
    }
    if (has_num) {
        return false;
    }
    return !ops.empty();
}

bool cigar_to_query_interval(
    const std::vector<CigarStringOp>& ops,
    int32_t read_len,
    QueryInterval& out,
    int32_t& ref_aligned_len,
    bool& has_large_indel) {
    out = QueryInterval{};
    ref_aligned_len = 0;
    has_large_indel = false;
    if (ops.empty() || read_len <= 0) {
        return false;
    }

    int32_t i = 0;
    while (i < static_cast<int32_t>(ops.size()) &&
           (ops[static_cast<size_t>(i)].op == 'H' || ops[static_cast<size_t>(i)].op == 'S')) {
        if (ops[static_cast<size_t>(i)].op == 'S') {
            out.leading_s += ops[static_cast<size_t>(i)].len;
        }
        ++i;
    }

    int32_t j = static_cast<int32_t>(ops.size()) - 1;
    while (j >= 0 &&
           (ops[static_cast<size_t>(j)].op == 'H' || ops[static_cast<size_t>(j)].op == 'S')) {
        if (ops[static_cast<size_t>(j)].op == 'S') {
            out.trailing_s += ops[static_cast<size_t>(j)].len;
        }
        --j;
    }

    int32_t q_aligned = 0;
    for (const auto& op : ops) {
        if (consumes_query_char(op.op) && op.op != 'S') {
            q_aligned += op.len;
        }
        if (consumes_ref_char(op.op)) {
            ref_aligned_len += op.len;
        }
        if ((op.op == 'I' || op.op == 'D') && op.len >= 20) {
            has_large_indel = true;
        }
    }
    if (q_aligned <= 0) {
        return false;
    }
    out.q_aligned = q_aligned;

    int32_t qstart = 0;
    if (out.leading_s > 0) {
        qstart = out.leading_s;
    } else if (out.trailing_s > 0) {
        qstart = read_len - out.trailing_s - q_aligned;
    } else if (q_aligned == read_len) {
        qstart = 0;
    } else {
        return false;
    }

    if (qstart < 0) {
        return false;
    }
    const int32_t qend = qstart + q_aligned;
    if (qend > read_len || qend <= qstart) {
        return false;
    }

    out.qstart = qstart;
    out.qend = qend;
    return true;
}

bool bam_to_query_interval(const ReadView& read, const ClipInfo& clip, QueryInterval& out) {
    out = QueryInterval{};
    const int32_t read_len = read.seq_len();
    if (read_len <= 0) {
        return false;
    }

    int32_t q_aligned = 0;
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return false;
    }

    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CINS) {
            q_aligned += len;
        }
    }
    if (q_aligned <= 0) {
        return false;
    }

    out.leading_s = clip.leading;
    out.trailing_s = clip.trailing;
    out.qstart = clip.leading;
    out.qend = read_len - clip.trailing;
    if (out.qend <= out.qstart || out.qend > read_len) {
        out.qstart = clip.leading;
        out.qend = out.qstart + q_aligned;
    }
    if (out.qend <= out.qstart || out.qend > read_len) {
        return false;
    }
    out.q_aligned = out.qend - out.qstart;
    return true;
}

double interval_iou(int32_t a0, int32_t a1, int32_t b0, int32_t b1) {
    const int32_t inter = std::max(0, std::min(a1, b1) - std::max(a0, b0));
    const int32_t uni = std::max(a1, b1) - std::min(a0, b0);
    if (uni <= 0) {
        return 0.0;
    }
    return static_cast<double>(inter) / static_cast<double>(uni);
}

int32_t overlap_len(int32_t a0, int32_t a1, int32_t b0, int32_t b1) {
    return std::max(0, std::min(a1, b1) - std::max(a0, b0));
}

int32_t nonoverlap_len(int32_t a0, int32_t a1, int32_t b0, int32_t b1) {
    const int32_t a_len = std::max(0, a1 - a0);
    return std::max(0, a_len - overlap_len(a0, a1, b0, b1));
}

int32_t distance_to_breakpoint(int32_t ref_start, int32_t ref_end, int32_t bp) {
    if (ref_end < ref_start) {
        std::swap(ref_end, ref_start);
    }
    if (bp < ref_start) {
        return ref_start - bp;
    }
    if (bp > ref_end) {
        return bp - ref_end;
    }
    return 0;
}

int32_t overlap_with_window(int32_t seg_start, int32_t seg_end, int32_t win_start, int32_t win_end) {
    return std::max(0, std::min(seg_end, win_end) - std::max(seg_start, win_start));
}

int32_t anchor_len_from_bam(const ReadView& read, int32_t bp, int32_t w_anchor, bool& has_large_indel_near_bp) {
    has_large_indel_near_bp = false;
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return 0;
    }

    const int32_t win_start = bp - w_anchor;
    const int32_t win_end = bp + w_anchor + 1;

    int32_t rpos = read.pos();
    int32_t anchored = 0;
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));

        if (is_match_like(op)) {
            anchored += overlap_with_window(rpos, rpos + len, win_start, win_end);
            rpos += len;
            continue;
        }
        if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            const int32_t ov = overlap_with_window(rpos, rpos + len, win_start, win_end);
            anchored += ov;
            if (len >= 20 && ov > 0) {
                has_large_indel_near_bp = true;
            }
            rpos += len;
            continue;
        }
        if (op == BAM_CINS) {
            if (rpos >= win_start && rpos < win_end) {
                anchored += len;
                if (len >= 20) {
                    has_large_indel_near_bp = true;
                }
            }
            continue;
        }
    }

    return anchored;
}

int32_t anchor_len_from_ops(
    const std::vector<CigarStringOp>& ops,
    int32_t ref_start,
    int32_t bp,
    int32_t w_anchor,
    bool& has_large_indel_near_bp) {
    has_large_indel_near_bp = false;
    if (ops.empty()) {
        return 0;
    }

    const int32_t win_start = bp - w_anchor;
    const int32_t win_end = bp + w_anchor + 1;

    int32_t rpos = ref_start;
    int32_t anchored = 0;
    for (const auto& op : ops) {
        if (is_match_like_char(op.op)) {
            anchored += overlap_with_window(rpos, rpos + op.len, win_start, win_end);
            rpos += op.len;
            continue;
        }
        if (op.op == 'D' || op.op == 'N') {
            const int32_t ov = overlap_with_window(rpos, rpos + op.len, win_start, win_end);
            anchored += ov;
            if (op.len >= 20 && ov > 0) {
                has_large_indel_near_bp = true;
            }
            rpos += op.len;
            continue;
        }
        if (op.op == 'I') {
            if (rpos >= win_start && rpos < win_end) {
                anchored += op.len;
                if (op.len >= 20) {
                    has_large_indel_near_bp = true;
                }
            }
        }
    }
    return anchored;
}

std::string sanitize_fragment_token(std::string s) {
    for (char& c : s) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            c = '_';
        }
        if (c == '|') {
            c = '/';
        }
    }
    return s;
}

std::string fasta_safe_header(std::string s) {
    for (char& c : s) {
        if (c == '\n' || c == '\r' || c == '\t') {
            c = '_';
            continue;
        }
        if (std::isspace(static_cast<unsigned char>(c))) {
            c = '_';
        }
    }
    return s;
}

const char* source_tag(InsertionFragmentSource source) {
    switch (source) {
        case InsertionFragmentSource::kClipRefLeft: return "clipRefL";
        case InsertionFragmentSource::kClipRefRight: return "clipRefR";
        case InsertionFragmentSource::kCigarInsertion: return "ins";
        case InsertionFragmentSource::kSplitSa: return "splitSA";
        default: return "unknown";
    }
}

void write_wrapped_sequence(std::ostream& out, const std::string& seq) {
    constexpr size_t kLineWidth = 80;
    for (size_t i = 0; i < seq.size(); i += kLineWidth) {
        out << seq.substr(i, std::min(kLineWidth, seq.size() - i)) << "\n";
    }
}

int compare_flank_candidate(const NormalizedAln& a, const NormalizedAln& b, int32_t bp, double alpha) {
    const double penalty_a =
        static_cast<double>(std::max(0, a.nm)) + 5.0 * static_cast<double>(a.has_large_indel_near_bp ? 1 : 0);
    const double penalty_b =
        static_cast<double>(std::max(0, b.nm)) + 5.0 * static_cast<double>(b.has_large_indel_near_bp ? 1 : 0);
    const double score_a = static_cast<double>(a.anchor_len) - alpha * penalty_a;
    const double score_b = static_cast<double>(b.anchor_len) - alpha * penalty_b;
    if (score_a != score_b) {
        return (score_a > score_b) ? -1 : 1;
    }
    if (a.anchor_len != b.anchor_len) {
        return (a.anchor_len > b.anchor_len) ? -1 : 1;
    }
    if (a.mapq != b.mapq) {
        return (a.mapq > b.mapq) ? -1 : 1;
    }
    if (a.nm != b.nm) {
        return (a.nm < b.nm) ? -1 : 1;
    }
    if (a.is_primary != b.is_primary) {
        return a.is_primary ? -1 : 1;
    }
    if (a.ref_start != b.ref_start) {
        return (a.ref_start < b.ref_start) ? -1 : 1;
    }
    const int32_t dist_a = distance_to_breakpoint(a.ref_start, a.ref_end, bp);
    const int32_t dist_b = distance_to_breakpoint(b.ref_start, b.ref_end, bp);
    if (dist_a != dist_b) {
        return (dist_a < dist_b) ? -1 : 1;
    }
    return 0;
}

std::mutex g_ins_fragments_fasta_mutex;

}  // namespace

SplitSAFragmentModule::SplitSAFragmentModule(PipelineConfig config)
    : config_(std::move(config)) {}

std::vector<InsertionFragment> SplitSAFragmentModule::extract(
    const ComponentCall& component,
    const std::vector<BamRecordPtr>& bin_records) const {
    std::vector<InsertionFragment> out;

    const int32_t min_len = config_.min_sa_aln_len_for_seq_extract;
    if (config_.max_sa_per_read <= 0) {
        return out;
    }
    const int32_t max_sa = config_.max_sa_per_read;

    const int32_t w_ref = std::max(200, config_.bin_size / 2);
    const int32_t w_anchor = 150;
    const int32_t mh_max = 20;
    const int32_t l_left = std::max(min_len, 100);
    const int32_t l_right = std::max(min_len, 100);
    const int32_t trim_q = 12;
    const int32_t eps_bp = 5;
    const double alpha = 1.0;

    std::unordered_map<std::string, std::vector<NormalizedAln>> supplementary_by_read;
    supplementary_by_read.reserve(component.read_indices.size());
    for (size_t idx_in_bin : component.read_indices) {
        if (idx_in_bin >= bin_records.size() || !bin_records[idx_in_bin]) {
            continue;
        }
        ReadView read(bin_records[idx_in_bin].get());
        const uint16_t flag = read.flag();
        if ((flag & BAM_FSUPPLEMENTARY) == 0 || (flag & BAM_FSECONDARY) != 0) {
            continue;
        }

        ClipInfo clip = analyze_clip_info(read);
        QueryInterval q_interval;
        if (!bam_to_query_interval(read, clip, q_interval)) {
            continue;
        }
        bool has_large_indel = false;
        const int32_t anchor_len = anchor_len_from_bam(read, component.anchor_pos, w_anchor, has_large_indel);

        int64_t nm = -1;
        read.get_int_tag("NM", nm);

        NormalizedAln aln;
        aln.chrom = component.chrom;
        aln.ref_start = read.pos();
        aln.ref_end = clip.ref_end;
        aln.qstart = q_interval.qstart;
        aln.qend = q_interval.qend;
        aln.is_reverse = (flag & BAM_FREVERSE) != 0;
        aln.mapq = read.mapq();
        aln.nm = static_cast<int32_t>(nm);
        aln.is_primary = false;
        aln.reliable = true;
        aln.anchor_len = anchor_len;
        aln.has_large_indel_near_bp = has_large_indel;

        const std::string read_id = sanitize_fragment_token(std::string(read.qname()));
        supplementary_by_read[read_id].push_back(std::move(aln));
    }

    for (size_t idx_in_bin : component.read_indices) {
        if (idx_in_bin >= bin_records.size() || !bin_records[idx_in_bin]) {
            continue;
        }

        ReadView read(bin_records[idx_in_bin].get());
        const uint16_t flag = read.flag();

        if ((flag & BAM_FSUPPLEMENTARY) != 0 || (flag & BAM_FSECONDARY) != 0) {
            continue;
        }
        if (!read.has_sa_tag()) {
            continue;
        }

        std::string sa_z;
        if (!read.get_string_tag("SA", sa_z) || sa_z.empty()) {
            continue;
        }
        const auto entries = parse_sa_tag_z(sa_z);
        if (entries.empty()) {
            continue;
        }

        const std::string read_id = sanitize_fragment_token(std::string(read.qname()));
        const int32_t read_len = read.seq_len();
        const bool is_reverse = (flag & BAM_FREVERSE) != 0;

        ClipInfo clip = analyze_clip_info(read);
        QueryInterval primary_q;
        if (!bam_to_query_interval(read, clip, primary_q)) {
            continue;
        }
        bool primary_large_indel = false;
        const int32_t primary_anchor =
            anchor_len_from_bam(read, component.anchor_pos, w_anchor, primary_large_indel);
        int64_t primary_nm = -1;
        read.get_int_tag("NM", primary_nm);

        std::vector<NormalizedAln> core_records;
        core_records.reserve(entries.size() + 1);
        NormalizedAln primary;
        primary.chrom = component.chrom;
        primary.ref_start = read.pos();
        primary.ref_end = clip.ref_end;
        primary.qstart = primary_q.qstart;
        primary.qend = primary_q.qend;
        primary.is_reverse = is_reverse;
        primary.mapq = read.mapq();
        primary.nm = static_cast<int32_t>(primary_nm);
        primary.is_primary = true;
        primary.reliable = true;
        primary.anchor_len = primary_anchor;
        primary.has_large_indel_near_bp = primary_large_indel;
        core_records.push_back(primary);

        std::vector<std::pair<SAEntry, NormalizedAln>> sa_all;
        sa_all.reserve(entries.size());

        const auto supp_it = supplementary_by_read.find(read_id);
        const bool has_explicit_supp =
            (supp_it != supplementary_by_read.end() && !supp_it->second.empty());

        for (const auto& entry : entries) {
            std::vector<CigarStringOp> ops;
            if (!parse_cigar_ops(entry.cigar, ops)) {
                continue;
            }

            QueryInterval q_interval;
            int32_t ref_aligned = 0;
            bool has_large_indel = false;
            if (!cigar_to_query_interval(ops, read_len, q_interval, ref_aligned, has_large_indel)) {
                continue;
            }
            if (q_interval.qend - q_interval.qstart < min_len) {
                continue;
            }

            NormalizedAln sa;
            sa.chrom = entry.rname;
            sa.ref_start = std::max(0, entry.pos - 1);
            sa.ref_end = sa.ref_start + ref_aligned;
            sa.qstart = q_interval.qstart;
            sa.qend = q_interval.qend;
            sa.is_reverse = (entry.strand == '-');
            sa.mapq = entry.mapq;
            sa.nm = entry.nm;
            sa.is_primary = false;
            sa.reliable = !has_explicit_supp;
            sa.anchor_len = anchor_len_from_ops(ops, sa.ref_start, component.anchor_pos, w_anchor, has_large_indel);
            sa.has_large_indel_near_bp = has_large_indel;

            if (has_explicit_supp) {
                double best_iou = 0.0;
                int32_t best_q_start_diff = std::numeric_limits<int32_t>::max();
                int32_t best_q_end_diff = std::numeric_limits<int32_t>::max();
                for (const auto& supp : supp_it->second) {
                    if (supp.is_reverse != sa.is_reverse) {
                        continue;
                    }
                    const double iou = interval_iou(supp.qstart, supp.qend, sa.qstart, sa.qend);
                    if (iou > best_iou) {
                        best_iou = iou;
                        best_q_start_diff = std::abs(supp.qstart - sa.qstart);
                        best_q_end_diff = std::abs(supp.qend - sa.qend);
                    }
                }
                sa.reliable =
                    best_iou >= 0.9 &&
                    best_q_start_diff <= 5 &&
                    best_q_end_diff <= 5;
            }

            if (sa.reliable) {
                core_records.push_back(sa);
            }
            sa_all.push_back({entry, sa});
        }

        int32_t flank_idx = -1;
        for (size_t i = 0; i < core_records.size(); ++i) {
            const auto& row = core_records[i];
            if (row.chrom != component.chrom) {
                continue;
            }
            if (distance_to_breakpoint(row.ref_start, row.ref_end, component.anchor_pos) > w_ref) {
                continue;
            }
            if (flank_idx < 0 ||
                compare_flank_candidate(row, core_records[static_cast<size_t>(flank_idx)], component.anchor_pos, alpha) < 0) {
                flank_idx = static_cast<int32_t>(i);
            }
        }

        bool emitted_robust = false;
        if (flank_idx >= 0) {
            int32_t mate_idx = -1;
            int32_t best_nonoverlap = -1;
            for (size_t i = 0; i < core_records.size(); ++i) {
                if (static_cast<int32_t>(i) == flank_idx) {
                    continue;
                }
                const auto& mate = core_records[i];
                const auto& flank = core_records[static_cast<size_t>(flank_idx)];
                const int32_t nonov = nonoverlap_len(
                    mate.qstart, mate.qend, flank.qstart, flank.qend);
                if (nonov > best_nonoverlap) {
                    best_nonoverlap = nonov;
                    mate_idx = static_cast<int32_t>(i);
                }
            }

            if (mate_idx >= 0 && best_nonoverlap >= min_len) {
                const auto& flank = core_records[static_cast<size_t>(flank_idx)];
                const auto& mate = core_records[static_cast<size_t>(mate_idx)];

                int32_t q_junc = -1;
                bool overlap_too_large = false;
                if (flank.qend <= mate.qstart + mh_max) {
                    q_junc = flank.qend;
                } else if (mate.qend <= flank.qstart + mh_max) {
                    q_junc = flank.qstart;
                } else {
                    const int32_t ov_start = std::max(flank.qstart, mate.qstart);
                    const int32_t ov_end = std::min(flank.qend, mate.qend);
                    const int32_t ov_len = std::max(0, ov_end - ov_start);
                    if (ov_len <= mh_max) {
                        q_junc = (ov_start + ov_end) / 2;
                    } else {
                        overlap_too_large = true;
                    }
                }

                if (!overlap_too_large && q_junc >= 0) {
                    int32_t opp_start = std::max(0, q_junc - l_left);
                    int32_t opp_end = std::min(read_len, q_junc + l_right);

                    int32_t core_start = flank.qstart + trim_q;
                    int32_t core_end = flank.qend - trim_q;
                    if (core_end <= core_start) {
                        core_start = flank.qstart;
                        core_end = flank.qend;
                    }

                    if (opp_start < core_end && opp_end > core_start) {
                        if (q_junc <= core_start) {
                            opp_end = std::min(opp_end, core_start);
                        } else if (q_junc >= core_end) {
                            opp_start = std::max(opp_start, core_end);
                        } else {
                            const int32_t left_len = std::max(0, core_start - opp_start);
                            const int32_t right_len = std::max(0, opp_end - core_end);
                            if (left_len >= right_len) {
                                opp_end = core_start;
                            } else {
                                opp_start = core_end;
                            }
                        }
                    }

                    const int32_t opp_len = opp_end - opp_start;
                    if (opp_len >= min_len) {
                        int32_t ref_junc_pos = -1;
                        const bool flank_on_left =
                            std::abs(q_junc - flank.qend) <= std::abs(q_junc - flank.qstart);
                        if (flank_on_left) {
                            ref_junc_pos = flank.is_reverse ? flank.ref_start : flank.ref_end;
                        } else {
                            ref_junc_pos = flank.is_reverse ? flank.ref_end : flank.ref_start;
                        }

                        const bool left = (ref_junc_pos <= component.anchor_pos + eps_bp);
                        const bool right = (ref_junc_pos >= component.anchor_pos - eps_bp);
                        ReferenceSide ref_side = ReferenceSide::kUnknown;
                        if (!(left && right)) {
                            if (left) {
                                ref_side = ReferenceSide::kRefLeft;
                            } else if (right) {
                                ref_side = ReferenceSide::kRefRight;
                            }
                        }

                        InsertionFragment frag;
                        frag.chrom = component.chrom;
                        frag.anchor_pos = component.anchor_pos;
                        frag.read_id = read_id;
                        frag.read_index = idx_in_bin;
                        frag.class_mask = kCandidateSplitSaSupplementary;
                        frag.is_reverse = is_reverse;
                        frag.source = InsertionFragmentSource::kSplitSa;
                        frag.start = opp_start;
                        frag.length = opp_len;
                        frag.read_len = read_len;
                        frag.anchor_len = flank.anchor_len;
                        frag.ref_side = ref_side;
                        frag.ref_junc_pos = ref_junc_pos;
                        frag.nm = flank.nm;
                        frag.split_sa_reliable = flank.reliable && mate.reliable;
                        frag.sequence = read.decode_subsequence(frag.start, frag.length);
                        frag.fragment_id =
                            component.chrom + ":" + std::to_string(component.anchor_pos) +
                            "|read=" + read_id +
                            "|idx=" + std::to_string(idx_in_bin) +
                            "|strand=" + std::string(is_reverse ? "-" : "+") +
                            "|src=" + source_tag(frag.source) +
                            "|qj=" + std::to_string(q_junc) +
                            "|rj=" + std::to_string(ref_junc_pos) +
                            "|len=" + std::to_string(frag.length);
                        if (!frag.sequence.empty()) {
                            out.push_back(std::move(frag));
                            emitted_robust = true;
                        }
                    }
                }
            }
        }

        int32_t emitted_diag = emitted_robust ? 1 : 0;
        int32_t sa_index = 0;
        for (const auto& row : sa_all) {
            if (emitted_diag >= max_sa) {
                break;
            }
            ++sa_index;
            const auto& entry = row.first;
            const auto& sa = row.second;
            if (sa.qstart < 0 || sa.qend <= sa.qstart || sa.qend > read_len) {
                continue;
            }

            InsertionFragment frag;
            frag.chrom = component.chrom;
            frag.anchor_pos = component.anchor_pos;
            frag.read_id = read_id;
            frag.read_index = idx_in_bin;
            frag.class_mask = kCandidateSplitSaSupplementary;
            frag.is_reverse = is_reverse;
            frag.source = InsertionFragmentSource::kSplitSa;
            frag.start = sa.qstart;
            frag.length = sa.qend - sa.qstart;
            frag.read_len = read_len;
            frag.anchor_len = sa.anchor_len;
            frag.ref_side = ReferenceSide::kUnknown;
            frag.ref_junc_pos = -1;
            frag.nm = sa.nm;
            frag.split_sa_reliable = sa.reliable;
            frag.sequence = read.decode_subsequence(frag.start, frag.length);
            frag.fragment_id =
                component.chrom + ":" + std::to_string(component.anchor_pos) +
                "|read=" + read_id +
                "|idx=" + std::to_string(idx_in_bin) +
                "|strand=" + std::string(is_reverse ? "-" : "+") +
                "|src=" + source_tag(frag.source) +
                "|sa_idx=" + std::to_string(sa_index) +
                "|sa=" + sanitize_fragment_token(entry.rname) + ":" + std::to_string(entry.pos) + entry.strand +
                "|qlen=" + std::to_string(frag.length);
            if (frag.length < min_len || frag.sequence.empty()) {
                continue;
            }
            out.push_back(std::move(frag));
            ++emitted_diag;
        }
    }

    return out;
}

CigarInsertionFragmentModule::CigarInsertionFragmentModule(PipelineConfig config)
    : config_(std::move(config)) {
    if (!config_.ins_fragments_fasta_path.empty()) {
        std::lock_guard<std::mutex> lock(g_ins_fragments_fasta_mutex);
        std::ofstream out(config_.ins_fragments_fasta_path, std::ios::out | std::ios::trunc);
        if (!out.is_open()) {
            std::cerr << "[InsertFragments] cannot open for writing: "
                      << config_.ins_fragments_fasta_path << "\n";
        }
    }
}

std::vector<InsertionFragment> CigarInsertionFragmentModule::extract(
    const ComponentCall& component,
    const std::vector<BamRecordPtr>& bin_records) const {
    std::vector<InsertionFragment> out;

    const int32_t min_clip = config_.min_soft_clip_for_seq_extract;
    const int32_t min_ins = config_.min_long_ins_for_seq_extract;

    const bool write_fasta = !config_.ins_fragments_fasta_path.empty();
    std::ostringstream fasta_buffer;

    auto emit_fragment = [&](InsertionFragment&& frag) {
        if (frag.sequence.empty()) {
            return;
        }
        out.push_back(frag);
        if (write_fasta) {
            fasta_buffer << ">" << fasta_safe_header(out.back().fragment_id) << "\n";
            write_wrapped_sequence(fasta_buffer, out.back().sequence);
        }
    };

    const auto& indices = component.read_indices;
    for (size_t idx_in_bin : indices) {
        if (idx_in_bin >= bin_records.size() || !bin_records[idx_in_bin]) {
            continue;
        }

        ReadView read(bin_records[idx_in_bin].get());
        const std::string read_id = sanitize_fragment_token(std::string(read.qname()));

        const uint16_t flag = read.flag();
        const bool is_reverse = (flag & BAM_FREVERSE) != 0;
        uint8_t class_mask = 0;
        if (read.has_sa_tag() || ((flag & BAM_FSUPPLEMENTARY) != 0)) {
            class_mask |= kCandidateSplitSaSupplementary;
        }

        const ClipInfo clip = analyze_clip_info(read);
        if (clip.leading >= min_clip || clip.trailing >= min_clip) {
            class_mask |= kCandidateSoftClip;
        }

        auto ins_ops = find_long_insertions(read, min_ins);
        if (!ins_ops.empty()) {
            class_mask |= kCandidateLongInsertion;
        }

        int64_t nm = -1;
        read.get_int_tag("NM", nm);
        const int32_t read_len = read.seq_len();

        if (clip.leading >= min_clip) {
            InsertionFragment frag;
            frag.chrom = component.chrom;
            frag.anchor_pos = component.anchor_pos;
            frag.read_id = read_id;
            frag.read_index = idx_in_bin;
            frag.class_mask = class_mask;
            frag.is_reverse = is_reverse;
            frag.source = InsertionFragmentSource::kClipRefLeft;
            frag.start = 0;
            frag.length = clip.leading;
            frag.read_len = read_len;
            frag.anchor_len = clip.right_anchor_after_leading;
            frag.ref_side = ReferenceSide::kRefLeft;
            frag.ref_junc_pos = read.pos();
            frag.nm = static_cast<int32_t>(nm);
            frag.sequence = read.decode_subsequence(frag.start, frag.length);
            frag.fragment_id =
                component.chrom + ":" + std::to_string(component.anchor_pos) +
                "|read=" + read_id +
                "|idx=" + std::to_string(idx_in_bin) +
                "|strand=" + std::string(is_reverse ? "-" : "+") +
                "|src=" + source_tag(frag.source) +
                "|rj=" + std::to_string(frag.ref_junc_pos) +
                "|len=" + std::to_string(frag.length);
            emit_fragment(std::move(frag));
        }

        if (clip.trailing >= min_clip) {
            InsertionFragment frag;
            frag.chrom = component.chrom;
            frag.anchor_pos = component.anchor_pos;
            frag.read_id = read_id;
            frag.read_index = idx_in_bin;
            frag.class_mask = class_mask;
            frag.is_reverse = is_reverse;
            frag.source = InsertionFragmentSource::kClipRefRight;
            frag.length = clip.trailing;
            frag.start = std::max(0, read_len - frag.length);
            frag.read_len = read_len;
            frag.anchor_len = clip.left_anchor_before_trailing;
            frag.ref_side = ReferenceSide::kRefRight;
            frag.ref_junc_pos = clip.ref_end;
            frag.nm = static_cast<int32_t>(nm);
            frag.sequence = read.decode_subsequence(frag.start, frag.length);
            frag.fragment_id =
                component.chrom + ":" + std::to_string(component.anchor_pos) +
                "|read=" + read_id +
                "|idx=" + std::to_string(idx_in_bin) +
                "|strand=" + std::string(is_reverse ? "-" : "+") +
                "|src=" + source_tag(frag.source) +
                "|rj=" + std::to_string(frag.ref_junc_pos) +
                "|len=" + std::to_string(frag.length);
            emit_fragment(std::move(frag));
        }

        if (!ins_ops.empty()) {
            std::sort(ins_ops.begin(), ins_ops.end(), [](const InsOp& a, const InsOp& b) {
                return a.len > b.len;
            });
            if (ins_ops.size() > 2) {
                ins_ops.resize(2);
            }

            for (const auto& op : ins_ops) {
                InsertionFragment frag;
                frag.chrom = component.chrom;
                frag.anchor_pos = component.anchor_pos;
                frag.read_id = read_id;
                frag.read_index = idx_in_bin;
                frag.class_mask = class_mask;
                frag.is_reverse = is_reverse;
                frag.source = InsertionFragmentSource::kCigarInsertion;
                frag.start = op.start;
                frag.length = op.len;
                frag.read_len = read_len;
                frag.anchor_len = std::min(op.left_anchor, op.right_anchor);
                frag.ref_side = ReferenceSide::kUnknown;
                frag.ref_junc_pos = op.ref_pos;
                frag.nm = static_cast<int32_t>(nm);
                frag.sequence = read.decode_subsequence(frag.start, frag.length);
                frag.fragment_id =
                    component.chrom + ":" + std::to_string(component.anchor_pos) +
                    "|read=" + read_id +
                    "|idx=" + std::to_string(idx_in_bin) +
                    "|strand=" + std::string(is_reverse ? "-" : "+") +
                    "|src=" + source_tag(frag.source) +
                    "|start=" + std::to_string(frag.start) +
                    "|rj=" + std::to_string(frag.ref_junc_pos) +
                    "|len=" + std::to_string(frag.length);
                emit_fragment(std::move(frag));
            }
        }
    }

    SplitSAFragmentModule split_sa_module(config_);
    auto sa_fragments = split_sa_module.extract(component, bin_records);
    for (auto& frag : sa_fragments) {
        emit_fragment(std::move(frag));
    }

    if (write_fasta) {
        const std::string rows = fasta_buffer.str();
        if (!rows.empty()) {
            std::lock_guard<std::mutex> lock(g_ins_fragments_fasta_mutex);
            std::ofstream fasta(config_.ins_fragments_fasta_path, std::ios::out | std::ios::app);
            if (fasta.is_open()) {
                fasta << rows;
            } else {
                std::cerr << "[InsertFragments] cannot open for appending: "
                          << config_.ins_fragments_fasta_path << "\n";
            }
        }
    }

    return out;
}

}  // namespace placer
