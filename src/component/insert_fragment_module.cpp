#include "pipeline.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <htslib/sam.h>

namespace placer {
namespace {

struct ClipInfo {
    int32_t leading = 0;
    int32_t trailing = 0;
};

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

ClipInfo get_soft_clip_ends(const ReadView& read) {
    ClipInfo info;

    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return info;
    }

    const int first = find_first_non_hard_clip(cigar, n_cigar);
    const int last = find_last_non_hard_clip(cigar, n_cigar);
    if (first < 0 || last < 0 || first > last) {
        return info;
    }

    if (bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
        info.leading = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
    }
    if (bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
        info.trailing = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
    }

    return info;
}

struct InsOp {
    int32_t start = 0;
    int32_t len = 0;
};

struct SAEntry {
    std::string rname;
    int32_t pos = 0;
    char strand = '+';
    std::string cigar;
    int32_t mapq = 0;
    int32_t nm = 0;
};

struct QueryInterval {
    int32_t qstart = 0;
    int32_t qlen = 0;
    int32_t leading_s = 0;
    int32_t trailing_s = 0;
};

std::vector<InsOp> find_long_insertions(const ReadView& read, int32_t min_long_ins) {
    std::vector<InsOp> ops;

    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return ops;
    }

    int32_t qpos = 0;
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));

        if (op == BAM_CINS && len >= min_long_ins) {
            ops.push_back({qpos, len});
        }

        if ((bam_cigar_type(op) & 1) != 0) {
            qpos += len;
        }
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
        const size_t rec_len = (rec_end == std::string::npos) ? (sa_z.size() - start) : (rec_end - start);
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

bool cigar_to_query_interval_mvp(const std::string& cigar, QueryInterval& out) {
    out = QueryInterval{};
    if (cigar.empty()) {
        return false;
    }

    struct Op {
        int32_t len = 0;
        char op = '\0';
    };

    std::vector<Op> ops;
    ops.reserve(16);
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
    if (ops.empty()) {
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

    // Without any soft clip in SA cigar, query start cannot be inferred in this MVP.
    if (out.leading_s == 0 && out.trailing_s == 0) {
        return false;
    }

    int32_t q_aligned = 0;
    for (const auto& op : ops) {
        switch (op.op) {
            case 'M':
            case '=':
            case 'X':
            case 'I':
                q_aligned += op.len;
                break;
            default:
                break;
        }
    }

    out.qstart = out.leading_s;
    out.qlen = q_aligned;
    return out.qlen > 0;
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

    const auto& indices = component.read_indices;
    for (size_t idx_in_bin : indices) {
        if (idx_in_bin >= bin_records.size() || !bin_records[idx_in_bin]) {
            continue;
        }

        ReadView read(bin_records[idx_in_bin].get());
        const uint16_t flag = read.flag();

        // Parse SA only on primary records to avoid duplicating fragments.
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

        int32_t emitted = 0;
        int32_t sa_index = 0;
        for (const auto& entry : entries) {
            if (emitted >= max_sa) {
                break;
            }
            ++sa_index;

            QueryInterval q_interval;
            if (!cigar_to_query_interval_mvp(entry.cigar, q_interval)) {
                continue;
            }
            if (q_interval.qlen < min_len) {
                continue;
            }
            if (q_interval.qstart < 0 || q_interval.qstart >= read_len) {
                continue;
            }
            if (q_interval.qstart + q_interval.qlen > read_len) {
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
            frag.start = q_interval.qstart;
            frag.length = q_interval.qlen;
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

            if (frag.sequence.empty()) {
                continue;
            }
            out.push_back(std::move(frag));
            ++emitted;
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

        const ClipInfo clips = get_soft_clip_ends(read);
        if (clips.leading >= min_clip || clips.trailing >= min_clip) {
            class_mask |= kCandidateSoftClip;
        }

        auto ins_ops = find_long_insertions(read, min_ins);
        if (!ins_ops.empty()) {
            class_mask |= kCandidateLongInsertion;
        }

        // Soft-clip fragments
        if (clips.leading >= min_clip) {
            InsertionFragment frag;
            frag.chrom = component.chrom;
            frag.anchor_pos = component.anchor_pos;
            frag.read_id = read_id;
            frag.read_index = idx_in_bin;
            frag.class_mask = class_mask;
            frag.is_reverse = is_reverse;
            frag.source = InsertionFragmentSource::kClipRefLeft;
            frag.start = 0;
            frag.length = clips.leading;
            frag.sequence = read.decode_subsequence(frag.start, frag.length);
            frag.fragment_id =
                component.chrom + ":" + std::to_string(component.anchor_pos) +
                "|read=" + read_id +
                "|idx=" + std::to_string(idx_in_bin) +
                "|strand=" + std::string(is_reverse ? "-" : "+") +
                "|src=" + source_tag(frag.source) +
                "|len=" + std::to_string(frag.length);
            emit_fragment(std::move(frag));
        }

        if (clips.trailing >= min_clip) {
            InsertionFragment frag;
            frag.chrom = component.chrom;
            frag.anchor_pos = component.anchor_pos;
            frag.read_id = read_id;
            frag.read_index = idx_in_bin;
            frag.class_mask = class_mask;
            frag.is_reverse = is_reverse;
            frag.source = InsertionFragmentSource::kClipRefRight;
            frag.length = clips.trailing;
            frag.start = std::max(0, read.seq_len() - frag.length);
            frag.sequence = read.decode_subsequence(frag.start, frag.length);
            frag.fragment_id =
                component.chrom + ":" + std::to_string(component.anchor_pos) +
                "|read=" + read_id +
                "|idx=" + std::to_string(idx_in_bin) +
                "|strand=" + std::string(is_reverse ? "-" : "+") +
                "|src=" + source_tag(frag.source) +
                "|len=" + std::to_string(frag.length);
            emit_fragment(std::move(frag));
        }

        // Long insertion fragments - keep top 2 by length
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
                frag.sequence = read.decode_subsequence(frag.start, frag.length);
                frag.fragment_id =
                    component.chrom + ":" + std::to_string(component.anchor_pos) +
                    "|read=" + read_id +
                    "|idx=" + std::to_string(idx_in_bin) +
                    "|strand=" + std::string(is_reverse ? "-" : "+") +
                    "|src=" + source_tag(frag.source) +
                    "|start=" + std::to_string(frag.start) +
                    "|len=" + std::to_string(frag.length);
                emit_fragment(std::move(frag));
            }
        }
    }

    // Include split-SA-derived query segments as additional fragments.
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
