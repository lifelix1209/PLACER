#include "pipeline.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/sam.h>

namespace placer {
namespace {

constexpr int32_t kComponentSoftClipSignalMin = 20;
constexpr int32_t kComponentLongInsertionSignalMin = 50;
constexpr int32_t kCompatibilityWindowBp = 4000;
constexpr int32_t kStrongDbscanMinPts = 5;
constexpr int32_t kWeakDbscanMinPts = 1;
constexpr double kIntraSignatureEpsilon = 250.0;
constexpr double kInterSignatureEpsilon = 500.0;

enum class SignatureSource : uint8_t {
    kCigarInsertion = 0,
    kSplitInsertion = 1,
    kClipHint = 2
};

struct InsertionSignature {
    int32_t pos = -1;
    int32_t end = -1;
    int32_t length = 0;
    size_t read_index = 0;
    SignatureSource source = SignatureSource::kCigarInsertion;
    uint8_t class_mask = 0;
    bool is_reverse = false;
    int32_t anchor_len = 0;
    std::string read_id;
};

struct CigarStringOp {
    int32_t len = 0;
    char op = '\0';
};

struct QueryInterval {
    int32_t qstart = 0;
    int32_t qend = 0;
    int32_t leading_s = 0;
    int32_t trailing_s = 0;
};

struct SAEntry {
    std::string rname;
    int32_t pos = 0;  // 1-based
    char strand = '+';
    std::string cigar;
};

struct NormalizedAlignment {
    std::string chrom;
    int32_t ref_start = -1;
    int32_t ref_end = -1;
    int32_t qstart = -1;
    int32_t qend = -1;
    bool is_reverse = false;
};

struct SignatureBlock {
    size_t begin = 0;
    size_t end = 0;
};

bool is_match_like(int op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}

bool consumes_query_char(char op) {
    return op == 'M' || op == '=' || op == 'X' || op == 'I' || op == 'S';
}

bool consumes_ref_char(char op) {
    return op == 'M' || op == '=' || op == 'X' || op == 'D' || op == 'N';
}

int find_first_non_hard_clip_local(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

int find_last_non_hard_clip_local(const uint32_t* cigar, int32_t n_cigar) {
    for (int32_t i = n_cigar - 1; i >= 0; --i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            return i;
        }
    }
    return -1;
}

int32_t compute_ref_end_local(const ReadView& read) {
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (cigar == nullptr || n_cigar <= 0) {
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

bool parse_int32(const char* text, int32_t& out) {
    if (text == nullptr || *text == '\0') {
        return false;
    }
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (end == text || (end != nullptr && *end != '\0')) {
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

        if (fields.size() >= 4) {
            SAEntry entry;
            entry.rname = fields[0];
            parse_int32(fields[1].c_str(), entry.pos);
            entry.strand = fields[2].empty() ? '+' : fields[2][0];
            entry.cigar = fields[3];
            if (!entry.rname.empty() && entry.pos > 0 && !entry.cigar.empty()) {
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
    bool have_num = false;
    for (char c : cigar) {
        if (c >= '0' && c <= '9') {
            have_num = true;
            num = (num * 10) + (c - '0');
            continue;
        }
        if (!have_num || num <= 0) {
            return false;
        }
        ops.push_back(CigarStringOp{num, c});
        num = 0;
        have_num = false;
    }
    return !have_num && !ops.empty();
}

bool cigar_to_query_interval(
    const std::vector<CigarStringOp>& ops,
    int32_t read_len,
    QueryInterval& out,
    int32_t& ref_aligned_len) {
    out = QueryInterval{};
    ref_aligned_len = 0;
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
    }
    if (q_aligned <= 0) {
        return false;
    }

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
    out.qstart = qstart;
    out.qend = qstart + q_aligned;
    return out.qend > out.qstart && out.qend <= read_len;
}

bool bam_to_query_interval(const ReadView& read, QueryInterval& out) {
    out = QueryInterval{};
    const int32_t read_len = read.seq_len();
    if (read_len <= 0) {
        return false;
    }
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (cigar == nullptr || n_cigar <= 0) {
        return false;
    }

    const int first = find_first_non_hard_clip_local(cigar, n_cigar);
    const int last = find_last_non_hard_clip_local(cigar, n_cigar);
    if (first >= 0 && bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
        out.leading_s = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
    }
    if (last >= 0 && bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
        out.trailing_s = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
    }

    int32_t q_aligned = 0;
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

    out.qstart = out.leading_s;
    out.qend = read_len - out.trailing_s;
    if (out.qend <= out.qstart || out.qend > read_len) {
        out.qend = out.qstart + q_aligned;
    }
    return out.qend > out.qstart && out.qend <= read_len;
}

bool normalized_primary_alignment(const ReadView& read, NormalizedAlignment& out) {
    out = NormalizedAlignment{};
    QueryInterval q_interval;
    if (!bam_to_query_interval(read, q_interval)) {
        return false;
    }
    out.chrom = std::string(read.qname().empty() ? "" : "");
    out.ref_start = read.pos();
    out.ref_end = compute_ref_end_local(read);
    out.qstart = q_interval.qstart;
    out.qend = q_interval.qend;
    out.is_reverse = (read.flag() & BAM_FREVERSE) != 0;
    return out.ref_end > out.ref_start && out.qend > out.qstart;
}

bool normalized_sa_alignment(
    const SAEntry& entry,
    int32_t read_len,
    NormalizedAlignment& out) {
    out = NormalizedAlignment{};
    std::vector<CigarStringOp> ops;
    QueryInterval q_interval;
    int32_t ref_aligned_len = 0;
    if (!parse_cigar_ops(entry.cigar, ops) ||
        !cigar_to_query_interval(ops, read_len, q_interval, ref_aligned_len)) {
        return false;
    }
    out.chrom = entry.rname;
    out.ref_start = entry.pos - 1;
    out.ref_end = out.ref_start + ref_aligned_len;
    out.qstart = q_interval.qstart;
    out.qend = q_interval.qend;
    out.is_reverse = (entry.strand == '-');
    return out.ref_end > out.ref_start && out.qend > out.qstart;
}

void append_insertion_signatures_from_cigar(
    const ReadView& read,
    size_t read_index,
    const std::string& read_id,
    bool is_reverse,
    std::vector<InsertionSignature>& out) {
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (cigar == nullptr || n_cigar <= 0) {
        return;
    }

    int32_t ref_pos = read.pos();
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));
        if (op == BAM_CINS && len >= kComponentLongInsertionSignalMin) {
            InsertionSignature sign;
            sign.pos = ref_pos;
            sign.end = ref_pos + 1;
            sign.length = len;
            sign.read_index = read_index;
            sign.source = SignatureSource::kCigarInsertion;
            sign.class_mask = kCandidateLongInsertion;
            sign.is_reverse = is_reverse;
            sign.read_id = read_id;
            sign.anchor_len = len;
            out.push_back(std::move(sign));
        }
        if ((bam_cigar_type(op) & 2) != 0) {
            ref_pos += len;
        }
    }
}

bool append_best_split_signature(
    const ReadView& read,
    const std::string& chrom,
    size_t read_index,
    const std::string& read_id,
    bool is_reverse,
    std::vector<InsertionSignature>& out) {
    std::string sa_z;
    if (!read.get_string_tag("SA", sa_z) || sa_z.empty()) {
        return false;
    }

    NormalizedAlignment primary;
    if (!normalized_primary_alignment(read, primary)) {
        return false;
    }
    primary.chrom = chrom;

    int32_t best_length = -1;
    int32_t best_anchor = -1;
    for (const SAEntry& entry : parse_sa_tag_z(sa_z)) {
        if (entry.rname != chrom) {
            continue;
        }
        NormalizedAlignment mate;
        if (!normalized_sa_alignment(entry, read.seq_len(), mate)) {
            continue;
        }
        if (mate.is_reverse != primary.is_reverse) {
            continue;
        }
        const NormalizedAlignment* left = &primary;
        const NormalizedAlignment* right = &mate;
        if (mate.qstart < primary.qstart) {
            left = &mate;
            right = &primary;
        }
        const int32_t query_gap = right->qstart - left->qend;
        const int32_t ref_gap = std::max(0, right->ref_start - left->ref_end);
        const int32_t insertion_len = query_gap - ref_gap;
        if (insertion_len >= kComponentLongInsertionSignalMin &&
            insertion_len > best_length) {
            best_length = insertion_len;
            best_anchor = left->ref_end;
        }
    }

    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    const int first = (cigar != nullptr && n_cigar > 0) ? find_first_non_hard_clip_local(cigar, n_cigar) : -1;
    const int last = (cigar != nullptr && n_cigar > 0) ? find_last_non_hard_clip_local(cigar, n_cigar) : -1;
    const int32_t leading_soft =
        (first >= 0 && cigar != nullptr && bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP)
        ? static_cast<int32_t>(bam_cigar_oplen(cigar[first]))
        : 0;
    const int32_t trailing_soft =
        (last >= 0 && cigar != nullptr && bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP)
        ? static_cast<int32_t>(bam_cigar_oplen(cigar[last]))
        : 0;

    if (best_length < 0) {
        if (leading_soft >= trailing_soft && leading_soft >= kComponentSoftClipSignalMin) {
            best_length = leading_soft;
            best_anchor = read.pos();
        } else if (trailing_soft >= kComponentSoftClipSignalMin) {
            best_length = trailing_soft;
            best_anchor = compute_ref_end_local(read);
        }
    }

    if (best_length < 0 || best_anchor < 0) {
        return false;
    }

    InsertionSignature sign;
    sign.pos = best_anchor;
    sign.end = best_anchor + 1;
    sign.length = best_length;
    sign.read_index = read_index;
    sign.source = SignatureSource::kSplitInsertion;
    sign.class_mask = kCandidateSplitSaSupplementary;
    sign.is_reverse = is_reverse;
    sign.read_id = read_id;
    sign.anchor_len = best_length;
    out.push_back(std::move(sign));
    return true;
}

void append_clip_hint_signatures(
    const ReadView& read,
    size_t read_index,
    const std::string& read_id,
    bool is_reverse,
    std::vector<InsertionSignature>& out) {
    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (cigar == nullptr || n_cigar <= 0) {
        return;
    }

    const int first = find_first_non_hard_clip_local(cigar, n_cigar);
    const int last = find_last_non_hard_clip_local(cigar, n_cigar);
    if (first >= 0 && bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
        if (len >= kComponentSoftClipSignalMin) {
            InsertionSignature sign;
            sign.pos = read.pos();
            sign.end = sign.pos + 1;
            sign.length = len;
            sign.read_index = read_index;
            sign.source = SignatureSource::kClipHint;
            sign.class_mask = kCandidateSoftClip;
            sign.is_reverse = is_reverse;
            sign.read_id = read_id;
            sign.anchor_len = len;
            out.push_back(std::move(sign));
        }
    }
    if (last >= 0 && bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
        if (len >= kComponentSoftClipSignalMin) {
            InsertionSignature sign;
            sign.pos = compute_ref_end_local(read);
            sign.end = sign.pos + 1;
            sign.length = len;
            sign.read_index = read_index;
            sign.source = SignatureSource::kClipHint;
            sign.class_mask = kCandidateSoftClip;
            sign.is_reverse = is_reverse;
            sign.read_id = read_id;
            sign.anchor_len = len;
            out.push_back(std::move(sign));
        }
    }
}

std::vector<InsertionSignature> extract_signatures(
    const std::vector<const bam1_t*>& records,
    const std::string& chrom,
    int32_t tid) {
    std::vector<InsertionSignature> out;
    for (size_t read_idx = 0; read_idx < records.size(); ++read_idx) {
        const bam1_t* record = records[read_idx];
        if (record == nullptr) {
            continue;
        }

        ReadView read(record);
        if (read.tid() != tid || (read.flag() & BAM_FSUPPLEMENTARY) != 0) {
            continue;
        }
        const std::string read_id(read.qname());
        const bool is_reverse = (read.flag() & BAM_FREVERSE) != 0;

        append_insertion_signatures_from_cigar(read, read_idx, read_id, is_reverse, out);
        const bool emitted_split =
            read.has_sa_tag() && append_best_split_signature(read, chrom, read_idx, read_id, is_reverse, out);
        if (!emitted_split) {
            append_clip_hint_signatures(read, read_idx, read_id, is_reverse, out);
        }
    }
    std::sort(out.begin(), out.end(), [](const InsertionSignature& lhs, const InsertionSignature& rhs) {
        if (lhs.pos != rhs.pos) {
            return lhs.pos < rhs.pos;
        }
        if (lhs.end != rhs.end) {
            return lhs.end < rhs.end;
        }
        if (lhs.length != rhs.length) {
            return lhs.length < rhs.length;
        }
        return lhs.read_index < rhs.read_index;
    });
    return out;
}

std::vector<SignatureBlock> build_blocks(const std::vector<InsertionSignature>& signs) {
    std::vector<SignatureBlock> blocks;
    if (signs.empty()) {
        return blocks;
    }

    size_t begin = 0;
    for (size_t i = 1; i < signs.size(); ++i) {
        if ((signs[i].pos - signs[i - 1].pos) > kCompatibilityWindowBp) {
            blocks.push_back(SignatureBlock{begin, i});
            begin = i;
        }
    }
    blocks.push_back(SignatureBlock{begin, signs.size()});
    return blocks;
}

double signature_distance(const InsertionSignature& lhs, const InsertionSignature& rhs) {
    const double dx = static_cast<double>(lhs.pos - rhs.pos);
    const double dy = static_cast<double>(lhs.end - rhs.end);
    const double dz = static_cast<double>(lhs.length - rhs.length);
    return std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

double signature_epsilon(const InsertionSignature& lhs, const InsertionSignature& rhs) {
    const bool both_intra =
        lhs.source == SignatureSource::kCigarInsertion &&
        rhs.source == SignatureSource::kCigarInsertion;
    return both_intra ? kIntraSignatureEpsilon : kInterSignatureEpsilon;
}

std::vector<std::vector<size_t>> dbscan_cluster_indices(
    const std::vector<InsertionSignature>& signs,
    const std::vector<size_t>& candidate_indices,
    int32_t min_pts) {
    std::vector<std::vector<size_t>> adjacency(candidate_indices.size());
    for (size_t i = 0; i < candidate_indices.size(); ++i) {
        for (size_t j = i + 1; j < candidate_indices.size(); ++j) {
            const auto& lhs = signs[candidate_indices[i]];
            const auto& rhs = signs[candidate_indices[j]];
            if (signature_distance(lhs, rhs) < signature_epsilon(lhs, rhs)) {
                adjacency[i].push_back(j);
                adjacency[j].push_back(i);
            }
        }
    }

    std::vector<int32_t> labels(candidate_indices.size(), -1);
    int32_t next_label = 0;

    for (size_t i = 0; i < candidate_indices.size(); ++i) {
        if (labels[i] != -1) {
            continue;
        }
        if (static_cast<int32_t>(adjacency[i].size()) + 1 < min_pts) {
            labels[i] = 0;
            continue;
        }

        ++next_label;
        labels[i] = next_label;
        std::vector<size_t> stack = adjacency[i];
        while (!stack.empty()) {
            const size_t j = stack.back();
            stack.pop_back();
            if (labels[j] == next_label) {
                continue;
            }
            const bool was_unvisited = (labels[j] == -1);
            labels[j] = next_label;
            if (was_unvisited &&
                static_cast<int32_t>(adjacency[j].size()) + 1 >= min_pts) {
                stack.insert(stack.end(), adjacency[j].begin(), adjacency[j].end());
            }
        }
    }

    std::vector<std::vector<size_t>> clusters(static_cast<size_t>(next_label));
    for (size_t i = 0; i < candidate_indices.size(); ++i) {
        if (labels[i] > 0) {
            clusters[static_cast<size_t>(labels[i] - 1)].push_back(candidate_indices[i]);
        }
    }
    return clusters;
}

std::vector<size_t> collect_noise_indices(
    const std::vector<InsertionSignature>& signs,
    const std::vector<size_t>& candidate_indices,
    const std::vector<std::vector<size_t>>& clusters) {
    (void)signs;
    std::unordered_set<size_t> clustered;
    for (const auto& cluster : clusters) {
        for (size_t idx : cluster) {
            clustered.insert(idx);
        }
    }
    std::vector<size_t> noise;
    for (size_t idx : candidate_indices) {
        if (clustered.find(idx) == clustered.end()) {
            noise.push_back(idx);
        }
    }
    return noise;
}

int32_t median_i32(std::vector<int32_t> values) {
    if (values.empty()) {
        return -1;
    }
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

ComponentCall project_cluster(
    const std::vector<InsertionSignature>& signs,
    const std::vector<size_t>& cluster_members,
    const std::string& chrom,
    int32_t tid) {
    ComponentCall call;
    call.chrom = chrom;
    call.tid = tid;
    call.bin_start = std::numeric_limits<int32_t>::max();
    call.bin_end = std::numeric_limits<int32_t>::min();

    std::vector<int32_t> anchors;
    std::unordered_set<size_t> all_reads;
    std::unordered_set<size_t> clip_reads;
    std::unordered_set<size_t> split_reads;
    std::unordered_set<size_t> insertion_reads;

    for (size_t idx : cluster_members) {
        const auto& sign = signs[idx];
        call.bin_start = std::min(call.bin_start, sign.pos);
        call.bin_end = std::max(call.bin_end, sign.end);
        anchors.push_back(sign.pos);
        all_reads.insert(sign.read_index);

        BreakpointCandidate bp;
        bp.chrom = chrom;
        bp.pos = sign.pos;
        bp.read_id = sign.read_id;
        bp.read_index = sign.read_index;
        bp.class_mask = sign.class_mask;
        bp.is_reverse = sign.is_reverse;
        bp.anchor_len = sign.anchor_len;
        switch (sign.source) {
            case SignatureSource::kCigarInsertion:
                insertion_reads.insert(sign.read_index);
                call.evidence_indel_count += 1;
                bp.ins_len = sign.length;
                break;
            case SignatureSource::kSplitInsertion:
                split_reads.insert(sign.read_index);
                call.evidence_sa_hint_count += 1;
                bp.ins_len = sign.length;
                break;
            case SignatureSource::kClipHint:
                clip_reads.insert(sign.read_index);
                call.evidence_soft_clip_count += 1;
                bp.clip_len = sign.length;
                break;
        }
        call.breakpoint_candidates.push_back(std::move(bp));
    }

    call.anchor_pos = median_i32(std::move(anchors));
    call.peak_weight = static_cast<double>(all_reads.size());

    call.read_indices.assign(all_reads.begin(), all_reads.end());
    call.soft_clip_read_indices.assign(clip_reads.begin(), clip_reads.end());
    call.split_sa_read_indices.assign(split_reads.begin(), split_reads.end());
    call.insertion_read_indices.assign(insertion_reads.begin(), insertion_reads.end());
    std::sort(call.read_indices.begin(), call.read_indices.end());
    std::sort(call.soft_clip_read_indices.begin(), call.soft_clip_read_indices.end());
    std::sort(call.split_sa_read_indices.begin(), call.split_sa_read_indices.end());
    std::sort(call.insertion_read_indices.begin(), call.insertion_read_indices.end());

    if (call.bin_start == std::numeric_limits<int32_t>::max()) {
        call.bin_start = call.anchor_pos;
    }
    if (call.bin_end == std::numeric_limits<int32_t>::min()) {
        call.bin_end = call.anchor_pos + 1;
    }
    return call;
}

}  // namespace

std::vector<ComponentCall> DbscanComponentModule::build(
    const std::vector<const bam1_t*>& records,
    const std::string& chrom,
    int32_t tid) const {
    const std::vector<InsertionSignature> signs = extract_signatures(records, chrom, tid);
    std::vector<ComponentCall> out;
    for (const SignatureBlock& block : build_blocks(signs)) {
        std::vector<size_t> indices;
        indices.reserve(block.end - block.begin);
        for (size_t idx = block.begin; idx < block.end; ++idx) {
            indices.push_back(idx);
        }

        std::vector<std::vector<size_t>> clusters =
            dbscan_cluster_indices(signs, indices, kStrongDbscanMinPts);
        std::vector<size_t> noise = collect_noise_indices(signs, indices, clusters);
        if (!noise.empty()) {
            auto weak_clusters = dbscan_cluster_indices(signs, noise, kWeakDbscanMinPts);
            clusters.insert(clusters.end(), weak_clusters.begin(), weak_clusters.end());
        }

        for (const auto& cluster : clusters) {
            if (cluster.empty()) {
                continue;
            }
            out.push_back(project_cluster(signs, cluster, chrom, tid));
        }
    }

    std::sort(out.begin(), out.end(), [](const ComponentCall& lhs, const ComponentCall& rhs) {
        if (lhs.anchor_pos != rhs.anchor_pos) {
            return lhs.anchor_pos < rhs.anchor_pos;
        }
        if (lhs.peak_weight != rhs.peak_weight) {
            return lhs.peak_weight > rhs.peak_weight;
        }
        return lhs.read_indices.size() > rhs.read_indices.size();
    });
    return out;
}

}  // namespace placer
