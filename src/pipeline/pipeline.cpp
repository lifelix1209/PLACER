#include "pipeline.h"
#include "decision_policy.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/sam.h>
#include <abpoa.h>

namespace placer {
namespace {

constexpr int32_t kSoftClipSignalMin = 20;
constexpr int32_t kLongInsertionSignalMin = 50;
constexpr int32_t kLargeIndelEvidenceMin = 40;
constexpr double kSaHintWeight = 0.35;
constexpr int32_t kHistBinSize = 20;
constexpr int32_t kHistPadding = 120;
constexpr double kPeakMinWeight = 1.2;
constexpr int32_t kPeakMergeDistance = 100;
constexpr int32_t kWindowMergeGap = 60;
// Allow component windows to cross bin boundaries so insertion evidence carried
// by long reads (whose starts may fall in adjacent bins) is not dropped.
constexpr int32_t kWindowBinSlackBp = 5000;
// Merge current bin with a small number of previous bins so split support from
// long reads can be assembled in one component.
constexpr int32_t kCrossBinContextBins = 2;
constexpr int32_t kWindowMaxExpandBins = 15;
constexpr double kWindowDropRatio = 0.50;
constexpr int32_t kMinWindowSpan = 120;
constexpr double kReadAssignMinScore = 0.8;
constexpr double kReadAmbiguousRatio = 0.75;
constexpr double kReadAssignDecayBp = 150.0;
constexpr int32_t kFinalCallDedupDistanceBp = 100;
constexpr int32_t kLocalEventFetchSlackBp = 1000;
constexpr int32_t kEventConsensusFlankBp = 80;
constexpr int32_t kEventConsensusMinAnchorBp = 50;
constexpr int32_t kEventSegmentationMinFlankAlignBp = 50;
constexpr double kEventSegmentationMinFlankIdentity = 0.90;
constexpr int32_t kEventSegmentationBreakpointSlackBp = 200;
constexpr int32_t kEventSegmentationEndpointSlackBp = 32;
constexpr int32_t kEventSegmentationMinInsertBp = 1;
constexpr int32_t kEventSegmentationMaxFlankQueryBp = 120;
constexpr double kEventSegmentationUniquenessMargin = 0.02;

enum class EvidenceKind : uint8_t {
    kSoftClip = 0,
    kIndel = 1,
    kSAHint = 2
};

struct EvidencePoint {
    size_t read_index = 0;
    int32_t pos = -1;
    double weight = 0.0;
    EvidenceKind kind = EvidenceKind::kSoftClip;
    int32_t signal_len = 0;
    uint8_t class_mask = 0;
};

struct ReadSignalSummary {
    std::string read_id;
    bool is_reverse = false;
    bool has_sa_or_supp = false;
    int32_t max_soft_clip = 0;
    int32_t max_ins = 0;
    uint8_t class_mask = 0;
};

struct EvidenceBundle {
    std::vector<EvidencePoint> points;
    std::vector<ReadSignalSummary> read_summaries;
};

struct ReadWindowAssignmentScore {
    bool valid = false;
    int32_t specificity_rank = std::numeric_limits<int32_t>::max();
    double specificity_score = 0.0;
    double total_score = 0.0;
};

struct LocalFetchedReads {
    std::vector<BamRecordPtr> owned_records;
    std::vector<const bam1_t*> records;
    std::vector<ReadReferenceSpan> read_spans;
};

ComponentCall build_local_fragment_component(
    const PipelineConfig& config,
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records) {
    ComponentCall local_component = component;
    local_component.read_indices.clear();
    local_component.soft_clip_read_indices.clear();
    local_component.split_sa_read_indices.clear();
    local_component.insertion_read_indices.clear();

    local_component.read_indices.reserve(local_records.size());
    local_component.soft_clip_read_indices.reserve(local_records.size());
    local_component.split_sa_read_indices.reserve(local_records.size());
    local_component.insertion_read_indices.reserve(local_records.size());

    for (size_t idx = 0; idx < local_records.size(); ++idx) {
        const bam1_t* record = local_records[idx];
        if (!record) {
            continue;
        }

        ReadView read(record);
        if (read.tid() != component.tid) {
            continue;
        }

        local_component.read_indices.push_back(idx);

        const uint16_t flag = read.flag();
        if (read.has_sa_tag() || ((flag & BAM_FSUPPLEMENTARY) != 0)) {
            local_component.split_sa_read_indices.push_back(idx);
        }

        const uint32_t* cigar = read.cigar();
        const int32_t n_cigar = read.n_cigar();
        if (!cigar || n_cigar <= 0) {
            continue;
        }

        int32_t max_soft_clip = 0;
        int32_t max_insertion = 0;
        for (int32_t ci = 0; ci < n_cigar; ++ci) {
            const int op = bam_cigar_op(cigar[ci]);
            const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[ci]));
            if (op == BAM_CSOFT_CLIP) {
                max_soft_clip = std::max(max_soft_clip, len);
            } else if (op == BAM_CINS) {
                max_insertion = std::max(max_insertion, len);
            }
        }

        if (max_soft_clip >= config.min_soft_clip_for_seq_extract) {
            local_component.soft_clip_read_indices.push_back(idx);
        }
        if (max_insertion >= config.min_long_ins_for_seq_extract) {
            local_component.insertion_read_indices.push_back(idx);
        }
    }

    return local_component;
}

struct DensityPeak {
    int32_t pos = -1;
    double weight = 0.0;
};

struct CandidateWindow {
    int32_t start = -1;
    int32_t end = -1;
    int32_t center = -1;
    double peak_weight = 0.0;
};

struct EventFlankPlacement {
    int32_t query_start = -1;
    int32_t query_end = -1;
    int32_t ref_start = -1;
    int32_t ref_end = -1;
    int32_t align_len = 0;
    int32_t breakpoint_delta = std::numeric_limits<int32_t>::max();
    int32_t endpoint_offset = std::numeric_limits<int32_t>::max();
    double identity = 0.0;
};

struct EventFlankSearchResult {
    std::vector<EventFlankPlacement> candidates;
};

enum class SaHintSide : uint8_t {
    kNone = 0,
    kLeft = 1,
    kRight = 2
};

SaHintSide choose_sa_hint_side(int32_t leading_soft, int32_t trailing_soft) {
    const bool left = leading_soft >= kSoftClipSignalMin;
    const bool right = trailing_soft >= kSoftClipSignalMin;
    if (left && !right) {
        return SaHintSide::kLeft;
    }
    if (right && !left) {
        return SaHintSide::kRight;
    }
    if (left && right) {
        if (leading_soft > trailing_soft) {
            return SaHintSide::kLeft;
        }
        if (trailing_soft > leading_soft) {
            return SaHintSide::kRight;
        }
    }
    return SaHintSide::kNone;
}

void append_record_ptrs(
    const std::vector<BufferedRecord>& src,
    std::vector<const bam1_t*>& dst,
    int32_t min_ref_end) {
    dst.reserve(dst.size() + src.size());
    for (const auto& rec : src) {
        if (!rec.record) {
            continue;
        }
        if (rec.ref_end < min_ref_end) {
            continue;
        }
        dst.push_back(rec.record.get());
    }
}

int bool_as_int(bool value) {
    return value ? 1 : 0;
}

const char* final_hypothesis_kind_name(FinalHypothesisKind kind) {
    switch (kind) {
        case FinalHypothesisKind::kReference:
            return "REFERENCE";
        case FinalHypothesisKind::kInsertionNonTe:
            return "NON_TE_INSERTION";
        case FinalHypothesisKind::kTeUnknown:
            return "TE_UNKNOWN";
        case FinalHypothesisKind::kTeResolved:
            return "TE_RESOLVED";
    }
    return "UNKNOWN";
}

int32_t evidence_specificity_rank(EvidenceKind kind) {
    switch (kind) {
        case EvidenceKind::kIndel:
            return 0;
        case EvidenceKind::kSoftClip:
            return 1;
        case EvidenceKind::kSAHint:
            return 2;
    }
    return std::numeric_limits<int32_t>::max();
}

void accumulate_read_window_assignment_score(
    const EvidencePoint& point,
    double point_score,
    ReadWindowAssignmentScore& score) {
    score.valid = true;
    score.total_score += point_score;

    const int32_t specificity_rank = evidence_specificity_rank(point.kind);
    if (specificity_rank < score.specificity_rank) {
        score.specificity_rank = specificity_rank;
        score.specificity_score = point_score;
        return;
    }
    if (specificity_rank == score.specificity_rank &&
        point_score > score.specificity_score) {
        score.specificity_score = point_score;
    }
}

bool better_read_window_assignment_score(
    const ReadWindowAssignmentScore& lhs,
    const ReadWindowAssignmentScore& rhs) {
    if (!lhs.valid) {
        return false;
    }
    if (!rhs.valid) {
        return true;
    }
    if (lhs.specificity_rank != rhs.specificity_rank) {
        return lhs.specificity_rank < rhs.specificity_rank;
    }
    if (lhs.specificity_score != rhs.specificity_score) {
        return lhs.specificity_score > rhs.specificity_score;
    }
    return lhs.total_score > rhs.total_score;
}

bool better_event_flank_placement(
    const EventFlankPlacement& lhs,
    const EventFlankPlacement& rhs) {
    if (lhs.identity != rhs.identity) {
        return lhs.identity > rhs.identity;
    }
    if (lhs.breakpoint_delta != rhs.breakpoint_delta) {
        return lhs.breakpoint_delta < rhs.breakpoint_delta;
    }
    if (lhs.endpoint_offset != rhs.endpoint_offset) {
        return lhs.endpoint_offset < rhs.endpoint_offset;
    }
    if (lhs.align_len != rhs.align_len) {
        return lhs.align_len > rhs.align_len;
    }
    if (lhs.ref_start != rhs.ref_start) {
        return lhs.ref_start < rhs.ref_start;
    }
    return lhs.ref_end < rhs.ref_end;
}

std::string display_or_na(const std::string& value) {
    return value.empty() ? "NA" : value;
}

std::string summarize_seq_ends(const std::string& seq, size_t keep = 24) {
    if (seq.empty()) {
        return "NA";
    }
    if (seq.size() <= (2 * keep)) {
        return seq;
    }
    return seq.substr(0, keep) + "..." + seq.substr(seq.size() - keep);
}

bool is_primary_alignment(const bam1_t* record) {
    if (!record) {
        return false;
    }
    const uint16_t flag = record->core.flag;
    return (flag & BAM_FUNMAP) == 0 &&
           (flag & BAM_FSECONDARY) == 0 &&
           (flag & BAM_FSUPPLEMENTARY) == 0;
}

std::string record_qname(const bam1_t* record) {
    if (!record) {
        return {};
    }
    ReadView view(record);
    return std::string(view.qname());
}

struct LocalEventSignal {
    bool split = false;
    bool indel = false;
    bool left_clip = false;
    bool right_clip = false;
    int32_t split_left_pos = -1;
    int32_t split_right_pos = -1;
    int32_t indel_pos = -1;
    int32_t left_clip_pos = -1;
    int32_t right_clip_pos = -1;

    bool any() const {
        return split || indel || left_clip || right_clip;
    }
};

LocalEventSignal classify_local_event_signal(
    const bam1_t* record,
    int32_t window_start,
    int32_t window_end) {
    LocalEventSignal signal;
    if (!record) {
        return signal;
    }

    ReadView read(record);
    const bool has_sa_or_supp =
        read.has_sa_tag() || ((read.flag() & BAM_FSUPPLEMENTARY) != 0);

    const uint32_t* cigar = read.cigar();
    const int32_t n_cigar = read.n_cigar();
    if (!cigar || n_cigar <= 0) {
        return signal;
    }

    int first = -1;
    int last = -1;
    for (int32_t i = 0; i < n_cigar; ++i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            first = i;
            break;
        }
    }
    for (int32_t i = n_cigar - 1; i >= 0; --i) {
        if (bam_cigar_op(cigar[i]) != BAM_CHARD_CLIP) {
            last = i;
            break;
        }
    }

    int32_t ref_pos = read.pos();
    const int32_t window_center = window_start + ((window_end - window_start) / 2);
    int32_t best_indel_dist = std::numeric_limits<int32_t>::max();
    for (int32_t i = 0; i < n_cigar; ++i) {
        const int op = bam_cigar_op(cigar[i]);
        const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[i]));

        if (i == first && op == BAM_CSOFT_CLIP && len >= kSoftClipSignalMin &&
            ref_pos >= window_start && ref_pos <= window_end) {
            signal.left_clip = true;
            signal.left_clip_pos = ref_pos;
        }
        if (i == last && op == BAM_CSOFT_CLIP && len >= kSoftClipSignalMin &&
            ref_pos >= window_start && ref_pos <= window_end) {
            signal.right_clip = true;
            signal.right_clip_pos = ref_pos;
        }
        if (op == BAM_CINS && len >= kLongInsertionSignalMin &&
            ref_pos >= window_start && ref_pos <= window_end) {
            signal.indel = true;
            const int32_t dist = std::abs(ref_pos - window_center);
            if (dist < best_indel_dist) {
                best_indel_dist = dist;
                signal.indel_pos = ref_pos;
            }
        }

        if ((bam_cigar_type(op) & 2) != 0) {
            ref_pos += len;
        }
    }

    signal.split = has_sa_or_supp && (signal.left_clip || signal.right_clip);
    if (signal.split) {
        signal.split_left_pos = signal.left_clip_pos;
        signal.split_right_pos = signal.right_clip_pos;
    }
    return signal;
}

bool read_has_local_event_signal(
    const bam1_t* record,
    int32_t window_start,
    int32_t window_end) {
    return classify_local_event_signal(record, window_start, window_end).any();
}

std::vector<std::string> sorted_string_set(const std::unordered_set<std::string>& values) {
    std::vector<std::string> out;
    out.reserve(values.size());
    for (const auto& value : values) {
        out.push_back(value);
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::string upper_acgt(const std::string& s);

std::string collect_aligned_query_bases_before(
    const ReadView& read,
    int32_t query_limit,
    int32_t max_bases) {
    if (query_limit <= 0 || max_bases <= 0) {
        return {};
    }
    const int32_t take_start = std::max(0, query_limit - max_bases);
    const int32_t take_len = query_limit - take_start;
    return upper_acgt(read.decode_subsequence(take_start, take_len));
}

std::string collect_aligned_query_bases_after(
    const ReadView& read,
    int32_t query_start,
    int32_t max_bases) {
    if (query_start < 0 || max_bases <= 0) {
        return {};
    }
    return upper_acgt(read.decode_subsequence(query_start, max_bases));
}

std::string build_event_string_from_fragment(
    const bam1_t* record,
    const InsertionFragment& fragment) {
    if (!record || fragment.sequence.empty() || fragment.start < 0 || fragment.length <= 0) {
        return {};
    }

    ReadView read(record);
    const int32_t insert_start = fragment.start;
    const int32_t insert_end = fragment.start + fragment.length;
    if (insert_end > read.seq_len()) {
        return {};
    }

    const std::string insert_seq = upper_acgt(fragment.sequence);
    if (insert_seq.empty()) {
        return {};
    }

    const auto fetch_left_flank = [&](int32_t flank_len) {
        if (flank_len <= 0) {
            return std::string();
        }
        return collect_aligned_query_bases_before(read, insert_start, flank_len);
    };
    const auto fetch_right_flank = [&](int32_t flank_len) {
        if (flank_len <= 0) {
            return std::string();
        }
        return collect_aligned_query_bases_after(read, insert_end, flank_len);
    };

    switch (fragment.source) {
        case InsertionFragmentSource::kClipRefLeft: {
            const std::string right_flank = fetch_right_flank(kEventConsensusFlankBp);
            if (static_cast<int32_t>(right_flank.size()) < kEventConsensusMinAnchorBp) {
                return {};
            }
            return insert_seq + right_flank;
        }
        case InsertionFragmentSource::kClipRefRight: {
            const std::string left_flank = fetch_left_flank(kEventConsensusFlankBp);
            if (static_cast<int32_t>(left_flank.size()) < kEventConsensusMinAnchorBp) {
                return {};
            }
            return left_flank + insert_seq;
        }
        case InsertionFragmentSource::kCigarInsertion:
        case InsertionFragmentSource::kSplitSa: {
            const std::string left_flank = fetch_left_flank(kEventConsensusFlankBp);
            const std::string right_flank = fetch_right_flank(kEventConsensusFlankBp);
            if (static_cast<int32_t>(left_flank.size()) < kEventConsensusMinAnchorBp ||
                static_cast<int32_t>(right_flank.size()) < kEventConsensusMinAnchorBp) {
                return {};
            }
            return left_flank + insert_seq + right_flank;
        }
        case InsertionFragmentSource::kUnknown:
        default:
            return {};
    }
}

bool fragment_supports_event_consensus(
    const InsertionFragment& fragment) {
    switch (fragment.source) {
        case InsertionFragmentSource::kClipRefLeft:
        case InsertionFragmentSource::kClipRefRight:
        case InsertionFragmentSource::kCigarInsertion:
        case InsertionFragmentSource::kSplitSa:
            return true;
        case InsertionFragmentSource::kUnknown:
        default:
            return false;
    }
}

bool fragment_has_full_event_context(
    const InsertionFragment& fragment) {
    switch (fragment.source) {
        case InsertionFragmentSource::kCigarInsertion:
        case InsertionFragmentSource::kSplitSa:
            return true;
        case InsertionFragmentSource::kClipRefLeft:
        case InsertionFragmentSource::kClipRefRight:
        case InsertionFragmentSource::kUnknown:
        default:
            return false;
    }
}

bool fragment_is_local_to_event(
    const InsertionFragment& fragment,
    int32_t bp_left,
    int32_t bp_right) {
    if (fragment.ref_junc_pos < 0 || bp_left < 0 || bp_right < 0) {
        return false;
    }
    constexpr int32_t kEventFragmentJunctionSlackBp = 25;
    const int32_t span_start = std::max(0, std::min(bp_left, bp_right) - kEventFragmentJunctionSlackBp);
    const int32_t span_end = std::max(span_start, std::max(bp_left, bp_right) + kEventFragmentJunctionSlackBp);
    return fragment.ref_junc_pos >= span_start && fragment.ref_junc_pos <= span_end;
}

void emit_pipeline_log_line(const std::string& line) {
    static std::mutex log_mutex;
    const std::lock_guard<std::mutex> lock(log_mutex);
    std::cerr << line << '\n';
}

using SharedRecordBatch = std::shared_ptr<const std::vector<BufferedRecord>>;

bool final_call_sort_less(const FinalCall& a, const FinalCall& b) {
    if (a.tid != b.tid) {
        return a.tid < b.tid;
    }
    if (a.pos != b.pos) {
        return a.pos < b.pos;
    }
    if (a.window_start != b.window_start) {
        return a.window_start < b.window_start;
    }
    if (a.window_end != b.window_end) {
        return a.window_end < b.window_end;
    }
    if (a.chrom != b.chrom) {
        return a.chrom < b.chrom;
    }
    return a.te_name < b.te_name;
}

bool same_call_locus(const FinalCall& a, const FinalCall& b) {
    if (a.tid != b.tid) {
        return false;
    }
    return std::abs(a.pos - b.pos) <= kFinalCallDedupDistanceBp;
}

bool prefer_new_call(const FinalCall& cur, const FinalCall& incumbent) {
    if (cur.support_reads != incumbent.support_reads) {
        return cur.support_reads > incumbent.support_reads;
    }
    if (cur.gq != incumbent.gq) {
        return cur.gq > incumbent.gq;
    }
    if (cur.cross_family_margin != incumbent.cross_family_margin) {
        return cur.cross_family_margin > incumbent.cross_family_margin;
    }
    if (cur.best_te_identity != incumbent.best_te_identity) {
        return cur.best_te_identity > incumbent.best_te_identity;
    }
    if (cur.te_name != incumbent.te_name) {
        if (incumbent.te_name.empty() || incumbent.te_name == "UNK") {
            return true;
        }
        if (cur.te_name.empty() || cur.te_name == "UNK") {
            return false;
        }
    }
    if (cur.event_consensus_len != incumbent.event_consensus_len) {
        return cur.event_consensus_len > incumbent.event_consensus_len;
    }
    if (cur.best_te_query_coverage != incumbent.best_te_query_coverage) {
        return cur.best_te_query_coverage > incumbent.best_te_query_coverage;
    }
    return cur.pos < incumbent.pos;
}

double sigmoid(double x) {
    if (x >= 0.0) {
        const double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(x);
    return z / (1.0 + z);
}

double log_binomial_likelihood(int32_t alt_count, int32_t ref_count, double alt_rate) {
    const double p = std::clamp(alt_rate, 1e-6, 1.0 - 1e-6);
    return (static_cast<double>(alt_count) * std::log(p)) +
           (static_cast<double>(ref_count) * std::log(1.0 - p));
}

int32_t median_i32(std::vector<int32_t> values) {
    if (values.empty()) {
        return -1;
    }
    std::sort(values.begin(), values.end());
    const size_t mid = values.size() / 2;
    if (values.size() % 2 == 1) {
        return values[mid];
    }
    return static_cast<int32_t>((static_cast<int64_t>(values[mid - 1]) + values[mid]) / 2);
}

std::vector<int32_t> anchor_local_breakpoint_cluster(
    std::vector<int32_t> positions,
    int32_t anchor_pos) {
    constexpr int32_t kBreakpointClusterGapBp = 75;

    if (positions.empty()) {
        return {};
    }
    std::sort(positions.begin(), positions.end());

    struct ClusterBounds {
        size_t begin = 0;
        size_t end = 0;
        int32_t center = -1;
    };
    std::vector<ClusterBounds> clusters;
    size_t cluster_begin = 0;
    for (size_t i = 1; i <= positions.size(); ++i) {
        const bool split_cluster =
            i == positions.size() ||
            (positions[i] - positions[i - 1]) > kBreakpointClusterGapBp;
        if (!split_cluster) {
            continue;
        }

        ClusterBounds cluster;
        cluster.begin = cluster_begin;
        cluster.end = i;
        cluster.center = median_i32(std::vector<int32_t>(
            positions.begin() + static_cast<std::ptrdiff_t>(cluster_begin),
            positions.begin() + static_cast<std::ptrdiff_t>(i)));
        clusters.push_back(cluster);
        cluster_begin = i;
    }

    const ClusterBounds* best = nullptr;
    for (const auto& cluster : clusters) {
        if (!best) {
            best = &cluster;
            continue;
        }
        const int32_t cur_dist = std::abs(cluster.center - anchor_pos);
        const int32_t best_dist = std::abs(best->center - anchor_pos);
        if (cur_dist != best_dist) {
            if (cur_dist < best_dist) {
                best = &cluster;
            }
            continue;
        }
        const size_t cur_size = cluster.end - cluster.begin;
        const size_t best_size = best->end - best->begin;
        if (cur_size != best_size) {
            if (cur_size > best_size) {
                best = &cluster;
            }
            continue;
        }
        if (cluster.center < best->center) {
            best = &cluster;
        }
    }

    if (!best) {
        return positions;
    }
    return std::vector<int32_t>(
        positions.begin() + static_cast<std::ptrdiff_t>(best->begin),
        positions.begin() + static_cast<std::ptrdiff_t>(best->end));
}

int32_t anchor_local_median_i32(
    std::vector<int32_t> positions,
    int32_t anchor_pos) {
    return median_i32(anchor_local_breakpoint_cluster(std::move(positions), anchor_pos));
}

struct BreakpointHypothesis {
    bool valid = false;
    int32_t left = -1;
    int32_t right = -1;
    int32_t center = -1;
    int32_t support = 0;
    int32_t priority = std::numeric_limits<int32_t>::max();
};

int32_t breakpoint_hypothesis_support_weight(int32_t priority) {
    switch (priority) {
        case 0: return 8;  // fragment split
        case 1: return 8;  // fragment indel
        case 2: return 7;  // raw split
        case 3: return 1;  // fragment clip pair
        case 4: return 6;  // raw indel
        case 5: return 1;  // raw clip pair
        default: return 1;
    }
}

int32_t breakpoint_hypothesis_score(const BreakpointHypothesis& hypothesis) {
    return hypothesis.support * breakpoint_hypothesis_support_weight(hypothesis.priority);
}

BreakpointHypothesis make_single_breakpoint_hypothesis(
    std::vector<int32_t> positions,
    int32_t anchor_pos,
    int32_t priority) {
    BreakpointHypothesis hypothesis;
    std::vector<int32_t> cluster = anchor_local_breakpoint_cluster(std::move(positions), anchor_pos);
    if (cluster.empty()) {
        return hypothesis;
    }
    const int32_t pos = median_i32(cluster);
    hypothesis.valid = pos >= 0;
    hypothesis.left = pos;
    hypothesis.right = pos;
    hypothesis.center = pos;
    hypothesis.support = static_cast<int32_t>(cluster.size());
    hypothesis.priority = priority;
    return hypothesis;
}

BreakpointHypothesis make_paired_breakpoint_hypothesis(
    std::vector<int32_t> left_positions,
    std::vector<int32_t> right_positions,
    int32_t anchor_pos,
    int32_t priority) {
    BreakpointHypothesis hypothesis;
    std::vector<int32_t> left_cluster = anchor_local_breakpoint_cluster(
        std::move(left_positions),
        anchor_pos);
    std::vector<int32_t> right_cluster = anchor_local_breakpoint_cluster(
        std::move(right_positions),
        anchor_pos);
    if (left_cluster.empty() && right_cluster.empty()) {
        return hypothesis;
    }

    const int32_t left = !left_cluster.empty()
        ? median_i32(left_cluster)
        : median_i32(right_cluster);
    const int32_t right = !right_cluster.empty()
        ? median_i32(right_cluster)
        : median_i32(left_cluster);
    if (left < 0 || right < 0) {
        return hypothesis;
    }

    hypothesis.valid = true;
    hypothesis.left = std::min(left, right);
    hypothesis.right = std::max(left, right);
    hypothesis.center =
        static_cast<int32_t>(
            (static_cast<int64_t>(hypothesis.left) + hypothesis.right) / 2);
    hypothesis.support =
        static_cast<int32_t>(left_cluster.size() + right_cluster.size());
    hypothesis.priority = priority;
    return hypothesis;
}

bool better_breakpoint_hypothesis(
    const BreakpointHypothesis& candidate,
    const BreakpointHypothesis& incumbent,
    int32_t anchor_pos) {
    if (!candidate.valid) {
        return false;
    }
    if (!incumbent.valid) {
        return true;
    }

    const int32_t candidate_score = breakpoint_hypothesis_score(candidate);
    const int32_t incumbent_score = breakpoint_hypothesis_score(incumbent);
    if (candidate_score != incumbent_score) {
        return candidate_score > incumbent_score;
    }
    const int32_t candidate_dist = std::abs(candidate.center - anchor_pos);
    const int32_t incumbent_dist = std::abs(incumbent.center - anchor_pos);
    if (candidate_dist != incumbent_dist) {
        return candidate_dist < incumbent_dist;
    }
    if (candidate.support != incumbent.support) {
        return candidate.support > incumbent.support;
    }
    if (candidate.priority != incumbent.priority) {
        return candidate.priority < incumbent.priority;
    }
    if (candidate.left != incumbent.left) {
        return candidate.left < incumbent.left;
    }
    return candidate.right < incumbent.right;
}

double mad_from_positions(const std::vector<int32_t>& positions) {
    if (positions.empty()) {
        return 0.0;
    }
    const int32_t center = median_i32(positions);
    std::vector<int32_t> abs_dev;
    abs_dev.reserve(positions.size());
    for (int32_t pos : positions) {
        abs_dev.push_back(std::abs(pos - center));
    }
    return static_cast<double>(median_i32(std::move(abs_dev)));
}

bool is_softclip_source(InsertionFragmentSource source) {
    return source == InsertionFragmentSource::kClipRefLeft ||
           source == InsertionFragmentSource::kClipRefRight;
}

int32_t max_homopolymer_run(const std::string& seq) {
    if (seq.empty()) {
        return 0;
    }
    int32_t best = 1;
    int32_t run = 1;
    for (size_t i = 1; i < seq.size(); ++i) {
        if (seq[i] == seq[i - 1]) {
            ++run;
            best = std::max(best, run);
        } else {
            run = 1;
        }
    }
    return best;
}

double at_fraction(const std::string& seq) {
    int32_t at_count = 0;
    int32_t total = 0;
    for (char c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            continue;
        }
        ++total;
        if (c == 'A' || c == 'T') {
            ++at_count;
        }
    }
    if (total <= 0) {
        return 0.0;
    }
    return static_cast<double>(at_count) / static_cast<double>(total);
}

double shannon_entropy_acgt(const std::string& seq) {
    std::array<int32_t, 4> counts = {0, 0, 0, 0};
    int32_t total = 0;
    for (char c : seq) {
        switch (c) {
            case 'A':
                counts[0] += 1;
                total += 1;
                break;
            case 'C':
                counts[1] += 1;
                total += 1;
                break;
            case 'G':
                counts[2] += 1;
                total += 1;
                break;
            case 'T':
                counts[3] += 1;
                total += 1;
                break;
            default:
                break;
        }
    }
    if (total <= 0) {
        return 0.0;
    }

    double entropy = 0.0;
    for (int32_t count : counts) {
        if (count <= 0) {
            continue;
        }
        const double p = static_cast<double>(count) / static_cast<double>(total);
        entropy -= p * std::log2(p);
    }
    return entropy;
}

uint8_t char_to_2bit(char c);
bool build_kmer(const std::string& s, int32_t start, int32_t k, uint64_t& out);
double kmer_uniqueness_ratio(const std::string& seq, int32_t k);

bool is_low_complexity_softclip(
    const InsertionFragment& fragment,
    const PipelineConfig& config) {
    if (!is_softclip_source(fragment.source)) {
        return false;
    }
    if (fragment.sequence.empty()) {
        return true;
    }
    const double at_frac = at_fraction(fragment.sequence);
    const int32_t homopolymer = max_homopolymer_run(fragment.sequence);
    const double entropy = shannon_entropy_acgt(fragment.sequence);
    const double kmer_uniqueness = kmer_uniqueness_ratio(fragment.sequence, 5);
    return at_frac >= config.te_softclip_low_complexity_at_frac_min ||
           homopolymer >= config.te_softclip_low_complexity_homopolymer_min ||
           entropy < std::max(0.0, config.te_softclip_entropy_min) ||
           kmer_uniqueness < std::clamp(config.te_softclip_kmer_uniqueness_min, 0.0, 1.0);
}

bool is_low_quality_softclip_anchor(
    const InsertionFragment& fragment,
    const PipelineConfig& config) {
    if (!is_softclip_source(fragment.source)) {
        return false;
    }
    if (fragment.anchor_len > 0 &&
        fragment.anchor_len < std::max(1, config.te_softclip_min_anchor_len)) {
        return true;
    }
    if (fragment.anchor_len > 0 && fragment.nm >= 0) {
        const double nm_per_bp =
            static_cast<double>(fragment.nm) / static_cast<double>(std::max(1, fragment.anchor_len));
        if (nm_per_bp > std::max(0.0, config.te_softclip_max_nm_per_bp)) {
            return true;
        }
    }
    return false;
}

uint8_t char_to_2bit(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

bool build_kmer(const std::string& s, int32_t start, int32_t k, uint64_t& out) {
    out = 0;
    for (int32_t i = 0; i < k; ++i) {
        const uint8_t code = char_to_2bit(s[static_cast<size_t>(start + i)]);
        if (code > 3) {
            return false;
        }
        out = (out << 2) | code;
    }
    return true;
}

std::unordered_set<uint64_t> kmer_set(const std::string& seq, int32_t k) {
    std::unordered_set<uint64_t> out;
    if (k <= 0 || static_cast<int32_t>(seq.size()) < k) {
        return out;
    }
    out.reserve(static_cast<size_t>(seq.size()));
    for (int32_t i = 0; i + k <= static_cast<int32_t>(seq.size()); ++i) {
        uint64_t key = 0;
        if (!build_kmer(seq, i, k, key)) {
            continue;
        }
        out.insert(key);
    }
    return out;
}

double kmer_uniqueness_ratio(const std::string& seq, int32_t k) {
    if (k <= 0 || static_cast<int32_t>(seq.size()) < k) {
        return 0.0;
    }
    int32_t total = 0;
    std::unordered_set<uint64_t> uniq;
    uniq.reserve(static_cast<size_t>(seq.size()));
    for (int32_t i = 0; i + k <= static_cast<int32_t>(seq.size()); ++i) {
        uint64_t key = 0;
        if (!build_kmer(seq, i, k, key)) {
            continue;
        }
        uniq.insert(key);
        total += 1;
    }
    if (total <= 0) {
        return 0.0;
    }
    return static_cast<double>(uniq.size()) / static_cast<double>(total);
}

std::string upper_acgt(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    return out;
}

struct EditDistanceWorkspace {
    std::vector<int32_t> prev;
    std::vector<int32_t> curr;

    void ensure_columns(int32_t m) {
        const size_t columns = static_cast<size_t>(m + 1);
        if (prev.size() < columns) {
            prev.resize(columns, 0);
        }
        if (curr.size() < columns) {
            curr.resize(columns, 0);
        }
    }
};

int32_t max_edits_for_identity_threshold(
    int32_t lhs_len,
    int32_t rhs_len,
    double min_identity) {
    const int32_t denom = std::max(lhs_len, rhs_len);
    if (denom <= 0) {
        return 0;
    }
    const double clamped = std::clamp(min_identity, 0.0, 1.0);
    return std::max(
        0,
        static_cast<int32_t>(
            std::floor(((1.0 - clamped) * static_cast<double>(denom)) + 1e-9)));
}

bool edit_identity_if_at_least(
    std::string_view lhs,
    std::string_view rhs,
    int32_t max_edits,
    EditDistanceWorkspace& workspace,
    double& identity_out) {
    identity_out = 0.0;
    const int32_t n = static_cast<int32_t>(lhs.size());
    const int32_t m = static_cast<int32_t>(rhs.size());
    if (n <= 0 || m <= 0 || max_edits < 0) {
        return false;
    }
    if (std::abs(n - m) > max_edits) {
        return false;
    }

    workspace.ensure_columns(m);
    auto& prev = workspace.prev;
    auto& curr = workspace.curr;
    const size_t columns = static_cast<size_t>(m + 1);
    const int32_t inf = max_edits + 1;
    std::fill(prev.begin(), prev.begin() + columns, inf);
    std::fill(curr.begin(), curr.begin() + columns, inf);
    for (int32_t j = 0; j <= std::min(m, max_edits); ++j) {
        prev[static_cast<size_t>(j)] = j;
    }

    for (int32_t i = 1; i <= n; ++i) {
        std::fill(curr.begin(), curr.begin() + columns, inf);
        if (i <= max_edits) {
            curr[0] = i;
        }

        const int32_t j_lo = std::max(1, i - max_edits);
        const int32_t j_hi = std::min(m, i + max_edits);
        if (j_lo > j_hi) {
            return false;
        }

        for (int32_t j = j_lo; j <= j_hi; ++j) {
            const int32_t sub_cost =
                prev[static_cast<size_t>(j - 1)] +
                ((lhs[static_cast<size_t>(i - 1)] == rhs[static_cast<size_t>(j - 1)]) ? 0 : 1);
            const int32_t del_cost = prev[static_cast<size_t>(j)] + 1;
            const int32_t ins_cost = curr[static_cast<size_t>(j - 1)] + 1;
            curr[static_cast<size_t>(j)] = std::min({sub_cost, del_cost, ins_cost});
        }
        std::swap(prev, curr);
    }

    const int32_t dist = prev[static_cast<size_t>(m)];
    if (dist > max_edits) {
        return false;
    }

    const double denom = static_cast<double>(std::max(n, m));
    identity_out = std::clamp(1.0 - (static_cast<double>(dist) / denom), 0.0, 1.0);
    return true;
}

uint8_t nt_to_abpoa_code(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

char abpoa_code_to_nt(uint8_t c) {
    switch (c) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

std::string abpoa_consensus(const std::vector<std::string>& seqs) {
    if (seqs.empty()) {
        return "";
    }

    abpoa_t* ab = abpoa_init();
    abpoa_para_t* abpt = abpoa_init_para();
    if (!ab || !abpt) {
        if (ab) {
            abpoa_free(ab);
        }
        if (abpt) {
            abpoa_free_para(abpt);
        }
        return "";
    }

    abpt->out_cons = 1;
    abpt->out_msa = 0;
    abpt->out_gfa = 0;
    abpt->max_n_cons = 1;
    abpoa_post_set_para(abpt);

    const int32_t n_seq = static_cast<int32_t>(seqs.size());
    std::vector<int> seq_lens(static_cast<size_t>(n_seq), 0);
    std::vector<std::vector<uint8_t>> packed(static_cast<size_t>(n_seq));
    std::vector<uint8_t*> seq_ptrs(static_cast<size_t>(n_seq), nullptr);

    for (int32_t i = 0; i < n_seq; ++i) {
        const std::string seq = upper_acgt(seqs[static_cast<size_t>(i)]);
        seq_lens[static_cast<size_t>(i)] = static_cast<int>(seq.size());
        packed[static_cast<size_t>(i)].resize(seq.size(), 4);
        for (size_t j = 0; j < seq.size(); ++j) {
            packed[static_cast<size_t>(i)][j] = nt_to_abpoa_code(seq[j]);
        }
        seq_ptrs[static_cast<size_t>(i)] = packed[static_cast<size_t>(i)].data();
    }

    const int ret = abpoa_msa(
        ab,
        abpt,
        n_seq,
        nullptr,
        seq_lens.data(),
        seq_ptrs.data(),
        nullptr,
        nullptr);

    std::string consensus;
    if (ret == 0 && ab->abc && ab->abc->n_cons > 0) {
        const int cons_len = ab->abc->cons_len[0];
        if (cons_len > 0 && ab->abc->cons_base && ab->abc->cons_base[0]) {
            consensus.reserve(static_cast<size_t>(cons_len));
            for (int i = 0; i < cons_len; ++i) {
                consensus.push_back(abpoa_code_to_nt(ab->abc->cons_base[0][i]));
            }
        }
    }

    abpoa_free(ab);
    abpoa_free_para(abpt);
    return consensus;
}

bool has_prefix(const std::string& value, const char* prefix) {
    return value.rfind(prefix, 0) == 0;
}

bool has_qc_token(const std::string& qc, const char* token) {
    if (!token || *token == '\0') {
        return false;
    }
    size_t start = 0;
    while (start <= qc.size()) {
        const size_t end = qc.find('|', start);
        const size_t len = (end == std::string::npos) ? (qc.size() - start) : (end - start);
        if (len == std::strlen(token) && qc.compare(start, len, token) == 0) {
            return true;
        }
        if (end == std::string::npos) {
            break;
        }
        start = end + 1;
    }
    return false;
}

void append_qc_token(std::string& qc, const char* token) {
    if (!token || *token == '\0' || has_qc_token(qc, token)) {
        return;
    }
    if (!qc.empty()) {
        qc += "|";
    }
    qc += token;
}

void update_best_bp_by_read(
    std::unordered_map<size_t, int32_t>& bp_by_read,
    size_t read_index,
    int32_t pos,
    int32_t anchor_pos) {
    if (pos < 0) {
        return;
    }
    const auto it = bp_by_read.find(read_index);
    if (it == bp_by_read.end()) {
        bp_by_read.emplace(read_index, pos);
        return;
    }
    const int32_t cur_dist = std::abs(it->second - anchor_pos);
    const int32_t new_dist = std::abs(pos - anchor_pos);
    if (new_dist < cur_dist || (new_dist == cur_dist && pos < it->second)) {
        it->second = pos;
    }
}

std::pair<int32_t, int32_t> infer_component_breakpoint_bounds(
    const ComponentCall& component) {
    int32_t bp_min = std::numeric_limits<int32_t>::max();
    int32_t bp_max = std::numeric_limits<int32_t>::min();

    for (const auto& bp : component.breakpoint_candidates) {
        if (bp.pos < 0) {
            continue;
        }
        bp_min = std::min(bp_min, bp.pos);
        bp_max = std::max(bp_max, bp.pos);
    }

    if (bp_min <= bp_max) {
        return {bp_min, bp_max};
    }

    const int32_t fallback = std::max(0, component.anchor_pos);
    return {fallback, fallback};
}

std::vector<BreakpointHypothesis> enumerate_breakpoint_hypotheses(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right,
    size_t top_k) {
    constexpr int32_t kEventBreakpointSearchSlackBp = 200;

    const int32_t search_start = std::max(
        0,
        std::min(seed_left, seed_right) - kEventBreakpointSearchSlackBp);
    const int32_t search_end = std::max(
        search_start + 1,
        std::max(seed_left, seed_right) + kEventBreakpointSearchSlackBp);

    std::vector<int32_t> fragment_split_left_positions;
    std::vector<int32_t> fragment_split_right_positions;
    std::vector<int32_t> fragment_indel_positions;
    std::vector<int32_t> fragment_clip_left_positions;
    std::vector<int32_t> fragment_clip_right_positions;
    std::vector<int32_t> raw_split_left_positions;
    std::vector<int32_t> raw_split_right_positions;
    std::vector<int32_t> raw_indel_positions;
    std::vector<int32_t> raw_clip_left_positions;
    std::vector<int32_t> raw_clip_right_positions;

    auto add_position = [&](std::vector<int32_t>& out, int32_t pos) {
        if (pos >= search_start && pos <= search_end) {
            out.push_back(pos);
        }
    };

    for (const bam1_t* record : local_records) {
        if (!record) {
            continue;
        }
        ReadView read(record);
        if (read.tid() != component.tid) {
            continue;
        }

        const LocalEventSignal signal = classify_local_event_signal(
            record,
            search_start,
            search_end);
        const bool has_raw_split_left = signal.split_left_pos >= 0;
        const bool has_raw_split_right = signal.split_right_pos >= 0;
        if (has_raw_split_left && has_raw_split_right) {
            const int32_t left_dist = std::abs(signal.split_left_pos - component.anchor_pos);
            const int32_t right_dist = std::abs(signal.split_right_pos - component.anchor_pos);
            if (left_dist < right_dist ||
                (left_dist == right_dist && signal.split_left_pos <= signal.split_right_pos)) {
                add_position(raw_split_left_positions, signal.split_left_pos);
            } else {
                add_position(raw_split_right_positions, signal.split_right_pos);
            }
        } else {
            add_position(raw_split_left_positions, signal.split_left_pos);
            add_position(raw_split_right_positions, signal.split_right_pos);
        }
        add_position(raw_clip_left_positions, signal.left_clip_pos);
        add_position(raw_clip_right_positions, signal.right_clip_pos);
        add_position(raw_indel_positions, signal.indel_pos);
    }

    for (const auto& fragment : fragments) {
        if (fragment.ref_junc_pos < search_start || fragment.ref_junc_pos > search_end) {
            continue;
        }
        switch (fragment.source) {
            case InsertionFragmentSource::kSplitSa:
                if (fragment.ref_side == ReferenceSide::kRefLeft) {
                    fragment_split_left_positions.push_back(fragment.ref_junc_pos);
                } else if (fragment.ref_side == ReferenceSide::kRefRight) {
                    fragment_split_right_positions.push_back(fragment.ref_junc_pos);
                }
                break;
            case InsertionFragmentSource::kCigarInsertion:
                fragment_indel_positions.push_back(fragment.ref_junc_pos);
                break;
            case InsertionFragmentSource::kClipRefLeft:
                fragment_clip_left_positions.push_back(fragment.ref_junc_pos);
                break;
            case InsertionFragmentSource::kClipRefRight:
                fragment_clip_right_positions.push_back(fragment.ref_junc_pos);
                break;
            default:
                break;
        }
    }

    std::vector<BreakpointHypothesis> hypotheses;
    hypotheses.reserve(6);
    const auto consider = [&](const BreakpointHypothesis& hypothesis) {
        if (!hypothesis.valid) {
            return;
        }
        for (const auto& existing : hypotheses) {
            if (existing.left == hypothesis.left &&
                existing.right == hypothesis.right) {
                return;
            }
        }
        hypotheses.push_back(hypothesis);
    };

    consider(make_paired_breakpoint_hypothesis(
        fragment_split_left_positions,
        fragment_split_right_positions,
        component.anchor_pos,
        0));
    consider(make_single_breakpoint_hypothesis(
        fragment_indel_positions,
        component.anchor_pos,
        1));
    consider(make_paired_breakpoint_hypothesis(
        raw_split_left_positions,
        raw_split_right_positions,
        component.anchor_pos,
        2));
    consider(make_paired_breakpoint_hypothesis(
        fragment_clip_left_positions,
        fragment_clip_right_positions,
        component.anchor_pos,
        3));
    consider(make_single_breakpoint_hypothesis(
        raw_indel_positions,
        component.anchor_pos,
        4));
    consider(make_paired_breakpoint_hypothesis(
        raw_clip_left_positions,
        raw_clip_right_positions,
        component.anchor_pos,
        5));

    if (hypotheses.empty()) {
        BreakpointHypothesis fallback;
        fallback.valid = true;
        fallback.left = std::max(0, component.anchor_pos);
        fallback.right = fallback.left;
        fallback.center = fallback.left;
        fallback.support = 1;
        fallback.priority = std::numeric_limits<int32_t>::max();
        hypotheses.push_back(fallback);
    }

    std::sort(
        hypotheses.begin(),
        hypotheses.end(),
        [&](const BreakpointHypothesis& lhs, const BreakpointHypothesis& rhs) {
            return better_breakpoint_hypothesis(lhs, rhs, component.anchor_pos);
        });
    if (top_k > 0 && hypotheses.size() > top_k) {
        hypotheses.resize(top_k);
    }
    return hypotheses;
}

std::pair<int32_t, int32_t> resolve_event_breakpoint_bounds(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right) {
    const auto hypotheses = enumerate_breakpoint_hypotheses(
        component,
        local_records,
        fragments,
        seed_left,
        seed_right,
        1);
    if (!hypotheses.empty()) {
        return {hypotheses.front().left, hypotheses.front().right};
    }
    return {std::max(0, component.anchor_pos), std::max(0, component.anchor_pos)};
}
template <typename T>
class SafeQueue {
public:
    explicit SafeQueue(size_t max_size = 0) : max_size_(max_size) {}

    void push(T&& value) {
        std::unique_lock<std::mutex> lock(mu_);
        if (max_size_ > 0) {
            not_full_cv_.wait(lock, [this]() { return finished_ || queue_.size() < max_size_; });
        }
        if (finished_) {
            return;
        }
        queue_.push(std::move(value));
        lock.unlock();
        cv_.notify_one();
    }

    bool pop(T& out) {
        std::unique_lock<std::mutex> lock(mu_);
        cv_.wait(lock, [this]() { return finished_ || !queue_.empty(); });
        if (queue_.empty()) {
            return false;
        }
        out = std::move(queue_.front());
        queue_.pop();
        lock.unlock();
        if (max_size_ > 0) {
            not_full_cv_.notify_one();
        }
        return true;
    }

    void close() {
        {
            std::lock_guard<std::mutex> lock(mu_);
            finished_ = true;
        }
        cv_.notify_all();
        if (max_size_ > 0) {
            not_full_cv_.notify_all();
        }
    }

private:
    std::queue<T> queue_;
    std::mutex mu_;
    std::condition_variable cv_;
    std::condition_variable not_full_cv_;
    bool finished_ = false;
    size_t max_size_ = 0;
};

struct BinTask {
    std::vector<const bam1_t*> records;
    std::array<SharedRecordBatch, static_cast<size_t>(kCrossBinContextBins) + 1> owners;
    size_t owner_count = 0;
    int32_t tid = -1;
    int32_t bin_index = -1;
};

bool is_match_like(int op) {
    return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
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

int32_t weighted_median_position(std::vector<std::pair<int32_t, double>> pos_weights) {
    pos_weights.erase(
        std::remove_if(pos_weights.begin(), pos_weights.end(), [](const auto& row) {
            return row.second <= 0.0;
        }),
        pos_weights.end());
    if (pos_weights.empty()) {
        return -1;
    }

    std::sort(pos_weights.begin(), pos_weights.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    double total_weight = 0.0;
    for (const auto& row : pos_weights) {
        total_weight += row.second;
    }
    if (total_weight <= 0.0) {
        return pos_weights[pos_weights.size() / 2].first;
    }

    const double half_weight = 0.5 * total_weight;
    double acc_weight = 0.0;
    for (const auto& row : pos_weights) {
        acc_weight += row.second;
        if (acc_weight >= half_weight) {
            return row.first;
        }
    }
    return pos_weights.back().first;
}

std::vector<double> smooth_histogram(const std::vector<double>& hist) {
    std::vector<double> smooth(hist.size(), 0.0);
    for (size_t i = 0; i < hist.size(); ++i) {
        double value = hist[i];
        if (i >= 1) {
            value += 0.60 * hist[i - 1];
        }
        if (i + 1 < hist.size()) {
            value += 0.60 * hist[i + 1];
        }
        if (i >= 2) {
            value += 0.25 * hist[i - 2];
        }
        if (i + 2 < hist.size()) {
            value += 0.25 * hist[i + 2];
        }
        smooth[i] = value;
    }
    return smooth;
}

EvidenceBundle extract_evidence_points(
    const std::vector<const bam1_t*>& bin_records,
    int32_t expected_tid) {
    EvidenceBundle bundle;
    bundle.read_summaries.resize(bin_records.size());

    for (size_t idx = 0; idx < bin_records.size(); ++idx) {
        if (!bin_records[idx]) {
            continue;
        }

        ReadView view(bin_records[idx]);
        if (view.tid() != expected_tid) {
            continue;
        }
        if ((view.flag() & BAM_FSUPPLEMENTARY) != 0) {
            continue;
        }

        ReadSignalSummary summary;
        summary.read_id = std::string(view.qname());
        summary.is_reverse = (view.flag() & BAM_FREVERSE) != 0;
        summary.has_sa_or_supp = view.has_sa_tag() || ((view.flag() & BAM_FSUPPLEMENTARY) != 0);

        const uint32_t* cigar = view.cigar();
        const int32_t n_cigar = view.n_cigar();
        if (!cigar || n_cigar <= 0) {
            if (summary.has_sa_or_supp) {
                summary.class_mask |= kCandidateSplitSaSupplementary;
            }
            bundle.read_summaries[idx] = std::move(summary);
            continue;
        }

        int32_t leading_soft = 0;
        int32_t trailing_soft = 0;
        const int first = find_first_non_hard_clip(cigar, n_cigar);
        const int last = find_last_non_hard_clip(cigar, n_cigar);
        if (first >= 0 && bam_cigar_op(cigar[first]) == BAM_CSOFT_CLIP) {
            leading_soft = static_cast<int32_t>(bam_cigar_oplen(cigar[first]));
        }
        if (last >= 0 && bam_cigar_op(cigar[last]) == BAM_CSOFT_CLIP) {
            trailing_soft = static_cast<int32_t>(bam_cigar_oplen(cigar[last]));
        }

        std::vector<EvidencePoint> local_points;
        local_points.reserve(8);

        int32_t ref_pos = view.pos();
        for (int32_t ci = 0; ci < n_cigar; ++ci) {
            const int op = bam_cigar_op(cigar[ci]);
            const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[ci]));

            if (op == BAM_CSOFT_CLIP) {
                summary.max_soft_clip = std::max(summary.max_soft_clip, len);
            } else if (op == BAM_CINS) {
                summary.max_ins = std::max(summary.max_ins, len);
            }

            if (op == BAM_CINS && len >= kLargeIndelEvidenceMin) {
                EvidencePoint point;
                point.read_index = idx;
                point.pos = ref_pos;
                point.weight = 0.80 + (static_cast<double>(std::min(len, 400)) / 200.0);
                point.kind = EvidenceKind::kIndel;
                point.signal_len = len;
                point.class_mask = kCandidateLongInsertion;
                local_points.push_back(point);
            }

            if ((bam_cigar_type(op) & 2) != 0) {
                ref_pos += len;
            }
        }
        const int32_t ref_end = ref_pos;

        if (leading_soft >= kSoftClipSignalMin) {
            EvidencePoint point;
            point.read_index = idx;
            point.pos = view.pos();
            point.weight = 1.00 + (static_cast<double>(std::min(leading_soft, 500)) / 180.0);
            point.kind = EvidenceKind::kSoftClip;
            point.signal_len = leading_soft;
            point.class_mask = kCandidateSoftClip;
            local_points.push_back(point);
        }
        if (trailing_soft >= kSoftClipSignalMin) {
            EvidencePoint point;
            point.read_index = idx;
            point.pos = ref_end;
            point.weight = 1.00 + (static_cast<double>(std::min(trailing_soft, 500)) / 180.0);
            point.kind = EvidenceKind::kSoftClip;
            point.signal_len = trailing_soft;
            point.class_mask = kCandidateSoftClip;
            local_points.push_back(point);
        }

        if (summary.has_sa_or_supp) {
            const SaHintSide hint_side = choose_sa_hint_side(leading_soft, trailing_soft);
            if (hint_side != SaHintSide::kNone) {
                EvidencePoint hint;
                hint.read_index = idx;
                hint.pos = (hint_side == SaHintSide::kLeft) ? view.pos() : ref_end;
                hint.weight = kSaHintWeight;
                hint.kind = EvidenceKind::kSAHint;
                hint.class_mask = kCandidateSplitSaSupplementary;
                local_points.push_back(hint);
            }
        }

        if (summary.has_sa_or_supp) {
            summary.class_mask |= kCandidateSplitSaSupplementary;
        }
        if (summary.max_soft_clip >= kSoftClipSignalMin) {
            summary.class_mask |= kCandidateSoftClip;
        }
        if (summary.max_ins >= kLongInsertionSignalMin) {
            summary.class_mask |= kCandidateLongInsertion;
        }

        for (auto& point : local_points) {
            bundle.points.push_back(point);
        }
        bundle.read_summaries[idx] = std::move(summary);
    }

    return bundle;
}

std::vector<CandidateWindow> build_density_windows(
    const std::vector<EvidencePoint>& evidence,
    int32_t bin_start,
    int32_t bin_end) {
    if (evidence.empty()) {
        return {};
    }

    int32_t min_pos = evidence.front().pos;
    int32_t max_pos = evidence.front().pos;
    for (const auto& point : evidence) {
        min_pos = std::min(min_pos, point.pos);
        max_pos = std::max(max_pos, point.pos);
    }

    const int32_t hist_start = std::min(min_pos, bin_start) - kHistPadding;
    const int32_t hist_end = std::max(max_pos, bin_end) + kHistPadding;
    const int32_t num_bins = std::max(1, ((hist_end - hist_start) / kHistBinSize) + 1);

    std::vector<double> hist(static_cast<size_t>(num_bins), 0.0);
    for (const auto& point : evidence) {
        int32_t idx = (point.pos - hist_start) / kHistBinSize;
        idx = std::max(0, std::min(num_bins - 1, idx));
        hist[static_cast<size_t>(idx)] += point.weight;
    }
    const std::vector<double> smooth = smooth_histogram(hist);

    std::vector<DensityPeak> peaks;
    for (int32_t i = 0; i < num_bins; ++i) {
        const double center = smooth[static_cast<size_t>(i)];
        if (center < kPeakMinWeight) {
            continue;
        }
        const double left = (i > 0) ? smooth[static_cast<size_t>(i - 1)] : center;
        const double right = (i + 1 < num_bins) ? smooth[static_cast<size_t>(i + 1)] : center;
        if (center >= left && center >= right) {
            DensityPeak peak;
            peak.pos = hist_start + (i * kHistBinSize) + (kHistBinSize / 2);
            peak.weight = center;
            peaks.push_back(peak);
        }
    }

    // Preserve direct long-indel loci even when nearby clip-driven density
    // dominates the smoothed histogram. These points are breakpoint-specific
    // enough to seed their own candidate windows.
    std::vector<const EvidencePoint*> indel_points;
    indel_points.reserve(evidence.size());
    for (const auto& point : evidence) {
        if (point.kind == EvidenceKind::kIndel) {
            indel_points.push_back(&point);
        }
    }
    std::sort(indel_points.begin(), indel_points.end(), [](const EvidencePoint* a, const EvidencePoint* b) {
        if (a == nullptr || b == nullptr) {
            return a < b;
        }
        return a->pos < b->pos;
    });
    size_t indel_cluster_begin = 0;
    while (indel_cluster_begin < indel_points.size()) {
        size_t indel_cluster_end = indel_cluster_begin + 1;
        while (indel_cluster_end < indel_points.size() &&
               std::abs(indel_points[indel_cluster_end]->pos -
                        indel_points[indel_cluster_end - 1]->pos) <= kPeakMergeDistance) {
            ++indel_cluster_end;
        }

        double total_weight = 0.0;
        double weighted_pos = 0.0;
        for (size_t i = indel_cluster_begin; i < indel_cluster_end; ++i) {
            total_weight += indel_points[i]->weight;
            weighted_pos += indel_points[i]->weight * static_cast<double>(indel_points[i]->pos);
        }
        if (total_weight >= kPeakMinWeight) {
            DensityPeak peak;
            peak.pos = static_cast<int32_t>(std::llround(weighted_pos / total_weight));
            peak.weight = total_weight;
            peaks.push_back(peak);
        }
        indel_cluster_begin = indel_cluster_end;
    }

    if (peaks.empty()) {
        std::vector<std::pair<int32_t, double>> pos_weights;
        pos_weights.reserve(evidence.size());
        for (const auto& point : evidence) {
            pos_weights.push_back({point.pos, point.weight});
        }
        const int32_t fallback_center = weighted_median_position(std::move(pos_weights));
        if (fallback_center >= 0) {
            DensityPeak fallback_peak;
            fallback_peak.pos = fallback_center;
            fallback_peak.weight = 1.0;
            peaks.push_back(fallback_peak);
        }
    }
    if (peaks.empty()) {
        return {};
    }

    std::sort(peaks.begin(), peaks.end(), [](const DensityPeak& a, const DensityPeak& b) {
        return a.pos < b.pos;
    });

    std::vector<DensityPeak> merged_peaks;
    for (const auto& peak : peaks) {
        if (merged_peaks.empty() ||
            (peak.pos - merged_peaks.back().pos) > kPeakMergeDistance) {
            merged_peaks.push_back(peak);
            continue;
        }

        DensityPeak& back = merged_peaks.back();
        const double total_weight = back.weight + peak.weight;
        if (total_weight > 0.0) {
            const double weighted_center =
                ((back.weight * static_cast<double>(back.pos)) +
                 (peak.weight * static_cast<double>(peak.pos))) / total_weight;
            back.pos = static_cast<int32_t>(std::llround(weighted_center));
        }
        back.weight = total_weight;
    }

    std::vector<CandidateWindow> windows;
    windows.reserve(merged_peaks.size());

    for (const auto& peak : merged_peaks) {
        int32_t center_idx = (peak.pos - hist_start) / kHistBinSize;
        center_idx = std::max(0, std::min(num_bins - 1, center_idx));
        const double center_weight = smooth[static_cast<size_t>(center_idx)];
        const double expand_threshold = std::max(0.50, center_weight * kWindowDropRatio);

        int32_t left = center_idx;
        int32_t right = center_idx;

        while (left > 0 &&
               (center_idx - left) < kWindowMaxExpandBins &&
               smooth[static_cast<size_t>(left - 1)] >= expand_threshold) {
            --left;
        }
        while ((right + 1) < num_bins &&
               (right - center_idx) < kWindowMaxExpandBins &&
               smooth[static_cast<size_t>(right + 1)] >= expand_threshold) {
            ++right;
        }

        int32_t start = hist_start + (left * kHistBinSize);
        int32_t end = hist_start + ((right + 1) * kHistBinSize);
        if ((end - start) < kMinWindowSpan) {
            start = peak.pos - (kMinWindowSpan / 2);
            end = start + kMinWindowSpan;
        }

        const int32_t relaxed_bin_start = bin_start - kWindowBinSlackBp;
        const int32_t relaxed_bin_end = bin_end + kWindowBinSlackBp;
        start = std::max(relaxed_bin_start, start);
        end = std::min(relaxed_bin_end, end);
        if (end <= start) {
            continue;
        }

        CandidateWindow window;
        window.start = start;
        window.end = end;
        window.center = std::max(window.start, std::min(window.end - 1, peak.pos));
        window.peak_weight = std::max(peak.weight, center_weight);
        windows.push_back(window);
    }

    if (windows.empty()) {
        return windows;
    }

    std::sort(windows.begin(), windows.end(), [](const CandidateWindow& a, const CandidateWindow& b) {
        if (a.start != b.start) {
            return a.start < b.start;
        }
        return a.end < b.end;
    });

    std::vector<CandidateWindow> merged_windows;
    for (const auto& window : windows) {
        if (merged_windows.empty() ||
            window.start > (merged_windows.back().end + kWindowMergeGap)) {
            merged_windows.push_back(window);
            continue;
        }

        CandidateWindow& back = merged_windows.back();
        const double total_weight = back.peak_weight + window.peak_weight;
        if (total_weight > 0.0) {
            const double center =
                ((back.peak_weight * static_cast<double>(back.center)) +
                 (window.peak_weight * static_cast<double>(window.center))) / total_weight;
            back.center = static_cast<int32_t>(std::llround(center));
        }
        back.start = std::min(back.start, window.start);
        back.end = std::max(back.end, window.end);
        back.peak_weight = std::max(back.peak_weight, window.peak_weight);
    }

    return merged_windows;
}

}  // namespace

EventReadEvidence Pipeline::collect_event_read_evidence(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<ReadReferenceSpan>& read_spans,
    const std::vector<InsertionFragment>& fragments) const {
    const auto seed_bounds = infer_component_breakpoint_bounds(component);
    const auto bp_bounds = resolve_event_breakpoint_bounds(
        component,
        local_records,
        fragments,
        seed_bounds.first,
        seed_bounds.second);
    return collect_event_read_evidence_for_bounds(
        component,
        local_records,
        read_spans,
        fragments,
        seed_bounds.first,
        seed_bounds.second,
        bp_bounds.first,
        bp_bounds.second);
}

EventReadEvidence Pipeline::collect_event_read_evidence_for_bounds(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<ReadReferenceSpan>& read_spans,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right,
    int32_t bp_left,
    int32_t bp_right) const {
    EventReadEvidence evidence;
    (void)seed_left;
    (void)seed_right;
    constexpr int32_t kAltSignalSlackBp = 25;
    constexpr int32_t kRefSignalSlackBp = 75;
    constexpr int32_t kCandidateClipSlackBp = 25;
    constexpr int32_t kCandidateClipSupplementMaxOffsetBp = 75;
    const int32_t left = std::min(bp_left, bp_right);
    const int32_t right = std::max(bp_left, bp_right);
    evidence.bp_left = left;
    evidence.bp_right = right;

    std::unordered_set<std::string> split_qnames;
    std::unordered_set<std::string> indel_qnames;
    std::unordered_set<std::string> left_clip_qnames;
    std::unordered_set<std::string> right_clip_qnames;
    std::unordered_set<std::string> nearby_candidate_left_clip_qnames;
    std::unordered_set<std::string> nearby_candidate_right_clip_qnames;

    // Alt support should stay hypothesis-specific, while clean ref-spanning
    // reads need a slightly wider exclusion band so nearby clipped noise does
    // not masquerade as reference support.
    const int32_t alt_signal_start = std::max(0, evidence.bp_left - kAltSignalSlackBp);
    const int32_t alt_signal_end = std::max(alt_signal_start, evidence.bp_right + kAltSignalSlackBp);
    const int32_t ref_signal_start = std::max(0, evidence.bp_left - kRefSignalSlackBp);
    const int32_t ref_signal_end = std::max(ref_signal_start, evidence.bp_right + kRefSignalSlackBp);
    const int32_t ref_span_start = alt_signal_start;
    const int32_t ref_span_end = alt_signal_end;

    std::vector<int32_t> nearby_component_clip_candidates;
    nearby_component_clip_candidates.reserve(component.breakpoint_candidates.size());
    for (const auto& candidate : component.breakpoint_candidates) {
        if (candidate.pos < 0) {
            continue;
        }
        const int32_t dist_to_hypothesis = (candidate.pos < left)
            ? (left - candidate.pos)
            : ((candidate.pos > right) ? (candidate.pos - right) : 0);
        if (dist_to_hypothesis <= kCandidateClipSupplementMaxOffsetBp) {
            nearby_component_clip_candidates.push_back(candidate.pos);
        }
    }
    std::sort(
        nearby_component_clip_candidates.begin(),
        nearby_component_clip_candidates.end());
    nearby_component_clip_candidates.erase(
        std::unique(
            nearby_component_clip_candidates.begin(),
            nearby_component_clip_candidates.end()),
        nearby_component_clip_candidates.end());

    const bool have_nearby_component_clip_candidates =
        !nearby_component_clip_candidates.empty();
    int32_t candidate_clip_start = 0;
    int32_t candidate_clip_end = 0;
    if (have_nearby_component_clip_candidates) {
        candidate_clip_start = std::max(
            0,
            nearby_component_clip_candidates.front() - kCandidateClipSlackBp);
        candidate_clip_end = std::max(
            candidate_clip_start,
            nearby_component_clip_candidates.back() + kCandidateClipSlackBp);
    }
    const auto matches_nearby_component_clip_candidate = [&](int32_t pos) {
        if (pos < 0) {
            return false;
        }
        for (const int32_t candidate_pos : nearby_component_clip_candidates) {
            if (std::abs(pos - candidate_pos) <= kCandidateClipSlackBp) {
                return true;
            }
        }
        return false;
    };

    for (size_t read_idx = 0; read_idx < local_records.size(); ++read_idx) {
        const bam1_t* record = local_records[read_idx];
        if (!record) {
            continue;
        }
        ReadView read(record);
        if (read.tid() != component.tid) {
            continue;
        }

        const std::string qname = record_qname(record);
        if (qname.empty()) {
            continue;
        }

        const LocalEventSignal signal = classify_local_event_signal(
            record,
            alt_signal_start,
            alt_signal_end);
        if (signal.split) {
            split_qnames.insert(qname);
        }
        if (signal.indel) {
            indel_qnames.insert(qname);
        }
        if (signal.left_clip) {
            left_clip_qnames.insert(qname);
        }
        if (signal.right_clip) {
            right_clip_qnames.insert(qname);
        }

        if (have_nearby_component_clip_candidates) {
            const LocalEventSignal nearby_signal = classify_local_event_signal(
                record,
                candidate_clip_start,
                candidate_clip_end);
            if (nearby_signal.left_clip &&
                matches_nearby_component_clip_candidate(nearby_signal.left_clip_pos)) {
                nearby_candidate_left_clip_qnames.insert(qname);
            }
            if (nearby_signal.right_clip &&
                matches_nearby_component_clip_candidate(nearby_signal.right_clip_pos)) {
                nearby_candidate_right_clip_qnames.insert(qname);
            }
        }
    }

    for (const auto& fragment : fragments) {
        if (fragment.read_id.empty() ||
            fragment.ref_junc_pos < alt_signal_start ||
            fragment.ref_junc_pos > alt_signal_end) {
            continue;
        }
        switch (fragment.source) {
            case InsertionFragmentSource::kSplitSa:
                split_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kCigarInsertion:
                indel_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefLeft:
                left_clip_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefRight:
                right_clip_qnames.insert(fragment.read_id);
                break;
            default:
                break;
        }

        if (have_nearby_component_clip_candidates &&
            matches_nearby_component_clip_candidate(fragment.ref_junc_pos)) {
            switch (fragment.source) {
                case InsertionFragmentSource::kClipRefLeft:
                    nearby_candidate_left_clip_qnames.insert(fragment.read_id);
                    break;
                case InsertionFragmentSource::kClipRefRight:
                    nearby_candidate_right_clip_qnames.insert(fragment.read_id);
                    break;
                default:
                    break;
            }
        }
    }

    const bool has_precise_nonclip_support =
        !split_qnames.empty() || !indel_qnames.empty();
    if (has_precise_nonclip_support &&
        !nearby_candidate_left_clip_qnames.empty() &&
        !nearby_candidate_right_clip_qnames.empty()) {
        left_clip_qnames.insert(
            nearby_candidate_left_clip_qnames.begin(),
            nearby_candidate_left_clip_qnames.end());
        right_clip_qnames.insert(
            nearby_candidate_right_clip_qnames.begin(),
            nearby_candidate_right_clip_qnames.end());
    }

    evidence.alt_split_reads = static_cast<int32_t>(split_qnames.size());
    evidence.alt_indel_reads = static_cast<int32_t>(indel_qnames.size());
    evidence.alt_left_clip_reads = static_cast<int32_t>(left_clip_qnames.size());
    evidence.alt_right_clip_reads = static_cast<int32_t>(right_clip_qnames.size());

    std::unordered_set<std::string> alt_qnames = split_qnames;
    alt_qnames.insert(indel_qnames.begin(), indel_qnames.end());
    const bool clip_support_has_partner =
        !alt_qnames.empty() ||
        (!left_clip_qnames.empty() && !right_clip_qnames.empty());
    if (clip_support_has_partner) {
        alt_qnames.insert(left_clip_qnames.begin(), left_clip_qnames.end());
        alt_qnames.insert(right_clip_qnames.begin(), right_clip_qnames.end());
    }

    evidence.alt_struct_reads = static_cast<int32_t>(alt_qnames.size());
    evidence.support_qnames = sorted_string_set(alt_qnames);

    std::unordered_set<std::string> ref_qnames;
    std::unordered_set<std::string> low_mapq_ref_qnames;
    for (size_t read_idx = 0; read_idx < read_spans.size(); ++read_idx) {
        const auto& span = read_spans[read_idx];
        if (!span.valid || span.tid != component.tid) {
            continue;
        }
        if (span.start > ref_span_start || span.end < ref_span_end) {
            continue;
        }
        if (read_idx >= local_records.size() || !local_records[read_idx]) {
            continue;
        }
        const bam1_t* record = local_records[read_idx];
        if (!is_primary_alignment(record)) {
            continue;
        }
        const std::string qname = record_qname(record);
        if (qname.empty() || alt_qnames.find(qname) != alt_qnames.end()) {
            continue;
        }
        if (read_has_local_event_signal(record, ref_signal_start, ref_signal_end)) {
            continue;
        }
        if (record->core.qual >= 20) {
            ref_qnames.insert(qname);
        } else {
            low_mapq_ref_qnames.insert(qname);
        }
    }

    evidence.ref_span_reads = static_cast<int32_t>(ref_qnames.size());
    evidence.low_mapq_ref_span_reads = static_cast<int32_t>(low_mapq_ref_qnames.size());
    evidence.ref_span_qnames = sorted_string_set(ref_qnames);
    return evidence;
}

EventConsensus Pipeline::build_event_consensus(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    const EventReadEvidence& event_evidence) const {
    EventConsensus consensus;
    (void)component;

    std::unordered_set<std::string> support_qnames(
        event_evidence.support_qnames.begin(),
        event_evidence.support_qnames.end());
    std::unordered_map<std::string, const bam1_t*> records_by_qname;
    records_by_qname.reserve(local_records.size());
    for (const bam1_t* record : local_records) {
        if (!record) {
            continue;
        }
        const std::string qname = record_qname(record);
        if (qname.empty()) {
            continue;
        }

        auto it = records_by_qname.find(qname);
        if (it == records_by_qname.end()) {
            records_by_qname.emplace(qname, record);
            continue;
        }

        ReadView candidate(record);
        ReadView incumbent(it->second);
        const int32_t candidate_seq_len = candidate.seq_len();
        const int32_t incumbent_seq_len = incumbent.seq_len();
        if (candidate_seq_len > incumbent_seq_len ||
            (candidate_seq_len == incumbent_seq_len &&
             is_primary_alignment(record) &&
             !is_primary_alignment(it->second))) {
            it->second = record;
        }
    }

    std::unordered_map<std::string, std::string> full_event_by_qname;
    std::unordered_map<std::string, std::string> partial_event_by_qname;
    full_event_by_qname.reserve(fragments.size());
    partial_event_by_qname.reserve(fragments.size());

    for (const auto& fragment : fragments) {
        if (fragment.read_id.empty() ||
            support_qnames.find(fragment.read_id) == support_qnames.end() ||
            !fragment_supports_event_consensus(fragment) ||
            !fragment_is_local_to_event(fragment, event_evidence.bp_left, event_evidence.bp_right)) {
            continue;
        }
        const auto record_it = records_by_qname.find(fragment.read_id);
        if (record_it == records_by_qname.end()) {
            continue;
        }
        const std::string event_string = build_event_string_from_fragment(
            record_it->second,
            fragment);
        if (event_string.empty()) {
            continue;
        }
        auto& event_by_qname = fragment_has_full_event_context(fragment)
            ? full_event_by_qname
            : partial_event_by_qname;
        auto it = event_by_qname.find(fragment.read_id);
        if (it == event_by_qname.end() || event_string.size() > it->second.size()) {
            event_by_qname[fragment.read_id] = event_string;
        }
    }

    const bool use_full_context = !full_event_by_qname.empty();
    const auto& event_by_qname = use_full_context
        ? full_event_by_qname
        : partial_event_by_qname;
    std::vector<std::string> event_strings;
    event_strings.reserve(event_by_qname.size());
    for (const auto& row : event_by_qname) {
        event_strings.push_back(row.second);
    }

    consensus.input_event_reads = static_cast<int32_t>(event_strings.size());
    if (consensus.input_event_reads <= 0) {
        consensus.qc_reason = "NO_EVENT_STRING_READS";
        return consensus;
    }
    const int32_t min_event_reads = use_full_context
        ? 1
        : std::max(1, config_.event_consensus_poa_min_reads);
    if (consensus.input_event_reads < min_event_reads) {
        consensus.qc_reason = "INSUFFICIENT_EVENT_READS";
        return consensus;
    }

    std::sort(event_strings.begin(), event_strings.end(), [](const std::string& a, const std::string& b) {
        if (a.size() != b.size()) {
            return a.size() > b.size();
        }
        return a < b;
    });
    if (static_cast<int32_t>(event_strings.size()) > config_.event_consensus_poa_max_reads) {
        event_strings.resize(static_cast<size_t>(config_.event_consensus_poa_max_reads));
    }

    consensus.consensus_seq = upper_acgt(abpoa_consensus(event_strings));
    consensus.consensus_len = static_cast<int32_t>(consensus.consensus_seq.size());
    if (consensus.consensus_seq.empty()) {
        consensus.qc_reason = "EMPTY_EVENT_CONSENSUS";
        return consensus;
    }

    consensus.qc_pass = true;
    consensus.qc_reason = "PASS_EVENT_CONSENSUS";
    return consensus;
}

EventSegmentation Pipeline::segment_event_consensus(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence,
    const EventConsensus& event_consensus) const {
    EventSegmentation segmentation;

    if (!event_consensus.qc_pass || event_consensus.consensus_seq.empty()) {
        segmentation.qc_reason = "NO_EVENT_CONSENSUS_TO_SEGMENT";
        return segmentation;
    }
    if (!tsd_detector_.can_fetch_reference()) {
        segmentation.qc_reason = "REFERENCE_FETCH_UNAVAILABLE";
        return segmentation;
    }

    const std::string& consensus = event_consensus.consensus_seq;
    const int32_t consensus_len = static_cast<int32_t>(consensus.size());
    const int32_t min_total_len =
        (2 * kEventSegmentationMinFlankAlignBp) + kEventSegmentationMinInsertBp;
    if (consensus_len < min_total_len) {
        segmentation.qc_reason = "EVENT_CONSENSUS_TOO_SHORT";
        return segmentation;
    }

    const int32_t max_flank_query_len = std::min(
        kEventSegmentationMaxFlankQueryBp,
        consensus_len - kEventSegmentationMinFlankAlignBp - kEventSegmentationMinInsertBp);
    if (max_flank_query_len < kEventSegmentationMinFlankAlignBp) {
        segmentation.qc_reason = "EVENT_CONSENSUS_TOO_SHORT";
        return segmentation;
    }

    const int32_t bp_left = std::min(event_evidence.bp_left, event_evidence.bp_right);
    const int32_t bp_right = std::max(event_evidence.bp_left, event_evidence.bp_right);
    if (component.chrom.empty() || bp_left < 0 || bp_right < 0) {
        segmentation.qc_reason = "INVALID_EVENT_BREAKPOINTS";
        return segmentation;
    }

    const int32_t left_window_start = std::max(
        0,
        bp_left - kEventSegmentationBreakpointSlackBp - max_flank_query_len);
    const int32_t left_window_end = std::max(
        left_window_start + 1,
        bp_left + kEventSegmentationBreakpointSlackBp);
    const int32_t right_window_start = std::max(
        0,
        bp_right - kEventSegmentationBreakpointSlackBp);
    const int32_t right_window_end = std::max(
        right_window_start + 1,
        bp_right + kEventSegmentationBreakpointSlackBp + max_flank_query_len);

    const std::string left_ref_window = tsd_detector_.fetch_window(
        component.chrom,
        left_window_start,
        left_window_end);
    const std::string right_ref_window = tsd_detector_.fetch_window(
        component.chrom,
        right_window_start,
        right_window_end);
    if (left_ref_window.empty() || right_ref_window.empty()) {
        segmentation.qc_reason = "REFERENCE_WINDOW_FETCH_FAILED";
        return segmentation;
    }

    EditDistanceWorkspace edit_workspace;
    const auto collect_candidates = [&](
        bool is_left,
        int32_t breakpoint,
        int32_t ref_window_start,
        const std::string& ref_window,
        int32_t query_endpoint_slack) {
        EventFlankSearchResult result;
        if (ref_window.empty()) {
            return result;
        }

        const int32_t ref_window_end =
            ref_window_start + static_cast<int32_t>(ref_window.size());
        const int32_t search_lo = std::max(0, breakpoint - kEventSegmentationBreakpointSlackBp);
        const int32_t search_hi = breakpoint + kEventSegmentationBreakpointSlackBp;

        for (int32_t align_len = max_flank_query_len;
             align_len >= kEventSegmentationMinFlankAlignBp;
             --align_len) {
            const int32_t max_query_start = consensus_len - align_len;
            if (max_query_start < 0) {
                continue;
            }
            const int32_t max_edits = max_edits_for_identity_threshold(
                align_len,
                align_len,
                kEventSegmentationMinFlankIdentity);

            const int32_t endpoint_slack = std::max(0, std::min(
                query_endpoint_slack,
                max_query_start));
            const int32_t query_lo = is_left
                ? 0
                : std::max(0, max_query_start - endpoint_slack);
            const int32_t query_hi = is_left
                ? std::min(max_query_start, endpoint_slack)
                : max_query_start;

            std::vector<EventFlankPlacement> placements;
            placements.reserve(
                static_cast<size_t>(std::max(0, search_hi - search_lo + 1)) *
                static_cast<size_t>(std::max(1, query_hi - query_lo + 1)));

            for (int32_t query_start = query_lo; query_start <= query_hi; ++query_start) {
                const std::string_view query(
                    consensus.data() + query_start,
                    static_cast<size_t>(align_len));
                const int32_t endpoint_offset = is_left
                    ? query_start
                    : (consensus_len - (query_start + align_len));

                for (int32_t anchor_pos = search_lo; anchor_pos <= search_hi; ++anchor_pos) {
                    const int32_t ref_start = is_left ? (anchor_pos - align_len) : anchor_pos;
                    const int32_t ref_end = is_left ? anchor_pos : (anchor_pos + align_len);
                    if (ref_start < ref_window_start || ref_end > ref_window_end) {
                        continue;
                    }

                    const size_t offset = static_cast<size_t>(ref_start - ref_window_start);
                    const std::string_view ref_seq(
                        ref_window.data() + offset,
                        static_cast<size_t>(align_len));
                    double identity = 0.0;
                    if (!edit_identity_if_at_least(
                            query,
                            ref_seq,
                            max_edits,
                            edit_workspace,
                            identity)) {
                        continue;
                    }

                    EventFlankPlacement placement;
                    placement.query_start = query_start;
                    placement.query_end = query_start + align_len;
                    placement.ref_start = ref_start;
                    placement.ref_end = ref_end;
                    placement.align_len = align_len;
                    placement.identity = identity;
                    placement.breakpoint_delta = is_left
                        ? std::abs(ref_end - breakpoint)
                        : std::abs(ref_start - breakpoint);
                    placement.endpoint_offset = endpoint_offset;
                    placements.push_back(placement);
                }
            }

            if (placements.empty()) {
                continue;
            }

            std::sort(
                placements.begin(),
                placements.end(),
                better_event_flank_placement);
            const double best_identity = placements.front().identity;
            for (const auto& placement : placements) {
                if ((placement.identity + kEventSegmentationUniquenessMargin) < best_identity) {
                    break;
                }
                result.candidates.push_back(placement);
            }
        }
        std::sort(
            result.candidates.begin(),
            result.candidates.end(),
            better_event_flank_placement);
        return result;
    };

    EventFlankSearchResult left_search = collect_candidates(
        true,
        bp_left,
        left_window_start,
        left_ref_window,
        0);
    if (left_search.candidates.empty() && kEventSegmentationEndpointSlackBp > 0) {
        left_search = collect_candidates(
            true,
            bp_left,
            left_window_start,
            left_ref_window,
            kEventSegmentationEndpointSlackBp);
    }
    if (left_search.candidates.empty()) {
        segmentation.qc_reason = "NO_LEFT_FLANK_MATCH";
        return segmentation;
    }

    EventFlankSearchResult right_search = collect_candidates(
        false,
        bp_right,
        right_window_start,
        right_ref_window,
        0);
    if (right_search.candidates.empty() && kEventSegmentationEndpointSlackBp > 0) {
        right_search = collect_candidates(
            false,
            bp_right,
            right_window_start,
            right_ref_window,
            kEventSegmentationEndpointSlackBp);
    }
    if (right_search.candidates.empty()) {
        segmentation.qc_reason = "NO_RIGHT_FLANK_MATCH";
        return segmentation;
    }

    EventFlankPlacement best_left;
    EventFlankPlacement best_right;
    bool have_pair = false;
    auto better_pair = [&](const EventFlankPlacement& lhs_left,
                           const EventFlankPlacement& lhs_right,
                           const EventFlankPlacement& rhs_left,
                           const EventFlankPlacement& rhs_right) {
        const double lhs_min_identity = std::min(lhs_left.identity, lhs_right.identity);
        const double rhs_min_identity = std::min(rhs_left.identity, rhs_right.identity);
        if (lhs_min_identity != rhs_min_identity) {
            return lhs_min_identity > rhs_min_identity;
        }

        const double lhs_mean_identity = 0.5 * (lhs_left.identity + lhs_right.identity);
        const double rhs_mean_identity = 0.5 * (rhs_left.identity + rhs_right.identity);
        if (lhs_mean_identity != rhs_mean_identity) {
            return lhs_mean_identity > rhs_mean_identity;
        }

        const int32_t lhs_delta = lhs_left.breakpoint_delta + lhs_right.breakpoint_delta;
        const int32_t rhs_delta = rhs_left.breakpoint_delta + rhs_right.breakpoint_delta;
        if (lhs_delta != rhs_delta) {
            return lhs_delta < rhs_delta;
        }

        const int32_t lhs_total = lhs_left.align_len + lhs_right.align_len;
        const int32_t rhs_total = rhs_left.align_len + rhs_right.align_len;
        if (lhs_total != rhs_total) {
            return lhs_total > rhs_total;
        }

        if (lhs_left.ref_start != rhs_left.ref_start) {
            return lhs_left.ref_start < rhs_left.ref_start;
        }
        return lhs_right.ref_start < rhs_right.ref_start;
    };

    for (const auto& left_candidate : left_search.candidates) {
        for (const auto& right_candidate : right_search.candidates) {
            if (left_candidate.query_end + kEventSegmentationMinInsertBp >
                right_candidate.query_start) {
                continue;
            }

            if (!have_pair ||
                better_pair(left_candidate, right_candidate, best_left, best_right)) {
                best_left = left_candidate;
                best_right = right_candidate;
                have_pair = true;
            }
        }
    }

    if (!have_pair) {
        segmentation.qc_reason = "NO_TRIPARTITE_EVENT_SEGMENTATION";
        return segmentation;
    }

    const int32_t insert_start = best_left.query_end;
    const int32_t insert_len = best_right.query_start - best_left.query_end;
    if (insert_len < kEventSegmentationMinInsertBp) {
        segmentation.qc_reason = "EMPTY_EVENT_INSERT_SEGMENT";
        return segmentation;
    }

    segmentation.left_flank_seq = consensus.substr(
        static_cast<size_t>(best_left.query_start),
        static_cast<size_t>(best_left.align_len));
    segmentation.insert_seq = consensus.substr(
        static_cast<size_t>(insert_start),
        static_cast<size_t>(insert_len));
    segmentation.right_flank_seq = consensus.substr(
        static_cast<size_t>(best_right.query_start),
        static_cast<size_t>(best_right.align_len));
    segmentation.left_ref_start = best_left.ref_start;
    segmentation.left_ref_end = best_left.ref_end;
    segmentation.right_ref_start = best_right.ref_start;
    segmentation.right_ref_end = best_right.ref_end;
    segmentation.left_flank_align_len = best_left.align_len;
    segmentation.right_flank_align_len = best_right.align_len;
    segmentation.left_flank_identity = best_left.identity;
    segmentation.right_flank_identity = best_right.identity;
    segmentation.pass = true;
    segmentation.qc_reason = "PASS_EVENT_SEGMENTATION";
    return segmentation;
}

EventSegmentationEvidence Pipeline::analyze_event_segmentation(
    const EventConsensus& event_consensus,
    const EventSegmentation& event_segmentation) const {
    return analyze_event_segmentation_for_test(
        event_consensus.qc_pass,
        event_segmentation);
}

TEAlignmentEvidence Pipeline::align_insert_seq_to_te(
    const EventSegmentation& event_segmentation) const {
    if (!event_segmentation.pass) {
        TEAlignmentEvidence evidence;
        evidence.qc_reason = "NO_EVENT_SEGMENTATION_FOR_TE_ALIGNMENT";
        return evidence;
    }
    return te_classifier_module_.align_insert_sequence(event_segmentation.insert_seq);
}

void finalize_final_calls(PipelineResult& result) {
    std::sort(result.final_calls.begin(), result.final_calls.end(), final_call_sort_less);

    std::vector<FinalCall> deduped;
    deduped.reserve(result.final_calls.size());

    size_t cluster_start = 0;
    while (cluster_start < result.final_calls.size()) {
        const FinalCall* best = &result.final_calls[cluster_start];
        size_t cluster_end = cluster_start + 1;
        while (cluster_end < result.final_calls.size() &&
               same_call_locus(result.final_calls[cluster_start], result.final_calls[cluster_end])) {
            if (prefer_new_call(result.final_calls[cluster_end], *best)) {
                best = &result.final_calls[cluster_end];
            }
            ++cluster_end;
        }
        deduped.push_back(*best);
        cluster_start = cluster_end;
    }

    result.final_calls = std::move(deduped);

    result.final_pass_calls = static_cast<int64_t>(result.final_calls.size());
}

std::vector<ComponentCall> LinearBinComponentModule::build(
    const std::vector<const bam1_t*>& bin_records,
    const std::string& chrom,
    int32_t tid,
    int32_t bin_start,
    int32_t bin_end) const {
    if (bin_records.empty()) {
        return {};
    }

    const EvidenceBundle evidence_bundle = extract_evidence_points(bin_records, tid);
    const auto windows = build_density_windows(evidence_bundle.points, bin_start, bin_end);
    if (windows.empty()) {
        return {};
    }

    std::vector<std::vector<const EvidencePoint*>> evidence_by_read(bin_records.size());
    for (const auto& point : evidence_bundle.points) {
        if (point.read_index < evidence_by_read.size()) {
            evidence_by_read[point.read_index].push_back(&point);
        }
    }

    std::vector<int32_t> assigned_window(bin_records.size(), -1);
    std::vector<int32_t> assigned_count(windows.size(), 0);
    std::vector<int32_t> ambiguous_count(windows.size(), 0);

    for (size_t read_idx = 0; read_idx < evidence_by_read.size(); ++read_idx) {
        const auto& read_points = evidence_by_read[read_idx];
        if (read_points.empty()) {
            continue;
        }

        int32_t best_window = -1;
        ReadWindowAssignmentScore best_score;
        ReadWindowAssignmentScore second_score;

        for (size_t wi = 0; wi < windows.size(); ++wi) {
            const auto& window = windows[wi];
            ReadWindowAssignmentScore score;
            for (const EvidencePoint* point : read_points) {
                if (!point) {
                    continue;
                }
                const int32_t dist = std::abs(point->pos - window.center);
                const double proximity =
                    (point->pos >= window.start && point->pos <= window.end)
                    ? 1.0
                    : std::exp(-static_cast<double>(dist) / kReadAssignDecayBp);
                accumulate_read_window_assignment_score(
                    *point,
                    point->weight * proximity,
                    score);
            }

            if (better_read_window_assignment_score(score, best_score)) {
                second_score = best_score;
                best_score = score;
                best_window = static_cast<int32_t>(wi);
            } else if (better_read_window_assignment_score(score, second_score)) {
                second_score = score;
            }
        }

        if (best_window < 0 ||
            !best_score.valid ||
            best_score.specificity_score < kReadAssignMinScore) {
            continue;
        }
        if (second_score.valid &&
            best_score.specificity_rank == second_score.specificity_rank &&
            second_score.specificity_score >=
                best_score.specificity_score * kReadAmbiguousRatio) {
            ambiguous_count[static_cast<size_t>(best_window)] += 1;
            continue;
        }

        assigned_window[read_idx] = best_window;
        assigned_count[static_cast<size_t>(best_window)] += 1;
    }

    std::vector<ComponentCall> components;
    components.reserve(windows.size());

    for (size_t wi = 0; wi < windows.size(); ++wi) {
        if (assigned_count[wi] <= 0) {
            continue;
        }

        ComponentCall call;
        call.chrom = chrom;
        call.tid = tid;
        call.bin_start = windows[wi].start;
        call.bin_end = windows[wi].end;
        call.anchor_pos = windows[wi].center;
        call.peak_weight = windows[wi].peak_weight;

        std::vector<std::pair<int32_t, double>> anchor_points;

        const auto better_point = [&](const EvidencePoint* candidate, const EvidencePoint* incumbent) {
            if (candidate == nullptr) {
                return false;
            }
            if (incumbent == nullptr) {
                return true;
            }
            const int32_t candidate_dist = std::abs(candidate->pos - windows[wi].center);
            const int32_t incumbent_dist = std::abs(incumbent->pos - windows[wi].center);
            if (candidate_dist != incumbent_dist) {
                return candidate_dist < incumbent_dist;
            }
            if (candidate->weight != incumbent->weight) {
                return candidate->weight > incumbent->weight;
            }
            return candidate->pos < incumbent->pos;
        };

        const auto append_point = [&](const EvidencePoint& point) {
            anchor_points.push_back({point.pos, point.weight});
            switch (point.kind) {
                case EvidenceKind::kSoftClip:
                    call.evidence_soft_clip_count += 1;
                    break;
                case EvidenceKind::kIndel:
                    call.evidence_indel_count += 1;
                    break;
                case EvidenceKind::kSAHint:
                    call.evidence_sa_hint_count += 1;
                    break;
            }

            BreakpointCandidate candidate;
            candidate.chrom = chrom;
            candidate.pos = point.pos;
            candidate.read_index = point.read_index;
            candidate.class_mask = point.class_mask;
            if (point.read_index < evidence_bundle.read_summaries.size()) {
                const auto& summary = evidence_bundle.read_summaries[point.read_index];
                candidate.read_id = summary.read_id;
                candidate.is_reverse = summary.is_reverse;
                candidate.anchor_len = summary.max_soft_clip;
            }
            if (point.kind == EvidenceKind::kSoftClip) {
                candidate.clip_len = point.signal_len;
            } else if (point.kind == EvidenceKind::kIndel) {
                candidate.ins_len = point.signal_len;
            }
            call.breakpoint_candidates.push_back(std::move(candidate));
        };

        const auto point_is_local_to_window = [&](const EvidencePoint* point) {
            if (point == nullptr) {
                return false;
            }
            return point->pos >= windows[wi].start && point->pos <= windows[wi].end;
        };

        for (size_t read_idx = 0; read_idx < assigned_window.size(); ++read_idx) {
            if (assigned_window[read_idx] != static_cast<int32_t>(wi)) {
                continue;
            }
            if (read_idx >= evidence_by_read.size()) {
                continue;
            }
            const auto& read_points = evidence_by_read[read_idx];
            const EvidencePoint* best_softclip = nullptr;
            const EvidencePoint* best_indel = nullptr;
            const EvidencePoint* best_sa_hint = nullptr;
            for (const EvidencePoint* point : read_points) {
                if (!point || point->read_index >= assigned_window.size()) {
                    continue;
                }
                if (assigned_window[point->read_index] != static_cast<int32_t>(wi)) {
                    continue;
                }
                if (!point_is_local_to_window(point)) {
                    continue;
                }
                switch (point->kind) {
                    case EvidenceKind::kSoftClip:
                        if (better_point(point, best_softclip)) {
                            best_softclip = point;
                        }
                        break;
                    case EvidenceKind::kIndel:
                        if (better_point(point, best_indel)) {
                            best_indel = point;
                        }
                        break;
                    case EvidenceKind::kSAHint:
                        if (better_point(point, best_sa_hint)) {
                            best_sa_hint = point;
                        }
                        break;
                }
            }

            if (best_softclip == nullptr && best_indel == nullptr && best_sa_hint == nullptr) {
                continue;
            }

            call.read_indices.push_back(read_idx);
            if (best_softclip != nullptr) {
                call.soft_clip_read_indices.push_back(read_idx);
                append_point(*best_softclip);
            }
            if (best_indel != nullptr) {
                call.insertion_read_indices.push_back(read_idx);
                append_point(*best_indel);
            }
            if (best_sa_hint != nullptr) {
                call.split_sa_read_indices.push_back(read_idx);
                append_point(*best_sa_hint);
            }
        }

        if (!anchor_points.empty()) {
            const int32_t anchor = weighted_median_position(std::move(anchor_points));
            if (anchor >= 0) {
                call.anchor_pos = anchor;
            }
        }

        if (!call.read_indices.empty()) {
            components.push_back(std::move(call));
        }
    }

    return components;
}

Pipeline::Pipeline(PipelineConfig config, std::unique_ptr<BamStreamReader> bam_reader)
    : config_(std::move(config)),
      bam_reader_(std::move(bam_reader)),
      gate1_module_(),
      component_module_(),
      ins_fragment_module_(config_),
      te_classifier_module_(config_),
      tsd_detector_(config_) {}

PipelineResult Pipeline::run() const {
    if (!bam_reader_ || !bam_reader_->is_valid()) {
        throw std::runtime_error("BAM reader is not valid");
    }
    if (!bam_reader_->can_fetch()) {
        throw std::runtime_error(
            "BAM reader does not support indexed local fetch; a BAM index is required for event-level recollection");
    }
    if (!tsd_detector_.can_fetch_reference()) {
        throw std::runtime_error(
            "Reference FASTA cannot be fetched; event consensus segmentation requires a readable indexed reference");
    }
    if (!te_classifier_module_.is_enabled()) {
        throw std::runtime_error(
            "TE library alignment is unavailable; Stage 5 requires a readable TE FASTA and a usable shortlist index");
    }

    {
        std::ostringstream oss;
        oss << "[Pipeline] starting mode="
            << (config_.enable_parallel ? "parallel" : "streaming")
            << " bam_threads=" << config_.bam_threads
            << " progress_interval=" << config_.progress_interval
            << " log_stage_bins=" << bool_as_int(config_.log_stage_bins)
            << " log_stage_components=" << bool_as_int(config_.log_stage_components);
        emit_pipeline_log_line(oss.str());
    }

    if (config_.enable_parallel) {
        return run_parallel();
    }
    return run_streaming();
}

PipelineResult Pipeline::run_streaming() const {
    PipelineResult result;
    StreamingState state;

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    result.total_reads = bam_reader_->stream(
        [this, &state, &result](BamRecordPtr&& record) {
            consume_record(std::move(record), state, result);
        },
        progress_cb,
        config_.progress_interval);

    flush_current_bin(state, result);
    finalize_final_calls(result);
    return result;
}

PipelineResult Pipeline::run_parallel() const {
    PipelineResult result;
    const int32_t worker_count = std::max(
        1,
        (config_.parallel_workers > 0) ? config_.parallel_workers : config_.bam_threads);
    const int32_t queue_cap_cfg = config_.parallel_queue_max_tasks;
    const size_t queue_cap = (queue_cap_cfg > 0) ? static_cast<size_t>(queue_cap_cfg) : size_t{0};
    SafeQueue<BinTask> bin_queue(queue_cap);
    std::vector<PipelineResult> worker_results(static_cast<size_t>(worker_count));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(worker_count));
    {
        std::ostringstream oss;
        oss << "[Pipeline] parallel workers=" << worker_count
            << " queue_cap=" << queue_cap;
        emit_pipeline_log_line(oss.str());
    }

    for (int32_t wi = 0; wi < worker_count; ++wi) {
        workers.emplace_back([this, wi, &bin_queue, &worker_results]() {
            BinTask task;
            while (bin_queue.pop(task)) {
                process_bin_records(
                    std::move(task.records),
                    task.tid,
                    task.bin_index,
                    worker_results[static_cast<size_t>(wi)]);
                for (size_t oi = 0; oi < task.owner_count; ++oi) {
                    task.owners[oi].reset();
                }
                task.owner_count = 0;
            }
        });
    }

    int64_t gate1_passed = 0;
    int32_t current_tid = -1;
    int32_t current_bin_index = -1;
    std::vector<BufferedRecord> current_bin_records;
    current_bin_records.reserve(config_.batch_size);
    struct BinSnapshot {
        int32_t tid = -1;
        int32_t bin_index = -1;
        SharedRecordBatch records;
    };
    std::deque<BinSnapshot> recent_bin_snapshots;

    const auto trim_recent_snapshots = [&](int32_t tid, int32_t bin_index) {
        while (!recent_bin_snapshots.empty()) {
            const auto& front = recent_bin_snapshots.front();
            if (front.tid != tid) {
                recent_bin_snapshots.pop_front();
                continue;
            }
            if ((bin_index - front.bin_index) > kCrossBinContextBins) {
                recent_bin_snapshots.pop_front();
                continue;
            }
            break;
        }
    };

    const auto flush_bin_task = [&]() {
        if (current_bin_records.empty()) {
            return;
        }
        SharedRecordBatch current_batch = std::make_shared<const std::vector<BufferedRecord>>(
            std::move(current_bin_records));
        BinTask task;
        task.tid = current_tid;
        task.bin_index = current_bin_index;
        task.owner_count = 0;
        const int32_t bin_start = current_bin_index * config_.bin_size;
        const int32_t snapshot_min_ref_end = bin_start - kWindowBinSlackBp;
        for (const auto& snapshot : recent_bin_snapshots) {
            if (snapshot.tid != current_tid) {
                continue;
            }
            if ((current_bin_index - snapshot.bin_index) > kCrossBinContextBins) {
                continue;
            }
            if (snapshot.records && task.owner_count < task.owners.size()) {
                task.owners[task.owner_count++] = snapshot.records;
            }
        }
        if (current_batch && task.owner_count < task.owners.size()) {
            task.owners[task.owner_count++] = current_batch;
        }
        size_t total_records = 0;
        for (size_t oi = 0; oi < task.owner_count; ++oi) {
            const auto& owner = task.owners[oi];
            if (owner) {
                total_records += owner->size();
            }
        }
        task.records.reserve(total_records);
        for (size_t oi = 0; oi < task.owner_count; ++oi) {
            const auto& owner = task.owners[oi];
            if (!owner) {
                continue;
            }
            const int32_t min_ref_end = (oi + 1 == task.owner_count)
                ? std::numeric_limits<int32_t>::min()
                : snapshot_min_ref_end;
            append_record_ptrs(*owner, task.records, min_ref_end);
        }
        if (!task.records.empty()) {
            bin_queue.push(std::move(task));
        }

        BinSnapshot snapshot;
        snapshot.tid = current_tid;
        snapshot.bin_index = current_bin_index;
        snapshot.records = std::move(current_batch);
        recent_bin_snapshots.push_back(std::move(snapshot));
        trim_recent_snapshots(current_tid, current_bin_index);

        current_bin_records.clear();
        current_bin_records.reserve(config_.batch_size);
    };

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    const int64_t total_reads = bam_reader_->stream(
        [this, &gate1_passed, &current_tid, &current_bin_index, &current_bin_records, &flush_bin_task](
            BamRecordPtr&& record) {
            if (!record) {
                return;
            }

            ReadView view(record.get());
            if (!gate1_module_.pass_preliminary(view)) {
                return;
            }
            ++gate1_passed;

            const int32_t bin_index = (config_.bin_size > 0) ? (view.pos() / config_.bin_size) : 0;
            if (current_tid < 0) {
                current_tid = view.tid();
                current_bin_index = bin_index;
            }

            if (view.tid() != current_tid || bin_index != current_bin_index) {
                flush_bin_task();
                current_tid = view.tid();
                current_bin_index = bin_index;
            }
            BufferedRecord buffered;
            buffered.ref_end = compute_ref_end(view);
            buffered.record = std::move(record);
            current_bin_records.push_back(std::move(buffered));
        },
        progress_cb,
        config_.progress_interval);

    flush_bin_task();

    bin_queue.close();
    for (auto& worker : workers) {
        worker.join();
    }

    result.total_reads = total_reads;
    result.gate1_passed = gate1_passed;

    for (auto& worker_result : worker_results) {
        result.processed_bins += worker_result.processed_bins;
        result.built_components += worker_result.built_components;
        result.event_consensus_calls += worker_result.event_consensus_calls;
        result.genotype_calls += worker_result.genotype_calls;
        for (auto& call : worker_result.final_calls) {
            result.final_calls.push_back(std::move(call));
        }
    }

    finalize_final_calls(result);
    return result;
}

void Pipeline::consume_record(
    BamRecordPtr&& record,
    StreamingState& state,
    PipelineResult& result) const {
    if (!record) {
        return;
    }

    ReadView view(record.get());
    if (!gate1_module_.pass_preliminary(view)) {
        return;
    }

    result.gate1_passed += 1;

    while (!state.active_window.empty()) {
        const auto& front = state.active_window.front();
        const bool different_tid = (view.tid() != front.tid);
        const bool beyond_window = (!different_tid && (view.pos() - front.pos > config_.window_size));
        if (!different_tid && !beyond_window) {
            break;
        }
        state.active_window.pop_front();
    }
    state.active_window.push_back({view.tid(), view.pos()});

    const int32_t bin_index = (config_.bin_size > 0) ? (view.pos() / config_.bin_size) : 0;

    if (state.current_tid < 0) {
        state.current_tid = view.tid();
        state.current_bin_index = bin_index;
    }

    if (view.tid() != state.current_tid || bin_index != state.current_bin_index) {
        flush_current_bin(state, result);
        state.current_tid = view.tid();
        state.current_bin_index = bin_index;
    }

    BufferedRecord buffered;
    buffered.ref_end = compute_ref_end(view);
    buffered.record = std::move(record);
    state.current_bin_records.push_back(std::move(buffered));
}

void Pipeline::flush_current_bin(
    StreamingState& state,
    PipelineResult& result) const {
    if (state.current_bin_records.empty()) {
        return;
    }

    std::vector<const bam1_t*> merged_records;
    const int32_t bin_start = state.current_bin_index * config_.bin_size;
    const int32_t snapshot_min_ref_end = bin_start - kWindowBinSlackBp;
    for (const auto& snapshot : state.recent_bin_snapshots) {
        if (snapshot.tid != state.current_tid) {
            continue;
        }
        if ((state.current_bin_index - snapshot.bin_index) > kCrossBinContextBins) {
            continue;
        }
        append_record_ptrs(snapshot.records, merged_records, snapshot_min_ref_end);
    }
    append_record_ptrs(
        state.current_bin_records,
        merged_records,
        std::numeric_limits<int32_t>::min());

    if (!merged_records.empty()) {
        process_bin_records(
            std::move(merged_records),
            state.current_tid,
            state.current_bin_index,
            result);
    }

    StreamingState::BinSnapshot snapshot;
    snapshot.tid = state.current_tid;
    snapshot.bin_index = state.current_bin_index;
    snapshot.records = std::move(state.current_bin_records);
    state.recent_bin_snapshots.push_back(std::move(snapshot));
    while (!state.recent_bin_snapshots.empty()) {
        const auto& front = state.recent_bin_snapshots.front();
        if (front.tid != state.current_tid) {
            state.recent_bin_snapshots.pop_front();
            continue;
        }
        if ((state.current_bin_index - front.bin_index) > kCrossBinContextBins) {
            state.recent_bin_snapshots.pop_front();
            continue;
        }
        break;
    }

    state.current_bin_records.clear();
}

GenotypeCall Pipeline::genotype_call(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence) const {
    GenotypeCall call;
    call.tid = component.tid;
    call.pos = component.anchor_pos;

    EventGenotypeInput input;
    input.alt_struct_reads = event_evidence.alt_struct_reads;
    input.ref_span_reads = event_evidence.ref_span_reads;
    input.min_depth = config_.genotype_min_depth;
    input.error_rate = config_.genotype_error_rate;
    const EventGenotypeDecision decision = genotype_event_from_alt_vs_ref(input);

    call.genotype = decision.best_gt;
    call.af = decision.allele_fraction;
    call.gq = decision.gq;
    return call;
}

BoundaryEvidence Pipeline::analyze_boundary(
    const EventSegmentation& event_segmentation,
    int32_t breakpoint_envelope_width) const {
    FinalBoundaryInput boundary_input;
    boundary_input.left_ref_start = event_segmentation.left_ref_start;
    boundary_input.left_ref_end = event_segmentation.left_ref_end;
    boundary_input.right_ref_start = event_segmentation.right_ref_start;
    boundary_input.right_ref_end = event_segmentation.right_ref_end;
    boundary_input.tsd_min_len = std::max(1, config_.tsd_min_len);
    boundary_input.tsd_max_len = std::max(boundary_input.tsd_min_len, config_.tsd_max_len);
    return evaluate_boundary_evidence(boundary_input, breakpoint_envelope_width);
}

FinalCall Pipeline::emit_final_te_call(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence,
    const EventConsensus& event_consensus,
    const EventSegmentation& event_segmentation,
    const TEAlignmentEvidence& te_alignment,
    const GenotypeCall& genotype,
    const FinalBoundaryDecision& boundary,
    const FinalTeAcceptanceDecision& acceptance) const {
    FinalCall call;
    call.chrom = component.chrom;
    call.tid = component.tid;
    call.bp_left = std::min(
        event_segmentation.left_ref_end,
        event_segmentation.right_ref_start);
    call.bp_right = std::max(
        event_segmentation.left_ref_end,
        event_segmentation.right_ref_start);
    call.pos = (call.bp_left >= 0 && call.bp_right >= 0)
        ? (call.bp_left + ((call.bp_right - call.bp_left) / 2))
        : component.anchor_pos;
    call.window_start = component.bin_start;
    call.window_end = component.bin_end;

    call.family = te_alignment.best_family.empty() ? "NA" : te_alignment.best_family;
    call.subfamily = te_alignment.best_subfamily.empty() ? "NA" : te_alignment.best_subfamily;
    call.te_name = (call.subfamily != "NA") ? call.subfamily : call.family;
    call.strand = "NA";
    call.insert_len = static_cast<int32_t>(event_segmentation.insert_seq.size());
    call.best_te_identity = te_alignment.best_identity;
    call.best_te_query_coverage = te_alignment.best_query_coverage;
    call.cross_family_margin = te_alignment.cross_family_margin;
    call.left_flank_align_len = event_segmentation.left_flank_align_len;
    call.right_flank_align_len = event_segmentation.right_flank_align_len;
    call.event_consensus_len = event_consensus.consensus_len;

    call.te_qc = te_alignment.qc_reason;
    call.final_qc = acceptance.qc;
    append_qc_token(call.final_qc, te_alignment.qc_reason.c_str());
    append_qc_token(call.final_qc, boundary.qc.c_str());
    append_qc_token(call.final_qc, event_consensus.qc_reason.c_str());
    append_qc_token(call.final_qc, event_segmentation.qc_reason.c_str());

    call.tsd_type = boundary.boundary_type;
    call.tsd_len = boundary.boundary_len;
    call.tsd_seq = "NA";
    call.tsd_bg_p = 1.0;

    call.support_reads = event_evidence.alt_struct_reads;
    call.alt_struct_reads = event_evidence.alt_struct_reads;
    call.ref_span_reads = event_evidence.ref_span_reads;
    call.low_mapq_ref_span_reads = event_evidence.low_mapq_ref_span_reads;

    call.genotype = genotype.genotype;
    call.af = genotype.af;
    call.gq = genotype.gq;
    return call;
}

void Pipeline::process_bin_records(
    std::vector<const bam1_t*>&& bin_records,
    int32_t tid,
    int32_t bin_index,
    PipelineResult& result) const {
    if (bin_records.empty()) {
        return;
    }

    result.processed_bins += 1;

    const int32_t bin_start = bin_index * config_.bin_size;
    const int32_t bin_end = bin_start + config_.bin_size;
    const std::string chrom = bam_reader_->chromosome_name(tid);

    auto components = component_module_.build(bin_records, chrom, tid, bin_start, bin_end);
    result.built_components += static_cast<int64_t>(components.size());

    struct BinStageStats {
        int64_t event_consensus_calls = 0;
        int64_t event_consensus_rejected = 0;
        int64_t genotype_calls = 0;
        int64_t final_calls = 0;
    } bin_stats;

    if (config_.log_stage_bins) {
        std::ostringstream oss;
        oss << "[Pipeline][bin][start]"
            << " chrom=" << chrom
            << " tid=" << tid
            << " start=" << bin_start
            << " end=" << bin_end
            << " records=" << bin_records.size()
            << " components=" << components.size();
        emit_pipeline_log_line(oss.str());
    }

    const auto log_component = [&](
        const char* outcome,
        const ComponentCall& component,
        const std::vector<InsertionFragment>& fragments,
        const GenotypeCall* genotype,
        const FinalCall* final_call) {
        if (!config_.log_stage_components) {
            return;
        }

        std::ostringstream oss;
        oss << "[Pipeline][component]"
            << " outcome=" << outcome
            << " chrom=" << component.chrom
            << " anchor=" << component.anchor_pos
            << " bin=" << component.bin_start << "-" << component.bin_end
            << " reads=" << component.read_indices.size()
            << " softclip_reads=" << component.soft_clip_read_indices.size()
            << " split_reads=" << component.split_sa_read_indices.size()
            << " indel_reads=" << component.insertion_read_indices.size()
            << " breakpoints=" << component.breakpoint_candidates.size()
            << " fragments=" << fragments.size();

        if (genotype != nullptr) {
            oss << " gt=" << display_or_na(genotype->genotype)
                << " af=" << genotype->af
                << " gq=" << genotype->gq;
        }

        if (final_call != nullptr) {
            oss << " final_te=" << display_or_na(final_call->te_name)
                << " final_support=" << final_call->support_reads
                << " final_qc=" << display_or_na(final_call->final_qc)
                << " te_qc=" << display_or_na(final_call->te_qc);
        }

        emit_pipeline_log_line(oss.str());
    };

    for (const auto& component : components) {
        const std::vector<InsertionFragment> empty_fragments;

        const auto seed_bounds = infer_component_breakpoint_bounds(component);
        const int32_t fetch_start = std::max(
            0,
            std::min(seed_bounds.first, seed_bounds.second) - kLocalEventFetchSlackBp);
        const int32_t fetch_end = std::max(
            fetch_start + 1,
            std::max(seed_bounds.first, seed_bounds.second) + kLocalEventFetchSlackBp);
        LocalFetchedReads local_reads;
        const bool fetch_ok = bam_reader_->fetch(
            component.chrom,
            fetch_start,
            fetch_end,
            [&](BamRecordPtr&& record) {
                if (!record) {
                    return true;
                }
                local_reads.owned_records.push_back(std::move(record));
                const bam1_t* ptr = local_reads.owned_records.back().get();
                local_reads.records.push_back(ptr);
                ReadView view(ptr);
                ReadReferenceSpan span;
                span.valid = true;
                span.tid = view.tid();
                span.start = view.pos();
                span.end = compute_ref_end(view);
                local_reads.read_spans.push_back(span);
                return true;
            });
        if (!fetch_ok || local_reads.records.empty()) {
            log_component(
                "LOCAL_EVENT_RECOLLECTION_FAILED",
                component,
                empty_fragments,
                nullptr,
                nullptr);
            continue;
        }

        const ComponentCall local_component = build_local_fragment_component(
            config_,
            component,
            local_reads.records);
        const std::vector<InsertionFragment> fragments = ins_fragment_module_.extract(
            local_component,
            local_reads.records);

        const auto breakpoint_hypotheses = enumerate_breakpoint_hypotheses(
            component,
            local_reads.records,
            fragments,
            seed_bounds.first,
            seed_bounds.second,
            3);

        bool have_best = false;
        double best_total = -1e9;
        JointDecisionResult best_joint;
        EventReadEvidence best_event_evidence;
        EventConsensus best_event_consensus;
        EventSegmentation best_event_segmentation;
        TEAlignmentEvidence best_te_alignment;
        BoundaryEvidence best_boundary_evidence;
        GenotypeCall best_genotype;

        for (const auto& hypothesis : breakpoint_hypotheses) {
            const EventReadEvidence event_evidence = collect_event_read_evidence_for_bounds(
                component,
                local_reads.records,
                local_reads.read_spans,
                fragments,
                hypothesis.left,
                hypothesis.right,
                hypothesis.left,
                hypothesis.right);
            const EventConsensus event_consensus = build_event_consensus(
                component,
                local_reads.records,
                fragments,
                event_evidence);
            result.event_consensus_calls += 1;
            bin_stats.event_consensus_calls += 1;
            if (!event_consensus.qc_pass) {
                bin_stats.event_consensus_rejected += 1;
            }

            const GenotypeCall genotype = genotype_call(component, event_evidence);
            result.genotype_calls += 1;
            bin_stats.genotype_calls += 1;

            EventGenotypeInput existence_input;
            existence_input.alt_struct_reads = event_evidence.alt_struct_reads;
            existence_input.ref_span_reads = event_evidence.ref_span_reads;
            existence_input.min_depth = config_.genotype_min_depth;
            existence_input.min_gq = 20;
            existence_input.error_rate = config_.genotype_error_rate;
            const EventExistenceEvidence existence = build_event_existence_evidence(
                existence_input);

            const EventSegmentation event_segmentation = segment_event_consensus(
                component,
                event_evidence,
                event_consensus);
            const EventSegmentationEvidence seg_evidence = analyze_event_segmentation(
                event_consensus,
                event_segmentation);

            TEAlignmentEvidence te_alignment;
            if (seg_evidence.has_insert_seq) {
                te_alignment = align_insert_seq_to_te(event_segmentation);
            }

            const int32_t breakpoint_envelope_width = std::max(
                1,
                event_evidence.bp_right - event_evidence.bp_left);
            BoundaryEvidence boundary_evidence;
            if (event_segmentation.pass) {
                boundary_evidence = analyze_boundary(
                    event_segmentation,
                    breakpoint_envelope_width);
            }

            const JointDecisionResult joint = evaluate_joint_hypotheses(
                existence,
                seg_evidence,
                te_alignment,
                boundary_evidence);
            const double candidate_total = joint.best.hard_veto ? -1e9 : joint.best.total;
            if (!have_best || candidate_total > best_total) {
                have_best = true;
                best_total = candidate_total;
                best_joint = joint;
                best_event_evidence = event_evidence;
                best_event_consensus = event_consensus;
                best_event_segmentation = event_segmentation;
                best_te_alignment = te_alignment;
                best_boundary_evidence = boundary_evidence;
                best_genotype = genotype;
            }
        }

        if (config_.log_stage_components && have_best) {
            std::ostringstream oss;
            oss << "[Pipeline][joint]"
                << " chrom=" << component.chrom
                << " anchor=" << component.anchor_pos
                << " bp=" << best_event_evidence.bp_left << "-" << best_event_evidence.bp_right
                << " gt=" << display_or_na(best_genotype.genotype)
                << " af=" << best_genotype.af
                << " gq=" << best_genotype.gq
                << " alt_reads=" << best_event_evidence.alt_struct_reads
                << " ref_reads=" << best_event_evidence.ref_span_reads
                << " consensus_qc=" << display_or_na(best_event_consensus.qc_reason)
                << " consensus_len=" << best_event_consensus.consensus_len
                << " seg_qc=" << display_or_na(best_event_segmentation.qc_reason)
                << " seg_insert_len=" << best_event_segmentation.insert_seq.size()
                << " seg_left=" << best_event_segmentation.left_flank_align_len
                << " seg_right=" << best_event_segmentation.right_flank_align_len
                << " te_qc=" << display_or_na(best_te_alignment.qc_reason)
                << " te_pass=" << bool_as_int(best_te_alignment.pass)
                << " te_family=" << display_or_na(best_te_alignment.best_family)
                << " te_subfamily=" << display_or_na(best_te_alignment.best_subfamily)
                << " te_identity=" << best_te_alignment.best_identity
                << " te_cov=" << best_te_alignment.best_query_coverage
                << " te_margin=" << best_te_alignment.cross_family_margin
                << " boundary_qc=" << display_or_na(best_boundary_evidence.qc)
                << " boundary_defined=" << bool_as_int(best_boundary_evidence.geometry_defined)
                << " boundary_canonical=" << bool_as_int(best_boundary_evidence.canonical_pass)
                << " boundary_consistent=" << bool_as_int(best_boundary_evidence.evidence_consistent)
                << " boundary_type=" << display_or_na(best_boundary_evidence.boundary_type)
                << " boundary_len=" << best_boundary_evidence.boundary_len
                << " best_kind=" << final_hypothesis_kind_name(best_joint.best.kind)
                << " best_total=" << best_joint.best.total
                << " best_exist=" << best_joint.best.existence
                << " best_seg=" << best_joint.best.segmentation
                << " best_te=" << best_joint.best.te
                << " best_boundary=" << best_joint.best.boundary
                << " runner_kind=" << final_hypothesis_kind_name(best_joint.runner_up.kind)
                << " runner_total=" << best_joint.runner_up.total
                << " final_qc=" << display_or_na(best_joint.final_qc)
                << " emit_te=" << bool_as_int(best_joint.emit_te_call)
                << " emit_unknown=" << bool_as_int(best_joint.emit_unknown_te);
            emit_pipeline_log_line(oss.str());
        }

        if (!have_best || !best_joint.emit_te_call) {
            log_component(
                have_best ? best_joint.final_qc.c_str() : "NO_BREAKPOINT_HYPOTHESIS",
                component,
                fragments,
                have_best ? &best_genotype : nullptr,
                nullptr);
            continue;
        }

        FinalBoundaryDecision boundary;
        boundary.pass =
            best_boundary_evidence.canonical_pass ||
            best_boundary_evidence.evidence_consistent;
        boundary.boundary_type = best_boundary_evidence.boundary_type;
        boundary.boundary_len = best_boundary_evidence.boundary_len;
        boundary.qc = best_boundary_evidence.qc;

        FinalTeAcceptanceDecision acceptance;
        acceptance.pass = true;
        acceptance.qc = best_joint.final_qc;

        FinalCall final_call = emit_final_te_call(
            component,
            best_event_evidence,
            best_event_consensus,
            best_event_segmentation,
            best_te_alignment,
            best_genotype,
            boundary,
            acceptance);
        result.final_calls.push_back(final_call);
        bin_stats.final_calls += 1;

        log_component(
            acceptance.qc.c_str(),
            component,
            fragments,
            &best_genotype,
            &result.final_calls.back());
    }

    if (config_.log_stage_bins) {
        std::ostringstream oss;
        oss << "[Pipeline][bin][done]"
            << " chrom=" << chrom
            << " tid=" << tid
            << " start=" << bin_start
            << " end=" << bin_end
            << " records=" << bin_records.size()
            << " components=" << components.size()
            << " event_consensus_calls=" << bin_stats.event_consensus_calls
            << " event_consensus_rejected=" << bin_stats.event_consensus_rejected
            << " genotype_calls=" << bin_stats.genotype_calls
            << " final_calls=" << bin_stats.final_calls;
        emit_pipeline_log_line(oss.str());
    }
}

std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config) {
    auto reader = make_bam_reader(config.bam_path, config.bam_threads);
    return std::make_unique<Pipeline>(config, std::move(reader));
}

}  // namespace placer
