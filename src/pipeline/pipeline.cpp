#include "pipeline.h"
#include "decision_policy.h"
#include "local_interval_cache.h"
#include "parallel_executor_internal.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <exception>
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

constexpr int32_t kFinalCallDedupDistanceBp = 100;

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
constexpr int32_t kEventSegmentationSeedK = 11;
constexpr int32_t kEventSegmentationSeedBinBp = 8;
constexpr int32_t kEventSegmentationSeedTopBins = 8;

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

std::unordered_map<std::string, const bam1_t*> best_records_by_qname(
    const std::vector<const bam1_t*>& local_records) {
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
    return records_by_qname;
}

bool can_build_event_string_from_fragment(
    const bam1_t* record,
    const InsertionFragment& fragment) {
    if (!record || fragment.sequence.empty() || fragment.start < 0 || fragment.length <= 0) {
        return false;
    }

    ReadView read(record);
    const int32_t insert_start = fragment.start;
    const int32_t insert_end = fragment.start + fragment.length;
    if (insert_end > read.seq_len()) {
        return false;
    }

    const int32_t left_flank_len = insert_start - std::max(0, insert_start - kEventConsensusFlankBp);
    const int32_t right_flank_len = std::max(0, std::min(read.seq_len(), insert_end + kEventConsensusFlankBp) - insert_end);

    switch (fragment.source) {
        case InsertionFragmentSource::kClipRefLeft:
            return right_flank_len >= kEventConsensusMinAnchorBp;
        case InsertionFragmentSource::kClipRefRight:
            return left_flank_len >= kEventConsensusMinAnchorBp;
        case InsertionFragmentSource::kCigarInsertion:
        case InsertionFragmentSource::kSplitSa:
            return left_flank_len >= kEventConsensusMinAnchorBp &&
                right_flank_len >= kEventConsensusMinAnchorBp;
        case InsertionFragmentSource::kUnknown:
        default:
            return false;
    }
}

std::string build_event_string_from_fragment(
    const bam1_t* record,
    const InsertionFragment& fragment) {
    if (!can_build_event_string_from_fragment(record, fragment)) {
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

int32_t exact_bin_index_for_record(const bam1_t* record, int32_t bin_size) {
    if (!record || bin_size <= 0) {
        return 0;
    }
    ReadView view(record);
    if (view.pos() < 0) {
        return 0;
    }
    return view.pos() / bin_size;
}

void append_partial_pipeline_result(PipelineResult& dst, PipelineResult&& src) {
    dst.processed_bins += src.processed_bins;
    dst.built_components += src.built_components;
    dst.event_consensus_calls += src.event_consensus_calls;
    dst.genotype_calls += src.genotype_calls;
    if (!src.final_calls.empty()) {
        dst.final_calls.insert(
            dst.final_calls.end(),
            std::make_move_iterator(src.final_calls.begin()),
            std::make_move_iterator(src.final_calls.end()));
    }
}

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

constexpr int32_t kBreakpointClusterGapBp = 75;
constexpr size_t kBreakpointHypothesisMaxClustersPerSide = 6;
constexpr size_t kBreakpointHypothesisMaxSingleClusters = 2;
constexpr int32_t kBreakpointPairCompatibilityBp = 250;

struct BreakpointPositionCluster {
    int32_t center = -1;
    int32_t start = -1;
    int32_t end = -1;
    int32_t support = 0;
};

std::vector<BreakpointPositionCluster> collect_breakpoint_position_clusters(
    std::vector<int32_t> positions) {
    std::vector<BreakpointPositionCluster> clusters;
    if (positions.empty()) {
        return clusters;
    }
    std::sort(positions.begin(), positions.end());

    size_t cluster_begin = 0;
    for (size_t i = 1; i <= positions.size(); ++i) {
        const bool split_cluster =
            i == positions.size() ||
            (positions[i] - positions[i - 1]) > kBreakpointClusterGapBp;
        if (!split_cluster) {
            continue;
        }

        BreakpointPositionCluster cluster;
        cluster.start = positions[cluster_begin];
        cluster.end = positions[i - 1];
        cluster.support = static_cast<int32_t>(i - cluster_begin);
        cluster.center = median_i32(std::vector<int32_t>(
            positions.begin() + static_cast<std::ptrdiff_t>(cluster_begin),
            positions.begin() + static_cast<std::ptrdiff_t>(i)));
        clusters.push_back(cluster);
        cluster_begin = i;
    }
    return clusters;
}

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

struct LocalBreakpointHypothesis {
    bool valid = false;
    int32_t left = -1;
    int32_t right = -1;
    int32_t center = -1;
    int32_t support = 0;
    int32_t priority = std::numeric_limits<int32_t>::max();
};

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

template <typename Fn>
void for_each_valid_kmer(const std::string& seq, int32_t k, Fn&& fn) {
    if (k <= 0 || static_cast<int32_t>(seq.size()) < k) {
        return;
    }

    const uint64_t mask = (k >= 32)
        ? std::numeric_limits<uint64_t>::max()
        : ((uint64_t{1} << (2 * k)) - 1);
    uint64_t key = 0;
    int32_t valid_bases = 0;
    for (int32_t i = 0; i < static_cast<int32_t>(seq.size()); ++i) {
        const uint8_t code = char_to_2bit(seq[static_cast<size_t>(i)]);
        if (code > 3) {
            key = 0;
            valid_bases = 0;
            continue;
        }
        key = ((key << 2) | code) & mask;
        if (valid_bases < k) {
            ++valid_bases;
        }
        if (valid_bases >= k) {
            fn(i - k + 1, key);
        }
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

std::vector<LocalBreakpointHypothesis> collect_breakpoint_hypotheses(
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

    const auto sort_clusters = [&](std::vector<BreakpointPositionCluster>& clusters) {
        std::sort(clusters.begin(), clusters.end(), [&](const auto& lhs, const auto& rhs) {
            if (lhs.support != rhs.support) {
                return lhs.support > rhs.support;
            }
            const int32_t lhs_dist = std::abs(lhs.center - component.anchor_pos);
            const int32_t rhs_dist = std::abs(rhs.center - component.anchor_pos);
            if (lhs_dist != rhs_dist) {
                return lhs_dist < rhs_dist;
            }
            return lhs.center < rhs.center;
        });
    };

    const auto make_single_breakpoint_hypotheses = [&](std::vector<int32_t> positions,
                                                       int32_t priority) {
        std::vector<LocalBreakpointHypothesis> out;
        auto clusters = collect_breakpoint_position_clusters(std::move(positions));
        sort_clusters(clusters);
        if (clusters.size() > kBreakpointHypothesisMaxSingleClusters) {
            clusters.resize(kBreakpointHypothesisMaxSingleClusters);
        }

        out.reserve(clusters.size());
        for (const auto& cluster : clusters) {
            LocalBreakpointHypothesis hypothesis;
            hypothesis.valid = cluster.center >= 0;
            hypothesis.left = cluster.center;
            hypothesis.right = cluster.center;
            hypothesis.center = cluster.center;
            hypothesis.support = cluster.support;
            hypothesis.priority = priority;
            out.push_back(hypothesis);
        }
        return out;
    };

    const auto make_paired_breakpoint_hypotheses = [&](std::vector<int32_t> left_positions,
                                                       std::vector<int32_t> right_positions,
                                                       int32_t priority) {
        std::vector<LocalBreakpointHypothesis> out;
        auto left_clusters = collect_breakpoint_position_clusters(std::move(left_positions));
        auto right_clusters = collect_breakpoint_position_clusters(std::move(right_positions));
        sort_clusters(left_clusters);
        sort_clusters(right_clusters);

        if (left_clusters.size() > kBreakpointHypothesisMaxClustersPerSide) {
            left_clusters.resize(kBreakpointHypothesisMaxClustersPerSide);
        }
        if (right_clusters.size() > kBreakpointHypothesisMaxClustersPerSide) {
            right_clusters.resize(kBreakpointHypothesisMaxClustersPerSide);
        }

        if (left_clusters.empty() && right_clusters.empty()) {
            return out;
        }

        if (left_clusters.empty()) {
            out.reserve(right_clusters.size());
            for (const auto& right_cluster : right_clusters) {
                LocalBreakpointHypothesis hypothesis;
                hypothesis.valid = right_cluster.center >= 0;
                hypothesis.left = right_cluster.center;
                hypothesis.right = right_cluster.center;
                hypothesis.center = right_cluster.center;
                hypothesis.support = right_cluster.support;
                hypothesis.priority = priority;
                out.push_back(hypothesis);
            }
            return out;
        }

        if (right_clusters.empty()) {
            out.reserve(left_clusters.size());
            for (const auto& left_cluster : left_clusters) {
                LocalBreakpointHypothesis hypothesis;
                hypothesis.valid = left_cluster.center >= 0;
                hypothesis.left = left_cluster.center;
                hypothesis.right = left_cluster.center;
                hypothesis.center = left_cluster.center;
                hypothesis.support = left_cluster.support;
                hypothesis.priority = priority;
                out.push_back(hypothesis);
            }
            return out;
        }

        out.reserve(left_clusters.size() * right_clusters.size());
        for (const auto& left_cluster : left_clusters) {
            for (const auto& right_cluster : right_clusters) {
                const int32_t left = std::min(left_cluster.center, right_cluster.center);
                const int32_t right = std::max(left_cluster.center, right_cluster.center);
                if ((right - left) > kBreakpointPairCompatibilityBp) {
                    continue;
                }

                LocalBreakpointHypothesis hypothesis;
                hypothesis.valid = true;
                hypothesis.left = left;
                hypothesis.right = right;
                hypothesis.center =
                    static_cast<int32_t>((static_cast<int64_t>(left) + right) / 2);
                hypothesis.support = left_cluster.support + right_cluster.support;
                hypothesis.priority = priority;
                out.push_back(hypothesis);
            }
        }
        return out;
    };

    std::vector<LocalBreakpointHypothesis> hypotheses;
    hypotheses.reserve(24);
    const auto append_unique = [&](const std::vector<LocalBreakpointHypothesis>& batch) {
        for (const auto& hypothesis : batch) {
            if (!hypothesis.valid) {
                continue;
            }
            bool duplicate = false;
            for (const auto& existing : hypotheses) {
                if (existing.left == hypothesis.left &&
                    existing.right == hypothesis.right) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) {
                hypotheses.push_back(hypothesis);
            }
        }
    };

    append_unique(make_paired_breakpoint_hypotheses(
        fragment_split_left_positions,
        fragment_split_right_positions,
        0));
    append_unique(make_single_breakpoint_hypotheses(
        fragment_indel_positions,
        1));
    append_unique(make_paired_breakpoint_hypotheses(
        raw_split_left_positions,
        raw_split_right_positions,
        2));
    append_unique(make_paired_breakpoint_hypotheses(
        fragment_clip_left_positions,
        fragment_clip_right_positions,
        3));
    append_unique(make_single_breakpoint_hypotheses(
        raw_indel_positions,
        4));
    append_unique(make_paired_breakpoint_hypotheses(
        raw_clip_left_positions,
        raw_clip_right_positions,
        5));

    if (hypotheses.empty()) {
        LocalBreakpointHypothesis fallback;
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
        [&](const LocalBreakpointHypothesis& lhs, const LocalBreakpointHypothesis& rhs) {
            const int32_t lhs_score =
                lhs.support * breakpoint_hypothesis_support_weight(lhs.priority);
            const int32_t rhs_score =
                rhs.support * breakpoint_hypothesis_support_weight(rhs.priority);
            if (lhs_score != rhs_score) {
                return lhs_score > rhs_score;
            }
            const int32_t lhs_dist = std::abs(lhs.center - component.anchor_pos);
            const int32_t rhs_dist = std::abs(rhs.center - component.anchor_pos);
            if (lhs_dist != rhs_dist) {
                return lhs_dist < rhs_dist;
            }
            if (lhs.support != rhs.support) {
                return lhs.support > rhs.support;
            }
            if (lhs.priority != rhs.priority) {
                return lhs.priority < rhs.priority;
            }
            if (lhs.left != rhs.left) {
                return lhs.left < rhs.left;
            }
            return lhs.right < rhs.right;
        });
    if (top_k > 0 && hypotheses.size() > top_k) {
        hypotheses.resize(top_k);
    }
    return hypotheses;
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

std::vector<size_t> select_component_final_call_indices(
    const std::vector<ComponentFinalCallCandidate>& candidates) {
    std::vector<size_t> order;
    order.reserve(candidates.size());
    for (size_t i = 0; i < candidates.size(); ++i) {
        if (candidates[i].emit_te) {
            order.push_back(i);
        }
    }
    if (order.empty()) {
        return {};
    }

    std::sort(order.begin(), order.end(), [&](size_t lhs, size_t rhs) {
        if (candidates[lhs].pos != candidates[rhs].pos) {
            return candidates[lhs].pos < candidates[rhs].pos;
        }
        return lhs < rhs;
    });

    std::vector<size_t> selected;
    selected.reserve(order.size());
    size_t cluster_start = 0;
    while (cluster_start < order.size()) {
        const int32_t anchor_pos = candidates[order[cluster_start]].pos;
        size_t best = order[cluster_start];
        size_t cluster_end = cluster_start;
        while (cluster_end < order.size() &&
               std::abs(candidates[order[cluster_end]].pos - anchor_pos) <=
                   kFinalCallDedupDistanceBp) {
            const size_t idx = order[cluster_end];
            if (candidates[idx].score > candidates[best].score ||
                (candidates[idx].score == candidates[best].score && idx < best)) {
                best = idx;
            }
            ++cluster_end;
        }
        selected.push_back(best);
        cluster_start = cluster_end;
    }

    std::sort(selected.begin(), selected.end(), [&](size_t lhs, size_t rhs) {
        if (candidates[lhs].pos != candidates[rhs].pos) {
            return candidates[lhs].pos < candidates[rhs].pos;
        }
        return lhs < rhs;
    });
    return selected;
}

std::vector<Pipeline::BreakpointHypothesis> Pipeline::enumerate_breakpoint_hypotheses(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right,
    size_t top_k) const {
    const auto local_hypotheses = collect_breakpoint_hypotheses(
        component,
        local_records,
        fragments,
        seed_left,
        seed_right,
        0);
    std::vector<BreakpointHypothesis> hypotheses;
    hypotheses.reserve(local_hypotheses.size());
    for (const auto& local : local_hypotheses) {
        BreakpointHypothesis hypothesis;
        hypothesis.valid = local.valid;
        hypothesis.left = local.left;
        hypothesis.right = local.right;
        hypothesis.center = local.center;
        hypothesis.support = local.support;
        hypothesis.priority = local.priority;
        hypotheses.push_back(hypothesis);
    }
    return select_diverse_breakpoint_hypotheses(hypotheses, top_k, component.anchor_pos);
}

std::vector<Pipeline::BreakpointHypothesis> Pipeline::select_diverse_breakpoint_hypotheses(
    const std::vector<BreakpointHypothesis>& hypotheses,
    size_t top_k,
    int32_t anchor_pos) const {
    constexpr int32_t kBreakpointDiversitySlackBp = 30;
    constexpr int32_t kAnchorProximalSlackBp = 100;
    constexpr int32_t kAnchorProximalMinScore = 8;
    if (top_k == 0 || hypotheses.size() <= top_k) {
        return hypotheses;
    }

    const auto same_locus = [](const BreakpointHypothesis& lhs, const BreakpointHypothesis& rhs) {
        return std::abs(lhs.left - rhs.left) <= kBreakpointDiversitySlackBp &&
            std::abs(lhs.right - rhs.right) <= kBreakpointDiversitySlackBp;
    };
    const auto hypothesis_score = [](const BreakpointHypothesis& hypothesis) {
        return hypothesis.support * breakpoint_hypothesis_support_weight(hypothesis.priority);
    };
    const auto selected_contains = [&](const std::vector<BreakpointHypothesis>& selected,
                                       const BreakpointHypothesis& hypothesis) {
        return std::any_of(
            selected.begin(),
            selected.end(),
            [&](const auto& incumbent) {
                return incumbent.left == hypothesis.left && incumbent.right == hypothesis.right;
            });
    };

    std::vector<BreakpointHypothesis> selected;
    selected.reserve(std::min(top_k, hypotheses.size()));
    for (const auto& hypothesis : hypotheses) {
        bool redundant = false;
        for (const auto& incumbent : selected) {
            if (same_locus(hypothesis, incumbent)) {
                redundant = true;
                break;
            }
        }
        if (redundant) {
            continue;
        }
        selected.push_back(hypothesis);
        if (selected.size() >= top_k) {
            break;
        }
    }

    const bool has_anchor_proximal = std::any_of(
        selected.begin(),
        selected.end(),
        [&](const auto& hypothesis) {
            return std::abs(hypothesis.center - anchor_pos) <= kAnchorProximalSlackBp;
        });
    if (has_anchor_proximal) {
        return selected;
    }

    const BreakpointHypothesis* best_anchor_candidate = nullptr;
    for (const auto& hypothesis : hypotheses) {
        if (selected_contains(selected, hypothesis)) {
            continue;
        }
        if (hypothesis_score(hypothesis) < kAnchorProximalMinScore) {
            continue;
        }
        if (best_anchor_candidate == nullptr) {
            best_anchor_candidate = &hypothesis;
            continue;
        }
        const int32_t dist = std::abs(hypothesis.center - anchor_pos);
        const int32_t best_dist = std::abs(best_anchor_candidate->center - anchor_pos);
        if (dist != best_dist) {
            if (dist < best_dist) {
                best_anchor_candidate = &hypothesis;
            }
            continue;
        }
        const int32_t score = hypothesis_score(hypothesis);
        const int32_t best_score = hypothesis_score(*best_anchor_candidate);
        if (score != best_score) {
            if (score > best_score) {
                best_anchor_candidate = &hypothesis;
            }
            continue;
        }
        if (hypothesis.support != best_anchor_candidate->support) {
            if (hypothesis.support > best_anchor_candidate->support) {
                best_anchor_candidate = &hypothesis;
            }
            continue;
        }
        if (hypothesis.priority != best_anchor_candidate->priority) {
            if (hypothesis.priority < best_anchor_candidate->priority) {
                best_anchor_candidate = &hypothesis;
            }
            continue;
        }
        if (hypothesis.left < best_anchor_candidate->left ||
            (hypothesis.left == best_anchor_candidate->left &&
             hypothesis.right < best_anchor_candidate->right)) {
            best_anchor_candidate = &hypothesis;
        }
    }
    if (best_anchor_candidate == nullptr) {
        return selected;
    }

    if (selected.size() < top_k) {
        selected.push_back(*best_anchor_candidate);
    } else if (!selected.empty()) {
        selected.back() = *best_anchor_candidate;
    }
    return selected;
}

std::pair<int32_t, int32_t> Pipeline::resolve_event_breakpoint_bounds(
    const ComponentCall& component,
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    int32_t seed_left,
    int32_t seed_right) const {
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

Pipeline::ComponentSignalCache Pipeline::build_component_signal_cache(
    const std::vector<const bam1_t*>& local_records,
    const std::vector<ReadReferenceSpan>& read_spans) const {
    ComponentSignalCache cache;
    cache.reads.reserve(local_records.size());
    for (size_t read_idx = 0; read_idx < local_records.size(); ++read_idx) {
        const bam1_t* record = local_records[read_idx];
        if (!record) {
            continue;
        }

        ReadView read(record);
        const uint32_t* cigar = read.cigar();
        const int32_t n_cigar = read.n_cigar();

        CachedLocalSignal cached;
        cached.qname = record_qname(record);
        cached.tid = read.tid();
        cached.mapq = read.mapq();
        cached.is_primary = is_primary_alignment(record);
        cached.has_sa_or_supp =
            read.has_sa_tag() || ((read.flag() & BAM_FSUPPLEMENTARY) != 0);
        if (read_idx < read_spans.size() && read_spans[read_idx].valid) {
            cached.span_valid = true;
            cached.span_start = read_spans[read_idx].start;
            cached.span_end = read_spans[read_idx].end;
        }

        if (cigar && n_cigar > 0) {
            const int first = find_first_non_hard_clip(cigar, n_cigar);
            const int last = find_last_non_hard_clip(cigar, n_cigar);
            int32_t ref_pos = read.pos();
            for (int32_t ci = 0; ci < n_cigar; ++ci) {
                const int op = bam_cigar_op(cigar[ci]);
                const int32_t len = static_cast<int32_t>(bam_cigar_oplen(cigar[ci]));

                if (ci == first && op == BAM_CSOFT_CLIP && len >= kSoftClipSignalMin) {
                    cached.left_clip_pos = ref_pos;
                }
                if (ci == last && op == BAM_CSOFT_CLIP && len >= kSoftClipSignalMin) {
                    cached.right_clip_pos = ref_pos;
                }
                if (op == BAM_CINS && len >= kLongInsertionSignalMin) {
                    cached.indel_positions.push_back(ref_pos);
                }
                if ((bam_cigar_type(op) & 2) != 0) {
                    ref_pos += len;
                }
            }
        }

        cache.reads.push_back(std::move(cached));
    }
    return cache;
}

EventReadEvidence Pipeline::collect_event_read_evidence_for_bounds_cached(
    const ComponentCall& component,
    const ComponentSignalCache& cache,
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

    const auto classify_cached_signal = [&](
        const CachedLocalSignal& cached,
        int32_t window_start,
        int32_t window_end) {
        LocalEventSignal signal;
        if (window_end < window_start) {
            return signal;
        }

        if (cached.left_clip_pos >= window_start && cached.left_clip_pos <= window_end) {
            signal.left_clip = true;
            signal.left_clip_pos = cached.left_clip_pos;
        }
        if (cached.right_clip_pos >= window_start && cached.right_clip_pos <= window_end) {
            signal.right_clip = true;
            signal.right_clip_pos = cached.right_clip_pos;
        }

        const int32_t window_center = window_start + ((window_end - window_start) / 2);
        int32_t best_indel_dist = std::numeric_limits<int32_t>::max();
        for (int32_t pos : cached.indel_positions) {
            if (pos < window_start || pos > window_end) {
                continue;
            }
            signal.indel = true;
            const int32_t dist = std::abs(pos - window_center);
            if (dist < best_indel_dist) {
                best_indel_dist = dist;
                signal.indel_pos = pos;
            }
        }

        signal.split = cached.has_sa_or_supp && (signal.left_clip || signal.right_clip);
        if (signal.split) {
            signal.split_left_pos = signal.left_clip_pos;
            signal.split_right_pos = signal.right_clip_pos;
        }
        return signal;
    };

    std::unordered_set<std::string> split_qnames;
    std::unordered_set<std::string> indel_qnames;
    std::unordered_set<std::string> left_clip_qnames;
    std::unordered_set<std::string> right_clip_qnames;
    std::unordered_set<std::string> nearby_candidate_left_clip_qnames;
    std::unordered_set<std::string> nearby_candidate_right_clip_qnames;

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

    for (const auto& cached : cache.reads) {
        if (cached.tid != component.tid || cached.qname.empty()) {
            continue;
        }

        const LocalEventSignal signal = classify_cached_signal(
            cached,
            alt_signal_start,
            alt_signal_end);
        if (signal.split) {
            split_qnames.insert(cached.qname);
        }
        if (signal.indel) {
            indel_qnames.insert(cached.qname);
        }
        if (signal.left_clip) {
            left_clip_qnames.insert(cached.qname);
        }
        if (signal.right_clip) {
            right_clip_qnames.insert(cached.qname);
        }

        if (have_nearby_component_clip_candidates) {
            const LocalEventSignal nearby_signal = classify_cached_signal(
                cached,
                candidate_clip_start,
                candidate_clip_end);
            if (nearby_signal.left_clip &&
                matches_nearby_component_clip_candidate(nearby_signal.left_clip_pos)) {
                nearby_candidate_left_clip_qnames.insert(cached.qname);
            }
            if (nearby_signal.right_clip &&
                matches_nearby_component_clip_candidate(nearby_signal.right_clip_pos)) {
                nearby_candidate_right_clip_qnames.insert(cached.qname);
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
    for (const auto& cached : cache.reads) {
        if (!cached.span_valid || cached.tid != component.tid) {
            continue;
        }
        if (cached.span_start > ref_span_start || cached.span_end < ref_span_end) {
            continue;
        }
        if (!cached.is_primary) {
            continue;
        }
        if (cached.qname.empty() || alt_qnames.find(cached.qname) != alt_qnames.end()) {
            continue;
        }
        if (classify_cached_signal(cached, ref_signal_start, ref_signal_end).any()) {
            continue;
        }
        if (cached.mapq >= 20) {
            ref_qnames.insert(cached.qname);
        } else {
            low_mapq_ref_qnames.insert(cached.qname);
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

    const ConsensusInputSummary inputs = collect_event_consensus_inputs(
        local_records,
        fragments,
        event_evidence);
    const bool use_full_context = !inputs.full_event_by_qname.empty();
    consensus.full_context_input_reads = inputs.full_context_input_reads;
    consensus.partial_context_input_reads = inputs.partial_context_input_reads;
    consensus.left_anchor_input_reads = inputs.left_anchor_input_reads;
    consensus.right_anchor_input_reads = inputs.right_anchor_input_reads;
    consensus.used_full_context = use_full_context;
    const auto& event_by_qname = use_full_context
        ? inputs.full_event_by_qname
        : inputs.partial_event_by_qname;
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

Pipeline::ConsensusInputSummary Pipeline::collect_event_consensus_inputs(
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    const EventReadEvidence& event_evidence) const {
    ConsensusInputSummary summary;

    std::unordered_set<std::string> support_qnames(
        event_evidence.support_qnames.begin(),
        event_evidence.support_qnames.end());
    const std::unordered_map<std::string, const bam1_t*> records_by_qname =
        best_records_by_qname(local_records);

    std::unordered_set<std::string> full_context_qnames;
    std::unordered_set<std::string> partial_context_qnames;
    std::unordered_set<std::string> left_anchor_qnames;
    std::unordered_set<std::string> right_anchor_qnames;
    summary.full_event_by_qname.reserve(fragments.size());
    summary.partial_event_by_qname.reserve(fragments.size());
    full_context_qnames.reserve(fragments.size());
    partial_context_qnames.reserve(fragments.size());
    left_anchor_qnames.reserve(fragments.size());
    right_anchor_qnames.reserve(fragments.size());

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
        if (!can_build_event_string_from_fragment(record_it->second, fragment)) {
            continue;
        }
        const std::string event_string = build_event_string_from_fragment(
            record_it->second,
            fragment);
        if (event_string.empty()) {
            continue;
        }

        switch (fragment.source) {
            case InsertionFragmentSource::kCigarInsertion:
            case InsertionFragmentSource::kSplitSa:
                full_context_qnames.insert(fragment.read_id);
                left_anchor_qnames.insert(fragment.read_id);
                right_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefLeft:
                partial_context_qnames.insert(fragment.read_id);
                right_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefRight:
                partial_context_qnames.insert(fragment.read_id);
                left_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kUnknown:
            default:
                break;
        }

        auto& event_by_qname = fragment_has_full_event_context(fragment)
            ? summary.full_event_by_qname
            : summary.partial_event_by_qname;
        auto it = event_by_qname.find(fragment.read_id);
        if (it == event_by_qname.end() || event_string.size() > it->second.size()) {
            event_by_qname[fragment.read_id] = event_string;
        }
    }

    summary.full_context_input_reads = static_cast<int32_t>(full_context_qnames.size());
    summary.partial_context_input_reads = static_cast<int32_t>(partial_context_qnames.size());
    summary.left_anchor_input_reads = static_cast<int32_t>(left_anchor_qnames.size());
    summary.right_anchor_input_reads = static_cast<int32_t>(right_anchor_qnames.size());
    return summary;
}

Pipeline::ConsensusInputCounts Pipeline::collect_event_consensus_input_counts(
    const std::vector<const bam1_t*>& local_records,
    const std::vector<InsertionFragment>& fragments,
    const EventReadEvidence& event_evidence) const {
    ConsensusInputCounts counts;

    std::unordered_set<std::string> support_qnames(
        event_evidence.support_qnames.begin(),
        event_evidence.support_qnames.end());
    const std::unordered_map<std::string, const bam1_t*> records_by_qname =
        best_records_by_qname(local_records);

    std::unordered_set<std::string> full_context_qnames;
    std::unordered_set<std::string> partial_context_qnames;
    std::unordered_set<std::string> left_anchor_qnames;
    std::unordered_set<std::string> right_anchor_qnames;
    full_context_qnames.reserve(fragments.size());
    partial_context_qnames.reserve(fragments.size());
    left_anchor_qnames.reserve(fragments.size());
    right_anchor_qnames.reserve(fragments.size());

    for (const auto& fragment : fragments) {
        if (fragment.read_id.empty() ||
            support_qnames.find(fragment.read_id) == support_qnames.end() ||
            !fragment_supports_event_consensus(fragment) ||
            !fragment_is_local_to_event(fragment, event_evidence.bp_left, event_evidence.bp_right)) {
            continue;
        }
        const auto record_it = records_by_qname.find(fragment.read_id);
        if (record_it == records_by_qname.end() ||
            !can_build_event_string_from_fragment(record_it->second, fragment)) {
            continue;
        }

        switch (fragment.source) {
            case InsertionFragmentSource::kCigarInsertion:
            case InsertionFragmentSource::kSplitSa:
                full_context_qnames.insert(fragment.read_id);
                left_anchor_qnames.insert(fragment.read_id);
                right_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefLeft:
                partial_context_qnames.insert(fragment.read_id);
                right_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kClipRefRight:
                partial_context_qnames.insert(fragment.read_id);
                left_anchor_qnames.insert(fragment.read_id);
                break;
            case InsertionFragmentSource::kUnknown:
            default:
                break;
        }
    }

    counts.full_context_input_reads = static_cast<int32_t>(full_context_qnames.size());
    counts.partial_context_input_reads = static_cast<int32_t>(partial_context_qnames.size());
    counts.left_anchor_input_reads = static_cast<int32_t>(left_anchor_qnames.size());
    counts.right_anchor_input_reads = static_cast<int32_t>(right_anchor_qnames.size());
    counts.input_event_reads =
        counts.full_context_input_reads > 0
            ? counts.full_context_input_reads
            : counts.partial_context_input_reads;
    return counts;
}

Pipeline::HypothesisValidatorEvidence Pipeline::collect_hypothesis_validator_evidence(
    const HypothesisSummary& summary,
    const ConsensusInputCounts& inputs,
    int32_t anchor_pos) const {
    HypothesisValidatorEvidence out;
    out.summary = summary;
    out.precise_support = summary.alt_split_reads + summary.alt_indel_reads;
    out.breakpoint_width = std::max(0, summary.bp_right - summary.bp_left);
    out.anchor_distance = std::min(
        std::abs(summary.bp_left - anchor_pos),
        std::abs(summary.bp_right - anchor_pos));
    out.full_context_input_reads = inputs.full_context_input_reads;
    out.partial_context_input_reads = inputs.partial_context_input_reads;
    out.left_anchor_input_reads = inputs.left_anchor_input_reads;
    out.right_anchor_input_reads = inputs.right_anchor_input_reads;
    out.input_event_reads = inputs.input_event_reads;

    const bool has_bilateral_anchor =
        out.left_anchor_input_reads > 0 && out.right_anchor_input_reads > 0;
    const bool has_precise_or_full_context =
        out.full_context_input_reads > 0 || out.precise_support > 0;
    const bool has_strong_bilateral_partial_context =
        has_strong_bilateral_partial_context_support(
            out.left_anchor_input_reads,
            out.right_anchor_input_reads,
            out.partial_context_input_reads,
            out.input_event_reads);
    const bool has_validator_support =
        has_precise_or_full_context || has_strong_bilateral_partial_context;

    if (!has_bilateral_anchor && !has_precise_or_full_context) {
        out.qc_reason = "VALIDATOR_NO_BILATERAL_ANCHOR";
        return out;
    }
    if (!has_validator_support) {
        out.qc_reason = "VALIDATOR_NO_PRECISE_OR_FULL_CONTEXT";
        return out;
    }

    out.feasible_for_expensive_stage = true;
    out.qc_reason = "VALIDATOR_PASS";
    return out;
}

bool Pipeline::has_strong_bilateral_partial_context_support(
    int32_t left_anchor_input_reads,
    int32_t right_anchor_input_reads,
    int32_t partial_context_input_reads,
    int32_t input_event_reads) const {
    const bool has_bilateral_anchor =
        left_anchor_input_reads > 0 && right_anchor_input_reads > 0;
    const int32_t min_partial_only_event_reads = std::max(
        4,
        config_.event_consensus_poa_min_reads * 2);
    return has_bilateral_anchor &&
        partial_context_input_reads >= min_partial_only_event_reads &&
        input_event_reads >= min_partial_only_event_reads;
}

bool Pipeline::compare_hypothesis_validator_priority(
    const HypothesisValidatorEvidence& lhs,
    const HypothesisValidatorEvidence& rhs) const {
    if (lhs.precise_support != rhs.precise_support) {
        return lhs.precise_support > rhs.precise_support;
    }
    if (lhs.full_context_input_reads != rhs.full_context_input_reads) {
        return lhs.full_context_input_reads > rhs.full_context_input_reads;
    }
    const int32_t lhs_anchor_balance =
        std::min(lhs.left_anchor_input_reads, lhs.right_anchor_input_reads);
    const int32_t rhs_anchor_balance =
        std::min(rhs.left_anchor_input_reads, rhs.right_anchor_input_reads);
    if (lhs_anchor_balance != rhs_anchor_balance) {
        return lhs_anchor_balance > rhs_anchor_balance;
    }
    if (lhs.breakpoint_width != rhs.breakpoint_width) {
        return lhs.breakpoint_width < rhs.breakpoint_width;
    }
    if (lhs.anchor_distance != rhs.anchor_distance) {
        return lhs.anchor_distance < rhs.anchor_distance;
    }
    if (lhs.summary.alt_struct_reads != rhs.summary.alt_struct_reads) {
        return lhs.summary.alt_struct_reads > rhs.summary.alt_struct_reads;
    }
    return lhs.summary.original_index < rhs.summary.original_index;
}

std::vector<Pipeline::ShortlistedHypothesis> Pipeline::build_expensive_stage_shortlist(
    const std::vector<HypothesisValidatorEvidence>& candidates) const {
    constexpr size_t kMaxExpensiveStageShortlist = 3;
    std::vector<HypothesisValidatorEvidence> feasible;
    feasible.reserve(candidates.size());
    for (const auto& candidate : candidates) {
        if (candidate.feasible_for_expensive_stage) {
            feasible.push_back(candidate);
        }
    }
    if (feasible.empty()) {
        return {};
    }

    std::sort(feasible.begin(), feasible.end(), [&](const auto& lhs, const auto& rhs) {
        return compare_hypothesis_validator_priority(lhs, rhs);
    });

    std::vector<ShortlistedHypothesis> out;
    out.push_back(ShortlistedHypothesis{feasible.front(), true});

    const auto support_jaccard = [](
                                     const std::vector<std::string>& lhs,
                                     const std::vector<std::string>& rhs) {
        size_t i = 0;
        size_t j = 0;
        size_t intersect = 0;
        while (i < lhs.size() && j < rhs.size()) {
            if (lhs[i] == rhs[j]) {
                ++intersect;
                ++i;
                ++j;
            } else if (lhs[i] < rhs[j]) {
                ++i;
            } else {
                ++j;
            }
        }
        const size_t union_size = lhs.size() + rhs.size() - intersect;
        return union_size == 0
            ? 0.0
            : static_cast<double>(intersect) / static_cast<double>(union_size);
    };

    for (size_t i = 1; i < feasible.size(); ++i) {
        const auto& candidate = feasible[i];
        bool redundant = false;
        for (const auto& selected : out) {
            const auto& incumbent = selected.validator;
            const bool spatially_distinct =
                std::abs(candidate.summary.bp_left - incumbent.summary.bp_left) > 60 ||
                std::abs(candidate.summary.bp_right - incumbent.summary.bp_right) > 60;
            const bool support_distinct =
                support_jaccard(
                    candidate.summary.support_qnames,
                    incumbent.summary.support_qnames) < 0.5;
            if (!spatially_distinct && !support_distinct) {
                redundant = true;
                break;
            }
        }
        if (redundant) {
            continue;
        }
        out.push_back(ShortlistedHypothesis{candidate, false});
        if (out.size() >= kMaxExpensiveStageShortlist) {
            break;
        }
    }
    return out;
}

std::string Pipeline::pre_segmentation_gate_reason(
    const EventReadEvidence& event_evidence,
    const EventConsensus& event_consensus) const {
    if (event_consensus.left_anchor_input_reads <= 0 ||
        event_consensus.right_anchor_input_reads <= 0) {
        return "PRESEG_NO_BILATERAL_ANCHOR";
    }
    const bool has_precise_or_full_context =
        event_consensus.full_context_input_reads > 0 ||
        (event_evidence.alt_split_reads + event_evidence.alt_indel_reads) > 0;
    const bool has_strong_bilateral_partial_context =
        has_strong_bilateral_partial_context_support(
            event_consensus.left_anchor_input_reads,
            event_consensus.right_anchor_input_reads,
            event_consensus.partial_context_input_reads,
            event_consensus.input_event_reads);
    if (!has_precise_or_full_context && !has_strong_bilateral_partial_context) {
        return "PRESEG_NO_PRECISE_OR_FULL_CONTEXT";
    }
    return {};
}

std::vector<Pipeline::AnchorSeedBin> Pipeline::collect_anchor_seed_bins(
    const std::string& query,
    int32_t breakpoint,
    int32_t ref_window_start,
    const std::string& ref_window,
    bool is_left) const {
    std::unordered_map<uint64_t, std::vector<int32_t>> ref_hits;
    ref_hits.reserve(ref_window.size());
    for_each_valid_kmer(ref_window, kEventSegmentationSeedK, [&](int32_t pos, uint64_t key) {
        ref_hits[key].push_back(pos);
    });

    std::unordered_map<int32_t, AnchorSeedBin> bins_by_start;
    bins_by_start.reserve(ref_window.size());
    for_each_valid_kmer(query, kEventSegmentationSeedK, [&](int32_t query_pos, uint64_t key) {
        const auto it = ref_hits.find(key);
        if (it == ref_hits.end()) {
            return;
        }
        for (int32_t ref_pos : it->second) {
            const int32_t ref_start = ref_window_start + ref_pos - query_pos;
            const int32_t ref_bin_start =
                (ref_start / kEventSegmentationSeedBinBp) * kEventSegmentationSeedBinBp;
            auto& bin = bins_by_start[ref_bin_start];
            bin.ref_bin_start = ref_bin_start;
            ++bin.support;

            const int32_t candidate_breakpoint = is_left
                ? (ref_start + static_cast<int32_t>(query.size()))
                : ref_start;
            bin.best_breakpoint_delta = std::min(
                bin.best_breakpoint_delta,
                std::abs(candidate_breakpoint - breakpoint));
        }
    });

    std::vector<AnchorSeedBin> bins;
    bins.reserve(bins_by_start.size());
    for (const auto& entry : bins_by_start) {
        bins.push_back(entry.second);
    }
    std::sort(bins.begin(), bins.end(), [](const AnchorSeedBin& lhs, const AnchorSeedBin& rhs) {
        if (lhs.support != rhs.support) {
            return lhs.support > rhs.support;
        }
        if (lhs.best_breakpoint_delta != rhs.best_breakpoint_delta) {
            return lhs.best_breakpoint_delta < rhs.best_breakpoint_delta;
        }
        return lhs.ref_bin_start < rhs.ref_bin_start;
    });
    if (static_cast<int32_t>(bins.size()) > kEventSegmentationSeedTopBins) {
        bins.resize(static_cast<size_t>(kEventSegmentationSeedTopBins));
    }
    return bins;
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
    const int32_t min_one_sided_len =
        kEventSegmentationMinFlankAlignBp + kEventSegmentationMinInsertBp;
    if (consensus_len < min_one_sided_len) {
        segmentation.qc_reason = "EVENT_CONSENSUS_TOO_SHORT";
        return segmentation;
    }

    const int32_t paired_max_flank_query_len = std::min(
        kEventSegmentationMaxFlankQueryBp,
        consensus_len - kEventSegmentationMinFlankAlignBp - kEventSegmentationMinInsertBp);
    const int32_t one_sided_max_flank_query_len = std::min(
        kEventSegmentationMaxFlankQueryBp,
        consensus_len - kEventSegmentationMinInsertBp);

    const int32_t bp_left = std::min(event_evidence.bp_left, event_evidence.bp_right);
    const int32_t bp_right = std::max(event_evidence.bp_left, event_evidence.bp_right);
    if (component.chrom.empty() || bp_left < 0 || bp_right < 0) {
        segmentation.qc_reason = "INVALID_EVENT_BREAKPOINTS";
        return segmentation;
    }

    const int32_t left_window_start = std::max(
        0,
        bp_left - kEventSegmentationBreakpointSlackBp - one_sided_max_flank_query_len);
    const int32_t left_window_end = std::max(
        left_window_start + 1,
        bp_left + kEventSegmentationBreakpointSlackBp);
    const int32_t right_window_start = std::max(
        0,
        bp_right - kEventSegmentationBreakpointSlackBp);
    const int32_t right_window_end = std::max(
        right_window_start + 1,
        bp_right + kEventSegmentationBreakpointSlackBp + one_sided_max_flank_query_len);

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
        int32_t flank_query_len,
        int32_t query_endpoint_slack) {
        EventFlankSearchResult result;
        if (ref_window.empty() || flank_query_len < kEventSegmentationMinFlankAlignBp) {
            return result;
        }

        const int32_t ref_window_end =
            ref_window_start + static_cast<int32_t>(ref_window.size());
        const int32_t search_lo = std::max(0, breakpoint - kEventSegmentationBreakpointSlackBp);
        const int32_t search_hi = breakpoint + kEventSegmentationBreakpointSlackBp;
        const int32_t query_seed_len = std::min(
            flank_query_len + query_endpoint_slack,
            consensus_len);
        const int32_t query_seed_start = is_left ? 0 : (consensus_len - query_seed_len);
        const std::string seed_query = is_left
            ? consensus.substr(0, static_cast<size_t>(query_seed_len))
            : consensus.substr(static_cast<size_t>(query_seed_start));
        const std::vector<AnchorSeedBin> seed_bins = collect_anchor_seed_bins(
            seed_query,
            breakpoint,
            ref_window_start,
            ref_window,
            is_left);
        if (seed_bins.empty()) {
            return result;
        }

        for (int32_t align_len = flank_query_len;
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
                seed_bins.size() *
                static_cast<size_t>(kEventSegmentationSeedBinBp) *
                static_cast<size_t>(std::max(1, query_hi - query_lo + 1)));

            for (int32_t query_start = query_lo; query_start <= query_hi; ++query_start) {
                const std::string_view query(
                    consensus.data() + query_start,
                    static_cast<size_t>(align_len));
                const int32_t endpoint_offset = is_left
                    ? query_start
                    : (consensus_len - (query_start + align_len));
                const int32_t seed_query_offset = query_start - query_seed_start;

                for (const auto& seed_bin : seed_bins) {
                    for (int32_t bin_offset = 0;
                         bin_offset < kEventSegmentationSeedBinBp;
                         ++bin_offset) {
                        const int32_t ref_start =
                            seed_bin.ref_bin_start + seed_query_offset + bin_offset;
                        const int32_t ref_end = ref_start + align_len;
                        const int32_t candidate_breakpoint = is_left ? ref_end : ref_start;
                        if (candidate_breakpoint < search_lo || candidate_breakpoint > search_hi) {
                            continue;
                        }
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
                        placement.breakpoint_delta = std::abs(candidate_breakpoint - breakpoint);
                        placement.endpoint_offset = endpoint_offset;
                        placements.push_back(placement);
                    }
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
        paired_max_flank_query_len,
        0);
    if (left_search.candidates.empty() && kEventSegmentationEndpointSlackBp > 0) {
        left_search = collect_candidates(
            true,
            bp_left,
            left_window_start,
            left_ref_window,
            paired_max_flank_query_len,
            kEventSegmentationEndpointSlackBp);
    }

    EventFlankSearchResult right_search = collect_candidates(
        false,
        bp_right,
        right_window_start,
        right_ref_window,
        paired_max_flank_query_len,
        0);
    if (right_search.candidates.empty() && kEventSegmentationEndpointSlackBp > 0) {
        right_search = collect_candidates(
            false,
            bp_right,
            right_window_start,
            right_ref_window,
            paired_max_flank_query_len,
            kEventSegmentationEndpointSlackBp);
    }

    const auto emit_left_sided_segmentation = [&](const EventFlankPlacement& best_left) {
        const int32_t insert_start = best_left.query_end;
        const int32_t insert_len = consensus_len - insert_start;
        if (insert_len < kEventSegmentationMinInsertBp) {
            segmentation.qc_reason = "EMPTY_EVENT_INSERT_SEGMENT";
            return;
        }

        segmentation.left_flank_seq = consensus.substr(
            static_cast<size_t>(best_left.query_start),
            static_cast<size_t>(best_left.align_len));
        segmentation.insert_seq = consensus.substr(
            static_cast<size_t>(insert_start),
            static_cast<size_t>(insert_len));
        segmentation.left_ref_start = best_left.ref_start;
        segmentation.left_ref_end = best_left.ref_end;
        segmentation.right_ref_start = bp_right;
        segmentation.right_ref_end = bp_right;
        segmentation.left_flank_align_len = best_left.align_len;
        segmentation.right_flank_align_len = 0;
        segmentation.left_flank_identity = best_left.identity;
        segmentation.right_flank_identity = 0.0;
        segmentation.pass = true;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";
    };

    const auto emit_right_sided_segmentation = [&](const EventFlankPlacement& best_right) {
        const int32_t insert_len = best_right.query_start;
        if (insert_len < kEventSegmentationMinInsertBp) {
            segmentation.qc_reason = "EMPTY_EVENT_INSERT_SEGMENT";
            return;
        }

        segmentation.insert_seq = consensus.substr(
            0,
            static_cast<size_t>(insert_len));
        segmentation.right_flank_seq = consensus.substr(
            static_cast<size_t>(best_right.query_start),
            static_cast<size_t>(best_right.align_len));
        segmentation.left_ref_start = bp_left;
        segmentation.left_ref_end = bp_left;
        segmentation.right_ref_start = best_right.ref_start;
        segmentation.right_ref_end = best_right.ref_end;
        segmentation.left_flank_align_len = 0;
        segmentation.right_flank_align_len = best_right.align_len;
        segmentation.left_flank_identity = 0.0;
        segmentation.right_flank_identity = best_right.identity;
        segmentation.pass = true;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION_ONE_SIDED_RIGHT";
    };

    const auto collect_one_sided_search = [&](bool is_left) {
        EventFlankSearchResult search = collect_candidates(
            is_left,
            is_left ? bp_left : bp_right,
            is_left ? left_window_start : right_window_start,
            is_left ? left_ref_window : right_ref_window,
            one_sided_max_flank_query_len,
            0);
        if (search.candidates.empty() && kEventSegmentationEndpointSlackBp > 0) {
            search = collect_candidates(
                is_left,
                is_left ? bp_left : bp_right,
                is_left ? left_window_start : right_window_start,
                is_left ? left_ref_window : right_ref_window,
                one_sided_max_flank_query_len,
                kEventSegmentationEndpointSlackBp);
        }
        return search;
    };

    const auto better_one_sided = [](
                                      const EventFlankPlacement& lhs,
                                      const EventFlankPlacement& rhs) {
        if (lhs.identity != rhs.identity) {
            return lhs.identity > rhs.identity;
        }
        if (lhs.breakpoint_delta != rhs.breakpoint_delta) {
            return lhs.breakpoint_delta < rhs.breakpoint_delta;
        }
        if (lhs.align_len != rhs.align_len) {
            return lhs.align_len > rhs.align_len;
        }
        if (lhs.endpoint_offset != rhs.endpoint_offset) {
            return lhs.endpoint_offset < rhs.endpoint_offset;
        }
        return lhs.ref_start < rhs.ref_start;
    };

    const auto emit_best_one_sided = [&](bool allow_left, bool allow_right) {
        EventFlankSearchResult one_sided_left;
        EventFlankSearchResult one_sided_right;
        if (allow_left) {
            one_sided_left = collect_one_sided_search(true);
        }
        if (allow_right) {
            one_sided_right = collect_one_sided_search(false);
        }

        const bool have_left = !one_sided_left.candidates.empty();
        const bool have_right = !one_sided_right.candidates.empty();
        if (!have_left && !have_right) {
            return false;
        }
        if (have_left && !have_right) {
            emit_left_sided_segmentation(one_sided_left.candidates.front());
            return segmentation.pass;
        }
        if (!have_left && have_right) {
            emit_right_sided_segmentation(one_sided_right.candidates.front());
            return segmentation.pass;
        }
        if (better_one_sided(
                one_sided_left.candidates.front(),
                one_sided_right.candidates.front())) {
            emit_left_sided_segmentation(one_sided_left.candidates.front());
        } else {
            emit_right_sided_segmentation(one_sided_right.candidates.front());
        }
        return segmentation.pass;
    };

    if (left_search.candidates.empty() && right_search.candidates.empty()) {
        if (emit_best_one_sided(true, true)) {
            return segmentation;
        }
        segmentation.qc_reason = "NO_LEFT_FLANK_MATCH";
        return segmentation;
    }
    if (left_search.candidates.empty()) {
        if (emit_best_one_sided(false, true)) {
            return segmentation;
        }
        segmentation.qc_reason = "NO_LEFT_FLANK_MATCH";
        return segmentation;
    }
    if (right_search.candidates.empty()) {
        if (emit_best_one_sided(true, false)) {
            return segmentation;
        }
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
        if (emit_best_one_sided(true, true)) {
            return segmentation;
        }
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
    (void)bin_start;
    (void)bin_end;
    DbscanComponentModule module;
    return module.build(bin_records, chrom, tid);
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
    const int32_t bin_size = std::max(1, config_.bin_size);
    int32_t worker_count = config_.parallel_workers;
    if (worker_count <= 0) {
        worker_count = static_cast<int32_t>(std::thread::hardware_concurrency());
    }
    worker_count = std::max(1, worker_count);

    struct ParallelScanState {
        int32_t current_tid = -1;
        int32_t current_bin_index = -1;
        std::vector<BufferedRecord> current_bin_records;
        std::deque<ExactBinSnapshot> recent_snapshots;
        size_t next_owner_offset = 0;
    } state;

    std::vector<ExactBinTask> tasks;
    std::vector<std::pair<int32_t, int32_t>> task_order;
    task_order.reserve(1024);

    const auto append_ready_parallel_tasks = [&](const std::vector<ExactBinTask>& ready_tasks) {
        for (const auto& task : ready_tasks) {
            task_order.push_back({task.tid, task.bin_index});
            tasks.push_back(task);
        }
    };

    const auto flush_parallel_bin = [&](
                                       bool flush_all_pending,
                                       int32_t latest_observed_bin_index) {
        if (state.current_tid < 0 || state.current_bin_index < 0) {
            return;
        }

        if (!state.current_bin_records.empty()) {
            SharedRecordBatch batch = std::make_shared<const std::vector<BufferedRecord>>(
                std::move(state.current_bin_records));
            state.current_bin_records.clear();
            state.recent_snapshots.push_back(ExactBinSnapshot{
                state.current_tid,
                state.current_bin_index,
                std::move(batch),
            });
        }

        if (state.recent_snapshots.empty()) {
            return;
        }

        append_ready_parallel_tasks(materialize_ready_exact_bin_tasks(
            state.recent_snapshots,
            state.next_owner_offset,
            latest_observed_bin_index >= 0
                ? latest_observed_bin_index
                : state.recent_snapshots.back().bin_index,
            flush_all_pending,
            bin_size,
            kCrossBinContextBins,
            kWindowBinSlackBp));

        if (flush_all_pending) {
            state.recent_snapshots.clear();
            state.next_owner_offset = 0;
        }
    };

    const auto progress_cb = [this](int64_t processed, int32_t tid) {
        std::cerr << "[Pipeline] processed=" << processed << " current_tid=" << tid << '\n';
        return true;
    };

    result.total_reads = bam_reader_->stream(
        [this, &result, &state, bin_size, &flush_parallel_bin](BamRecordPtr&& record) {
            if (!record) {
                return;
            }

            ReadView view(record.get());
            if (!gate1_module_.pass_preliminary(view)) {
                return;
            }

            result.gate1_passed += 1;
            const int32_t tid = view.tid();
            const int32_t bin_index = exact_bin_index_for_record(record.get(), bin_size);

            if (state.current_tid < 0) {
                state.current_tid = tid;
                state.current_bin_index = bin_index;
            }

            if (tid != state.current_tid || bin_index != state.current_bin_index) {
                flush_parallel_bin(
                    tid != state.current_tid,
                    tid == state.current_tid ? bin_index : state.current_bin_index);
                state.current_tid = tid;
                state.current_bin_index = bin_index;
            }

            BufferedRecord buffered;
            buffered.ref_end = compute_ref_end(view);
            buffered.record = std::move(record);
            state.current_bin_records.push_back(std::move(buffered));
        },
        progress_cb,
        config_.progress_interval);

    flush_parallel_bin(true, state.current_bin_index);

    if (tasks.empty()) {
        finalize_final_calls(result);
        return result;
    }
    if (tasks.size() == 1) {
        const ExactBinTask& task = tasks.front();
        std::vector<const bam1_t*> bin_records;
        bin_records.reserve(task.records.size());
        for (const auto& record_ref : task.records) {
            if (record_ref.record) {
                bin_records.push_back(record_ref.record);
            }
        }
        process_bin_records(
            std::move(bin_records),
            task.tid,
            task.bin_index,
            task.bin_start,
            task.bin_end,
            result);
        finalize_final_calls(result);
        return result;
    }

    worker_count = std::min<int32_t>(
        worker_count,
        static_cast<int32_t>(tasks.size()));
    std::vector<ExactBinTaskResult> task_results(tasks.size());
    std::atomic<size_t> next_task_index{0};
    std::exception_ptr worker_error;
    std::mutex worker_error_mutex;
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(worker_count));

    for (int32_t worker_id = 0; worker_id < worker_count; ++worker_id) {
        workers.emplace_back([&, worker_id]() {
            (void)worker_id;
            try {
                PipelineConfig worker_config = config_;
                worker_config.enable_parallel = false;
                worker_config.bam_threads = 1;
                auto worker_reader = bam_reader_->clone(1);
                Pipeline worker_pipeline(worker_config, std::move(worker_reader));

                while (true) {
                    const size_t task_index = next_task_index.fetch_add(1);
                    if (task_index >= tasks.size()) {
                        break;
                    }

                    const ExactBinTask& task = tasks[task_index];
                    std::vector<const bam1_t*> bin_records;
                    bin_records.reserve(task.records.size());
                    for (const auto& record_ref : task.records) {
                        if (record_ref.record) {
                            bin_records.push_back(record_ref.record);
                        }
                    }

                    PipelineResult partial;
                    const auto started = std::chrono::steady_clock::now();
                    worker_pipeline.process_bin_records(
                        std::move(bin_records),
                        task.tid,
                        task.bin_index,
                        task.bin_start,
                        task.bin_end,
                        partial);
                    ExactBinTaskResult task_result;
                    task_result.tid = task.tid;
                    task_result.bin_index = task.bin_index;
                    task_result.partial = std::move(partial);
                    task_result.elapsed_s = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - started).count();
                    task_results[task_index] = std::move(task_result);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(worker_error_mutex);
                if (!worker_error) {
                    worker_error = std::current_exception();
                }
            }
        });
    }

    for (auto& worker : workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
    if (worker_error) {
        std::rethrow_exception(worker_error);
    }

    for (size_t i = 0; i < task_results.size(); ++i) {
        const auto& expected = task_order[i];
        const auto& actual = task_results[i];
        if (actual.tid != expected.first || actual.bin_index != expected.second) {
            throw std::runtime_error("parallel exact-bin task result order mismatch");
        }
        append_partial_pipeline_result(result, std::move(task_results[i].partial));
    }

    if (config_.log_parallel_progress) {
        ParallelProgressSnapshot snapshot;
        snapshot.reads_scanned = result.total_reads;
        snapshot.raw_bins_discovered = static_cast<int64_t>(task_order.size());
        snapshot.tasks_materialized = static_cast<int64_t>(tasks.size());
        snapshot.tasks_completed = static_cast<int64_t>(task_results.size());
        snapshot.reducer_committed = static_cast<int64_t>(task_results.size());
        emit_pipeline_log_line(render_parallel_progress(snapshot));
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

    if (state.current_tid < 0) {
        state.current_tid = view.tid();
    }

    if (view.tid() != state.current_tid) {
        flush_current_bin(state, result);
        state.current_tid = view.tid();
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

    std::vector<const bam1_t*> contig_records;
    contig_records.reserve(state.current_bin_records.size());
    for (const auto& buffered : state.current_bin_records) {
        if (buffered.record) {
            contig_records.push_back(buffered.record.get());
        }
    }

    if (!contig_records.empty()) {
        process_bin_records(
            std::move(contig_records),
            state.current_tid,
            0,
            std::numeric_limits<int32_t>::min(),
            std::numeric_limits<int32_t>::max(),
            result);
    }

    state.current_bin_records.clear();
}

std::vector<int32_t> collect_alt_observed_lengths(
    const ComponentCall& component,
    const EventReadEvidence& event_evidence) {
    std::unordered_set<std::string> support_qnames(
        event_evidence.support_qnames.begin(),
        event_evidence.support_qnames.end());
    std::vector<int32_t> lengths;
    lengths.reserve(component.breakpoint_candidates.size());
    for (const auto& bp : component.breakpoint_candidates) {
        if (bp.read_id.empty() || support_qnames.find(bp.read_id) == support_qnames.end()) {
            continue;
        }
        if (bp.ins_len > 0) {
            lengths.push_back(bp.ins_len);
        } else if (bp.clip_len > 0) {
            lengths.push_back(bp.clip_len);
        }
    }
    return lengths;
}

int32_t infer_event_length_from_alt_support(const std::vector<int32_t>& lengths) {
    return median_i32(lengths);
}

double support_jaccard(
    const std::vector<std::string>& lhs,
    const std::vector<std::string>& rhs) {
    size_t i = 0;
    size_t j = 0;
    size_t intersect = 0;
    while (i < lhs.size() && j < rhs.size()) {
        if (lhs[i] == rhs[j]) {
            ++intersect;
            ++i;
            ++j;
        } else if (lhs[i] < rhs[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    const size_t union_size = lhs.size() + rhs.size() - intersect;
    return union_size == 0
        ? 0.0
        : static_cast<double>(intersect) / static_cast<double>(union_size);
}

Pipeline::HypothesisSummary Pipeline::build_hypothesis_summary(
    const ComponentCall& component,
    const ComponentSignalCache& signal_cache,
    const std::vector<InsertionFragment>& fragments,
    size_t original_index,
    int32_t bp_left,
    int32_t bp_right) const {
    const EventReadEvidence event_evidence = collect_event_read_evidence_for_bounds_cached(
        component,
        signal_cache,
        fragments,
        bp_left,
        bp_right,
        bp_left,
        bp_right);

    HypothesisSummary summary;
    summary.original_index = original_index;
    summary.bp_left = event_evidence.bp_left;
    summary.bp_right = event_evidence.bp_right;
    summary.alt_split_reads = event_evidence.alt_split_reads;
    summary.alt_indel_reads = event_evidence.alt_indel_reads;
    summary.alt_left_clip_reads = event_evidence.alt_left_clip_reads;
    summary.alt_right_clip_reads = event_evidence.alt_right_clip_reads;
    summary.alt_struct_reads = event_evidence.alt_struct_reads;
    summary.ref_span_reads = event_evidence.ref_span_reads;
    summary.inferred_event_length = infer_event_length_from_alt_support(
        collect_alt_observed_lengths(component, event_evidence));
    summary.support_qnames = event_evidence.support_qnames;
    return summary;
}

std::vector<Pipeline::HypothesisSummary> Pipeline::select_hypothesis_summaries_for_expensive_stage(
    const std::vector<HypothesisSummary>& summaries) const {
    const auto collapsed = collapse_hypothesis_summaries(summaries);
    std::vector<HypothesisSummary> selected;
    selected.reserve(collapsed.size());
    for (size_t i = 0; i < collapsed.size(); ++i) {
        if (should_keep_hypothesis_for_expensive_stage(collapsed[i], i == 0)) {
            selected.push_back(collapsed[i]);
        }
    }
    return selected;
}

std::vector<Pipeline::HypothesisSummary> Pipeline::collapse_hypothesis_summaries(
    const std::vector<HypothesisSummary>& summaries) const {
    constexpr int32_t kHypothesisCollapseSlackBp = 30;
    constexpr double kHypothesisCollapseMinJaccard = 0.8;

    const auto support_jaccard = [](
                                     const std::vector<std::string>& lhs,
                                     const std::vector<std::string>& rhs) {
        size_t i = 0;
        size_t j = 0;
        size_t intersect = 0;
        while (i < lhs.size() && j < rhs.size()) {
            if (lhs[i] == rhs[j]) {
                ++intersect;
                ++i;
                ++j;
            } else if (lhs[i] < rhs[j]) {
                ++i;
            } else {
                ++j;
            }
        }
        const size_t union_size = lhs.size() + rhs.size() - intersect;
        return union_size == 0
            ? 0.0
            : static_cast<double>(intersect) / static_cast<double>(union_size);
    };

    const auto better_summary_representative = [](
                                                   const HypothesisSummary& lhs,
                                                   const HypothesisSummary& rhs) {
        if (lhs.alt_struct_reads != rhs.alt_struct_reads) {
            return lhs.alt_struct_reads > rhs.alt_struct_reads;
        }
        const int32_t lhs_precise = lhs.alt_split_reads + lhs.alt_indel_reads;
        const int32_t rhs_precise = rhs.alt_split_reads + rhs.alt_indel_reads;
        if (lhs_precise != rhs_precise) {
            return lhs_precise > rhs_precise;
        }
        const int32_t lhs_width = std::abs(lhs.bp_right - lhs.bp_left);
        const int32_t rhs_width = std::abs(rhs.bp_right - rhs.bp_left);
        if (lhs_width != rhs_width) {
            return lhs_width < rhs_width;
        }
        return lhs.original_index < rhs.original_index;
    };

    std::vector<HypothesisSummary> kept;
    kept.reserve(summaries.size());
    for (const auto& summary : summaries) {
        bool merged = false;
        for (auto& incumbent : kept) {
            const bool close_left =
                std::abs(summary.bp_left - incumbent.bp_left) <= kHypothesisCollapseSlackBp;
            const bool close_right =
                std::abs(summary.bp_right - incumbent.bp_right) <= kHypothesisCollapseSlackBp;
            if (!close_left || !close_right) {
                continue;
            }
            if (support_jaccard(summary.support_qnames, incumbent.support_qnames) <
                kHypothesisCollapseMinJaccard) {
                continue;
            }
            if (better_summary_representative(summary, incumbent)) {
                incumbent = summary;
            }
            merged = true;
            break;
        }
        if (!merged) {
            kept.push_back(summary);
        }
    }
    return kept;
}

bool Pipeline::should_keep_hypothesis_for_expensive_stage(
    const HypothesisSummary& summary,
    bool is_top_ranked_survivor) const {
    if (is_top_ranked_survivor) {
        return true;
    }
    if (summary.alt_split_reads > 0 || summary.alt_indel_reads > 0) {
        return true;
    }
    const bool bilateral_clip =
        summary.alt_left_clip_reads > 0 && summary.alt_right_clip_reads > 0;
    return bilateral_clip && summary.alt_struct_reads >= 2;
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
    input.alt_observed_lengths = collect_alt_observed_lengths(component, event_evidence);
    input.event_length = infer_event_length_from_alt_support(input.alt_observed_lengths);
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

    const bool emit_unknown_te =
        acceptance.qc == "PASS_FINAL_TE_CALL_UNKNOWN" ||
        acceptance.qc.rfind("PASS_FINAL_TE_CALL_UNKNOWN|", 0) == 0;
    call.family = emit_unknown_te
        ? "UNKNOWN"
        : (te_alignment.best_family.empty() ? "NA" : te_alignment.best_family);
    call.subfamily = emit_unknown_te
        ? "UNKNOWN"
        : (te_alignment.best_subfamily.empty() ? "NA" : te_alignment.best_subfamily);
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
    int32_t owner_bin_start,
    int32_t owner_bin_end,
    PipelineResult& result) const {
    if (bin_records.empty()) {
        return;
    }

    result.processed_bins += 1;

    int32_t bin_start = std::numeric_limits<int32_t>::max();
    int32_t bin_end = std::numeric_limits<int32_t>::min();
    for (const bam1_t* record : bin_records) {
        if (!record) {
            continue;
        }
        ReadView read(record);
        bin_start = std::min(bin_start, read.pos());
        bin_end = std::max(bin_end, compute_ref_end(read));
    }
    if (bin_start == std::numeric_limits<int32_t>::max()) {
        bin_start = std::max(0, bin_index);
    }
    if (bin_end == std::numeric_limits<int32_t>::min()) {
        bin_end = bin_start;
    }
    const std::string chrom = bam_reader_->chromosome_name(tid);

    const auto component_build_started = std::chrono::steady_clock::now();
    auto components = component_module_.build(bin_records, chrom, tid, bin_start, bin_end);
    if (owner_bin_start < owner_bin_end) {
        components.erase(
            std::remove_if(
                components.begin(),
                components.end(),
                [owner_bin_start, owner_bin_end](const ComponentCall& component) {
                    return component.anchor_pos < owner_bin_start ||
                        component.anchor_pos >= owner_bin_end;
                }),
            components.end());
    }
    const double component_build_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - component_build_started).count();
    result.built_components += static_cast<int64_t>(components.size());
    std::vector<std::pair<int32_t, int32_t>> component_seed_bounds;
    component_seed_bounds.reserve(components.size());
    std::vector<LocalIntervalRequest> requests;
    requests.reserve(components.size());
    for (size_t ci = 0; ci < components.size(); ++ci) {
        const auto& component = components[ci];
        const auto seed_bounds = infer_component_breakpoint_bounds(component);
        component_seed_bounds.push_back(seed_bounds);
        requests.push_back(LocalIntervalRequest{
            component.chrom,
            std::max(0, std::min(seed_bounds.first, seed_bounds.second) - kLocalEventFetchSlackBp),
            std::max(
                std::max(0, std::min(seed_bounds.first, seed_bounds.second) - kLocalEventFetchSlackBp) + 1,
                std::max(seed_bounds.first, seed_bounds.second) + kLocalEventFetchSlackBp),
            ci,
        });
    }
    LocalIntervalReuseStats reuse_stats;
    reuse_stats.request_count = requests.size();
    const std::vector<CanonicalLocalInterval> intervals =
        build_canonical_local_intervals(requests, 128);
    reuse_stats.canonical_interval_count = intervals.size();
    std::vector<LocalIntervalCacheEntry> interval_cache;
    interval_cache.reserve(intervals.size());
    const auto local_fetch_started = std::chrono::steady_clock::now();
    for (const auto& interval : intervals) {
        LocalIntervalCacheEntry entry;
        entry.interval = interval;
        const bool fetch_ok = bam_reader_->fetch(
            interval.chrom,
            interval.start,
            interval.end,
            [&](BamRecordPtr&& record) {
                if (!record) {
                    return true;
                }
                entry.owned_records.push_back(std::move(record));
                const bam1_t* ptr = entry.owned_records.back().get();
                entry.records.push_back(ptr);
                ReadView view(ptr);
                ReadReferenceSpan span;
                span.valid = true;
                span.tid = view.tid();
                span.start = view.pos();
                span.end = compute_ref_end(view);
                entry.read_spans.push_back(span);
                return true;
            });
        if (!fetch_ok) {
            entry.records.clear();
            entry.read_spans.clear();
        }
        interval_cache.push_back(std::move(entry));
    }
    const double local_fetch_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - local_fetch_started).count();

    struct BinStageStats {
        int64_t breakpoint_hypotheses_total = 0;
        int64_t breakpoint_hypotheses_after_collapse = 0;
        int64_t breakpoint_hypotheses_after_gate = 0;
        int64_t validator_candidates_total = 0;
        int64_t validator_feasible_total = 0;
        int64_t shortlist_empty = 0;
        int64_t shortlist_primary_only = 0;
        int64_t shortlist_with_challenger = 0;
        int64_t expensive_stage_hypotheses_total = 0;
        int64_t event_consensus_calls = 0;
        int64_t event_consensus_rejected = 0;
        int64_t pre_segmentation_gate_rejected = 0;
        int64_t pre_segmentation_gate_no_bilateral_anchor = 0;
        int64_t pre_segmentation_gate_no_precise_or_full_context = 0;
        int64_t event_segmentation_calls = 0;
        int64_t event_segmentation_rejected = 0;
        int64_t te_alignment_calls = 0;
        int64_t te_alignment_rejected = 0;
        int64_t genotype_calls = 0;
        int64_t final_calls = 0;
        double component_build_seconds = 0.0;
        double local_fetch_seconds = 0.0;
        double projection_seconds = 0.0;
        double local_component_refresh_seconds = 0.0;
        double fragment_extract_seconds = 0.0;
        double signal_cache_build_seconds = 0.0;
        double hypothesis_summary_seconds = 0.0;
        double validator_input_count_seconds = 0.0;
        double event_evidence_seconds = 0.0;
        ExpensiveStageTiming timing;
    } bin_stats;
    bin_stats.component_build_seconds = component_build_seconds;
    bin_stats.local_fetch_seconds = local_fetch_seconds;

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

    for (size_t ci = 0; ci < components.size(); ++ci) {
        const auto& component = components[ci];
        const std::vector<InsertionFragment> empty_fragments;
        const auto seed_bounds = component_seed_bounds[ci];
        const auto projection_started = std::chrono::steady_clock::now();
        const LocalIntervalProjection projection =
            project_cached_interval_reads(requests[ci], interval_cache);
        bin_stats.projection_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - projection_started).count();
        if (projection.records.empty()) {
            log_component(
                "LOCAL_EVENT_RECOLLECTION_FAILED",
                component,
                empty_fragments,
                nullptr,
                nullptr);
            continue;
        }

        const auto local_component_refresh_started = std::chrono::steady_clock::now();
        const ComponentCall local_component = build_local_fragment_component(
            config_,
            component,
            projection.records);
        bin_stats.local_component_refresh_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - local_component_refresh_started).count();
        const auto fragment_extract_started = std::chrono::steady_clock::now();
        const std::vector<InsertionFragment> fragments = ins_fragment_module_.extract(
            local_component,
            projection.records);
        bin_stats.fragment_extract_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - fragment_extract_started).count();
        const auto signal_cache_build_started = std::chrono::steady_clock::now();
        const ComponentSignalCache signal_cache = build_component_signal_cache(
            projection.records,
            projection.read_spans);
        bin_stats.signal_cache_build_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - signal_cache_build_started).count();

        if (config_.log_stage_components) {
            const auto all_breakpoint_hypotheses = enumerate_breakpoint_hypotheses(
                component,
                projection.records,
                fragments,
                seed_bounds.first,
                seed_bounds.second,
                0);
            for (const auto& hypothesis : all_breakpoint_hypotheses) {
                std::ostringstream bp_oss;
                bp_oss << "[Pipeline][breakpoint_hypothesis]"
                       << " chrom=" << component.chrom
                       << " anchor=" << component.anchor_pos
                       << " bp=" << hypothesis.left << "-" << hypothesis.right
                       << " center=" << hypothesis.center
                       << " support=" << hypothesis.support
                       << " priority=" << hypothesis.priority
                       << " score="
                       << (hypothesis.support *
                           breakpoint_hypothesis_support_weight(hypothesis.priority));
                emit_pipeline_log_line(bp_oss.str());
            }
        }

        const auto breakpoint_hypotheses = enumerate_breakpoint_hypotheses(
            component,
            projection.records,
            fragments,
            seed_bounds.first,
            seed_bounds.second,
            3);
        bin_stats.breakpoint_hypotheses_total += static_cast<int64_t>(breakpoint_hypotheses.size());
        std::vector<HypothesisSummary> hypothesis_summaries;
        hypothesis_summaries.reserve(breakpoint_hypotheses.size());
        const auto hypothesis_summary_started = std::chrono::steady_clock::now();
        for (size_t hi = 0; hi < breakpoint_hypotheses.size(); ++hi) {
            const auto& hypothesis = breakpoint_hypotheses[hi];
            hypothesis_summaries.push_back(build_hypothesis_summary(
                component,
                signal_cache,
                fragments,
                hi,
                hypothesis.left,
                hypothesis.right));
        }
        bin_stats.hypothesis_summary_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - hypothesis_summary_started).count();
        const auto collapsed_summaries = collapse_hypothesis_summaries(hypothesis_summaries);
        bin_stats.breakpoint_hypotheses_after_collapse +=
            static_cast<int64_t>(collapsed_summaries.size());
        std::vector<HypothesisSummary> coarse_survivors;
        coarse_survivors.reserve(collapsed_summaries.size());
        for (size_t si = 0; si < collapsed_summaries.size(); ++si) {
            if (should_keep_hypothesis_for_expensive_stage(
                    collapsed_summaries[si],
                    si == 0)) {
                coarse_survivors.push_back(collapsed_summaries[si]);
            }
        }
        bin_stats.breakpoint_hypotheses_after_gate +=
            static_cast<int64_t>(coarse_survivors.size());

        std::vector<HypothesisValidatorEvidence> validator_candidates;
        validator_candidates.reserve(coarse_survivors.size());
        for (const auto& summary : coarse_survivors) {
            const EventReadEvidence event_evidence = collect_event_read_evidence_for_bounds_cached(
                component,
                signal_cache,
                fragments,
                summary.bp_left,
                summary.bp_right,
                summary.bp_left,
                summary.bp_right);
            const auto validator_input_started = std::chrono::steady_clock::now();
            const ConsensusInputCounts inputs = collect_event_consensus_input_counts(
                projection.records,
                fragments,
                event_evidence);
            bin_stats.validator_input_count_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - validator_input_started).count();
            validator_candidates.push_back(
                collect_hypothesis_validator_evidence(
                    summary,
                    inputs,
                    component.anchor_pos));
        }
        const auto shortlist = build_expensive_stage_shortlist(validator_candidates);
        bin_stats.validator_candidates_total += static_cast<int64_t>(validator_candidates.size());
        for (const auto& candidate : validator_candidates) {
            if (candidate.feasible_for_expensive_stage) {
                ++bin_stats.validator_feasible_total;
            }
        }
        if (shortlist.empty()) {
            ++bin_stats.shortlist_empty;
        } else if (shortlist.size() == 1) {
            ++bin_stats.shortlist_primary_only;
        } else {
            ++bin_stats.shortlist_with_challenger;
        }
        bin_stats.expensive_stage_hypotheses_total += static_cast<int64_t>(shortlist.size());

        if (config_.log_stage_components) {
            std::ostringstream shortlist_oss;
            shortlist_oss << "[Pipeline][shortlist]"
                          << " chrom=" << component.chrom
                          << " anchor=" << component.anchor_pos
                          << " collapsed_hypotheses=" << collapsed_summaries.size()
                          << " coarse_survivors=" << coarse_survivors.size()
                          << " validator_feasible=" << std::count_if(
                                 validator_candidates.begin(),
                                 validator_candidates.end(),
                                 [](const auto& row) { return row.feasible_for_expensive_stage; })
                          << " shortlist_size=" << shortlist.size();
            if (!shortlist.empty()) {
                shortlist_oss << " primary_bp="
                              << shortlist.front().validator.summary.bp_left
                              << "-" << shortlist.front().validator.summary.bp_right;
            }
            if (shortlist.size() > 1) {
                shortlist_oss << " challenger_bp="
                              << shortlist[1].validator.summary.bp_left
                              << "-" << shortlist[1].validator.summary.bp_right;
            }
            emit_pipeline_log_line(shortlist_oss.str());

            for (const auto& candidate : validator_candidates) {
                std::ostringstream candidate_oss;
                candidate_oss << "[Pipeline][shortlist_candidate]"
                              << " chrom=" << component.chrom
                              << " anchor=" << component.anchor_pos
                              << " bp=" << candidate.summary.bp_left
                              << "-" << candidate.summary.bp_right
                              << " precise=" << candidate.precise_support
                              << " alt_struct=" << candidate.summary.alt_struct_reads
                              << " ref_span=" << candidate.summary.ref_span_reads
                              << " full_context=" << candidate.full_context_input_reads
                              << " partial_context=" << candidate.partial_context_input_reads
                              << " input_event_reads=" << candidate.input_event_reads
                              << " left_anchor=" << candidate.left_anchor_input_reads
                              << " right_anchor=" << candidate.right_anchor_input_reads
                              << " feasible=" << bool_as_int(candidate.feasible_for_expensive_stage)
                              << " qc=" << candidate.qc_reason;
                emit_pipeline_log_line(candidate_oss.str());
            }
        }

        bool have_best = false;
        double best_total = -1e9;
        JointDecisionResult best_joint;
        EventReadEvidence best_event_evidence;
        EventConsensus best_event_consensus;
        EventSegmentation best_event_segmentation;
        TEAlignmentEvidence best_te_alignment;
        BoundaryEvidence best_boundary_evidence;
        GenotypeCall best_genotype;

        struct ComponentFinalCallEvaluation {
            FinalCall call;
            GenotypeCall genotype;
            double score = 0.0;
        };
        std::vector<ComponentFinalCallEvaluation> final_call_candidates;
        final_call_candidates.reserve(shortlist.size());

        for (const auto& shortlisted : shortlist) {
            const auto& summary = shortlisted.validator.summary;
            const auto event_evidence_started = std::chrono::steady_clock::now();
            const EventReadEvidence event_evidence = collect_event_read_evidence_for_bounds_cached(
                component,
                signal_cache,
                fragments,
                summary.bp_left,
                summary.bp_right,
                summary.bp_left,
                summary.bp_right);
            bin_stats.event_evidence_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - event_evidence_started).count();
            const auto consensus_started = std::chrono::steady_clock::now();
            const EventConsensus event_consensus = build_event_consensus(
                component,
                projection.records,
                fragments,
                event_evidence);
            bin_stats.timing.event_consensus_seconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - consensus_started).count();
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
            existence_input.alt_observed_lengths =
                collect_alt_observed_lengths(component, event_evidence);
            existence_input.event_length =
                infer_event_length_from_alt_support(existence_input.alt_observed_lengths);
            const EventExistenceEvidence existence = build_event_existence_evidence(
                existence_input);

            EventSegmentation event_segmentation;
            const std::string preseg_qc = pre_segmentation_gate_reason(
                event_evidence,
                event_consensus);
            if (!preseg_qc.empty()) {
                event_segmentation.qc_reason = preseg_qc;
                ++bin_stats.pre_segmentation_gate_rejected;
                if (preseg_qc == "PRESEG_NO_BILATERAL_ANCHOR") {
                    ++bin_stats.pre_segmentation_gate_no_bilateral_anchor;
                } else if (preseg_qc == "PRESEG_NO_PRECISE_OR_FULL_CONTEXT") {
                    ++bin_stats.pre_segmentation_gate_no_precise_or_full_context;
                }
            } else {
                const auto segmentation_started = std::chrono::steady_clock::now();
                event_segmentation = segment_event_consensus(
                    component,
                    event_evidence,
                    event_consensus);
                bin_stats.timing.event_segmentation_seconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - segmentation_started).count();
                bin_stats.event_segmentation_calls += 1;
                if (!event_segmentation.pass) {
                    bin_stats.event_segmentation_rejected += 1;
                }
            }
            const EventSegmentationEvidence seg_evidence = analyze_event_segmentation(
                event_consensus,
                event_segmentation);

            TEAlignmentEvidence te_alignment;
            if (seg_evidence.has_insert_seq) {
                const auto te_started = std::chrono::steady_clock::now();
                te_alignment = align_insert_seq_to_te(event_segmentation);
                bin_stats.timing.te_alignment_seconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - te_started).count();
                bin_stats.te_alignment_calls += 1;
                if (!te_alignment.pass) {
                    bin_stats.te_alignment_rejected += 1;
                }
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

            if (config_.log_stage_components) {
                std::ostringstream candidate_joint_oss;
                candidate_joint_oss << "[Pipeline][shortlist_result]"
                                    << " chrom=" << component.chrom
                                    << " anchor=" << component.anchor_pos
                                    << " bp=" << event_evidence.bp_left
                                    << "-" << event_evidence.bp_right
                                    << " primary=" << bool_as_int(shortlisted.is_primary)
                                    << " alt_reads=" << event_evidence.alt_struct_reads
                                    << " ref_reads=" << event_evidence.ref_span_reads
                                    << " consensus_qc=" << display_or_na(event_consensus.qc_reason)
                                    << " consensus_len=" << event_consensus.consensus_len
                                    << " seg_qc=" << display_or_na(event_segmentation.qc_reason)
                                    << " seg_insert_len=" << event_segmentation.insert_seq.size()
                                    << " te_qc=" << display_or_na(te_alignment.qc_reason)
                                    << " te_pass=" << bool_as_int(te_alignment.pass)
                                    << " best_kind=" << final_hypothesis_kind_name(joint.best.kind)
                                    << " best_total=" << joint.best.total
                                    << " runner_kind="
                                    << final_hypothesis_kind_name(joint.runner_up.kind)
                                    << " runner_total=" << joint.runner_up.total
                                    << " final_qc=" << display_or_na(joint.final_qc)
                                    << " emit_te=" << bool_as_int(joint.emit_te_call)
                                    << " emit_unknown=" << bool_as_int(joint.emit_unknown_te);
                emit_pipeline_log_line(candidate_joint_oss.str());
            }

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

            if (joint.emit_te_call) {
                FinalBoundaryDecision boundary;
                boundary.pass =
                    boundary_evidence.canonical_pass ||
                    boundary_evidence.evidence_consistent;
                boundary.boundary_type = boundary_evidence.boundary_type;
                boundary.boundary_len = boundary_evidence.boundary_len;
                boundary.qc = boundary_evidence.qc;

                FinalTeAcceptanceDecision acceptance;
                acceptance.pass = true;
                acceptance.qc = joint.final_qc;

                ComponentFinalCallEvaluation evaluation;
                evaluation.call = emit_final_te_call(
                    component,
                    event_evidence,
                    event_consensus,
                    event_segmentation,
                    te_alignment,
                    genotype,
                    boundary,
                    acceptance);
                evaluation.genotype = genotype;
                evaluation.score = candidate_total;
                final_call_candidates.push_back(std::move(evaluation));
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
                << " te_coarse_prefilter=" << best_te_alignment.coarse_prefilter_score
                << " te_coarse_cov=" << best_te_alignment.coarse_chain_coverage
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

        if (shortlist.empty()) {
            log_component(
                "NO_FEASIBLE_EVENT_HYPOTHESIS",
                component,
                fragments,
                nullptr,
                nullptr);
            continue;
        }

        if (!have_best || final_call_candidates.empty()) {
            log_component(
                have_best ? best_joint.final_qc.c_str() : "NO_BREAKPOINT_HYPOTHESIS",
                component,
                fragments,
                have_best ? &best_genotype : nullptr,
                nullptr);
            continue;
        }

        std::vector<ComponentFinalCallCandidate> selection_candidates;
        selection_candidates.reserve(final_call_candidates.size());
        for (const auto& candidate : final_call_candidates) {
            ComponentFinalCallCandidate selection_candidate;
            selection_candidate.pos = candidate.call.pos;
            selection_candidate.score = candidate.score;
            selection_candidate.emit_te = true;
            selection_candidates.push_back(selection_candidate);
        }
        const auto selected_indices = select_component_final_call_indices(selection_candidates);
        if (selected_indices.empty()) {
            log_component(
                best_joint.final_qc.c_str(),
                component,
                fragments,
                &best_genotype,
                nullptr);
            continue;
        }

        for (size_t selected_index : selected_indices) {
            const auto& selected = final_call_candidates[selected_index];
            result.final_calls.push_back(selected.call);
            bin_stats.final_calls += 1;

            log_component(
                selected.call.final_qc.c_str(),
                component,
                fragments,
                &selected.genotype,
                &result.final_calls.back());
        }
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
                << " breakpoint_hypotheses_total=" << bin_stats.breakpoint_hypotheses_total
                << " breakpoint_hypotheses_after_collapse="
                << bin_stats.breakpoint_hypotheses_after_collapse
                << " breakpoint_hypotheses_after_gate=" << bin_stats.breakpoint_hypotheses_after_gate
                << " validator_candidates_total=" << bin_stats.validator_candidates_total
                << " validator_feasible_total=" << bin_stats.validator_feasible_total
                << " shortlist_empty=" << bin_stats.shortlist_empty
                << " shortlist_primary_only=" << bin_stats.shortlist_primary_only
                << " shortlist_with_challenger=" << bin_stats.shortlist_with_challenger
                << " expensive_stage_hypotheses_total=" << bin_stats.expensive_stage_hypotheses_total
                << " event_consensus_calls=" << bin_stats.event_consensus_calls
                << " event_consensus_rejected=" << bin_stats.event_consensus_rejected
                << " pre_segmentation_gate_rejected=" << bin_stats.pre_segmentation_gate_rejected
                << " pre_segmentation_gate_no_bilateral_anchor="
                << bin_stats.pre_segmentation_gate_no_bilateral_anchor
                << " pre_segmentation_gate_no_precise_or_full_context="
                << bin_stats.pre_segmentation_gate_no_precise_or_full_context
                << " event_segmentation_calls=" << bin_stats.event_segmentation_calls
                << " event_segmentation_rejected=" << bin_stats.event_segmentation_rejected
                << " te_alignment_calls=" << bin_stats.te_alignment_calls
                << " te_alignment_rejected=" << bin_stats.te_alignment_rejected
                << " component_build_seconds=" << bin_stats.component_build_seconds
                << " local_fetch_seconds=" << bin_stats.local_fetch_seconds
                << " projection_seconds=" << bin_stats.projection_seconds
                << " local_component_refresh_seconds=" << bin_stats.local_component_refresh_seconds
                << " fragment_extract_seconds=" << bin_stats.fragment_extract_seconds
                << " signal_cache_build_seconds=" << bin_stats.signal_cache_build_seconds
                << " hypothesis_summary_seconds=" << bin_stats.hypothesis_summary_seconds
                << " validator_input_count_seconds=" << bin_stats.validator_input_count_seconds
                << " event_evidence_seconds=" << bin_stats.event_evidence_seconds
                << " event_consensus_seconds=" << bin_stats.timing.event_consensus_seconds
                << " event_segmentation_seconds=" << bin_stats.timing.event_segmentation_seconds
                << " te_alignment_seconds=" << bin_stats.timing.te_alignment_seconds
                << " genotype_calls=" << bin_stats.genotype_calls
                << " final_calls=" << bin_stats.final_calls
                << " local_reuse_ratio=" << local_interval_reuse_ratio(reuse_stats);
        emit_pipeline_log_line(oss.str());
    }
}

std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config) {
    auto reader = make_bam_reader(config.bam_path, config.bam_threads, config.bam_region_scope);
    return std::make_unique<Pipeline>(config, std::move(reader));
}

}  // namespace placer
