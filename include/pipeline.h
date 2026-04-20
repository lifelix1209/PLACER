#ifndef PLACER_PIPELINE_H
#define PLACER_PIPELINE_H

#include "bam_io.h"
#include "decision_policy.h"
#include "gate1_module.h"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace placer {

struct FinalBoundaryDecision;
struct FinalTeAcceptanceDecision;

enum CandidateClassMask : uint8_t {
    kCandidateSoftClip = 1 << 0,
    kCandidateSplitSaSupplementary = 1 << 1,
    kCandidateLongInsertion = 1 << 2
};

struct BreakpointCandidate {
    std::string chrom;
    int32_t pos = -1;
    bool is_reverse = false;

    int32_t anchor_len = 0;
    int32_t clip_len = 0;
    int32_t ins_len = 0;

    std::string read_id;
    size_t read_index = 0;
    uint8_t class_mask = 0;
};

struct ComponentCall {
    std::string chrom;
    int32_t tid = -1;
    int32_t bin_start = -1;
    int32_t bin_end = -1;
    int32_t anchor_pos = -1;
    double peak_weight = 0.0;
    int32_t evidence_soft_clip_count = 0;
    int32_t evidence_indel_count = 0;
    int32_t evidence_sa_hint_count = 0;

    std::vector<size_t> read_indices;
    std::vector<size_t> soft_clip_read_indices;
    std::vector<size_t> split_sa_read_indices;
    std::vector<size_t> insertion_read_indices;
    std::vector<BreakpointCandidate> breakpoint_candidates;
};

enum class InsertionFragmentSource : uint8_t {
    kUnknown = 0,
    kClipRefLeft = 1,   // soft-clip on reference-left side
    kClipRefRight = 2,  // soft-clip on reference-right side
    kCigarInsertion = 3,
    kSplitSa = 4
};

enum class ReferenceSide : uint8_t {
    kUnknown = 0,
    kRefLeft = 1,
    kRefRight = 2
};

struct InsertionFragment {
    std::string fragment_id;
    std::string chrom;
    int32_t anchor_pos = -1;

    std::string read_id;
    size_t read_index = 0;
    uint8_t class_mask = 0;
    bool is_reverse = false;  // BAM_FREVERSE

    InsertionFragmentSource source = InsertionFragmentSource::kUnknown;
    int32_t start = -1;  // read coordinate (0-based)
    int32_t length = 0;
    int32_t read_len = 0;

    // Anchor/junction diagnostics used by event-level breakpoint inference.
    int32_t anchor_len = 0;
    ReferenceSide ref_side = ReferenceSide::kUnknown;
    int32_t ref_junc_pos = -1;
    int32_t nm = -1;
    bool split_sa_reliable = false;

    std::string sequence;
};

struct FragmentTEHit {
    std::string fragment_id;
    std::string te_name;

    int32_t fragment_len = 0;
    int32_t aligned_len_est = 0;
    double kmer_support = 0.0;  // quick estimate [0,1], not sequence identity
    double coverage = 0.0;      // quick estimate [0,1]
    double multik_support = 0.0;
    bool rescue_used = false;

    int32_t hit_kmers = 0;
    int32_t total_kmers = 0;
};

struct GenotypeCall {
    int32_t tid = -1;
    int32_t pos = -1;
    std::string genotype = "./.";
    double af = 0.0;
    int32_t gq = 0;
};

struct ReadReferenceSpan {
    bool valid = false;
    int32_t tid = -1;
    int32_t start = -1;
    int32_t end = -1;
};

struct EventReadEvidence {
    int32_t bp_left = -1;
    int32_t bp_right = -1;
    int32_t alt_split_reads = 0;
    int32_t alt_indel_reads = 0;
    int32_t alt_left_clip_reads = 0;
    int32_t alt_right_clip_reads = 0;
    int32_t alt_struct_reads = 0;
    int32_t ref_span_reads = 0;
    int32_t low_mapq_ref_span_reads = 0;
    std::vector<std::string> support_qnames;
    std::vector<std::string> ref_span_qnames;
};

struct EventConsensus {
    std::string consensus_seq;
    int32_t input_event_reads = 0;
    int32_t consensus_len = 0;
    int32_t full_context_input_reads = 0;
    int32_t partial_context_input_reads = 0;
    int32_t left_anchor_input_reads = 0;
    int32_t right_anchor_input_reads = 0;
    bool used_full_context = false;
    bool qc_pass = false;
    std::string qc_reason = "NO_EVENT_CONSENSUS";
};

struct EventSegmentation {
    std::string left_flank_seq;
    std::string insert_seq;
    std::string right_flank_seq;
    int32_t left_ref_start = -1;
    int32_t left_ref_end = -1;
    int32_t right_ref_start = -1;
    int32_t right_ref_end = -1;
    int32_t left_flank_align_len = 0;
    int32_t right_flank_align_len = 0;
    double left_flank_identity = 0.0;
    double right_flank_identity = 0.0;
    bool pass = false;
    std::string qc_reason = "NO_EVENT_SEGMENTATION";
};

struct TEAlignmentEvidence {
    std::string best_family;
    std::string best_subfamily;
    double best_identity = 0.0;
    double best_query_coverage = 0.0;
    double best_score = 0.0;
    std::string second_family = "NA";
    double second_score = 0.0;
    double cross_family_margin = 0.0;
    bool pass = false;
    std::string qc_reason = "NO_TE_ALIGNMENT";
};

struct FinalCall {
    std::string chrom;
    int32_t tid = -1;
    int32_t pos = -1;
    int32_t bp_left = -1;
    int32_t bp_right = -1;

    int32_t window_start = -1;
    int32_t window_end = -1;

    std::string te_name;

    std::string tsd_type = "NONE";
    int32_t tsd_len = 0;
    std::string tsd_seq = "NA";
    double tsd_bg_p = 1.0;

    int32_t support_reads = 0;
    int32_t alt_struct_reads = 0;
    int32_t ref_span_reads = 0;
    int32_t low_mapq_ref_span_reads = 0;
    std::string genotype = "./.";
    double af = 0.0;
    int32_t gq = 0;

    std::string family = "NA";
    std::string subfamily = "NA";
    std::string strand = "NA";
    int32_t insert_len = 0;
    double best_te_identity = 0.0;
    double best_te_query_coverage = 0.0;
    double cross_family_margin = 0.0;
    int32_t left_flank_align_len = 0;
    int32_t right_flank_align_len = 0;
    int32_t event_consensus_len = 0;
    std::string te_qc = "NA";
    std::string final_qc = "NA";
};

struct PipelineConfig {
    std::string bam_path;
    std::string reference_fasta_path;
    std::string te_fasta_path;

    int32_t bam_threads = 2;
    int64_t progress_interval = 100000;
    bool log_stage_bins = false;
    bool log_stage_components = false;

    int32_t window_size = 10000;
    int32_t bin_size = 10000;

    // Module 2.1: insertion fragment extraction (optional).
    // Empty path disables fragment FASTA emission (recommended default for large runs).
    std::string ins_fragments_fasta_path;
    int32_t min_soft_clip_for_seq_extract = 50;
    int32_t min_long_ins_for_seq_extract = 50;
    int32_t min_sa_aln_len_for_seq_extract = 50;
    int32_t max_sa_per_read = 3;

    // Module 2.2: quick TE classification (optional, requires te_fasta_path).
    // Empty path disables per-fragment hit TSV emission (recommended default for large runs).
    std::string ins_fragment_hits_tsv_path;
    int32_t te_kmer_size = 13;
    std::string te_kmer_sizes_csv = "9,11,13";
    double te_low_kmer_support_trigger = 0.30;
    bool te_low_kmer_rescue_enable = true;
    int32_t te_low_kmer_rescue_topn = 3;
    int32_t te_low_kmer_rescue_min_frag_len = 40;
    double te_low_kmer_rescue_identity_min = 0.55;
    double te_low_kmer_rescue_margin_max = 0.08;
    double te_one_sided_breakpoint_mad_max = 40.0;
    bool short_ins_enable = true;
    int32_t short_ins_min_len = 35;
    int32_t short_ins_max_len = 300;
    int32_t short_ins_min_reads = 2;
    double short_ins_kmer_relax_identity = 0.15;
    double te_softclip_low_complexity_at_frac_min = 0.90;
    int32_t te_softclip_low_complexity_homopolymer_min = 80;
    double te_softclip_entropy_min = 1.25;
    double te_softclip_kmer_uniqueness_min = 0.35;
    int32_t te_softclip_min_anchor_len = 20;
    double te_softclip_max_nm_per_bp = 0.12;

    // Module 2.3: TSD detector (duplication / deletion).
    bool tsd_enable = true;
    int32_t tsd_min_len = 3;
    int32_t tsd_max_len = 50;
    int32_t tsd_flank_window = 150;
    double tsd_bg_p_max = 0.05;

    // Genotype likelihood model.
    int32_t genotype_min_depth = 3;
    double genotype_error_rate = 0.02;

    // Event consensus (abPOA on event strings).
    int32_t event_consensus_poa_min_reads = 2;
    int32_t event_consensus_poa_max_reads = 48;

    bool enable_parallel = false;
    size_t batch_size = 1000;
    int32_t parallel_workers = 0;
    int32_t parallel_queue_max_tasks = 0;  // <=0 means auto-size bounded queue
    int32_t parallel_result_buffer_max = 0;  // <=0 means auto-size bounded queue
    bool log_parallel_progress = false;
};

struct PipelineResult {
    int64_t total_reads = 0;
    int64_t gate1_passed = 0;
    int64_t processed_bins = 0;
    int64_t built_components = 0;
    int64_t event_consensus_calls = 0;
    int64_t genotype_calls = 0;
    int64_t final_pass_calls = 0;

    std::vector<FinalCall> final_calls;
};

void finalize_final_calls(PipelineResult& result);

class DbscanComponentModule final {
public:
    std::vector<ComponentCall> build(
        const std::vector<const bam1_t*>& records,
        const std::string& chrom,
        int32_t tid) const;
};

class LinearBinComponentModule final {
public:
    std::vector<ComponentCall> build(
        const std::vector<const bam1_t*>& bin_records,
        const std::string& chrom,
        int32_t tid,
        int32_t bin_start,
        int32_t bin_end) const;
};

class CigarInsertionFragmentModule final {
public:
    explicit CigarInsertionFragmentModule(PipelineConfig config);

    std::vector<InsertionFragment> extract(
        const ComponentCall& component,
        const std::vector<const bam1_t*>& bin_records) const;

private:
    PipelineConfig config_;
};

class SplitSAFragmentModule final {
public:
    explicit SplitSAFragmentModule(PipelineConfig config);

    std::vector<InsertionFragment> extract(
        const ComponentCall& component,
        const std::vector<const bam1_t*>& bin_records) const;

private:
    PipelineConfig config_;
};

class TEKmerQuickClassifierModule final {
public:
    explicit TEKmerQuickClassifierModule(PipelineConfig config);

    bool is_enabled() const;
    std::vector<FragmentTEHit> classify(
        const std::vector<InsertionFragment>& fragments) const;
    TEAlignmentEvidence align_insert_sequence(
        const std::string& insert_seq) const;

private:
    PipelineConfig config_;
    struct Index;
    struct AlignmentShortlistDb;
    std::shared_ptr<const Index> primary_index_;
    std::vector<std::shared_ptr<const Index>> indices_;
    std::shared_ptr<const AlignmentShortlistDb> alignment_shortlist_db_;
    std::shared_ptr<const std::vector<std::string>> te_names_;
    std::shared_ptr<const std::vector<std::string>> te_sequences_;
    std::shared_ptr<const std::vector<std::string>> te_reverse_complement_sequences_;
};

struct TsdDetection {
    std::string type = "NONE";  // DUP / DEL / NONE / UNCERTAIN
    int32_t length = 0;
    std::string sequence;
    double bg_p = 1.0;
    bool significant = false;
};

struct BufferedRecord {
    BamRecordPtr record;
    int32_t ref_end = -1;
};

class TSDDetector final {
public:
    explicit TSDDetector(PipelineConfig config);

    bool is_enabled() const;
    bool can_fetch_reference() const;
    std::string fetch_window(
        const std::string& chrom,
        int32_t start,
        int32_t end) const;
    TsdDetection detect(
        const std::string& chrom,
        int32_t left_bp,
        int32_t right_bp) const;

private:
    PipelineConfig config_;
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};

class Pipeline {
public:
    Pipeline(PipelineConfig config, std::unique_ptr<BamStreamReader> bam_reader);

    PipelineResult run() const;

private:
    struct WindowCoord {
        int32_t tid = -1;
        int32_t pos = -1;
    };

    struct StreamingState {
        struct BinSnapshot {
            int32_t tid = -1;
            int32_t bin_index = -1;
            std::vector<BufferedRecord> records;
        };

        std::deque<WindowCoord> active_window;
        std::deque<BinSnapshot> recent_bin_snapshots;
        std::vector<BufferedRecord> current_bin_records;
        int32_t current_tid = -1;
        int32_t current_bin_index = -1;
    };

    struct HypothesisSummary {
        size_t original_index = 0;
        int32_t bp_left = -1;
        int32_t bp_right = -1;
        int32_t alt_split_reads = 0;
        int32_t alt_indel_reads = 0;
        int32_t alt_left_clip_reads = 0;
        int32_t alt_right_clip_reads = 0;
        int32_t alt_struct_reads = 0;
        int32_t ref_span_reads = 0;
        int32_t inferred_event_length = 0;
        std::vector<std::string> support_qnames;
    };

    struct BreakpointHypothesis {
        bool valid = false;
        int32_t left = -1;
        int32_t right = -1;
        int32_t center = -1;
        int32_t support = 0;
        int32_t priority = std::numeric_limits<int32_t>::max();
    };

    struct ExpensiveStageTiming {
        double event_consensus_seconds = 0.0;
        double event_segmentation_seconds = 0.0;
        double te_alignment_seconds = 0.0;
    };

    struct ConsensusInputSummary {
        std::unordered_map<std::string, std::string> full_event_by_qname;
        std::unordered_map<std::string, std::string> partial_event_by_qname;
        int32_t full_context_input_reads = 0;
        int32_t partial_context_input_reads = 0;
        int32_t left_anchor_input_reads = 0;
        int32_t right_anchor_input_reads = 0;
    };

    struct ConsensusInputCounts {
        int32_t full_context_input_reads = 0;
        int32_t partial_context_input_reads = 0;
        int32_t left_anchor_input_reads = 0;
        int32_t right_anchor_input_reads = 0;
        int32_t input_event_reads = 0;
    };

    struct CachedLocalSignal {
        std::string qname;
        int32_t tid = -1;
        int32_t mapq = 0;
        bool is_primary = false;
        bool has_sa_or_supp = false;
        bool span_valid = false;
        int32_t span_start = -1;
        int32_t span_end = -1;
        int32_t left_clip_pos = -1;
        int32_t right_clip_pos = -1;
        std::vector<int32_t> indel_positions;
    };

    struct ComponentSignalCache {
        std::vector<CachedLocalSignal> reads;
    };

    struct HypothesisValidatorEvidence {
        HypothesisSummary summary;
        int32_t precise_support = 0;
        int32_t breakpoint_width = 0;
        int32_t anchor_distance = 0;
        int32_t full_context_input_reads = 0;
        int32_t partial_context_input_reads = 0;
        int32_t left_anchor_input_reads = 0;
        int32_t right_anchor_input_reads = 0;
        int32_t input_event_reads = 0;
        bool feasible_for_expensive_stage = false;
        std::string qc_reason = "VALIDATOR_UNSET";
    };

    struct ShortlistedHypothesis {
        HypothesisValidatorEvidence validator;
        bool is_primary = false;
    };

    struct AnchorSeedBin {
        int32_t ref_bin_start = -1;
        int32_t support = 0;
        int32_t best_breakpoint_delta = std::numeric_limits<int32_t>::max();
    };

    PipelineResult run_streaming() const;
    PipelineResult run_parallel() const;

    HypothesisSummary build_hypothesis_summary(
        const ComponentCall& component,
        const ComponentSignalCache& signal_cache,
        const std::vector<InsertionFragment>& fragments,
        size_t original_index,
        int32_t bp_left,
        int32_t bp_right) const;

    std::vector<HypothesisSummary> select_hypothesis_summaries_for_expensive_stage(
        const std::vector<HypothesisSummary>& summaries) const;

    std::vector<HypothesisSummary> collapse_hypothesis_summaries(
        const std::vector<HypothesisSummary>& summaries) const;

    bool should_keep_hypothesis_for_expensive_stage(
        const HypothesisSummary& summary,
        bool is_top_ranked_survivor) const;

    void consume_record(
        BamRecordPtr&& record,
        StreamingState& state,
        PipelineResult& result) const;

    void flush_current_bin(
        StreamingState& state,
        PipelineResult& result) const;

    void process_bin_records(
        std::vector<const bam1_t*>&& bin_records,
        int32_t tid,
        int32_t bin_index,
        int32_t owner_bin_start,
        int32_t owner_bin_end,
        PipelineResult& result) const;

    EventReadEvidence collect_event_read_evidence(
        const ComponentCall& component,
        const std::vector<const bam1_t*>& local_records,
        const std::vector<ReadReferenceSpan>& read_spans,
        const std::vector<InsertionFragment>& fragments) const;

    std::pair<int32_t, int32_t> resolve_event_breakpoint_bounds(
        const ComponentCall& component,
        const std::vector<const bam1_t*>& local_records,
        const std::vector<InsertionFragment>& fragments,
        int32_t seed_left,
        int32_t seed_right) const;

    std::vector<BreakpointHypothesis> enumerate_breakpoint_hypotheses(
        const ComponentCall& component,
        const std::vector<const bam1_t*>& local_records,
        const std::vector<InsertionFragment>& fragments,
        int32_t seed_left,
        int32_t seed_right,
        size_t top_k) const;

    std::vector<BreakpointHypothesis> select_diverse_breakpoint_hypotheses(
        const std::vector<BreakpointHypothesis>& hypotheses,
        size_t top_k) const;

    EventReadEvidence collect_event_read_evidence_for_bounds(
        const ComponentCall& component,
        const std::vector<const bam1_t*>& local_records,
        const std::vector<ReadReferenceSpan>& read_spans,
        const std::vector<InsertionFragment>& fragments,
        int32_t seed_left,
        int32_t seed_right,
        int32_t bp_left,
        int32_t bp_right) const;

    ComponentSignalCache build_component_signal_cache(
        const std::vector<const bam1_t*>& local_records,
        const std::vector<ReadReferenceSpan>& read_spans) const;

    EventReadEvidence collect_event_read_evidence_for_bounds_cached(
        const ComponentCall& component,
        const ComponentSignalCache& cache,
        const std::vector<InsertionFragment>& fragments,
        int32_t seed_left,
        int32_t seed_right,
        int32_t bp_left,
        int32_t bp_right) const;

    ConsensusInputSummary collect_event_consensus_inputs(
        const std::vector<const bam1_t*>& local_records,
        const std::vector<InsertionFragment>& fragments,
        const EventReadEvidence& event_evidence) const;

    ConsensusInputCounts collect_event_consensus_input_counts(
        const std::vector<const bam1_t*>& local_records,
        const std::vector<InsertionFragment>& fragments,
        const EventReadEvidence& event_evidence) const;

    HypothesisValidatorEvidence collect_hypothesis_validator_evidence(
        const HypothesisSummary& summary,
        const ConsensusInputCounts& inputs,
        int32_t anchor_pos) const;

    bool has_strong_bilateral_partial_context_support(
        int32_t left_anchor_input_reads,
        int32_t right_anchor_input_reads,
        int32_t partial_context_input_reads,
        int32_t input_event_reads) const;

    bool compare_hypothesis_validator_priority(
        const HypothesisValidatorEvidence& lhs,
        const HypothesisValidatorEvidence& rhs) const;

    std::vector<ShortlistedHypothesis> build_expensive_stage_shortlist(
        const std::vector<HypothesisValidatorEvidence>& candidates) const;

    EventConsensus build_event_consensus(
        const ComponentCall& component,
        const std::vector<const bam1_t*>& local_records,
        const std::vector<InsertionFragment>& fragments,
        const EventReadEvidence& event_evidence) const;

    std::string pre_segmentation_gate_reason(
        const EventReadEvidence& event_evidence,
        const EventConsensus& event_consensus) const;

    std::vector<AnchorSeedBin> collect_anchor_seed_bins(
        const std::string& query,
        int32_t breakpoint,
        int32_t ref_window_start,
        const std::string& ref_window,
        bool is_left) const;

    EventSegmentation segment_event_consensus(
        const ComponentCall& component,
        const EventReadEvidence& event_evidence,
        const EventConsensus& event_consensus) const;

    EventSegmentationEvidence analyze_event_segmentation(
        const EventConsensus& event_consensus,
        const EventSegmentation& event_segmentation) const;

    TEAlignmentEvidence align_insert_seq_to_te(
        const EventSegmentation& event_segmentation) const;

    GenotypeCall genotype_call(
        const ComponentCall& component,
        const EventReadEvidence& event_evidence) const;

    BoundaryEvidence analyze_boundary(
        const EventSegmentation& event_segmentation,
        int32_t breakpoint_envelope_width) const;

    FinalCall emit_final_te_call(
        const ComponentCall& component,
        const EventReadEvidence& event_evidence,
        const EventConsensus& event_consensus,
        const EventSegmentation& event_segmentation,
        const TEAlignmentEvidence& te_alignment,
        const GenotypeCall& genotype,
        const FinalBoundaryDecision& boundary,
        const FinalTeAcceptanceDecision& acceptance) const;

    PipelineConfig config_;
    std::unique_ptr<BamStreamReader> bam_reader_;
    SignalFirstGate1Module gate1_module_;
    LinearBinComponentModule component_module_;
    CigarInsertionFragmentModule ins_fragment_module_;
    TEKmerQuickClassifierModule te_classifier_module_;
    TSDDetector tsd_detector_;
};

std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config);

}  // namespace placer

#endif  // PLACER_PIPELINE_H
