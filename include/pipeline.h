#ifndef PLACER_PIPELINE_H
#define PLACER_PIPELINE_H

#include "bam_io.h"
#include "decision_policy.h"
#include "gate1_module.h"
#include "te_sequence_explainer.h"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef PLACER_BLASTN_EXECUTABLE
#define PLACER_BLASTN_EXECUTABLE "blastn"
#endif

#ifndef PLACER_MAKEBLASTDB_EXECUTABLE
#define PLACER_MAKEBLASTDB_EXECUTABLE "makeblastdb"
#endif

namespace placer {

struct FinalBoundaryDecision;
struct TeFamilyAlignmentIndex;
struct TeFamilyGroupCache;
class TSDDetector;

enum CandidateClassMask : uint8_t {
    kCandidateSoftClip = 1 << 0,
    kCandidateSplitSaSupplementary = 1 << 1,
    kCandidateLongInsertion = 1 << 2
};

constexpr int32_t kInsertionCandidateRequiredMapq = 60;

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
    int32_t raw_cigar_insert_reads = 0;
    int32_t max_raw_cigar_insert_len = 0;
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
    double best_identity = 0.0;        // Identity on the aligned TE core.
    double best_query_coverage = 0.0;  // Fraction of insert occupied by the aligned TE core.
    double best_score = 0.0;
    // L2 (profile depth): where the insert aligns on the TE consensus. The start
    // is the profile-depth signature of 5' truncation (a full-length insertion
    // starts near 0; a 5'-truncated L1 starts hundreds/thousands of bp in); the
    // end gives the 3' extent. -1 when no consensus alignment is available.
    int32_t te_consensus_start = -1;
    int32_t te_consensus_end = -1;
    double coarse_prefilter_score = 0.0;
    double coarse_chain_coverage = 0.0;
    std::string second_family = "NA";
    double second_score = 0.0;
    double cross_family_margin = 0.0;
    std::string sequence_model_label = "TE_MODEL_UNAVAILABLE";
    double sequence_model_score = 0.0;
    double sequence_model_gc = 0.0;
    double sequence_model_entropy = 0.0;
    double sequence_model_tandem_fraction = 0.0;
    double sequence_model_low_complexity_fraction = 0.0;
    double sequence_model_jsd_k5 = 0.0;
    double sequence_model_jsd_k6 = 0.0;
    double sequence_model_k9_containment = 0.0;
    std::string annotation_confidence = "NA";
    std::string annotation_class = "NA";
    std::string annotation_order = "NA";
    std::string annotation_intervals = "NA";
    double annotation_residual_fraction = 0.0;
    double annotation_masked_fraction = 0.0;
    TeSequenceExplanation te_sequence_explanation;
    bool pass = false;
    std::string qc_reason = "NO_TE_ALIGNMENT";
};

struct FinalCall {
    std::string chrom;
    int32_t tid = -1;
    int32_t pos = -1;
    int32_t bp_left = -1;
    int32_t bp_right = -1;

    // L1: breakpoint-position posterior summary. Each supporting signal (split/SA,
    // CIGAR insertion, soft clip) contributes a Gaussian kernel at its position
    // with a scale set by its precision; bp_ci_width is the 90% credible-interval
    // width (bp) of the combined posterior and bp_posterior_entropy its
    // normalised entropy. Small width / low entropy = a precise breakpoint; a call
    // can be TE-positive while breakpoint-imprecise.
    double bp_ci_width = 0.0;
    double bp_posterior_entropy = 0.0;

    int32_t window_start = -1;
    int32_t window_end = -1;

    std::string te_name;

    std::string tsd_type = "NONE";
    int32_t tsd_len = 0;
    std::string tsd_seq = "NA";
    double tsd_bg_p = 1.0;

    int32_t support_reads = 0;
    int32_t alt_struct_reads = 0;
    int32_t raw_cigar_insert_reads = 0;
    int32_t max_raw_cigar_insert_len = 0;
    int32_t ref_span_reads = 0;
    int32_t low_mapq_ref_span_reads = 0;
    std::vector<std::string> support_qnames;
    std::string genotype = "./.";
    double af = 0.0;
    int32_t gq = 0;

    std::string family = "NA";
    std::string subfamily = "NA";
    // Label value and abstention state are separate because a TE library may
    // legitimately contain a family named "Unknown".
    bool family_committed = false;
    // A resolved sequence annotation may be retained while the event-level
    // decision abstains on family. Finalization commits this label only after
    // detection selection is complete, so family coverage cannot influence
    // breakpoint choice, deduplication, or event emission.
    std::string sequence_family_candidate = "NA";
    std::string sequence_subfamily_candidate = "NA";
    bool sequence_family_commit_eligible = false;
    std::string strand = "NA";
    int32_t insert_len = 0;
    double best_te_identity = 0.0;
    double best_te_query_coverage = 0.0;
    double cross_family_margin = 0.0;
    // L2 profile depth: insert placement on the TE consensus (5'-truncation
    // signature = start; 3' extent = end). -1 when unavailable.
    int32_t te_consensus_start = -1;
    int32_t te_consensus_end = -1;
    std::string te_sequence_model_label = "TE_MODEL_UNAVAILABLE";
    double te_sequence_model_score = 0.0;
    double te_sequence_model_gc = 0.0;
    double te_sequence_model_entropy = 0.0;
    double te_sequence_model_tandem_fraction = 0.0;
    double te_sequence_model_low_complexity_fraction = 0.0;
    double te_sequence_model_jsd_k5 = 0.0;
    double te_sequence_model_jsd_k6 = 0.0;
    double te_sequence_model_k9_containment = 0.0;
    std::string te_annotation_confidence = "NA";
    std::string te_annotation_class = "NA";
    std::string te_annotation_order = "NA";
    std::string te_annotation_intervals = "NA";
    double te_annotation_residual_fraction = 0.0;
    double te_annotation_masked_fraction = 0.0;
    int32_t left_flank_align_len = 0;
    int32_t right_flank_align_len = 0;
    int32_t event_consensus_len = 0;
    std::string insert_seq;
    std::string te_qc = "NA";
    std::string final_qc = "NA";
    std::string best_explanation = "NA";
    std::string explanation_residual = "NA";
    std::string explanation_path = "NA";
    std::string te_structure_path = "NA";
    double te_structure_log_evidence = 0.0;
    double nonte_structure_log_evidence = 0.0;
    double artifact_structure_log_evidence = 0.0;
    double te_structure_path_confidence = 0.0;
    double polyA_posterior = 0.0;
    double transduction_posterior = 0.0;
    double te_core_coverage = 0.0;
    int32_t unexplained_high_complexity_bp = 0;
    double te_posterior = 0.0;
    double non_te_posterior = 0.0;
    double artifact_posterior = 0.0;
    double te_vs_artifact_log_odds = 0.0;
    double te_vs_non_te_log_odds = 0.0;
    std::string posterior_qc = "POSTERIOR_NOT_EVALUATED";
    std::string latent_mechanism = "NA";
    double family_activity_prior = 0.0;
    double lfdr = 1.0;
    double worst_case_lfdr = 1.0;
    std::string lfdr_qc = "LFDR_NOT_EVALUATED";
    double mechanistic_lower_log_bf_te_vs_artifact = 0.0;
    double mechanistic_lower_log_bf_te_vs_non_te = 0.0;
    double mechanistic_ref_conflict_signal = 0.0;
    double mechanistic_ambiguity_width = 0.0;
    std::string mechanistic_blocks = "NA";
    double robust_mechanistic_lfdr = 1.0;
    double robust_mechanistic_worst_case_lfdr = 1.0;
    std::string robust_mechanistic_qc = "TE_LFDR_HIGH";
    double conformal_null_p = 1.0;
    double conformal_by_threshold = 0.0;
    int32_t conformal_dominated_nulls = 0;
    int32_t conformal_null_count = 0;
    std::string conformal_qc = "CONFORMAL_NOT_EVALUATED";
};

struct EvidenceLedgerRow {
    std::string chrom;
    int32_t tid = -1;
    int32_t pos = -1;
    int32_t bp_left = -1;
    int32_t bp_right = -1;
    int32_t coverage_left = -1;
    int32_t coverage_right = -1;
    int32_t owner_context_left = -1;
    int32_t owner_context_right = -1;
    std::string family = "NA";
    std::string subfamily = "NA";
    bool family_alignment_resolved = false;
    std::string final_qc = "NA";
    std::string posterior_qc = "NA";
    std::string lfdr_qc = "NA";
    std::string candidate_retention_reason = "NA";
    int32_t alt_struct_reads = 0;
    int32_t alt_split_reads = 0;
    int32_t alt_indel_reads = 0;
    int32_t alt_left_clip_reads = 0;
    int32_t alt_right_clip_reads = 0;
    int32_t full_context_input_reads = 0;
    int32_t raw_cigar_insert_reads = 0;
    int32_t max_raw_cigar_insert_len = 0;
    int32_t partial_context_input_reads = 0;
    int32_t left_anchor_input_reads = 0;
    int32_t right_anchor_input_reads = 0;
    int32_t input_event_reads = 0;
    int32_t event_consensus_len = 0;
    int32_t left_flank_align_len = 0;
    int32_t right_flank_align_len = 0;
    std::string insert_seq;
    int32_t ref_span_reads = 0;
    std::vector<std::string> support_qnames;
    double best_te_identity = 0.0;
    double best_te_query_coverage = 0.0;
    double cross_family_margin = 0.0;
    std::string te_structure_path = "NA";
    double te_structure_log_evidence = 0.0;
    double nonte_structure_log_evidence = 0.0;
    double artifact_structure_log_evidence = 0.0;
    double te_structure_path_confidence = 0.0;
    double polyA_posterior = 0.0;
    double transduction_posterior = 0.0;
    double te_core_coverage = 0.0;
    int32_t unexplained_high_complexity_bp = 0;
    double te_posterior = 0.0;
    double non_te_posterior = 0.0;
    double artifact_posterior = 0.0;
    double lfdr = 1.0;
    double worst_case_lfdr = 1.0;
    double mechanistic_lower_log_bf_te_vs_artifact = 0.0;
    double mechanistic_lower_log_bf_te_vs_non_te = 0.0;
    double mechanistic_ref_conflict_signal = 0.0;
    double mechanistic_ambiguity_width = 0.0;
    std::string mechanistic_blocks = "NA";
    double robust_mechanistic_lfdr = 1.0;
    double robust_mechanistic_worst_case_lfdr = 1.0;
    std::string robust_mechanistic_qc = "TE_LFDR_HIGH";
    double conformal_null_p = 1.0;
    double conformal_by_threshold = 0.0;
    int32_t conformal_dominated_nulls = 0;
    int32_t conformal_null_count = 0;
    std::string conformal_qc = "CONFORMAL_NOT_EVALUATED";
};

enum class FinalReportMode {
    Legacy,
    TeCalibrated,
};

struct PipelineConfig {
    std::string bam_path;
    std::string reference_fasta_path;
    std::string te_fasta_path;
    std::string te_blastn_path = PLACER_BLASTN_EXECUTABLE;
    std::string te_makeblastdb_path = PLACER_MAKEBLASTDB_EXECUTABLE;
    BamRegionScope bam_region_scope;

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
    int32_t te_family_topn = 4;
    int32_t te_family_representatives = 3;
    int32_t te_template_refine_topn = 3;
    int32_t te_exact_align_topn = 2;
    double te_family_margin_min = 0.05;
    double te_subfamily_margin_min = 0.04;
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
    // Beta-binomial genotype overdispersion (intra-class correlation). Default is
    // the binomial-ish baseline; finalization re-estimates it per sample from the
    // run's own count distribution (see PipelineResult::estimated_overdispersion).
    double genotype_overdispersion = 0.02;

    // Final sample-local conformal FDR target. This is a statistical risk
    // level, not an evidence weight.
    double final_fdr_q = 0.10;
    int32_t min_final_raw_cigar_insert_len_bp = 50;
    FinalReportMode final_report_mode = FinalReportMode::TeCalibrated;

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

struct FinalCallFilterConfig {
    int32_t min_raw_cigar_insert_len_bp = 0;
    FinalReportMode report_mode = FinalReportMode::Legacy;
};

struct EventSegmentationSearchStats {
    int64_t edit_distance_calls = 0;
    int64_t edit_distance_cache_hits = 0;
    int64_t edit_distance_cache_misses = 0;
    int64_t seed_bins_total = 0;
    int64_t paired_searches = 0;
    int64_t one_sided_searches = 0;
    int64_t endpoint_slack_searches = 0;
    int64_t reverse_complement_retries = 0;
    int64_t segmentation_cache_hits = 0;
    int64_t segmentation_cache_misses = 0;
};

struct PipelinePerformanceSummary {
    double pipeline_wall_seconds = 0.0;
    double pipeline_construction_seconds = 0.0;
    double bam_stream_wall_seconds = 0.0;
    double bin_processing_wall_seconds = 0.0;
    double finalization_seconds = 0.0;
    double gate1_preliminary_seconds = 0.0;
    double exact_bin_scan_seconds = 0.0;
    double component_build_seconds = 0.0;
    double local_fetch_seconds = 0.0;
    double projection_seconds = 0.0;
    double local_component_refresh_seconds = 0.0;
    double fragment_extract_seconds = 0.0;
    double signal_cache_build_seconds = 0.0;
    double hypothesis_summary_seconds = 0.0;
    double validator_input_count_seconds = 0.0;
    double event_consensus_seconds = 0.0;
    double event_segmentation_seconds = 0.0;
    double te_alignment_seconds = 0.0;
    int64_t te_blast_batches = 0;
    int64_t te_blast_queries = 0;
    int64_t te_blast_cache_hits = 0;
    int64_t te_blast_cache_misses = 0;
    int64_t te_blast_deduplicated_queries = 0;
    int64_t segmentation_edit_distance_calls = 0;
    int64_t segmentation_edit_distance_cache_hits = 0;
    int64_t segmentation_edit_distance_cache_misses = 0;
    int64_t segmentation_seed_bins_total = 0;
    int64_t segmentation_paired_searches = 0;
    int64_t segmentation_one_sided_searches = 0;
    int64_t segmentation_endpoint_slack_searches = 0;
    int64_t segmentation_reverse_complement_retries = 0;
    int64_t segmentation_cache_hits = 0;
    int64_t segmentation_cache_misses = 0;
};

struct TEAlignmentBatchStats {
    int64_t blast_batches = 0;
    int64_t blast_queries = 0;
    int64_t cache_hits = 0;
    int64_t cache_misses = 0;
    int64_t deduplicated_queries = 0;
};

struct PipelineResult {
    int64_t total_reads = 0;
    int64_t gate1_passed = 0;
    int64_t processed_bins = 0;
    int64_t built_components = 0;
    int64_t event_consensus_calls = 0;
    int64_t genotype_calls = 0;
    int64_t final_pass_calls = 0;
    // Per-sample beta-binomial overdispersion (intra-class correlation)
    // re-estimated at finalization from the run's own alt/ref count distribution.
    double estimated_overdispersion = 0.02;
    PipelinePerformanceSummary performance;

    std::vector<FinalCall> final_calls;
    std::vector<EvidenceLedgerRow> evidence_ledger;
};

struct ComponentFinalCallCandidate {
    int32_t pos = -1;
    int32_t anchor_pos = -1;
    double score = 0.0;
    bool emit_te = false;
    bool evidence_te = false;
    bool resolved_te = false;
    bool one_sided_segmentation = false;
    int32_t anchor_support = 0;
    int32_t anchor_ref_span_reads = 0;
    int32_t anchor_priority = std::numeric_limits<int32_t>::max();
    double anchor_hypothesis_score = 0.0;
    int32_t component_anchor_pos = -1;
    int32_t retether_direction = 0;
};

void retether_evidence_supported_final_call_positions(
    std::vector<ComponentFinalCallCandidate>& candidates);

std::vector<size_t> select_component_final_call_indices(
    const std::vector<ComponentFinalCallCandidate>& candidates);

void finalize_final_calls(PipelineResult& result);
void finalize_final_calls(PipelineResult& result, double target_fdr);
void finalize_final_calls(
    PipelineResult& result,
    double target_fdr,
    const FinalCallFilterConfig& filter_config);
void finalize_final_calls(PipelineResult& result, const TSDDetector& reference);
void finalize_final_calls(PipelineResult& result, const TSDDetector& reference, double target_fdr);
void finalize_final_calls(
    PipelineResult& result,
    const TSDDetector& reference,
    double target_fdr,
    const FinalCallFilterConfig& filter_config);

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
    std::vector<TEAlignmentEvidence> align_insert_sequences(
        const std::vector<std::string>& insert_seqs) const;
    TEAlignmentBatchStats last_alignment_batch_stats() const;

private:
    PipelineConfig config_;
    struct Index;
    struct AlignmentShortlistDb;
    struct TemplateSeedDb;
    struct AlignmentEvidenceCache;
    std::shared_ptr<const Index> primary_index_;
    std::vector<std::shared_ptr<const Index>> indices_;
    std::shared_ptr<const AlignmentShortlistDb> alignment_shortlist_db_;
    std::shared_ptr<const TemplateSeedDb> template_seed_db_;
    std::shared_ptr<const TeFamilyAlignmentIndex> family_alignment_index_;
    std::shared_ptr<const TeFamilyGroupCache> family_rep_groups_;
    std::shared_ptr<const std::vector<std::string>> te_names_;
    std::shared_ptr<const std::vector<std::string>> te_sequences_;
    std::shared_ptr<const std::vector<std::string>> te_reverse_complement_sequences_;
    std::shared_ptr<const std::string> blast_db_prefix_;
    std::shared_ptr<AlignmentEvidenceCache> alignment_evidence_cache_;
    mutable TEAlignmentBatchStats last_alignment_batch_stats_;
    mutable int32_t last_exact_alignments_ = 0;
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
    bool reference_position_is_poly_n(
        const std::string& chrom,
        int32_t pos) const;
    bool reference_interval_is_n_rich(
        const std::string& chrom,
        int32_t start,
        int32_t end) const;
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
    EventSegmentationSearchStats last_segmentation_stats() const;

private:
    struct HypothesisSummary {
        size_t original_index = 0;
        EventReadEvidence event_evidence;
        int32_t bp_left = -1;
        int32_t bp_right = -1;
        int32_t alt_split_reads = 0;
        int32_t alt_indel_reads = 0;
        int32_t alt_left_clip_reads = 0;
        int32_t alt_right_clip_reads = 0;
        int32_t alt_struct_reads = 0;
        int32_t raw_cigar_insert_reads = 0;
        int32_t max_raw_cigar_insert_len = 0;
        int32_t ref_span_reads = 0;
        int32_t inferred_event_length = 0;
        int32_t hypothesis_support = 0;
        int32_t hypothesis_priority = std::numeric_limits<int32_t>::max();
        double hypothesis_score = 0.0;
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

    enum class ConsensusContextMode {
        kPreferFull,
        kPartialOnly,
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
        std::vector<int32_t> split_positions;
        std::vector<int32_t> indel_positions;
        std::vector<std::pair<int32_t, int32_t>> raw_cigar_insertions;
    };

    struct ComponentSignalCache {
        std::vector<CachedLocalSignal> reads;
    };

    struct HypothesisValidatorEvidence {
        HypothesisSummary summary;
        EventReadEvidence event_evidence;
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
        int32_t bp_right,
        int32_t hypothesis_support,
        int32_t hypothesis_priority) const;

    std::vector<HypothesisSummary> select_hypothesis_summaries_for_expensive_stage(
        const std::vector<HypothesisSummary>& summaries) const;

    std::vector<HypothesisSummary> collapse_hypothesis_summaries(
        const std::vector<HypothesisSummary>& summaries) const;

    bool should_keep_hypothesis_for_expensive_stage(
        const HypothesisSummary& summary,
        bool is_top_ranked_survivor) const;

    bool should_record_hypothesis_in_evidence_ledger(
        const HypothesisSummary& summary) const;

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
        size_t top_k,
        int32_t anchor_pos) const;

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
        const std::vector<ReadReferenceSpan>& read_spans,
        const std::string& chrom) const;

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

    bool has_identifiable_bilateral_partial_context_support(
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
        const EventReadEvidence& event_evidence,
        ConsensusContextMode mode = ConsensusContextMode::kPreferFull) const;

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

    static bool edit_identity_if_at_least_fast(
        std::string_view lhs,
        std::string_view rhs,
        int32_t max_edits,
        double& identity_out);

    EventSegmentation segment_event_consensus_cached(
        const ComponentCall& component,
        const EventReadEvidence& event_evidence,
        const EventConsensus& event_consensus,
        std::unordered_map<std::string, EventSegmentation>& segmentation_cache) const;

    EventSegmentationEvidence analyze_event_segmentation(
        const EventConsensus& event_consensus,
        const EventSegmentation& event_segmentation) const;

    TEAlignmentEvidence align_insert_seq_to_te(
        const EventSegmentation& event_segmentation) const;

    TEAlignmentEvidence align_insert_seq_to_te_cached(
        const EventSegmentation& event_segmentation,
        std::unordered_map<std::string, TEAlignmentEvidence>& te_alignment_cache) const;

    ClipInsertConcordanceEvidence analyze_clip_insert_concordance(
        const EventReadEvidence& event_evidence,
        const EventSegmentation& event_segmentation,
        const std::vector<InsertionFragment>& fragments) const;

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
        const JointDecisionResult& joint) const;

    PipelineConfig config_;
    std::unique_ptr<BamStreamReader> bam_reader_;
    SignalFirstGate1Module gate1_module_;
    LinearBinComponentModule component_module_;
    CigarInsertionFragmentModule ins_fragment_module_;
    TEKmerQuickClassifierModule te_classifier_module_;
    TSDDetector tsd_detector_;
    mutable EventSegmentationSearchStats last_segmentation_stats_;
};

std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config);

}  // namespace placer

#endif  // PLACER_PIPELINE_H
