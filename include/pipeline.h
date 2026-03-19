#ifndef PLACER_PIPELINE_H
#define PLACER_PIPELINE_H

#include "bam_io.h"
#include "gate1_module.h"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <string>
#include <vector>

namespace placer {

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
    double ambiguous_frac = 0.0;
    int32_t evidence_soft_clip_count = 0;
    int32_t evidence_indel_count = 0;
    int32_t evidence_sa_hint_count = 0;

    std::vector<size_t> read_indices;
    std::vector<size_t> soft_clip_read_indices;
    std::vector<size_t> split_sa_read_indices;
    std::vector<size_t> insertion_read_indices;
    std::vector<BreakpointCandidate> breakpoint_candidates;
};

struct EvidenceFeatures {
    int32_t tid = -1;
    int32_t pos = -1;

    int32_t global_cov_reads = 0;
    int32_t local_cov_reads = 0;
    int32_t depth_reads = 0;
    int32_t alt_support_reads = 0;
    int32_t ref_support_reads = 0;
    double effective_alt_support = 0.0;
    int32_t evidence_point_count = 0;

    double breakpoint_mad = 0.0;
    double split_breakpoint_mad = 0.0;
    double split_like_breakpoint_mad = 0.0;
    double clip_breakpoint_mad = 0.0;
    double indel_breakpoint_mad = 0.0;
    double split_clip_core_delta = 0.0;
    double breakpoint_inconsistency_penalty = 0.0;
    int32_t bp_source_split_count = 0;
    int32_t bp_source_clip_count = 0;
    int32_t bp_source_indel_count = 0;
    bool bp_fallback_used = false;
    double te_vote_fraction = 0.0;
    double te_median_identity = 0.0;
    double ambiguous_frac = 0.0;
    double low_complex_softclip_frac = 0.0;
    bool anchor_consensus_ok = false;

    int32_t min_support_required = 0;
    bool pass_min_support = true;
    bool pass_breakpoint_mad = true;
    bool pass_split_clip_consistency = true;
    bool pass_low_complexity = true;
    bool pass_te_consistency = true;
    bool hard_filtered = false;
};

struct AssemblyCall {
    int32_t tid = -1;
    int32_t pos = -1;
    std::string consensus;
    std::string assembly_mode = "NONE";
    int32_t input_fragments = 0;
    int32_t used_fragments = 0;
    int32_t consensus_len = 0;
    double identity_est = 0.0;
    bool qc_pass = false;
    std::string qc = "NO_ASSEMBLY";
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

enum class InsertionTheta : uint8_t {
    kUnknown = 0,
    kFwd = 1,
    kRev = 2
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

    // Anchor/split diagnostics used by TE consensus.
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

struct ClusterTECall {
    std::string te_name;
    double vote_fraction = 0.0;
    double median_identity = 0.0;  // median fragment k-mer support
    double multik_support = 0.0;
    double rescue_frac = 0.0;
    int32_t fragment_count = 0;
    std::string top1_te_name;
    std::string top2_te_name;
    double posterior_top1 = 0.0;    // proxy from vote fractions
    double posterior_top2 = 0.0;    // proxy from vote fractions
    double posterior_margin = 0.0;  // posterior_top1 - posterior_top2
    bool passed = false;
};

struct AnchorLockedReport {
    bool enabled = false;
    bool has_result = false;

    std::string te_name;
    InsertionTheta theta0 = InsertionTheta::kUnknown;
    bool fail_theta_uncertain = false;

    double mad_fwd = 0.0;
    double mad_rev = 0.0;
    double sum_w_fwd = 0.0;
    double sum_w_rev = 0.0;

    double center_c0 = 0.0;
    double center_c1 = 0.0;
    int32_t te_breakpoint_core = -1;
    int32_t te_breakpoint_window_start = -1;
    int32_t te_breakpoint_window_end = -1;
    double core_span = 0.0;

    int32_t total_fragments = 0;
    int32_t fragments_with_placements = 0;
    int32_t core_candidate_count = 0;
    int32_t core_set_count = 0;
    int32_t split_sa_core_count = 0;
    double split_sa_core_frac = 0.0;

    int32_t ref_junc_pos_min = -1;
    int32_t ref_junc_pos_max = -1;
};

struct PlaceabilityReport {
    int32_t tid = -1;
    int32_t pos = -1;
    int32_t depth_reads = 0;
    int32_t support_reads = 0;
    int32_t ref_support_reads = 0;
    double delta_score = 0.0;
    double evidence_p = 0.0;
    double evidence_q = 0.0;
    int32_t min_support_required = 0;
    double breakpoint_mad = 0.0;
    double low_complex_softclip_frac = 0.0;
    bool hard_filtered = false;
    int32_t tier = 3;
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

struct SequenceClosureEvidence {
    bool enabled = false;
    bool pass = false;
    bool certain = false;
    bool force_non_te = false;
    int32_t left_anchor_reads = 0;
    int32_t right_anchor_reads = 0;
    int32_t dual_anchor_reads = 0;
    int32_t total_anchor_reads = 0;
    int32_t empty_span_reads = 0;
    int32_t split_like_support_reads = 0;
    double mean_anchor_identity = 0.0;
    std::string qc = "SEQ_CLOSURE_DISABLED";
};

struct FinalCall {
    std::string chrom;
    int32_t tid = -1;
    int32_t pos = -1;

    int32_t window_start = -1;
    int32_t window_end = -1;

    std::string te_name;
    double te_vote_fraction = 0.0;
    double te_median_identity = 0.0;
    double te_multik_support = 0.0;
    double te_rescue_frac = 0.0;
    int32_t te_fragment_count = 0;
    std::string te_theta = "NA";
    double te_mad_fwd = 0.0;
    double te_mad_rev = 0.0;
    int32_t te_breakpoint_core = -1;
    int32_t te_breakpoint_window_start = -1;
    int32_t te_breakpoint_window_end = -1;
    int32_t te_core_candidates = 0;
    int32_t te_core_set = 0;
    double te_split_sa_core_frac = 0.0;
    int32_t te_ref_junc_pos_min = -1;
    int32_t te_ref_junc_pos_max = -1;
    std::string bp_source_counts = "split:0,clip:0,indel:0";
    bool bp_fallback_used = false;
    std::string insertion_qc = "NA";
    std::string te_qc = "NA";
    std::string tsd_type = "NONE";
    int32_t tsd_len = 0;
    std::string tsd_seq = "NA";
    double tsd_bg_p = 1.0;
    std::string te_status = "NON_TE";
    std::string te_top1_name;
    std::string te_top2_name;
    double te_posterior_top1 = 0.0;
    double te_posterior_top2 = 0.0;
    double te_posterior_margin = 0.0;
    double te_confidence_prob = 0.0;
    std::string confidence = "HIGH";

    int32_t tier = 3;
    int32_t support_reads = 0;
    int32_t softclip_support_reads = 0;
    int32_t split_sa_support_reads = 0;
    int32_t indel_support_reads = 0;
    int32_t split_like_support_reads = 0;
    int32_t support_low_mapq_reads = 0;
    double support_low_mapq_frac = 0.0;
    bool seq_closure_enabled = false;
    bool seq_closure_pass = false;
    bool seq_closure_certain = false;
    bool seq_closure_force_non_te = false;
    int32_t seq_closure_left_anchor_reads = 0;
    int32_t seq_closure_right_anchor_reads = 0;
    int32_t seq_closure_dual_anchor_reads = 0;
    int32_t seq_closure_total_anchor_reads = 0;
    int32_t seq_closure_empty_span_reads = 0;
    double seq_closure_mean_anchor_identity = 0.0;
    std::string seq_closure_qc = "SEQ_CLOSURE_DISABLED";
    std::string genotype = "./.";
    double af = 0.0;
    int32_t gq = 0;
    std::string asm_mode = "NONE";
    int32_t asm_input_fragments = 0;
    int32_t asm_used_fragments = 0;
    int32_t asm_consensus_len = 0;
    double asm_identity_est = 0.0;
    std::string asm_qc = "NO_ASSEMBLY";
};

struct PipelineConfig {
    std::string bam_path;
    std::string reference_fasta_path;
    std::string te_fasta_path;

    int32_t bam_threads = 2;
    int64_t progress_interval = 100000;

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
    double te_vote_fraction_min = 0.40;
    double te_median_identity_min = 0.30;
    int32_t te_min_fragments_for_vote = 2;
    double te_rescue_vote_fraction_min = 0.25;
    double te_rescue_median_identity_min = 0.20;
    bool te_low_kmer_rescue_enable = true;
    int32_t te_low_kmer_rescue_topn = 3;
    int32_t te_low_kmer_rescue_min_frag_len = 40;
    double te_low_kmer_rescue_identity_min = 0.55;
    double te_low_kmer_rescue_margin_max = 0.08;
    int32_t te_no_softclip_min_reads = 2;
    int32_t te_no_softclip_min_fragments = 2;
    double te_no_softclip_identity_min = 0.20;
    double te_one_sided_breakpoint_mad_max = 40.0;
    int32_t te_mixed_min_non_softclip_reads = 4;
    double te_mixed_min_non_softclip_frac = 0.35;
    double te_mixed_clip_dominant_ratio = 2.0;
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
    int32_t te_pure_softclip_min_reads = 6;
    int32_t te_pure_softclip_min_fragments = 6;
    double te_pure_softclip_min_identity = 0.35;
    // Open-set TE status gating using consensus proxy posterior.
    double te_proxy_posterior_top1_min = 0.90;
    double te_proxy_posterior_margin_min = 0.50;
    double te_proxy_identity_min = 0.60;
    // Proxy posterior calibration to continuous confidence probability.
    double te_confidence_bias = -3.0;
    double te_confidence_w_top1 = 2.4;
    double te_confidence_w_margin = 2.0;
    double te_confidence_w_asm_identity = 1.8;
    double te_confidence_w_support = 0.8;
    double te_confidence_w_vote = 1.0;
    double te_confidence_w_te_identity = 1.0;
    double te_confidence_w_breakpoint_mad = 0.8;
    double te_confidence_rescue_penalty = 0.35;
    double te_confidence_prob_certain_min = 0.85;
    double te_confidence_prob_uncertain_min = 0.35;
    double te_certain_posterior_margin_min = 0.25;
    double te_same_family_ambiguity_margin_max = 0.30;
    double te_confidence_anchor_fail_penalty = 0.18;
    double te_confidence_no_tsd_penalty = 0.10;
    double te_confidence_low_margin_penalty = 0.16;
    bool te_force_non_te_on_combined_weakness = true;
    bool te_force_non_te_on_anchor_weak_tsd = true;
    bool te_fail_on_tsd_inconsistent = false;
    bool te_sequence_closure_enable = true;
    int32_t te_sequence_closure_flank_bases = 48;
    int32_t te_sequence_closure_min_anchor_len = 24;
    double te_sequence_closure_min_anchor_identity = 0.72;
    int32_t te_sequence_closure_min_side_reads = 1;
    int32_t te_sequence_closure_min_total_reads = 2;
    int32_t te_sequence_closure_min_dual_reads = 1;
    int32_t te_sequence_closure_split_like_rescue_min_reads = 3;
    int32_t te_sequence_closure_empty_window = 25;
    double te_sequence_closure_max_empty_ratio_pass = 3.0;
    double te_sequence_closure_max_empty_ratio_certain = 1.5;
    // Optional low-confidence fallback acceptance for exploratory analysis.
    bool emit_low_confidence_calls = false;
    int32_t low_conf_min_support_reads = 2;
    int32_t low_conf_max_tier = 2;
    // Pass-1 bootstrap export for low-confidence/unknown TE calls.
    bool bootstrap_export_enable = false;
    bool bootstrap_export_include_non_te = true;
    int32_t bootstrap_export_min_consensus_len = 80;
    std::string bootstrap_consensus_fasta_path = "pass1_bootstrap_consensus.fasta";
    std::string bootstrap_metadata_tsv_path = "pass1_bootstrap_calls.tsv";

    // Module 2.3: deterministic TE consensus (anchor-locked + theta).
    bool te_consensus_enable = true;
    int32_t te_consensus_seed_k = 13;
    int32_t te_consensus_max_start_candidates = 5;
    double te_consensus_delta_tpl_min = 0.05;
    int32_t te_consensus_anchor_min = 80;
    double te_consensus_start_span_sqrt_scale = 2.0;
    double te_consensus_start_span_bias = 6.0;
    int32_t te_consensus_w_side_max = 300;
    double te_consensus_temp_t = 0.25;
    double te_consensus_beta = 0.20;
    double te_consensus_lambda = 0.35;
    double te_consensus_mu = 0.10;
    double te_consensus_w_core_gate = 80.0;
    double te_consensus_w_min = 25.0;
    double te_consensus_gamma = 2.0;
    double te_consensus_eps_theta = 0.05;
    double te_consensus_rel_sumw_eps = 0.10;
    int32_t te_consensus_min_core_set = 2;

    // Module 2.4: TSD detector (duplication / deletion).
    bool tsd_enable = true;
    int32_t tsd_min_len = 3;
    int32_t tsd_max_len = 50;
    int32_t tsd_flank_window = 150;
    double tsd_bg_p_max = 0.05;

    // Evidence scoring and hard filters.
    double evidence_min_support_alpha = 0.08;
    double evidence_min_support_lambda = 0.75;
    double evidence_breakpoint_mad_max = 80.0;
    double evidence_clip_breakpoint_mad_max = 120.0;
    double evidence_split_clip_core_delta_max = 150.0;
    double evidence_low_complex_softclip_frac_max = 0.80;
    double evidence_tier1_prob = 0.90;
    double evidence_tier2_prob = 0.60;
    double evidence_logit_bias = 3.0;

    // Genotype likelihood model.
    int32_t genotype_min_depth = 3;
    double genotype_error_rate = 0.02;

    // M2: local assembly (abPOA only).
    int32_t assembly_poa_min_reads = 2;
    int32_t assembly_poa_max_reads = 48;
    int32_t assembly_min_fragment_len = 80;
    int32_t assembly_min_consensus_len = 80;
    double assembly_min_identity_est = 0.55;
    int32_t assembly_kmer_size = 11;

    bool enable_parallel = false;
    size_t batch_size = 1000;
    int32_t parallel_workers = 0;
    int32_t parallel_queue_max_tasks = 0;  // <=0 means unbounded queue
};

struct PipelineResult {
    int64_t total_reads = 0;
    int64_t gate1_passed = 0;
    int64_t processed_bins = 0;
    int64_t built_components = 0;
    int64_t evidence_rows = 0;
    int64_t assembled_calls = 0;
    int64_t placeability_calls = 0;
    int64_t genotype_calls = 0;
    int64_t final_te_certain = 0;
    int64_t final_te_uncertain = 0;
    int64_t final_non_te = 0;
    int64_t final_high_confidence = 0;
    int64_t final_low_confidence = 0;
    int64_t bootstrap_exported_calls = 0;

    std::vector<FinalCall> final_calls;
};

void finalize_final_calls(PipelineResult& result);

struct TeFinalEvidenceDecision {
    bool force_te_uncertain = false;
    bool force_non_te = false;
};

TeFinalEvidenceDecision apply_te_evidence_gates(
    FinalCall& call,
    const PipelineConfig& config);

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
    ClusterTECall vote_cluster(
        const std::vector<FragmentTEHit>& hits) const;

private:
    PipelineConfig config_;
    struct Index;
    std::shared_ptr<const Index> primary_index_;
    std::vector<std::shared_ptr<const Index>> indices_;
    std::vector<std::string> te_names_;
    std::vector<std::string> te_sequences_;
};

class DeterministicAnchorLockedModule final {
public:
    explicit DeterministicAnchorLockedModule(PipelineConfig config);

    bool is_enabled() const;
    AnchorLockedReport resolve(
        const ComponentCall& component,
        const std::vector<InsertionFragment>& fragments,
        const std::vector<FragmentTEHit>& hits,
        const ClusterTECall& te_call) const;

private:
    PipelineConfig config_;
    struct TemplateDb;
    std::shared_ptr<const TemplateDb> template_db_;
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

    PipelineResult run_streaming() const;
    PipelineResult run_parallel() const;

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
        PipelineResult& result) const;

    EvidenceFeatures collect_evidence(
        const ComponentCall& component,
        const std::vector<const bam1_t*>& bin_records,
        const std::vector<ReadReferenceSpan>& read_spans,
        const std::vector<InsertionFragment>& fragments,
        const ClusterTECall& te_call,
        const AnchorLockedReport& anchor_report) const;

    AssemblyCall assemble_component(
        const ComponentCall& component,
        const std::vector<InsertionFragment>& fragments,
        const std::vector<FragmentTEHit>& hits,
        const ClusterTECall& te_call) const;

    SequenceClosureEvidence build_sequence_closure_evidence(
        const ComponentCall& component,
        const std::vector<const bam1_t*>& bin_records,
        const std::vector<ReadReferenceSpan>& read_spans,
        const std::vector<InsertionFragment>& fragments) const;

    PlaceabilityReport score_placeability(
        const AssemblyCall& assembly,
        const EvidenceFeatures& evidence) const;

    GenotypeCall genotype_call(
        const AssemblyCall& assembly,
        const PlaceabilityReport& placeability) const;

    PipelineConfig config_;
    std::unique_ptr<BamStreamReader> bam_reader_;
    SignalFirstGate1Module gate1_module_;
    LinearBinComponentModule component_module_;
    CigarInsertionFragmentModule ins_fragment_module_;
    TEKmerQuickClassifierModule te_classifier_module_;
    DeterministicAnchorLockedModule anchor_locked_module_;
    TSDDetector tsd_detector_;
};

std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config);

}  // namespace placer

#endif  // PLACER_PIPELINE_H
