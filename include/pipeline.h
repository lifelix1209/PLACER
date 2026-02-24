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
    double te_vote_fraction = 0.0;
    double te_median_identity = 0.0;
    double ambiguous_frac = 0.0;
    double low_complex_softclip_frac = 0.0;
    bool anchor_consensus_ok = false;

    int32_t min_support_required = 0;
    bool pass_min_support = true;
    bool pass_breakpoint_mad = true;
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

    int32_t hit_kmers = 0;
    int32_t total_kmers = 0;
};

struct ClusterTECall {
    std::string te_name;
    double vote_fraction = 0.0;
    double median_identity = 0.0;  // median fragment k-mer support
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

struct FinalCall {
    std::string chrom;
    int32_t tid = -1;
    int32_t pos = -1;

    int32_t window_start = -1;
    int32_t window_end = -1;

    std::string te_name;
    double te_vote_fraction = 0.0;
    double te_median_identity = 0.0;
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
    std::string te_qc = "NA";
    std::string te_status = "NON_TE";
    std::string te_top1_name;
    std::string te_top2_name;
    double te_posterior_top1 = 0.0;
    double te_posterior_top2 = 0.0;
    double te_posterior_margin = 0.0;

    int32_t tier = 3;
    int32_t support_reads = 0;
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
    std::string ins_fragments_fasta_path = "ins_fragments.fasta";
    int32_t min_soft_clip_for_seq_extract = 50;
    int32_t min_long_ins_for_seq_extract = 50;
    int32_t min_sa_aln_len_for_seq_extract = 50;
    int32_t max_sa_per_read = 3;

    // Module 2.2: quick TE classification (optional, requires te_fasta_path).
    std::string ins_fragment_hits_tsv_path = "ins_fragment_hits.tsv";
    int32_t te_kmer_size = 13;
    double te_vote_fraction_min = 0.40;
    double te_median_identity_min = 0.30;
    int32_t te_min_fragments_for_vote = 2;
    double te_rescue_vote_fraction_min = 0.25;
    double te_rescue_median_identity_min = 0.20;
    double te_softclip_low_complexity_at_frac_min = 0.90;
    int32_t te_softclip_low_complexity_homopolymer_min = 80;
    int32_t te_pure_softclip_min_reads = 6;
    int32_t te_pure_softclip_min_fragments = 6;
    double te_pure_softclip_min_identity = 0.35;

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

    // Evidence scoring and hard filters.
    double evidence_min_support_alpha = 0.08;
    double evidence_min_support_lambda = 0.75;
    double evidence_breakpoint_mad_max = 80.0;
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

    std::vector<FinalCall> final_calls;
};

class LinearBinComponentModule final {
public:
    std::vector<ComponentCall> build(
        const std::vector<BamRecordPtr>& bin_records,
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
        const std::vector<BamRecordPtr>& bin_records) const;

private:
    PipelineConfig config_;
};

class SplitSAFragmentModule final {
public:
    explicit SplitSAFragmentModule(PipelineConfig config);

    std::vector<InsertionFragment> extract(
        const ComponentCall& component,
        const std::vector<BamRecordPtr>& bin_records) const;

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
    std::shared_ptr<const Index> index_;
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
        std::deque<WindowCoord> active_window;
        std::vector<BamRecordPtr> current_bin_records;
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
        std::vector<BamRecordPtr>&& bin_records,
        int32_t tid,
        int32_t bin_index,
        PipelineResult& result) const;

    EvidenceFeatures collect_evidence(
        const ComponentCall& component,
        const std::vector<BamRecordPtr>& bin_records,
        const std::vector<InsertionFragment>& fragments,
        const ClusterTECall& te_call,
        const AnchorLockedReport& anchor_report) const;

    AssemblyCall assemble_component(
        const ComponentCall& component,
        const std::vector<InsertionFragment>& fragments,
        const std::vector<FragmentTEHit>& hits,
        const ClusterTECall& te_call) const;

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
};

std::unique_ptr<Pipeline> build_default_pipeline(const PipelineConfig& config);

}  // namespace placer

#endif  // PLACER_PIPELINE_H
