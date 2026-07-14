#ifndef PLACER_DECISION_POLICY_H
#define PLACER_DECISION_POLICY_H

#include "event_explanation.h"

#include <cstdint>
#include <string>
#include <vector>

namespace placer {

struct EventSegmentation;
struct TEAlignmentEvidence;

struct EventGenotypeInput {
    int32_t alt_struct_reads = 0;
    int32_t alt_split_reads = -1;
    int32_t alt_indel_reads = -1;
    int32_t alt_left_clip_reads = -1;
    int32_t alt_right_clip_reads = -1;
    int32_t ref_span_reads = 0;
    int32_t min_depth = 3;
    int32_t min_gq = 20;
    double error_rate = 0.02;
    // Beta-binomial overdispersion: the intra-class correlation rho in [0,1) of
    // the alt/ref read counts. rho -> 0 recovers the independent-read binomial
    // limit; larger rho models the count overdispersion that mapping bias and
    // local alignment ambiguity produce in repeats (a per-sample/context nuisance
    // the pipeline can estimate genome-wide and inject here).
    double overdispersion = 0.02;
    int32_t event_length = 0;
    std::vector<int32_t> alt_observed_lengths;
};

struct EventGenotypeDecision {
    std::string best_gt = "./.";
    double allele_fraction = 0.0;
    int32_t gq = 0;
    int32_t depth = 0;
    double best_nonref_minus_ref_ll = 0.0;
    bool pass = false;
};

EventGenotypeDecision genotype_event_from_alt_vs_ref(
    const EventGenotypeInput& input);

// A single (alt, depth) count observation used to estimate the sample-level
// beta-binomial overdispersion.
struct AltDepthObservation {
    int32_t alt = 0;
    int32_t depth = 0;
};

// Estimate the sample-level beta-binomial overdispersion (intra-class
// correlation rho in [0,1)) of alt/ref read counts by moments (a Fleiss-style
// ANOVA estimator over heterozygous-like sites). Returns `fallback` when there
// is too little information (few informative sites) to estimate it. This makes
// the genotyper's overdispersion a per-sample, data-driven quantity instead of a
// fixed default; rho -> 0 is the independent-read binomial limit.
double estimate_alt_depth_overdispersion(
    const std::vector<AltDepthObservation>& observations,
    double fallback = 0.02);

struct FinalBoundaryInput {
    int32_t left_ref_start = -1;
    int32_t left_ref_end = -1;
    int32_t right_ref_start = -1;
    int32_t right_ref_end = -1;
    int32_t tsd_min_len = 3;
    int32_t tsd_max_len = 50;
};

struct FinalBoundaryDecision {
    bool pass = false;
    std::string boundary_type = "REJECT";
    int32_t boundary_len = 0;
    std::string qc = "REJECT_BOUNDARY_UNSET";
};

FinalBoundaryDecision check_boundary_consistency(
    const FinalBoundaryInput& input);

enum class FinalHypothesisKind : uint8_t {
    kReference = 0,
    kInsertionNonTe = 1,
    kTeUnknown = 2,
    kTeResolved = 3
};

struct EventExistenceEvidence {
    std::string best_gt = "./.";
    double af = 0.0;
    int32_t gq = 0;
    int32_t alt_struct_reads = 0;
    int32_t alt_split_reads = -1;
    int32_t alt_indel_reads = -1;
    int32_t alt_left_clip_reads = -1;
    int32_t alt_right_clip_reads = -1;
    int32_t ref_span_reads = 0;
    int32_t depth = 0;
    double best_nonref_minus_ref_ll = 0.0;
    double score = -3.0;
};

struct EventSegmentationEvidence {
    bool has_consensus = false;
    bool has_left_flank = false;
    bool has_right_flank = false;
    bool has_insert_seq = false;
    bool pair_valid = false;
    int32_t left_align_len = 0;
    int32_t right_align_len = 0;
    double left_identity = 0.0;
    double right_identity = 0.0;
    int32_t insert_len = 0;
    double score = -3.0;
    std::string qc = "NO_EVENT_SEGMENTATION";
};

struct BoundaryEvidence {
    bool geometry_defined = false;
    bool canonical_pass = false;
    bool evidence_consistent = false;
    std::string boundary_type = "REJECT";
    int32_t boundary_len = 0;
    double score = -3.0;
    std::string qc = "REJECT_BOUNDARY_UNSET";
};

struct ClipInsertConcordanceEvidence {
    bool pass = false;
    int32_t full_insert_reads = 0;
    int32_t left_clip_reads = 0;
    int32_t right_clip_reads = 0;
    double max_left_identity = 0.0;
    double max_right_identity = 0.0;
    std::string qc = "NO_CLIP_INSERT_CONCORDANCE";
};

struct JointHypothesisScore {
    FinalHypothesisKind kind = FinalHypothesisKind::kReference;
    double total = -1e9;
    double existence = 0.0;
    double segmentation = 0.0;
    double te = 0.0;
    double boundary = 0.0;
    bool hard_veto = false;
    std::string reason;
};

struct JointDecisionResult {
    JointHypothesisScore best;
    JointHypothesisScore runner_up;
    bool emit_te_call = false;
    bool emit_structural_event_call = false;
    bool emit_unknown_te = false;
    bool emit_evidence_te_call = false;
    std::string final_qc = "REJECT_EVENT_EXISTENCE";
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
    ExplanationDecision explanation_decision;
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
};

struct LatentMechanismEvidence {
    std::string latent_mechanism = "artifact_reference";
    double family_activity_prior = 0.0;
    double te_posterior = 0.0;
    double non_te_posterior = 0.0;
    double artifact_posterior = 0.0;
    double lfdr = 1.0;
    double worst_case_lfdr = 1.0;
    std::string lfdr_qc = "TE_LFDR_HIGH";
};

EventExistenceEvidence build_event_existence_evidence(
    const EventGenotypeInput& input);

JointDecisionResult evaluate_joint_hypotheses(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    const ClipInsertConcordanceEvidence* clip_insert_concordance = nullptr);

ExplanationDecision evaluate_event_explanations_for_test(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    const ClipInsertConcordanceEvidence* clip_insert_concordance = nullptr);

// Structural-sanity predicates for the TE-unknown (h2) and TE-resolved (h3)
// hypotheses in evaluate_joint_hypotheses. They mark a hypothesis ineligible for
// ranking when the event cannot definitionally be a TE call (no insert sequence,
// no usable TE alignment, annotation quality too low, resolved call without a
// resolved qc_reason). They are necessary conditions only -- final emission is
// gated by the calibrated worst-case local FDR, not by these predicates.
// Exposed for direct unit testing; return true when the hypothesis is vetoed.
bool compute_te_unknown_hard_veto(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary);

bool compute_te_resolved_hard_veto(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary);

LatentMechanismEvidence evaluate_latent_mechanism_lfdr(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    const ClipInsertConcordanceEvidence* clip_insert_concordance = nullptr,
    double target_q = 0.10);

EventSegmentationEvidence analyze_event_segmentation_for_test(
    bool has_consensus,
    const EventSegmentation& segmentation);

BoundaryEvidence evaluate_boundary_evidence(
    const FinalBoundaryInput& input,
    int32_t breakpoint_envelope_width);

}  // namespace placer

#endif  // PLACER_DECISION_POLICY_H
