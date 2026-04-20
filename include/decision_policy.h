#ifndef PLACER_DECISION_POLICY_H
#define PLACER_DECISION_POLICY_H

#include <cstdint>
#include <string>
#include <vector>

namespace placer {

struct EventSegmentation;
struct TEAlignmentEvidence;

struct EventGenotypeInput {
    int32_t alt_struct_reads = 0;
    int32_t ref_span_reads = 0;
    int32_t min_depth = 3;
    int32_t min_gq = 20;
    double error_rate = 0.02;
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
    bool emit_unknown_te = false;
    std::string final_qc = "REJECT_EVENT_EXISTENCE";
};

EventExistenceEvidence build_event_existence_evidence(
    const EventGenotypeInput& input);

JointDecisionResult evaluate_joint_hypotheses(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary);

bool should_emit_te_call(
    const JointHypothesisScore& best_te,
    const JointHypothesisScore& best_non_te);

EventSegmentationEvidence analyze_event_segmentation_for_test(
    bool has_consensus,
    const EventSegmentation& segmentation);

BoundaryEvidence evaluate_boundary_evidence(
    const FinalBoundaryInput& input,
    int32_t breakpoint_envelope_width);

struct FinalTeAcceptanceInput {
    bool event_existence_pass = false;
    bool event_closure_pass = false;
    bool te_sequence_pass = false;
    bool boundary_pass = false;
};

struct FinalTeAcceptanceDecision {
    bool pass = false;
    std::string qc = "REJECT_EVENT_EXISTENCE";
};

FinalTeAcceptanceDecision evaluate_final_te_acceptance(
    const FinalTeAcceptanceInput& input);

}  // namespace placer

#endif  // PLACER_DECISION_POLICY_H
