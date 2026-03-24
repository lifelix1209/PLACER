#ifndef PLACER_DECISION_POLICY_H
#define PLACER_DECISION_POLICY_H

#include <cstdint>
#include <string>

namespace placer {

struct EventGenotypeInput {
    int32_t alt_struct_reads = 0;
    int32_t ref_span_reads = 0;
    int32_t min_depth = 3;
    int32_t min_gq = 20;
    double error_rate = 0.02;
};

struct EventGenotypeDecision {
    std::string best_gt = "./.";
    double allele_fraction = 0.0;
    int32_t gq = 0;
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
