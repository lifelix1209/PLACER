#ifndef PLACER_DECISION_POLICY_H
#define PLACER_DECISION_POLICY_H

#include <cstdint>
#include <string>

namespace placer {

struct InsertionEvidence {
    int32_t tier = 3;
    int32_t support_reads = 0;
};

struct InsertionAcceptanceParams {
    bool emit_low_confidence_calls = false;
    int32_t low_conf_min_support_reads = 2;
    int32_t low_conf_max_tier = 2;
};

struct InsertionAcceptanceDecision {
    bool pass = false;
    bool low_confidence = false;
};

InsertionAcceptanceDecision evaluate_insertion_acceptance(
    const InsertionEvidence& evidence,
    const InsertionAcceptanceParams& params);

struct TeOpenSetInput {
    bool te_gate_pass = false;
    bool te_gate_uncertain_path = false;
    bool has_proxy_signal = false;
    bool proxy_certain = false;
    double te_confidence_prob = 0.0;
};

struct TeOpenSetParams {
    double conf_certain_min = 0.85;
    double conf_uncertain_min = 0.35;
};

struct TeOpenSetDecision {
    std::string status = "NON_TE";
    bool promote_top1_name = false;
    bool add_te_gate_fail_tag = false;
    std::string qc_proxy_tag;
};

TeOpenSetDecision classify_te_open_set(
    const TeOpenSetInput& input,
    const TeOpenSetParams& params);

}  // namespace placer

#endif  // PLACER_DECISION_POLICY_H
