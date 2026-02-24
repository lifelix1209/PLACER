#include "decision_policy.h"

#include <algorithm>

namespace placer {

InsertionAcceptanceDecision evaluate_insertion_acceptance(
    const InsertionEvidence& evidence,
    const InsertionAcceptanceParams& params) {
    InsertionAcceptanceDecision decision;

    if (evidence.tier <= 2) {
        decision.pass = true;
        return decision;
    }

    if (!params.emit_low_confidence_calls) {
        return decision;
    }

    const int32_t min_support = std::max(1, params.low_conf_min_support_reads);
    const int32_t max_tier = std::clamp(params.low_conf_max_tier, 1, 3);
    const bool low_conf_ok =
        evidence.support_reads >= min_support &&
        evidence.tier <= max_tier;
    if (low_conf_ok) {
        decision.pass = true;
        decision.low_confidence = true;
    }
    return decision;
}

TeOpenSetDecision classify_te_open_set(
    const TeOpenSetInput& input,
    const TeOpenSetParams& params) {
    TeOpenSetDecision decision;

    const double conf_certain_min = std::clamp(params.conf_certain_min, 0.0, 1.0);
    const double conf_uncertain_min = std::clamp(
        std::min(params.conf_uncertain_min, conf_certain_min),
        0.0,
        1.0);
    const bool prob_certain = input.te_confidence_prob >= conf_certain_min;
    const bool prob_uncertain = input.te_confidence_prob >= conf_uncertain_min;

    const bool uncertain_path = !input.te_gate_pass || input.te_gate_uncertain_path;
    if (!uncertain_path) {
        decision.status = "TE_CERTAIN";
        return decision;
    }

    if (!input.te_gate_pass) {
        decision.add_te_gate_fail_tag = true;
    }
    if (input.has_proxy_signal && (input.proxy_certain || prob_certain)) {
        decision.status = "TE_CERTAIN";
        decision.promote_top1_name = true;
        decision.qc_proxy_tag = "TE_PROXY_CERTAIN";
    } else if (input.has_proxy_signal && prob_uncertain) {
        decision.status = "TE_UNCERTAIN";
        decision.qc_proxy_tag = "TE_PROXY_WEAK";
    } else {
        decision.status = "NON_TE";
        decision.qc_proxy_tag = "TE_PROXY_NONE";
    }

    return decision;
}

}  // namespace placer
