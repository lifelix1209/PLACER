#include "decision_policy.h"

#include <algorithm>
#include <cmath>

namespace placer {

namespace {

double log_binomial_likelihood(int32_t alt_count, int32_t ref_count, double alt_rate) {
    const double p = std::clamp(alt_rate, 1e-6, 1.0 - 1e-6);
    return (static_cast<double>(alt_count) * std::log(p)) +
           (static_cast<double>(ref_count) * std::log(1.0 - p));
}

}  // namespace

EventGenotypeDecision genotype_event_from_alt_vs_ref(
    const EventGenotypeInput& input) {
    EventGenotypeDecision decision;

    const int32_t alt = std::max(0, input.alt_struct_reads);
    const int32_t ref = std::max(0, input.ref_span_reads);
    const int32_t depth = alt + ref;
    if (depth <= 0) {
        return decision;
    }

    decision.allele_fraction = std::clamp(
        static_cast<double>(alt) / static_cast<double>(depth),
        0.0,
        1.0);

    if (depth < std::max(1, input.min_depth)) {
        return decision;
    }

    const double error_rate = std::clamp(input.error_rate, 1e-4, 0.25);
    const double ll_00 = log_binomial_likelihood(alt, ref, error_rate);
    const double ll_01 = log_binomial_likelihood(alt, ref, 0.5);
    const double ll_11 = log_binomial_likelihood(alt, ref, 1.0 - error_rate);

    const bool het_is_best_nonref = ll_01 >= ll_11;
    const double best_nonref_ll = het_is_best_nonref ? ll_01 : ll_11;
    decision.best_gt = het_is_best_nonref ? "0/1" : "1/1";

    if (best_nonref_ll <= ll_00) {
        decision.best_gt = "0/0";
        decision.gq = 0;
        decision.pass = false;
        return decision;
    }

    // Event existence is a non-reference vs reference decision. The caller
    // still emits the better non-reference genotype, but confidence should not
    // collapse just because 0/1 and 1/1 are both plausible.
    const double delta_ll = best_nonref_ll - ll_00;
    const double gq = 4.3429448190325175 * delta_ll;
    decision.gq = std::max(0, std::min(99, static_cast<int32_t>(std::lround(gq))));
    decision.pass = decision.gq >= std::max(0, input.min_gq);
    return decision;
}

FinalBoundaryDecision check_boundary_consistency(
    const FinalBoundaryInput& input) {
    FinalBoundaryDecision decision;

    if (input.left_ref_start < 0 ||
        input.left_ref_end < 0 ||
        input.right_ref_start < 0 ||
        input.right_ref_end < 0) {
        decision.qc = "REJECT_BOUNDARY_MISSING_REF_SEGMENTS";
        return decision;
    }
    if (input.left_ref_start >= input.left_ref_end ||
        input.right_ref_start >= input.right_ref_end) {
        decision.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";
        return decision;
    }

    const int32_t tsd_min_len = std::max(1, input.tsd_min_len);
    const int32_t tsd_max_len = std::max(tsd_min_len, input.tsd_max_len);
    const int32_t link_delta = input.right_ref_start - input.left_ref_end;

    if (link_delta < 0) {
        const int32_t overlap = -link_delta;
        if (overlap < tsd_min_len || overlap > tsd_max_len) {
            decision.qc = "REJECT_BOUNDARY_TSD_RANGE";
            return decision;
        }
        decision.pass = true;
        decision.boundary_type = "TSD";
        decision.boundary_len = overlap;
        decision.qc = "PASS_BOUNDARY_TSD";
        return decision;
    }

    if (link_delta == 0) {
        decision.pass = true;
        decision.boundary_type = "BLUNT";
        decision.boundary_len = 0;
        decision.qc = "PASS_BOUNDARY_BLUNT";
        return decision;
    }

    if (link_delta <= tsd_max_len) {
        decision.pass = true;
        decision.boundary_type = "SMALL_DEL";
        decision.boundary_len = link_delta;
        decision.qc = "PASS_BOUNDARY_SMALL_DEL";
        return decision;
    }

    decision.qc = "REJECT_BOUNDARY_DEL_RANGE";
    return decision;
}

FinalTeAcceptanceDecision evaluate_final_te_acceptance(
    const FinalTeAcceptanceInput& input) {
    FinalTeAcceptanceDecision decision;
    if (!input.event_existence_pass) {
        return decision;
    }
    if (!input.event_closure_pass) {
        decision.qc = "REJECT_EVENT_CLOSURE";
        return decision;
    }
    if (!input.te_sequence_pass) {
        decision.qc = "REJECT_TE_SEQUENCE";
        return decision;
    }
    if (!input.boundary_pass) {
        decision.qc = "REJECT_BOUNDARY";
        return decision;
    }
    decision.pass = true;
    decision.qc = "PASS_FINAL_TE_CALL";
    return decision;
}

}  // namespace placer
