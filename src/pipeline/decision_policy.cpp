#include "decision_policy.h"
#include "pipeline.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace placer {

namespace {

double safe_log(double p) {
    return std::log(std::max(1e-12, p));
}

double phred_from_error_probability(double p_error) {
    const double bounded = std::clamp(p_error, 1e-12, 1.0);
    return -10.0 * std::log10(bounded);
}

double logsumexp3(double a, double b, double c) {
    const double m = std::max(a, std::max(b, c));
    return m + std::log(std::exp(a - m) + std::exp(b - m) + std::exp(c - m));
}

int32_t median_i32(std::vector<int32_t> values) {
    if (values.empty()) {
        return 0;
    }
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

int32_t infer_event_length(const EventGenotypeInput& input) {
    if (input.event_length > 0) {
        return input.event_length;
    }
    return median_i32(input.alt_observed_lengths);
}

double alt_length_support_probability(
    int32_t observed_length,
    int32_t event_length,
    double error_rate) {
    if (event_length <= 0 || observed_length <= 0) {
        return std::clamp(1.0 - (error_rate * 2.0), 0.55, 0.95);
    }
    const double norm_factor = (error_rate <= 0.03) ? 20.0 : 30.0;
    const double sigma = std::max(1.0, static_cast<double>(event_length) / norm_factor);
    const double z =
        (static_cast<double>(observed_length) - static_cast<double>(event_length)) / sigma;
    return std::max(1e-6, std::exp(-0.5 * z * z));
}

double genotype_log_likelihood(
    const EventGenotypeInput& input,
    double alt_copy_fraction) {
    constexpr double kAltCallGivenRef = 1e-4;
    constexpr double kRefCallGivenAlt = 1e-3;
    constexpr double kRefCallGivenRef = 0.999;

    const int32_t event_length = infer_event_length(input);
    double ll = 0.0;

    for (int32_t observed_length : input.alt_observed_lengths) {
        const double p_alt = alt_length_support_probability(
            observed_length,
            event_length,
            input.error_rate);
        const double p =
            (alt_copy_fraction * p_alt) +
            ((1.0 - alt_copy_fraction) * kAltCallGivenRef);
        ll += safe_log(p);
    }

    const int32_t precise_alt_reads = static_cast<int32_t>(input.alt_observed_lengths.size());
    const int32_t coarse_alt_reads = std::max(0, input.alt_struct_reads - precise_alt_reads);
    const double coarse_alt_given_alt = std::clamp(1.0 - (input.error_rate * 3.0), 0.55, 0.95);
    for (int32_t i = 0; i < coarse_alt_reads; ++i) {
        const double p =
            (alt_copy_fraction * coarse_alt_given_alt) +
            ((1.0 - alt_copy_fraction) * kAltCallGivenRef);
        ll += safe_log(p);
    }

    const int32_t ref_reads = std::max(0, input.ref_span_reads);
    for (int32_t i = 0; i < ref_reads; ++i) {
        const double p =
            (alt_copy_fraction * kRefCallGivenAlt) +
            ((1.0 - alt_copy_fraction) * kRefCallGivenRef);
        ll += safe_log(p);
    }

    return ll;
}

double clamp_score(double value, double lo = -3.0, double hi = 3.0) {
    return std::clamp(value, lo, hi);
}

double joint_event_existence_score(const EventExistenceEvidence& existence) {
    constexpr int32_t kUnopposedAltMinReads = 3;
    if (existence.ref_span_reads == 0 &&
        existence.alt_struct_reads >= kUnopposedAltMinReads) {
        return std::max(0.0, existence.score);
    }
    return existence.score;
}

bool is_one_sided_segmentation_pass(const EventSegmentationEvidence& segmentation) {
    return segmentation.pair_valid &&
           segmentation.has_insert_seq &&
           (segmentation.has_left_flank != segmentation.has_right_flank);
}

double adjusted_nonref_existence_score(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    double te_presence) {
    double adjusted = joint_event_existence_score(existence);
    const bool has_strong_structural_te_support =
        segmentation.pair_valid &&
        segmentation.has_insert_seq &&
        te_presence >= 1.0;
    const bool alt_not_weaker_than_ref =
        existence.alt_struct_reads > 0 &&
        existence.alt_struct_reads >= existence.ref_span_reads;
    if (has_strong_structural_te_support && alt_not_weaker_than_ref) {
        adjusted = std::max(0.0, adjusted);
    }
    return adjusted;
}

double adjusted_segmentation_score(const EventSegmentationEvidence& segmentation) {
    if (is_one_sided_segmentation_pass(segmentation)) {
        return std::max(0.25, segmentation.score);
    }
    return segmentation.score;
}

double adjusted_te_boundary_score(
    const EventSegmentationEvidence& segmentation,
    const BoundaryEvidence& boundary) {
    if (!boundary.geometry_defined && is_one_sided_segmentation_pass(segmentation)) {
        return 0.0;
    }
    return boundary.score;
}

double te_resolved_score(const TEAlignmentEvidence& te_alignment) {
    if (te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY") {
        const double identity_term =
            0.45 * clamp_score((te_alignment.best_identity - 0.68) / 0.07, -2.0, 2.0);
        const double coverage_term =
            0.35 * clamp_score((te_alignment.best_query_coverage - 0.55) / 0.07, -2.0, 2.0);
        const double margin_term =
            0.20 * clamp_score((te_alignment.cross_family_margin - 0.05) / 0.03, -2.0, 2.0);
        return identity_term + coverage_term + margin_term;
    }

    const double identity_term =
        0.4 * clamp_score((te_alignment.best_identity - 0.80) / 0.05, -2.0, 2.0);
    const double coverage_term =
        0.4 * clamp_score((te_alignment.best_query_coverage - 0.80) / 0.05, -2.0, 2.0);
    const double margin_term =
        0.2 * clamp_score((te_alignment.cross_family_margin - 0.10) / 0.05, -2.0, 2.0);
    return identity_term + coverage_term + margin_term;
}

double te_existence_score(const TEAlignmentEvidence& te_alignment) {
    constexpr double kUnknownIdentityMin = 0.55;
    constexpr double kUnknownCoverageMin = 0.45;
    constexpr double kDivergentUnknownIdentityMin = 0.52;
    constexpr double kDivergentUnknownCoverageMin = 0.95;

    if (te_alignment.pass &&
        te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT_UNKNOWN") {
        const double identity_term =
            0.10 * clamp_score((te_alignment.best_identity - 0.45) / 0.10, -2.0, 2.0);
        const double coverage_term =
            0.40 * clamp_score((te_alignment.best_query_coverage - kUnknownCoverageMin) / 0.10, -2.0, 2.0);
        const double margin_term =
            0.20 * clamp_score((te_alignment.cross_family_margin - 0.05) / 0.05, -2.0, 2.0);
        const double raw = 0.35 + identity_term + coverage_term + margin_term;
        return std::max(0.0, raw);
    }

    if (te_alignment.qc_reason == "TE_ALIGNMENT_LOW_IDENTITY" &&
        te_alignment.best_identity + 1e-9 >= kDivergentUnknownIdentityMin &&
        te_alignment.best_query_coverage + 1e-9 >= kDivergentUnknownCoverageMin) {
        const double identity_term =
            0.15 * clamp_score((te_alignment.best_identity - kDivergentUnknownIdentityMin) / 0.03, -2.0, 2.0);
        const double coverage_term =
            0.35 * clamp_score((te_alignment.best_query_coverage - kDivergentUnknownCoverageMin) / 0.03, 0.0, 2.0);
        const double raw = 0.45 + identity_term + coverage_term;
        return std::max(0.0, raw);
    }

    const double identity_term =
        0.25 * clamp_score((te_alignment.best_identity - kUnknownIdentityMin) / 0.10, -2.0, 2.0);
    const double coverage_term =
        0.25 * clamp_score((te_alignment.best_query_coverage - kUnknownCoverageMin) / 0.10, -2.0, 2.0);
    const double raw = identity_term + coverage_term;

    const bool supports_unknown_te =
        te_alignment.best_identity + 1e-9 >= kUnknownIdentityMin &&
        te_alignment.best_query_coverage + 1e-9 >= kUnknownCoverageMin;
    return supports_unknown_te ? std::max(0.0, raw) : std::min(0.0, raw);
}

JointHypothesisScore make_hypothesis(FinalHypothesisKind kind) {
    JointHypothesisScore out;
    out.kind = kind;
    return out;
}

bool is_te_hypothesis(FinalHypothesisKind kind) {
    return kind == FinalHypothesisKind::kTeUnknown ||
           kind == FinalHypothesisKind::kTeResolved;
}

constexpr double kJointEmitMinTotalScore = 1.5;
constexpr double kJointEmitMinMargin = 0.75;

}  // namespace

EventGenotypeDecision genotype_event_from_alt_vs_ref(
    const EventGenotypeInput& input) {
    EventGenotypeDecision decision;

    const int32_t alt = std::max(0, input.alt_struct_reads);
    const int32_t ref = std::max(0, input.ref_span_reads);
    const int32_t depth = alt + ref;
    decision.depth = depth;
    if (depth <= 0) {
        return decision;
    }

    decision.allele_fraction = std::clamp(
        static_cast<double>(alt) / static_cast<double>(depth),
        0.0,
        1.0);

    const double ll_00 = genotype_log_likelihood(input, 0.0) + safe_log(0.25);
    const double ll_01 = genotype_log_likelihood(input, 0.5) + safe_log(0.50);
    const double ll_11 = genotype_log_likelihood(input, 1.0) + safe_log(0.25);

    const bool het_is_best_nonref = ll_01 >= ll_11;
    const double best_nonref_ll = het_is_best_nonref ? ll_01 : ll_11;
    decision.best_gt = het_is_best_nonref ? "0/1" : "1/1";

    if (best_nonref_ll <= ll_00) {
        decision.best_gt = "0/0";
        decision.gq = 0;
        decision.pass = false;
        return decision;
    }

    const double delta_ll = best_nonref_ll - ll_00;
    decision.best_nonref_minus_ref_ll = delta_ll;
    const double best_total_ll = std::max(ll_00, std::max(ll_01, ll_11));
    const double posterior_best = std::exp(best_total_ll - logsumexp3(ll_00, ll_01, ll_11));
    const double gq = phred_from_error_probability(1.0 - posterior_best);
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

EventExistenceEvidence build_event_existence_evidence(
    const EventGenotypeInput& input) {
    EventExistenceEvidence evidence;
    const EventGenotypeDecision decision = genotype_event_from_alt_vs_ref(input);
    evidence.best_gt = decision.best_gt;
    evidence.af = decision.allele_fraction;
    evidence.gq = decision.gq;
    evidence.alt_struct_reads = std::max(0, input.alt_struct_reads);
    evidence.ref_span_reads = std::max(0, input.ref_span_reads);
    evidence.depth = decision.depth;
    evidence.best_nonref_minus_ref_ll = decision.best_nonref_minus_ref_ll;
    evidence.score = clamp_score((static_cast<double>(decision.gq) - 20.0) / 20.0);
    return evidence;
}

bool should_emit_te_call(
    const JointHypothesisScore& best_te,
    const JointHypothesisScore& best_non_te) {
    if (!is_te_hypothesis(best_te.kind) || best_te.hard_veto) {
        return false;
    }
    const double margin = best_te.total - best_non_te.total;
    return best_te.total >= kJointEmitMinTotalScore &&
           margin >= kJointEmitMinMargin;
}

JointDecisionResult evaluate_joint_hypotheses(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary) {
    JointDecisionResult result;
    const double te_presence = te_existence_score(te_alignment);
    const double te_resolved = te_resolved_score(te_alignment);
    const double existence_support = adjusted_nonref_existence_score(
        existence,
        segmentation,
        te_presence);
    const double segmentation_support = adjusted_segmentation_score(segmentation);
    const double te_boundary_support = adjusted_te_boundary_score(segmentation, boundary);

    JointHypothesisScore h0 = make_hypothesis(FinalHypothesisKind::kReference);
    h0.existence = -1.0 * std::max(existence.score, 0.0);
    h0.segmentation = -0.5 * std::max(segmentation_support, 0.0);
    h0.te = -0.8 * std::max(te_presence, 0.0);
    h0.total = h0.existence + h0.segmentation + h0.te;
    h0.reason = "REFERENCE";

    JointHypothesisScore h1 = make_hypothesis(FinalHypothesisKind::kInsertionNonTe);
    h1.existence = existence_support;
    h1.segmentation = 0.7 * segmentation_support;
    h1.te = -1.0 * std::max(te_presence, 0.0);
    h1.total = h1.existence + h1.segmentation + h1.te;
    h1.reason = "NON_TE_INSERTION";

    JointHypothesisScore h2 = make_hypothesis(FinalHypothesisKind::kTeUnknown);
    h2.existence = existence_support;
    h2.segmentation = 0.8 * segmentation_support;
    h2.te = 0.6 * te_presence;
    h2.boundary = 0.2 * te_boundary_support;
    h2.total = h2.existence + h2.segmentation + h2.te + h2.boundary;
    h2.reason = "TE_UNKNOWN";
    h2.hard_veto = !segmentation.has_insert_seq || te_alignment.qc_reason == "NO_TE_ALIGNMENT";

    JointHypothesisScore h3 = make_hypothesis(FinalHypothesisKind::kTeResolved);
    h3.existence = existence_support;
    h3.segmentation = 0.8 * segmentation_support;
    h3.te = te_resolved;
    h3.boundary = 0.2 * te_boundary_support;
    h3.total = h3.existence + h3.segmentation + h3.te + h3.boundary;
    h3.reason = "TE_RESOLVED";
    h3.hard_veto = !segmentation.has_insert_seq || !te_alignment.pass ||
                   te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

    std::array<JointHypothesisScore, 4> all = {h0, h1, h2, h3};
    std::sort(all.begin(), all.end(), [](const JointHypothesisScore& lhs, const JointHypothesisScore& rhs) {
        const double lhs_total = lhs.hard_veto ? -1e9 : lhs.total;
        const double rhs_total = rhs.hard_veto ? -1e9 : rhs.total;
        if (lhs_total != rhs_total) {
            return lhs_total > rhs_total;
        }
        return static_cast<int>(lhs.kind) > static_cast<int>(rhs.kind);
    });

    result.best = all[0];
    result.runner_up = all[1];

    const JointHypothesisScore* best_te = nullptr;
    const JointHypothesisScore* best_non_te = nullptr;
    for (const JointHypothesisScore& score : all) {
        if (score.hard_veto) {
            continue;
        }
        if (is_te_hypothesis(score.kind)) {
            if (best_te == nullptr) {
                best_te = &score;
            }
        } else if (best_non_te == nullptr) {
            best_non_te = &score;
        }
        if (best_te != nullptr && best_non_te != nullptr) {
            break;
        }
    }

    // TE existence and TE label resolution are different questions.
    // Emit if the best TE hypothesis beats the best non-TE competitor,
    // regardless of how close H2 and H3 are to each other.
    result.emit_te_call =
        best_te != nullptr &&
        best_non_te != nullptr &&
        should_emit_te_call(*best_te, *best_non_te);
    result.emit_unknown_te =
        result.emit_te_call &&
        best_te->kind == FinalHypothesisKind::kTeUnknown;
    result.final_qc = result.emit_unknown_te
        ? "PASS_FINAL_TE_CALL_UNKNOWN"
        : (result.emit_te_call ? "PASS_FINAL_TE_CALL" : result.best.reason);
    return result;
}

EventSegmentationEvidence analyze_event_segmentation_for_test(
    bool has_consensus,
    const EventSegmentation& segmentation) {
    EventSegmentationEvidence evidence;
    evidence.has_consensus = has_consensus;
    evidence.has_left_flank = segmentation.left_flank_align_len > 0;
    evidence.has_right_flank = segmentation.right_flank_align_len > 0;
    evidence.has_insert_seq = !segmentation.insert_seq.empty();
    evidence.pair_valid = segmentation.pass;
    evidence.left_align_len = segmentation.left_flank_align_len;
    evidence.right_align_len = segmentation.right_flank_align_len;
    evidence.left_identity = segmentation.left_flank_identity;
    evidence.right_identity = segmentation.right_flank_identity;
    evidence.insert_len = static_cast<int32_t>(segmentation.insert_seq.size());
    evidence.qc = segmentation.qc_reason;

    int32_t scoring_flank_len = std::min(evidence.left_align_len, evidence.right_align_len);
    double scoring_flank_identity = std::min(evidence.left_identity, evidence.right_identity);
    const bool one_sided_pass =
        evidence.pair_valid &&
        evidence.has_insert_seq &&
        (evidence.has_left_flank != evidence.has_right_flank);
    if (one_sided_pass) {
        scoring_flank_len =
            evidence.has_left_flank ? evidence.left_align_len : evidence.right_align_len;
        scoring_flank_identity =
            evidence.has_left_flank ? evidence.left_identity : evidence.right_identity;
    }
    evidence.score =
        0.5 * clamp_score((static_cast<double>(scoring_flank_len) - 50.0) / 25.0, -2.0, 2.0) +
        0.5 * clamp_score((scoring_flank_identity - 0.90) / 0.05, -2.0, 2.0);
    if (one_sided_pass) {
        evidence.score -= 0.5;
    } else if (!evidence.has_left_flank || !evidence.has_right_flank) {
        evidence.score -= 1.0;
    }
    if (!evidence.pair_valid) {
        evidence.score -= 1.5;
    }
    return evidence;
}

BoundaryEvidence evaluate_boundary_evidence(
    const FinalBoundaryInput& input,
    int32_t breakpoint_envelope_width) {
    BoundaryEvidence evidence;
    const FinalBoundaryDecision canonical = check_boundary_consistency(input);
    evidence.geometry_defined =
        input.left_ref_start >= 0 && input.left_ref_end >= 0 &&
        input.right_ref_start >= 0 && input.right_ref_end >= 0 &&
        input.left_ref_start < input.left_ref_end &&
        input.right_ref_start < input.right_ref_end;
    evidence.canonical_pass = canonical.pass;
    evidence.boundary_type = canonical.boundary_type;
    evidence.boundary_len = canonical.boundary_len;
    evidence.qc = canonical.qc;

    if (!evidence.geometry_defined) {
        evidence.score = -2.0;
        return evidence;
    }

    if (canonical.pass) {
        evidence.evidence_consistent = true;
        evidence.score = 1.0;
        return evidence;
    }

    const int32_t link_delta = input.right_ref_start - input.left_ref_end;
    const int32_t noncanonical_span = std::abs(link_delta);
    evidence.evidence_consistent =
        breakpoint_envelope_width > 0 && noncanonical_span <= breakpoint_envelope_width;
    if (evidence.evidence_consistent) {
        evidence.boundary_type = "NONCANONICAL";
        evidence.boundary_len = noncanonical_span;
        evidence.qc = "PASS_BOUNDARY_NONCANONICAL_CONSISTENT";
        evidence.score = 0.25;
    } else {
        evidence.score = -2.0;
    }
    return evidence;
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
