#include "event_explanation.h"

#include <algorithm>
#include <sstream>

namespace placer {
namespace {

bool is_primary_no_worse(
    const ExplanationResidual& lhs,
    const ExplanationResidual& rhs) {
    return lhs.structural_conflicts <= rhs.structural_conflicts &&
           lhs.missing_required_components <= rhs.missing_required_components &&
           lhs.unexplained_high_complexity_bases <= rhs.unexplained_high_complexity_bases &&
           lhs.read_assignment_conflicts <= rhs.read_assignment_conflicts &&
           lhs.reference_counterevidence <= rhs.reference_counterevidence &&
           lhs.artifact_evidence <= rhs.artifact_evidence;
}

bool is_primary_strictly_better(
    const ExplanationResidual& lhs,
    const ExplanationResidual& rhs) {
    return lhs.structural_conflicts < rhs.structural_conflicts ||
           lhs.missing_required_components < rhs.missing_required_components ||
           lhs.unexplained_high_complexity_bases < rhs.unexplained_high_complexity_bases ||
           lhs.read_assignment_conflicts < rhs.read_assignment_conflicts ||
           lhs.reference_counterevidence < rhs.reference_counterevidence ||
           lhs.artifact_evidence < rhs.artifact_evidence;
}

int tie_break_tuple_sum(const ExplanationResidual& residual) {
    return residual.breakpoint_disagreement_bp +
           residual.label_ambiguity +
           residual.path_complexity;
}

bool explanation_sort_less(
    const EventExplanation& lhs,
    const EventExplanation& rhs) {
    const bool lhs_dominates = dominates_primary_residuals(lhs.residual, rhs.residual);
    const bool rhs_dominates = dominates_primary_residuals(rhs.residual, lhs.residual);
    if (lhs_dominates != rhs_dominates) {
        return lhs_dominates;
    }
    const int lhs_primary =
        lhs.residual.structural_conflicts +
        lhs.residual.missing_required_components +
        lhs.residual.unexplained_high_complexity_bases +
        lhs.residual.read_assignment_conflicts +
        lhs.residual.reference_counterevidence +
        lhs.residual.artifact_evidence;
    const int rhs_primary =
        rhs.residual.structural_conflicts +
        rhs.residual.missing_required_components +
        rhs.residual.unexplained_high_complexity_bases +
        rhs.residual.read_assignment_conflicts +
        rhs.residual.reference_counterevidence +
        rhs.residual.artifact_evidence;
    if (lhs_primary != rhs_primary) {
        return lhs_primary < rhs_primary;
    }
    const int lhs_tie = tie_break_tuple_sum(lhs.residual);
    const int rhs_tie = tie_break_tuple_sum(rhs.residual);
    if (lhs_tie != rhs_tie) {
        return lhs_tie < rhs_tie;
    }
    return static_cast<int>(lhs.kind) > static_cast<int>(rhs.kind);
}

bool has_missing_required_components(const EventExplanation& explanation) {
    return explanation.residual.missing_required_components > 0 ||
           explanation.residual.structural_conflicts > 0;
}

}  // namespace

const char* explanation_kind_name(ExplanationKind kind) {
    switch (kind) {
        case ExplanationKind::kReference: return "REFERENCE";
        case ExplanationKind::kInsertionNonTe: return "INSERTION_NON_TE";
        case ExplanationKind::kTe: return "TE";
        case ExplanationKind::kArtifact: return "ARTIFACT";
    }
    return "UNKNOWN";
}

std::string serialize_residual(const ExplanationResidual& residual) {
    std::ostringstream out;
    out << "structural=" << residual.structural_conflicts
        << ";missing=" << residual.missing_required_components
        << ";unexplained=" << residual.unexplained_high_complexity_bases
        << ";breakpoint=" << residual.breakpoint_disagreement_bp
        << ";read_conflicts=" << residual.read_assignment_conflicts
        << ";ref_counter=" << residual.reference_counterevidence
        << ";artifact=" << residual.artifact_evidence
        << ";label_ambiguity=" << residual.label_ambiguity
        << ";path_complexity=" << residual.path_complexity;
    return out.str();
}

std::string serialize_explanation_path(const EventExplanation& explanation) {
    if (explanation.path.empty()) {
        return "NA";
    }
    std::ostringstream out;
    for (size_t i = 0; i < explanation.path.size(); ++i) {
        if (i > 0) {
            out << ",";
        }
        const auto& segment = explanation.path[i];
        out << segment.query_start << "-" << segment.query_end
            << ":" << segment.label;
    }
    return out.str();
}

bool dominates_primary_residuals(
    const ExplanationResidual& lhs,
    const ExplanationResidual& rhs) {
    return is_primary_no_worse(lhs, rhs) &&
           is_primary_strictly_better(lhs, rhs);
}

ExplanationDecision compare_event_explanations(
    const std::vector<EventExplanation>& explanations,
    bool closed_breakpoints) {
    ExplanationDecision decision;
    if (explanations.empty()) {
        return decision;
    }

    std::vector<EventExplanation> ranked = explanations;
    std::sort(ranked.begin(), ranked.end(), explanation_sort_less);
    decision.best = ranked.front();
    decision.alternatives = ranked;

    bool uniquely_dominates_all = true;
    for (size_t i = 1; i < ranked.size(); ++i) {
        if (!dominates_primary_residuals(decision.best.residual, ranked[i].residual)) {
            uniquely_dominates_all = false;
            break;
        }
    }

    if (has_missing_required_components(decision.best)) {
        decision.final_qc = "NO_CALL_INCOMPLETE";
        return decision;
    }

    if (decision.best.kind == ExplanationKind::kTe) {
        if (!uniquely_dominates_all) {
            decision.final_qc = "TE_AMBIGUOUS";
            return decision;
        }
        decision.emit_te_call = true;
        decision.emit_evidence_te_call = !closed_breakpoints;
        decision.emit_unknown_te = !closed_breakpoints ||
                                   decision.best.subfamily.empty() ||
                                   decision.best.subfamily == "UNKNOWN" ||
                                   decision.best.family == "UNKNOWN";
        decision.final_qc = closed_breakpoints ? "PASS_TE_CLOSED" : "PASS_TE_IMPRECISE";
        return decision;
    }

    if (decision.best.kind == ExplanationKind::kInsertionNonTe && uniquely_dominates_all) {
        decision.final_qc = "PASS_NONTE_INSERTION";
        return decision;
    }

    decision.final_qc = "REFERENCE_OR_ARTIFACT";
    return decision;
}

}  // namespace placer
