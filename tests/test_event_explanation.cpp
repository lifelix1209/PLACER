#ifdef NDEBUG
#undef NDEBUG
#endif

#include "event_explanation.h"

#include <cassert>
#include <vector>

namespace {

placer::EventExplanation make_explanation(
    placer::ExplanationKind kind,
    int structural,
    int missing,
    int unexplained,
    int read_conflicts,
    int ref_counter,
    int artifact,
    int breakpoint,
    int ambiguity,
    int path_complexity) {
    placer::EventExplanation out;
    out.kind = kind;
    out.residual.structural_conflicts = structural;
    out.residual.missing_required_components = missing;
    out.residual.unexplained_high_complexity_bases = unexplained;
    out.residual.read_assignment_conflicts = read_conflicts;
    out.residual.reference_counterevidence = ref_counter;
    out.residual.artifact_evidence = artifact;
    out.residual.breakpoint_disagreement_bp = breakpoint;
    out.residual.label_ambiguity = ambiguity;
    out.residual.path_complexity = path_complexity;
    return out;
}

void test_dominance_requires_no_worse_primary_dimensions() {
    const auto te = make_explanation(
        placer::ExplanationKind::kTe,
        0, 0, 0, 1, 1, 0, 10, 0, 2);
    const auto non_te = make_explanation(
        placer::ExplanationKind::kInsertionNonTe,
        0, 0, 12, 1, 1, 0, 0, 0, 1);

    assert(placer::dominates_primary_residuals(te.residual, non_te.residual));
    assert(!placer::dominates_primary_residuals(non_te.residual, te.residual));
}

void test_non_dominated_te_becomes_ambiguous() {
    std::vector<placer::EventExplanation> explanations;
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kReference,
        1, 0, 80, 8, 0, 0, 0, 0, 0));
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kInsertionNonTe,
        0, 0, 5, 2, 1, 0, 0, 0, 1));
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kTe,
        0, 0, 0, 2, 1, 0, 0, 2, 2));
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kArtifact,
        0, 0, 0, 1, 1, 3, 0, 0, 1));

    const placer::ExplanationDecision decision =
        placer::compare_event_explanations(explanations, true);

    assert(decision.best.kind == placer::ExplanationKind::kTe);
    assert(decision.final_qc == "TE_AMBIGUOUS");
    assert(!decision.emit_te_call);
}

void test_unique_closed_te_dominance_emits_pass_te_closed() {
    std::vector<placer::EventExplanation> explanations;
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kReference,
        1, 0, 120, 10, 0, 0, 0, 0, 0));
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kInsertionNonTe,
        0, 0, 90, 6, 1, 0, 0, 0, 1));
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kTe,
        0, 0, 0, 0, 0, 0, 0, 0, 1));
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kArtifact,
        0, 0, 100, 8, 1, 3, 0, 0, 1));

    const placer::ExplanationDecision decision =
        placer::compare_event_explanations(explanations, true);

    assert(decision.best.kind == placer::ExplanationKind::kTe);
    assert(decision.final_qc == "PASS_TE_CLOSED");
    assert(decision.emit_te_call);
    assert(!decision.emit_unknown_te);
}

void test_unique_imprecise_te_dominance_emits_pass_te_imprecise() {
    std::vector<placer::EventExplanation> explanations;
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kReference,
        1, 0, 120, 10, 0, 0, 0, 0, 0));
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kInsertionNonTe,
        0, 0, 90, 6, 1, 0, 0, 0, 1));
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kTe,
        0, 0, 0, 0, 0, 0, 25, 1, 2));
    explanations.push_back(make_explanation(
        placer::ExplanationKind::kArtifact,
        0, 0, 100, 8, 1, 3, 0, 0, 1));

    const placer::ExplanationDecision decision =
        placer::compare_event_explanations(explanations, false);

    assert(decision.best.kind == placer::ExplanationKind::kTe);
    assert(decision.final_qc == "PASS_TE_IMPRECISE");
    assert(decision.emit_te_call);
    assert(decision.emit_unknown_te);
}

}  // namespace

int main() {
    test_dominance_requires_no_worse_primary_dimensions();
    test_non_dominated_te_becomes_ambiguous();
    test_unique_closed_te_dominance_emits_pass_te_closed();
    test_unique_imprecise_te_dominance_emits_pass_te_imprecise();
    return 0;
}
