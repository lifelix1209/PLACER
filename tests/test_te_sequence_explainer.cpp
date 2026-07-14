#ifdef NDEBUG
#undef NDEBUG
#endif

#include "pipeline.h"
#include "te_sequence_explainer.h"

#include <cassert>
#include <cmath>
#include <string>

namespace {

bool path_contains_kind(
    const placer::TeSequenceExplanation& explanation,
    const std::string& kind) {
    for (const auto& segment : explanation.path) {
        if (segment.kind == kind) {
            return true;
        }
    }
    return false;
}

void test_resolved_qc_maps_to_resolved_status() {
    const placer::TeSequenceExplanation explanation =
        placer::explain_te_alignment_shadow(
            200,
            "PASS_INSERT_TE_ALIGNMENT",
            "Gypsy",
            "Gypsy-1",
            0.91,
            0.86,
            0.14,
            0.02,
            0.01,
            0.80);

    assert(explanation.status == placer::TeAnnotationStatus::kResolved);
    assert(explanation.family == "Gypsy");
    assert(explanation.subfamily == "Gypsy-1");
    assert(!explanation.path.empty());
    assert(std::string(placer::te_annotation_status_name(explanation.status)) == "RESOLVED");
}

void test_family_only_qc_maps_to_family_only_status() {
    const placer::TeSequenceExplanation explanation =
        placer::explain_te_alignment_shadow(
            180,
            "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY",
            "Copia",
            "",
            0.75,
            0.70,
            0.30,
            0.0,
            0.90,
            0.10);

    assert(explanation.status == placer::TeAnnotationStatus::kFamilyOnly);
    assert(explanation.family == "Copia");
    assert(explanation.subfamily == "UNKNOWN");
}

void test_unknown_qc_maps_to_unknown_without_margin_gate() {
    const placer::TeSequenceExplanation low_margin =
        placer::explain_te_alignment_shadow(
            160,
            "PASS_INSERT_TE_ALIGNMENT_UNKNOWN",
            "UNKNOWN",
            "UNKNOWN",
            0.56,
            0.75,
            0.25,
            0.05,
            0.001,
            0.50);
    const placer::TeSequenceExplanation high_margin =
        placer::explain_te_alignment_shadow(
            160,
            "PASS_INSERT_TE_ALIGNMENT_UNKNOWN",
            "UNKNOWN",
            "UNKNOWN",
            0.56,
            0.75,
            0.25,
            0.05,
            0.80,
            0.50);

    assert(low_margin.status == placer::TeAnnotationStatus::kUnknownTe);
    assert(high_margin.status == placer::TeAnnotationStatus::kUnknownTe);
}

void test_low_identity_qc_is_non_te_like_with_residual() {
    const placer::TeSequenceExplanation explanation =
        placer::explain_te_alignment_shadow(
            100,
            "TE_ALIGNMENT_LOW_IDENTITY",
            "Gypsy",
            "Gypsy-1",
            0.30,
            0.60,
            0.40,
            0.10,
            0.20,
            0.10);

    assert(explanation.status == placer::TeAnnotationStatus::kNonTeLike);
    assert(explanation.residual.unexplained_high_complexity_bases > 0);
    assert(explanation.residual.edit_distance > 0);
}

void test_serialization_is_stable() {
    placer::TeSequenceExplanation explanation;
    explanation.status = placer::TeAnnotationStatus::kResolved;
    explanation.family = "Gypsy";
    explanation.subfamily = "Gypsy-1";
    explanation.residual.unexplained_high_complexity_bases = 7;
    explanation.residual.edit_distance = 3;
    explanation.residual.family_conflicts = 1;
    explanation.path.push_back({"TE", 0, 93, "Gypsy/Gypsy-1"});

    assert(placer::serialize_te_path_residual(explanation.residual) ==
           "unexplained=7;edit_distance=3;family_conflicts=1;"
           "subfamily_conflicts=0;orientation_conflicts=0;segment_breaks=0;"
           "low_complexity_only=0");
    assert(placer::serialize_te_sequence_path(explanation) ==
           "0-93:TE:Gypsy/Gypsy-1");
}

void test_te_alignment_evidence_carries_shadow_explanation() {
    placer::TEAlignmentEvidence evidence;
    assert(evidence.te_sequence_explanation.status ==
           placer::TeAnnotationStatus::kUnavailable);
}

void test_structure_model_splits_te_core_transduction_polya() {
    const std::string te_core(120, 'C');
    const std::string transduction =
        "ACGTGCACTGATCGTACGATGCTAGCTAGGTCAGTACGAT";
    const std::string poly_a(36, 'A');
    const std::string insert_seq = te_core + transduction + poly_a;

    const placer::TeSequenceExplanation explanation =
        placer::explain_te_sequence_structure(
            insert_seq,
            "PASS_INSERT_TE_ALIGNMENT",
            "Gypsy",
            "Gypsy-1",
            0.94,
            static_cast<double>(te_core.size()) /
                static_cast<double>(insert_seq.size()),
            1.0 - (static_cast<double>(te_core.size()) /
                   static_cast<double>(insert_seq.size())),
            0.0,
            0.20,
            0.10,
            "TE_MODEL_IN_DISTRIBUTION",
            0.45);

    assert(explanation.status == placer::TeAnnotationStatus::kResolved);
    assert(path_contains_kind(explanation, "TE_CORE"));
    assert(path_contains_kind(explanation, "TRANSDUCTION"));
    assert(path_contains_kind(explanation, "POLYA"));
    assert(!path_contains_kind(explanation, "TE"));
    assert(std::fabs(explanation.te_core_coverage -
                     (static_cast<double>(te_core.size()) /
                      static_cast<double>(insert_seq.size()))) < 0.05);
    assert(explanation.polyA_posterior > 0.80);
    assert(explanation.transduction_posterior > 0.50);
    assert(explanation.te_structure_log_evidence >
           explanation.nonte_structure_log_evidence);
    assert(explanation.te_structure_log_evidence >
           explanation.artifact_structure_log_evidence);
    assert(explanation.structure_path_confidence > 0.50);
}

void test_structure_model_keeps_high_complexity_residual_out_of_te_core() {
    const std::string te_core(90, 'G');
    const std::string residual =
        "GTCAGTCAGATCGATGCTAGCTAGGACCTAGTCAGATCGATGCTAGCTA";
    const std::string insert_seq = te_core + residual;

    const placer::TeSequenceExplanation explanation =
        placer::explain_te_sequence_structure(
            insert_seq,
            "PASS_INSERT_TE_ALIGNMENT",
            "L1",
            "L1-1",
            0.95,
            static_cast<double>(te_core.size()) /
                static_cast<double>(insert_seq.size()),
            static_cast<double>(residual.size()) /
                static_cast<double>(insert_seq.size()),
            0.0,
            0.25,
            0.0,
            "TE_MODEL_IN_DISTRIBUTION",
            0.40);

    assert(path_contains_kind(explanation, "TE_CORE"));
    assert(path_contains_kind(explanation, "UNEXPLAINED"));
    assert(explanation.unexplained_high_complexity_bp >=
           static_cast<int32_t>(residual.size() * 0.80));
    assert(explanation.te_core_coverage < 0.70);
}

}  // namespace

int main() {
    test_resolved_qc_maps_to_resolved_status();
    test_family_only_qc_maps_to_family_only_status();
    test_unknown_qc_maps_to_unknown_without_margin_gate();
    test_low_identity_qc_is_non_te_like_with_residual();
    test_serialization_is_stable();
    test_te_alignment_evidence_carries_shadow_explanation();
    test_structure_model_splits_te_core_transduction_polya();
    test_structure_model_keeps_high_complexity_residual_out_of_te_core();
    return 0;
}
