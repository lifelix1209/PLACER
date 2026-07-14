#ifdef NDEBUG
#undef NDEBUG
#endif

#include "decision_policy.h"
#include "mechanistic_evidence.h"
#include "pipeline.h"

#include <cassert>
#include <string>

namespace {

placer::EventExistenceEvidence strong_existence() {
    placer::EventExistenceEvidence e;
    e.alt_struct_reads = 18;
    e.alt_split_reads = 8;
    e.alt_indel_reads = 3;
    e.alt_left_clip_reads = 5;
    e.alt_right_clip_reads = 5;
    e.ref_span_reads = 0;
    e.af = 0.65;
    e.gq = 60;
    e.score = 2.0;
    return e;
}

placer::EventSegmentationEvidence closed_segmentation() {
    placer::EventSegmentationEvidence s;
    s.has_consensus = true;
    s.has_insert_seq = true;
    s.has_left_flank = true;
    s.has_right_flank = true;
    s.pair_valid = true;
    s.insert_len = 320;
    s.score = 1.5;
    s.qc = "PASS_EVENT_SEGMENTATION";
    return s;
}

placer::BoundaryEvidence tsd_boundary() {
    placer::BoundaryEvidence b;
    b.geometry_defined = true;
    b.canonical_pass = true;
    b.evidence_consistent = true;
    b.boundary_type = "TSD";
    b.boundary_len = 12;
    b.score = 1.0;
    b.qc = "PASS_BOUNDARY_TSD";
    return b;
}

placer::TEAlignmentEvidence resolved_te() {
    placer::TEAlignmentEvidence t;
    t.pass = true;
    t.best_family = "L1";
    t.best_subfamily = "L1HS";
    t.best_identity = 0.96;
    t.best_query_coverage = 0.90;
    t.cross_family_margin = 0.18;
    t.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
    t.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
    t.sequence_model_score = 0.35;
    t.annotation_confidence = "HIGH";
    return t;
}

placer::TEAlignmentEvidence structured_te(
    double coverage,
    double residual_fraction) {
    placer::TEAlignmentEvidence t = resolved_te();
    t.best_query_coverage = coverage;
    t.annotation_residual_fraction = residual_fraction;
    t.te_sequence_explanation = placer::explain_te_sequence_structure(
        std::string(220, 'G') + std::string(30, 'A'),
        t.qc_reason,
        t.best_family,
        t.best_subfamily,
        t.best_identity,
        coverage,
        residual_fraction,
        t.annotation_masked_fraction,
        t.cross_family_margin,
        t.second_score,
        t.sequence_model_label,
        t.sequence_model_score);
    return t;
}

void strong_mechanism_has_positive_lower_bfs() {
    const placer::MechanisticEvidenceCertificate cert =
        placer::build_mechanistic_evidence_certificate(
            strong_existence(),
            closed_segmentation(),
            resolved_te(),
            tsd_boundary(),
            nullptr);

    assert(cert.lower_log_bf_te_vs_artifact > 0.0);
    assert(cert.lower_log_bf_te_vs_non_te > 0.0);
    assert(cert.mechanistic_support_signal > 0.0);
    assert(!cert.blocks.empty());
}

void reference_conflict_decreases_te_vs_artifact_bf() {
    placer::EventExistenceEvidence clean = strong_existence();
    placer::EventExistenceEvidence conflicted = clean;
    conflicted.ref_span_reads = 18;

    const auto clean_cert = placer::build_mechanistic_evidence_certificate(
        clean, closed_segmentation(), resolved_te(), tsd_boundary(), nullptr);
    const auto conflicted_cert = placer::build_mechanistic_evidence_certificate(
        conflicted, closed_segmentation(), resolved_te(), tsd_boundary(), nullptr);

    assert(conflicted_cert.lower_log_bf_te_vs_artifact <
           clean_cert.lower_log_bf_te_vs_artifact);
    assert(conflicted_cert.ref_conflict_signal >
           clean_cert.ref_conflict_signal);
}

void outlier_sequence_decreases_sequence_block() {
    placer::TEAlignmentEvidence outlier = resolved_te();
    outlier.sequence_model_label = "TE_MODEL_OUTLIER";
    outlier.sequence_model_score = -0.80;

    const auto clean_cert = placer::build_mechanistic_evidence_certificate(
        strong_existence(), closed_segmentation(), resolved_te(), tsd_boundary(), nullptr);
    const auto outlier_cert = placer::build_mechanistic_evidence_certificate(
        strong_existence(), closed_segmentation(), outlier, tsd_boundary(), nullptr);

    assert(outlier_cert.sequence_lower_log_lr <
           clean_cert.sequence_lower_log_lr);
}

void certificate_exposes_structure_block() {
    const auto cert = placer::build_mechanistic_evidence_certificate(
        strong_existence(),
        closed_segmentation(),
        structured_te(0.88, 0.12),
        tsd_boundary(),
        nullptr);

    assert(cert.structure_te_log_evidence != 0.0);
    assert(cert.structure_lower_log_lr > 0.0);
    const std::string serialized = placer::serialize_mechanistic_blocks(cert);
    assert(serialized.find("structure") != std::string::npos);
}

void stronger_structure_evidence_decreases_robust_lfdr() {
    const auto weak_cert = placer::build_mechanistic_evidence_certificate(
        strong_existence(),
        closed_segmentation(),
        structured_te(0.55, 0.45),
        tsd_boundary(),
        nullptr);
    const auto strong_cert = placer::build_mechanistic_evidence_certificate(
        strong_existence(),
        closed_segmentation(),
        structured_te(0.92, 0.08),
        tsd_boundary(),
        nullptr);

    assert(strong_cert.structure_lower_log_lr > weak_cert.structure_lower_log_lr);
    const placer::PriorInterval prior;
    const auto weak_lfdr = placer::evaluate_robust_mechanistic_lfdr(weak_cert, prior, 0.10);
    const auto strong_lfdr = placer::evaluate_robust_mechanistic_lfdr(strong_cert, prior, 0.10);
    assert(strong_lfdr.worst_case_lfdr < weak_lfdr.worst_case_lfdr);
}

void robust_lfdr_decreases_with_stronger_te_bf() {
    placer::MechanisticEvidenceCertificate weak;
    weak.lower_log_bf_te_vs_artifact = 0.5;
    weak.lower_log_bf_te_vs_non_te = 0.5;

    placer::MechanisticEvidenceCertificate strong = weak;
    strong.lower_log_bf_te_vs_artifact = 5.0;
    strong.lower_log_bf_te_vs_non_te = 5.0;

    const placer::PriorInterval prior;
    const auto weak_eval = placer::evaluate_robust_mechanistic_lfdr(weak, prior, 0.10);
    const auto strong_eval = placer::evaluate_robust_mechanistic_lfdr(strong, prior, 0.10);

    assert(strong_eval.worst_case_lfdr < weak_eval.worst_case_lfdr);
}

void higher_ambiguity_increases_worst_case_lfdr() {
    placer::MechanisticEvidenceCertificate clean;
    clean.lower_log_bf_te_vs_artifact = 4.0;
    clean.lower_log_bf_te_vs_non_te = 4.0;
    clean.ambiguity_width = 0.1;

    placer::MechanisticEvidenceCertificate ambiguous = clean;
    ambiguous.ambiguity_width = 2.0;

    const placer::PriorInterval prior;
    const auto clean_eval = placer::evaluate_robust_mechanistic_lfdr(clean, prior, 0.10);
    const auto ambiguous_eval = placer::evaluate_robust_mechanistic_lfdr(ambiguous, prior, 0.10);

    assert(ambiguous_eval.worst_case_lfdr > clean_eval.worst_case_lfdr);
}

}  // namespace

int main() {
    strong_mechanism_has_positive_lower_bfs();
    reference_conflict_decreases_te_vs_artifact_bf();
    outlier_sequence_decreases_sequence_block();
    certificate_exposes_structure_block();
    stronger_structure_evidence_decreases_robust_lfdr();
    robust_lfdr_decreases_with_stronger_te_bf();
    higher_ambiguity_increases_worst_case_lfdr();
    return 0;
}
