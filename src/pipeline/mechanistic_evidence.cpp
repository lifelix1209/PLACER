#include "mechanistic_evidence.h"

#include "decision_policy.h"
#include "pipeline.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace placer {
namespace {

double clamp01(double value) {
    return std::clamp(value, 0.0, 1.0);
}

double count_signal(int32_t count, double scale) {
    if (count <= 0) {
        return 0.0;
    }
    return clamp01(1.0 - std::exp(-static_cast<double>(count) / std::max(scale, 1e-6)));
}

int32_t positive_or_zero(int32_t value) {
    return std::max(0, value);
}

bool is_one_sided_segmentation_pass(const EventSegmentationEvidence& segmentation) {
    return segmentation.has_insert_seq &&
           !segmentation.pair_valid &&
           (segmentation.has_left_flank != segmentation.has_right_flank);
}

double mechanistic_read_signal(
    const EventExistenceEvidence& existence,
    const ClipInsertConcordanceEvidence* clip_insert_concordance) {
    const int32_t alt = positive_or_zero(existence.alt_struct_reads);
    int32_t mechanistic_reads =
        positive_or_zero(existence.alt_split_reads) +
        positive_or_zero(existence.alt_indel_reads) +
        std::min(
            positive_or_zero(existence.alt_left_clip_reads),
            positive_or_zero(existence.alt_right_clip_reads));
    if (clip_insert_concordance != nullptr && clip_insert_concordance->pass) {
        mechanistic_reads += positive_or_zero(clip_insert_concordance->full_insert_reads);
        mechanistic_reads += std::min(
            positive_or_zero(clip_insert_concordance->left_clip_reads),
            positive_or_zero(clip_insert_concordance->right_clip_reads));
    }
    if (alt <= 0 || mechanistic_reads <= 0) {
        return 0.0;
    }
    const double count_part = count_signal(mechanistic_reads, 4.0);
    const double fraction_part = clamp01(
        static_cast<double>(mechanistic_reads) / static_cast<double>(alt));
    return std::sqrt(count_part * fraction_part);
}

double ref_conflict_signal(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation) {
    const int32_t alt = positive_or_zero(existence.alt_struct_reads);
    const int32_t ref = positive_or_zero(existence.ref_span_reads);
    double signal = (alt + ref) > 0
        ? static_cast<double>(ref) / static_cast<double>(alt + ref)
        : 0.0;
    signal = std::max(signal, 0.60 * count_signal(ref, 8.0));
    if (is_one_sided_segmentation_pass(segmentation) && ref > 0) {
        signal = std::max(signal, 0.35);
    }
    return clamp01(signal);
}

double event_signal(
    const EventExistenceEvidence& existence,
    double independent_signal) {
    const double support_signal = count_signal(positive_or_zero(existence.alt_struct_reads), 8.0);
    const double quality_signal = clamp01(static_cast<double>(existence.gq) / 60.0);
    return clamp01(
        (0.45 * support_signal) +
        (0.30 * independent_signal) +
        (0.15 * quality_signal) +
        (0.10 * clamp01(existence.af)));
}

double sequence_signal(const TEAlignmentEvidence& te_alignment) {
    if (!te_alignment.pass && te_alignment.qc_reason != "TE_ALIGNMENT_LOW_IDENTITY") {
        return 0.0;
    }
    double signal =
        (0.48 * clamp01(te_alignment.best_identity)) +
        (0.32 * clamp01(te_alignment.best_query_coverage)) +
        (0.20 * clamp01(te_alignment.cross_family_margin * 4.0));
    if (te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT") {
        signal += 0.15;
    } else if (te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY") {
        signal += 0.08;
    } else if (te_alignment.qc_reason == "PASS_INSERT_TE_ALIGNMENT_UNKNOWN") {
        signal -= 0.08;
    } else if (te_alignment.qc_reason == "TE_ALIGNMENT_LOW_IDENTITY") {
        signal -= 0.25;
    }
    if (te_alignment.sequence_model_label == "TE_MODEL_IN_DISTRIBUTION") {
        signal += 0.12 + (0.10 * clamp01(te_alignment.sequence_model_score));
    } else if (te_alignment.sequence_model_label == "TE_MODEL_EDGE") {
        signal -= 0.10;
    } else if (te_alignment.sequence_model_label == "TE_MODEL_OUTLIER") {
        signal -= 0.45;
    }
    if (te_alignment.annotation_confidence == "HIGH") {
        signal += 0.08;
    } else if (te_alignment.annotation_confidence == "LOW") {
        signal -= 0.12;
    }
    signal -= 0.20 * clamp01(te_alignment.annotation_residual_fraction);
    return clamp01(signal);
}

TeSequenceExplanation structure_explanation_for_certificate(
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment) {
    if (te_alignment.te_sequence_explanation.status !=
        TeAnnotationStatus::kUnavailable) {
        return te_alignment.te_sequence_explanation;
    }
    return explain_te_alignment_shadow(
        segmentation.insert_len,
        te_alignment.qc_reason,
        te_alignment.best_family,
        te_alignment.best_subfamily,
        te_alignment.best_identity,
        te_alignment.best_query_coverage,
        te_alignment.annotation_residual_fraction,
        te_alignment.annotation_masked_fraction,
        te_alignment.cross_family_margin,
        te_alignment.second_score);
}

double structure_ambiguity_width(
    const EventSegmentationEvidence& segmentation,
    const TeSequenceExplanation& explanation) {
    const double residual_fraction = segmentation.insert_len > 0
        ? std::clamp(
              static_cast<double>(std::max(0, explanation.unexplained_high_complexity_bp)) /
                  static_cast<double>(segmentation.insert_len),
              0.0,
              1.0)
        : 1.0;
    return 0.18 +
           (0.45 * (1.0 - std::clamp(explanation.structure_path_confidence, 0.0, 1.0))) +
           (0.30 * residual_fraction);
}

double boundary_signal(
    const EventSegmentationEvidence& segmentation,
    const BoundaryEvidence& boundary) {
    if (!segmentation.has_insert_seq) {
        return 0.0;
    }
    double signal = 0.20;
    if (segmentation.pair_valid && segmentation.has_left_flank && segmentation.has_right_flank) {
        signal += 0.35;
    } else if (segmentation.pair_valid) {
        signal += 0.20;
    } else if (is_one_sided_segmentation_pass(segmentation)) {
        signal += 0.10;
    }
    if (segmentation.has_left_flank != segmentation.has_right_flank) {
        signal += 0.05;
    }
    if (boundary.geometry_defined && boundary.canonical_pass) {
        signal += 0.35;
    } else if (boundary.geometry_defined && boundary.evidence_consistent) {
        signal += 0.22;
    } else if (!boundary.geometry_defined) {
        signal -= 0.10;
    }
    if (boundary.boundary_type == "TSD") {
        signal += 0.15;
    } else if (boundary.boundary_type == "BLUNT" ||
               boundary.boundary_type == "SMALL_DEL") {
        signal += 0.05;
    }
    return clamp01(signal);
}

double artifact_context_signal(
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    double ref_conflict) {
    double signal = 0.80 * ref_conflict;
    if (!segmentation.has_insert_seq) {
        signal = std::max(signal, 1.0);
    }
    if (te_alignment.sequence_model_label == "TE_MODEL_OUTLIER") {
        signal = std::max(signal, 0.90);
    } else if (te_alignment.sequence_model_label == "TE_MODEL_EDGE") {
        signal = std::max(signal, 0.45);
    }
    if (boundary.geometry_defined && !boundary.canonical_pass && !boundary.evidence_consistent) {
        signal = std::max(signal, 0.55);
    }
    if (te_alignment.annotation_confidence == "LOW") {
        signal = std::max(signal, 0.35);
    }
    return clamp01(signal);
}

void add_block(
    MechanisticEvidenceCertificate& certificate,
    const char* name,
    double raw_signal,
    double te_vs_artifact,
    double te_vs_non_te,
    double ambiguity_width) {
    MechanisticEvidenceBlock block;
    block.name = name;
    block.raw_signal = raw_signal;
    block.lower_log_lr_te_vs_artifact = te_vs_artifact;
    block.lower_log_lr_te_vs_non_te = te_vs_non_te;
    block.ambiguity_width = ambiguity_width;
    certificate.blocks.push_back(block);
}

}  // namespace

MechanisticEvidenceCertificate build_mechanistic_evidence_certificate(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    const ClipInsertConcordanceEvidence* clip_insert_concordance) {
    MechanisticEvidenceCertificate certificate;
    const double independent = mechanistic_read_signal(existence, clip_insert_concordance);
    const double event = event_signal(existence, independent);
    const double sequence = sequence_signal(te_alignment);
    const TeSequenceExplanation structure = structure_explanation_for_certificate(
        segmentation,
        te_alignment);
    const double boundary_mechanism = boundary_signal(segmentation, boundary);
    const double ref_conflict = ref_conflict_signal(existence, segmentation);
    const double artifact_context =
        artifact_context_signal(segmentation, te_alignment, boundary, ref_conflict);

    certificate.mechanistic_support_signal =
        clamp01((0.35 * independent) + (0.25 * event) +
                (0.25 * boundary_mechanism) + (0.15 * sequence));
    certificate.ref_conflict_signal = ref_conflict;
    certificate.artifact_context_signal = artifact_context;

    certificate.event_lower_log_lr =
        -0.75 + (3.00 * event) - (2.20 * ref_conflict);
    certificate.independent_lower_log_lr =
        -1.00 + (4.00 * independent) - (1.20 * ref_conflict);
    certificate.sequence_lower_log_lr =
        -1.40 + (4.20 * sequence) - (1.00 * artifact_context);
    certificate.structure_te_log_evidence = structure.te_structure_log_evidence;
    certificate.structure_nonte_log_evidence = structure.nonte_structure_log_evidence;
    certificate.structure_artifact_log_evidence =
        structure.artifact_structure_log_evidence;
    certificate.structure_lower_log_lr = std::clamp(
        certificate.structure_te_log_evidence -
            std::max(
                certificate.structure_nonte_log_evidence,
                certificate.structure_artifact_log_evidence),
        -4.0,
        5.0);
    certificate.boundary_lower_log_lr =
        -0.70 + (2.40 * boundary_mechanism) - (0.70 * ref_conflict);
    certificate.ref_conflict_lower_log_lr =
        -3.00 * ref_conflict - (0.80 * artifact_context);

    add_block(
        certificate,
        "event",
        event,
        certificate.event_lower_log_lr,
        0.45 * certificate.event_lower_log_lr,
        0.25 + (0.30 * (1.0 - event)));
    add_block(
        certificate,
        "independent",
        independent,
        certificate.independent_lower_log_lr,
        0.65 * certificate.independent_lower_log_lr,
        0.30 + (0.45 * (1.0 - independent)));
    add_block(
        certificate,
        "sequence",
        sequence,
        certificate.sequence_lower_log_lr,
        certificate.sequence_lower_log_lr,
        0.25 + (0.40 * (1.0 - sequence)));
    add_block(
        certificate,
        "structure",
        structure.structure_path_confidence,
        std::clamp(
            certificate.structure_te_log_evidence -
                certificate.structure_artifact_log_evidence,
            -4.0,
            5.0),
        std::clamp(
            certificate.structure_te_log_evidence -
                certificate.structure_nonte_log_evidence,
            -4.0,
            5.0),
        structure_ambiguity_width(segmentation, structure));
    add_block(
        certificate,
        "boundary",
        boundary_mechanism,
        certificate.boundary_lower_log_lr,
        0.80 * certificate.boundary_lower_log_lr,
        0.20 + (0.35 * (1.0 - boundary_mechanism)));
    add_block(
        certificate,
        "ref_conflict",
        ref_conflict,
        certificate.ref_conflict_lower_log_lr,
        0.50 * certificate.ref_conflict_lower_log_lr,
        0.20 + (0.50 * ref_conflict));

    certificate.ambiguity_width = 0.0;
    for (const auto& block : certificate.blocks) {
        certificate.ambiguity_width += block.ambiguity_width;
    }
    certificate.ambiguity_width /= std::max<size_t>(1, certificate.blocks.size());

    constexpr double kDependencyPenalty = 0.65;
    certificate.lower_log_bf_te_vs_artifact =
        certificate.event_lower_log_lr +
        certificate.independent_lower_log_lr +
        certificate.sequence_lower_log_lr +
        (0.90 * std::clamp(certificate.structure_lower_log_lr, -4.0, 5.0)) +
        certificate.boundary_lower_log_lr +
        certificate.ref_conflict_lower_log_lr -
        kDependencyPenalty;
    certificate.lower_log_bf_te_vs_non_te =
        (0.40 * certificate.event_lower_log_lr) +
        (0.70 * certificate.independent_lower_log_lr) +
        certificate.sequence_lower_log_lr +
        (0.90 * std::clamp(certificate.structure_lower_log_lr, -4.0, 5.0)) +
        (0.85 * certificate.boundary_lower_log_lr) +
        (0.35 * certificate.ref_conflict_lower_log_lr) -
        (0.60 * kDependencyPenalty);
    return certificate;
}

RobustMechanisticLfdrResult evaluate_robust_mechanistic_lfdr(
    const MechanisticEvidenceCertificate& certificate,
    const PriorInterval& prior,
    double target_q) {
    RobustMechanisticLfdrResult result;
    const double lower_log_bf = std::min(
        certificate.lower_log_bf_te_vs_artifact,
        certificate.lower_log_bf_te_vs_non_te);
    const double te_odds_lower =
        std::exp(std::clamp(lower_log_bf - certificate.ambiguity_width, -60.0, 60.0)) *
        std::max(prior.te_min, 1e-12);
    const double null_odds_upper =
        std::max(0.0, prior.artifact_max) + std::max(0.0, prior.non_te_max);
    result.worst_case_lfdr = null_odds_upper /
        std::max(1e-12, null_odds_upper + te_odds_lower);
    result.lfdr = result.worst_case_lfdr;
    result.qc = result.worst_case_lfdr <= std::clamp(target_q, 0.0, 1.0)
        ? "PASS_TE_LFDR"
        : "TE_LFDR_HIGH";
    return result;
}

std::string serialize_mechanistic_blocks(
    const MechanisticEvidenceCertificate& certificate) {
    if (certificate.blocks.empty()) {
        return "NA";
    }
    std::ostringstream out;
    for (size_t i = 0; i < certificate.blocks.size(); ++i) {
        const auto& block = certificate.blocks[i];
        if (i > 0) {
            out << ";";
        }
        out << block.name
            << ":raw=" << block.raw_signal
            << ",te_art=" << block.lower_log_lr_te_vs_artifact
            << ",te_non=" << block.lower_log_lr_te_vs_non_te
            << ",amb=" << block.ambiguity_width;
    }
    return out.str();
}

}  // namespace placer
