#ifndef PLACER_MECHANISTIC_EVIDENCE_H
#define PLACER_MECHANISTIC_EVIDENCE_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace placer {

struct BoundaryEvidence;
struct ClipInsertConcordanceEvidence;
struct EventExistenceEvidence;
struct EventSegmentationEvidence;
struct TEAlignmentEvidence;

struct MechanisticEvidenceBlock {
    std::string name = "NA";
    double raw_signal = 0.0;
    double lower_log_lr_te_vs_artifact = 0.0;
    double lower_log_lr_te_vs_non_te = 0.0;
    double ambiguity_width = 0.0;
};

struct PriorInterval {
    double te_min = 0.05;
    double te_max = 0.80;
    double non_te_min = 1e-6;
    double non_te_max = 0.30;
    double artifact_min = 0.50;
    double artifact_max = 0.95;
};

struct MechanisticEvidenceCertificate {
    std::vector<MechanisticEvidenceBlock> blocks;
    double event_lower_log_lr = 0.0;
    double sequence_lower_log_lr = 0.0;
    double structure_te_log_evidence = 0.0;
    double structure_nonte_log_evidence = 0.0;
    double structure_artifact_log_evidence = 0.0;
    double structure_lower_log_lr = 0.0;
    double boundary_lower_log_lr = 0.0;
    double independent_lower_log_lr = 0.0;
    double ref_conflict_lower_log_lr = 0.0;
    double lower_log_bf_te_vs_artifact = 0.0;
    double lower_log_bf_te_vs_non_te = 0.0;
    double mechanistic_support_signal = 0.0;
    double ref_conflict_signal = 0.0;
    double artifact_context_signal = 0.0;
    double ambiguity_width = 0.0;
};

struct RobustMechanisticLfdrResult {
    double lfdr = 1.0;
    double worst_case_lfdr = 1.0;
    std::string qc = "TE_LFDR_HIGH";
};

MechanisticEvidenceCertificate build_mechanistic_evidence_certificate(
    const EventExistenceEvidence& existence,
    const EventSegmentationEvidence& segmentation,
    const TEAlignmentEvidence& te_alignment,
    const BoundaryEvidence& boundary,
    const ClipInsertConcordanceEvidence* clip_insert_concordance);

RobustMechanisticLfdrResult evaluate_robust_mechanistic_lfdr(
    const MechanisticEvidenceCertificate& certificate,
    const PriorInterval& prior,
    double target_q);

std::string serialize_mechanistic_blocks(
    const MechanisticEvidenceCertificate& certificate);

}  // namespace placer

#endif  // PLACER_MECHANISTIC_EVIDENCE_H
