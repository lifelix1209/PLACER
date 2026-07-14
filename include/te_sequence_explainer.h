#ifndef PLACER_TE_SEQUENCE_EXPLAINER_H
#define PLACER_TE_SEQUENCE_EXPLAINER_H

#include <cstdint>
#include <string>
#include <vector>

namespace placer {

enum class TeAnnotationStatus : uint8_t {
    kUnavailable = 0,
    kEmpty = 1,
    kTooShort = 2,
    kNoCandidate = 3,
    kResolved = 4,
    kFamilyOnly = 5,
    kUnknownTe = 6,
    kAmbiguous = 7,
    kNonTeLike = 8
};

struct TePathResidual {
    int32_t unexplained_high_complexity_bases = 0;
    int32_t edit_distance = 0;
    int32_t family_conflicts = 0;
    int32_t subfamily_conflicts = 0;
    int32_t orientation_conflicts = 0;
    int32_t segment_breaks = 0;
    int32_t low_complexity_only_bases = 0;
};

struct TeSequenceSegment {
    std::string kind;
    int32_t query_start = 0;
    int32_t query_end = 0;
    std::string label = "NA";
};

struct TeSequenceExplanation {
    TeAnnotationStatus status = TeAnnotationStatus::kUnavailable;
    TePathResidual residual;
    std::string family = "NA";
    std::string subfamily = "NA";
    std::vector<TeSequenceSegment> path;
    double te_structure_log_evidence = 0.0;
    double nonte_structure_log_evidence = 0.0;
    double artifact_structure_log_evidence = 0.0;
    double structure_path_confidence = 0.0;
    double polyA_posterior = 0.0;
    double transduction_posterior = 0.0;
    double te_core_coverage = 0.0;
    int32_t unexplained_high_complexity_bp = 0;
};

TeSequenceExplanation explain_te_sequence_structure(
    const std::string& insert_seq,
    const std::string& qc_reason,
    const std::string& family,
    const std::string& subfamily,
    double best_identity,
    double effective_query_coverage,
    double annotation_residual_fraction,
    double annotation_masked_fraction,
    double cross_family_margin,
    double second_score,
    const std::string& sequence_model_label,
    double sequence_model_score);

TeSequenceExplanation explain_te_alignment_shadow(
    int32_t insert_len,
    const std::string& qc_reason,
    const std::string& family,
    const std::string& subfamily,
    double best_identity,
    double effective_query_coverage,
    double annotation_residual_fraction,
    double annotation_masked_fraction,
    double cross_family_margin,
    double second_score);

const char* te_annotation_status_name(TeAnnotationStatus status);
std::string serialize_te_path_residual(const TePathResidual& residual);
std::string serialize_te_sequence_path(const TeSequenceExplanation& explanation);

}  // namespace placer

#endif  // PLACER_TE_SEQUENCE_EXPLAINER_H
