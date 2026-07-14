#ifndef PLACER_EVENT_EXPLANATION_H
#define PLACER_EVENT_EXPLANATION_H

#include <cstdint>
#include <string>
#include <vector>

namespace placer {

enum class SequenceSegmentKind : uint8_t {
    kTeTemplate = 0,
    kLocalReference = 1,
    kTransduction = 2,
    kPolyA = 3,
    kTsd = 4,
    kLowComplexity = 5,
    kUnexplained = 6
};

enum class ExplanationKind : uint8_t {
    kReference = 0,
    kInsertionNonTe = 1,
    kTe = 2,
    kArtifact = 3
};

struct SequenceExplanationSegment {
    SequenceSegmentKind kind = SequenceSegmentKind::kUnexplained;
    int32_t query_start = 0;
    int32_t query_end = 0;
    std::string label;
    std::string family = "NA";
    std::string subfamily = "NA";
    std::string strand = "NA";
    int32_t edit_distance = 0;
    int32_t aligned_bases = 0;
    double identity = 0.0;
};

struct ExplanationResidual {
    int32_t structural_conflicts = 0;
    int32_t missing_required_components = 0;
    int32_t unexplained_high_complexity_bases = 0;
    int32_t breakpoint_disagreement_bp = 0;
    int32_t read_assignment_conflicts = 0;
    int32_t reference_counterevidence = 0;
    int32_t artifact_evidence = 0;
    int32_t label_ambiguity = 0;
    int32_t path_complexity = 0;
    std::vector<std::string> notes;
};

struct EventExplanation {
    ExplanationKind kind = ExplanationKind::kReference;
    ExplanationResidual residual;
    std::vector<SequenceExplanationSegment> path;
    std::string family = "NA";
    std::string subfamily = "NA";
    std::string status = "UNSET";
};

struct ExplanationDecision {
    EventExplanation best;
    std::vector<EventExplanation> alternatives;
    std::string final_qc = "NO_CALL_INCOMPLETE";
    bool emit_te_call = false;
    bool emit_unknown_te = false;
    bool emit_evidence_te_call = false;
};

const char* explanation_kind_name(ExplanationKind kind);
std::string serialize_residual(const ExplanationResidual& residual);
std::string serialize_explanation_path(const EventExplanation& explanation);

bool dominates_primary_residuals(
    const ExplanationResidual& lhs,
    const ExplanationResidual& rhs);

ExplanationDecision compare_event_explanations(
    const std::vector<EventExplanation>& explanations,
    bool closed_breakpoints);

}  // namespace placer

#endif  // PLACER_EVENT_EXPLANATION_H
