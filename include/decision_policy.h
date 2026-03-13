#ifndef PLACER_DECISION_POLICY_H
#define PLACER_DECISION_POLICY_H

#include <cstdint>
#include <string>

namespace placer {

struct ClusterTECall;
struct AssemblyCall;
struct ComponentCall;
struct PipelineConfig;

struct InsertionEvidence {
    int32_t tier = 3;
    int32_t support_reads = 0;
};

struct InsertionAcceptanceParams {
    bool emit_low_confidence_calls = false;
    int32_t low_conf_min_support_reads = 2;
    int32_t low_conf_max_tier = 2;
};

struct InsertionAcceptanceDecision {
    bool pass = false;
    bool low_confidence = false;
};

InsertionAcceptanceDecision evaluate_insertion_acceptance(
    const InsertionEvidence& evidence,
    const InsertionAcceptanceParams& params);

struct TeOpenSetInput {
    bool te_gate_pass = false;
    bool te_gate_uncertain_path = false;
    bool force_te_uncertain = false;
    bool force_non_te = false;
    bool has_proxy_signal = false;
    bool proxy_certain = false;
    double te_confidence_prob = 0.0;
};

struct TeOpenSetParams {
    double conf_certain_min = 0.85;
    double conf_uncertain_min = 0.35;
};

struct TeOpenSetDecision {
    std::string status = "NON_TE";
    bool promote_top1_name = false;
    bool add_te_gate_fail_tag = false;
    std::string qc_proxy_tag;
};

TeOpenSetDecision classify_te_open_set(
    const TeOpenSetInput& input,
    const TeOpenSetParams& params);

struct BreakpointConsistencyInput {
    int32_t split_count = 0;
    int32_t clip_count = 0;
    int32_t indel_count = 0;
    double all_breakpoint_mad = 0.0;
    double split_breakpoint_mad = 0.0;
    double split_like_breakpoint_mad = 0.0;
    double clip_breakpoint_mad = 0.0;
    double indel_breakpoint_mad = 0.0;
    double split_clip_core_delta = 0.0;
};

struct BreakpointConsistencyParams {
    double breakpoint_mad_max = 80.0;
    double clip_breakpoint_mad_max = 120.0;
    double split_clip_core_delta_max = 150.0;
    double tight_split_like_breakpoint_mad_max = 40.0;
    int32_t min_split_like_count = 1;
    int32_t min_clip_count = 2;
};

struct BreakpointConsistencyDecision {
    bool pass_breakpoint_mad = true;
    bool pass_split_clip_consistency = true;
    bool severe_inconsistency = false;
    double placeability_penalty = 0.0;
};

BreakpointConsistencyDecision evaluate_breakpoint_consistency(
    const BreakpointConsistencyInput& input,
    const BreakpointConsistencyParams& params);

struct PostAssemblyTeDecision {
    bool pass = false;
    bool force_te_uncertain = false;
    std::string qc = "FAIL_TE_CLASSIFICATION";
};

PostAssemblyTeDecision evaluate_post_assembly_te_decision(
    const ClusterTECall& te_call,
    const AssemblyCall& assembly,
    const ComponentCall& component,
    const PipelineConfig& config);

}  // namespace placer

#endif  // PLACER_DECISION_POLICY_H
