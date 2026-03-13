#include "decision_policy.h"
#include "pipeline.h"

#include <cassert>

int main() {
    using namespace placer;

    PipelineConfig config;
    config.te_no_softclip_min_reads = 2;
    config.te_no_softclip_identity_min = 0.20;
    config.te_one_sided_breakpoint_mad_max = 40.0;
    config.short_ins_enable = false;

    ComponentCall component_pass;
    component_pass.split_sa_read_indices = {0, 1};
    BreakpointCandidate bp1;
    bp1.pos = 100;
    bp1.read_index = 0;
    bp1.class_mask = kCandidateSplitSaSupplementary;
    BreakpointCandidate bp2;
    bp2.pos = 108;
    bp2.read_index = 1;
    bp2.class_mask = kCandidateSplitSaSupplementary;
    component_pass.breakpoint_candidates = {bp1, bp2};

    AssemblyCall assembly_pass;
    assembly_pass.qc_pass = true;
    assembly_pass.identity_est = 0.60;
    assembly_pass.consensus_len = 600;

    ClusterTECall te_empty;
    const auto decision_pass = evaluate_post_assembly_te_decision(
        te_empty,
        assembly_pass,
        component_pass,
        config);
    assert(decision_pass.pass);
    assert(decision_pass.qc == "PASS_SPLIT_INDEL_RESCUE");

    ComponentCall component_fail;
    component_fail.split_sa_read_indices = {0, 1};
    BreakpointCandidate bp3;
    bp3.pos = 100;
    bp3.read_index = 0;
    bp3.class_mask = kCandidateSplitSaSupplementary;
    BreakpointCandidate bp4;
    bp4.pos = 240;
    bp4.read_index = 1;
    bp4.class_mask = kCandidateSplitSaSupplementary;
    component_fail.breakpoint_candidates = {bp3, bp4};

    AssemblyCall assembly_fail;
    assembly_fail.qc_pass = true;
    assembly_fail.identity_est = 0.30;
    assembly_fail.consensus_len = 600;

    const auto decision_fail = evaluate_post_assembly_te_decision(
        te_empty,
        assembly_fail,
        component_fail,
        config);
    assert(!decision_fail.pass);
    assert(decision_fail.qc == "FAIL_SPLIT_INDEL_INCONSISTENT");

    ClusterTECall te_low_identity;
    te_low_identity.te_name = "L1:L1HS";
    te_low_identity.fragment_count = 2;
    te_low_identity.vote_fraction = 1.0;
    te_low_identity.median_identity = 0.05;

    AssemblyCall assembly_low_identity;
    assembly_low_identity.qc_pass = true;
    assembly_low_identity.identity_est = 0.98;
    assembly_low_identity.consensus_len = 900;

    const auto decision_low_identity = evaluate_post_assembly_te_decision(
        te_low_identity,
        assembly_low_identity,
        component_pass,
        config);
    assert(!decision_low_identity.pass);
    assert(decision_low_identity.qc == "FAIL_LOW_TE_SUPPORT");

    return 0;
}
