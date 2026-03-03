#include "decision_policy.h"
#include "pipeline.h"

#include <cassert>

int main() {
    using namespace placer;

    PipelineConfig config;
    config.short_ins_enable = true;
    config.short_ins_min_len = 35;
    config.short_ins_max_len = 300;
    config.short_ins_min_reads = 2;
    config.short_ins_kmer_relax_identity = 0.15;
    config.assembly_min_identity_est = 0.55;

    ComponentCall component;
    component.split_sa_read_indices = {0};
    component.insertion_read_indices = {1};
    BreakpointCandidate bp1;
    bp1.pos = 500;
    bp1.read_index = 0;
    bp1.class_mask = kCandidateSplitSaSupplementary;
    BreakpointCandidate bp2;
    bp2.pos = 506;
    bp2.read_index = 1;
    bp2.class_mask = kCandidateLongInsertion;
    component.breakpoint_candidates = {bp1, bp2};

    AssemblyCall assembly;
    assembly.qc_pass = true;
    assembly.identity_est = 0.42;  // below assembly_min_identity_est, above short-ins relaxed identity
    assembly.consensus_len = 120;

    ClusterTECall te_empty;
    const auto decision_short_ins = evaluate_post_assembly_te_decision(
        te_empty,
        assembly,
        component,
        config);
    assert(decision_short_ins.pass);
    assert(decision_short_ins.qc == "PASS_SHORT_INS_RESCUE");

    PipelineConfig config_disabled = config;
    config_disabled.short_ins_enable = false;
    const auto decision_disabled = evaluate_post_assembly_te_decision(
        te_empty,
        assembly,
        component,
        config_disabled);
    assert(!decision_disabled.pass);
    assert(decision_disabled.qc == "FAIL_NO_TE_LABEL");

    return 0;
}
