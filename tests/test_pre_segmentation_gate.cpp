#include <cassert>
#include <string>

#define private public
#include "pipeline.h"
#undef private

int main() {
    using namespace placer;

    PipelineConfig config;
    Pipeline pipeline(config, nullptr);

    {
        EventReadEvidence evidence;
        EventConsensus consensus;
        consensus.left_anchor_input_reads = 0;
        consensus.right_anchor_input_reads = 1;

        assert(
            pipeline.pre_segmentation_gate_reason(evidence, consensus) ==
            "PRESEG_NO_BILATERAL_ANCHOR");
    }

    {
        EventReadEvidence evidence;
        evidence.alt_split_reads = 0;
        evidence.alt_indel_reads = 0;

        EventConsensus consensus;
        consensus.left_anchor_input_reads = 1;
        consensus.right_anchor_input_reads = 1;
        consensus.full_context_input_reads = 0;

        assert(
            pipeline.pre_segmentation_gate_reason(evidence, consensus) ==
            "PRESEG_NO_PRECISE_OR_FULL_CONTEXT");
    }

    {
        EventReadEvidence evidence;
        evidence.alt_split_reads = 1;

        EventConsensus consensus;
        consensus.left_anchor_input_reads = 1;
        consensus.right_anchor_input_reads = 1;

        assert(pipeline.pre_segmentation_gate_reason(evidence, consensus).empty());
    }

    {
        EventReadEvidence evidence;
        evidence.alt_split_reads = 0;
        evidence.alt_indel_reads = 0;

        EventConsensus consensus;
        consensus.input_event_reads = 4;
        consensus.partial_context_input_reads = 4;
        consensus.left_anchor_input_reads = 2;
        consensus.right_anchor_input_reads = 2;

        assert(pipeline.pre_segmentation_gate_reason(evidence, consensus).empty());
    }

    {
        EventReadEvidence evidence;
        evidence.alt_split_reads = 0;
        evidence.alt_indel_reads = 0;

        EventConsensus consensus;
        consensus.input_event_reads = 3;
        consensus.partial_context_input_reads = 3;
        consensus.left_anchor_input_reads = 2;
        consensus.right_anchor_input_reads = 2;

        assert(
            pipeline.pre_segmentation_gate_reason(evidence, consensus) ==
            "PRESEG_NO_PRECISE_OR_FULL_CONTEXT");
    }

    {
        Pipeline::HypothesisSummary summary;
        summary.bp_left = 1000;
        summary.bp_right = 1005;
        summary.alt_split_reads = 0;
        summary.alt_indel_reads = 0;
        summary.alt_struct_reads = 2;

        Pipeline::ConsensusInputCounts inputs;
        inputs.full_context_input_reads = 0;
        inputs.left_anchor_input_reads = 1;
        inputs.right_anchor_input_reads = 0;

        EventReadEvidence evidence;
        evidence.bp_left = summary.bp_left;
        evidence.bp_right = summary.bp_right;
        evidence.alt_split_reads = summary.alt_split_reads;
        evidence.alt_indel_reads = summary.alt_indel_reads;
        evidence.alt_struct_reads = summary.alt_struct_reads;
        summary.event_evidence = evidence;

        const auto validator = pipeline.collect_hypothesis_validator_evidence(
            summary,
            inputs,
            1002);
        assert(!validator.feasible_for_expensive_stage);
        assert(validator.qc_reason == "VALIDATOR_NO_BILATERAL_ANCHOR");
    }

    {
        Pipeline::HypothesisSummary summary;
        summary.bp_left = 1000;
        summary.bp_right = 1005;
        summary.alt_split_reads = 1;
        summary.alt_indel_reads = 0;

        Pipeline::ConsensusInputCounts inputs;
        inputs.left_anchor_input_reads = 1;
        inputs.right_anchor_input_reads = 1;

        EventReadEvidence evidence;
        evidence.bp_left = summary.bp_left;
        evidence.bp_right = summary.bp_right;
        evidence.alt_split_reads = summary.alt_split_reads;
        evidence.alt_indel_reads = summary.alt_indel_reads;
        evidence.alt_struct_reads = summary.alt_struct_reads;
        summary.event_evidence = evidence;

        const auto validator = pipeline.collect_hypothesis_validator_evidence(
            summary,
            inputs,
            1002);
        assert(validator.feasible_for_expensive_stage);
        assert(validator.precise_support == 1);
        assert(validator.breakpoint_width == 5);
    }

    {
        Pipeline::HypothesisSummary summary;
        summary.bp_left = 1000;
        summary.bp_right = 1005;
        summary.alt_split_reads = 1;
        summary.alt_indel_reads = 0;

        Pipeline::ConsensusInputCounts inputs;
        inputs.left_anchor_input_reads = 1;
        inputs.right_anchor_input_reads = 0;

        EventReadEvidence evidence;
        evidence.bp_left = summary.bp_left;
        evidence.bp_right = summary.bp_right;
        evidence.alt_split_reads = summary.alt_split_reads;
        evidence.alt_indel_reads = summary.alt_indel_reads;
        evidence.alt_struct_reads = summary.alt_struct_reads;
        summary.event_evidence = evidence;

        const auto validator = pipeline.collect_hypothesis_validator_evidence(
            summary,
            inputs,
            1002);
        assert(validator.feasible_for_expensive_stage);
        assert(validator.qc_reason == "VALIDATOR_PASS");
    }

    {
        Pipeline::HypothesisSummary summary;
        summary.bp_left = 1000;
        summary.bp_right = 1020;
        summary.alt_struct_reads = 8;

        Pipeline::ConsensusInputCounts inputs;
        inputs.partial_context_input_reads = 4;
        inputs.left_anchor_input_reads = 2;
        inputs.right_anchor_input_reads = 2;
        inputs.input_event_reads = 4;

        EventReadEvidence evidence;
        evidence.bp_left = summary.bp_left;
        evidence.bp_right = summary.bp_right;
        evidence.alt_split_reads = summary.alt_split_reads;
        evidence.alt_indel_reads = summary.alt_indel_reads;
        evidence.alt_struct_reads = summary.alt_struct_reads;
        summary.event_evidence = evidence;

        const auto validator = pipeline.collect_hypothesis_validator_evidence(
            summary,
            inputs,
            1002);
        assert(validator.feasible_for_expensive_stage);
        assert(validator.qc_reason == "VALIDATOR_PASS");
    }

    {
        Pipeline::HypothesisSummary summary;
        summary.bp_left = 1000;
        summary.bp_right = 1020;
        summary.alt_struct_reads = 6;

        Pipeline::ConsensusInputCounts inputs;
        inputs.partial_context_input_reads = 3;
        inputs.left_anchor_input_reads = 2;
        inputs.right_anchor_input_reads = 2;
        inputs.input_event_reads = 3;

        EventReadEvidence evidence;
        evidence.bp_left = summary.bp_left;
        evidence.bp_right = summary.bp_right;
        evidence.alt_split_reads = summary.alt_split_reads;
        evidence.alt_indel_reads = summary.alt_indel_reads;
        evidence.alt_struct_reads = summary.alt_struct_reads;
        summary.event_evidence = evidence;

        const auto validator = pipeline.collect_hypothesis_validator_evidence(
            summary,
            inputs,
            1002);
        assert(!validator.feasible_for_expensive_stage);
        assert(validator.qc_reason == "VALIDATOR_NO_PRECISE_OR_FULL_CONTEXT");
    }

    return 0;
}
