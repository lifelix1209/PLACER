#ifdef NDEBUG
#undef NDEBUG
#endif
#include "decision_policy.h"
#define private public
#include "pipeline.h"
#undef private

#include <cassert>
#include <string>

int main() {
    using namespace placer;

    EventExistenceEvidence existence;
    existence.score = 2.0;

    EventSegmentationEvidence segmentation;
    segmentation.has_insert_seq = true;
    segmentation.score = 1.5;

    BoundaryEvidence boundary;
    boundary.score = 0.5;

    {
        TEAlignmentEvidence te;
        te.best_family = "Gypsy";
        te.best_identity = 0.72;
        te.best_query_coverage = 0.68;
        te.cross_family_margin = 0.08;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY";

        const JointDecisionResult result =
            evaluate_joint_hypotheses(existence, segmentation, te, boundary);

        assert(result.emit_te_call);
        assert(!result.emit_unknown_te);
        assert(result.final_qc == "PASS_FINAL_TE_CALL");
    }

    {
        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.72;
        te.best_query_coverage = 0.72;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        const JointDecisionResult result =
            evaluate_joint_hypotheses(existence, segmentation, te, boundary);

        assert(result.emit_te_call);
        assert(result.emit_unknown_te);
        assert(result.final_qc == "PASS_FINAL_TE_CALL_UNKNOWN");
    }

    {
        EventExistenceEvidence alt_only_existence;
        alt_only_existence.alt_struct_reads = 20;
        alt_only_existence.ref_span_reads = 0;
        alt_only_existence.score = 1.85;

        EventSegmentationEvidence long_insert_segmentation;
        long_insert_segmentation.has_insert_seq = true;
        long_insert_segmentation.score = 1.18;

        BoundaryEvidence canonical_boundary;
        canonical_boundary.score = 1.0;

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.540205;
        te.best_query_coverage = 0.516461;
        te.cross_family_margin = 0.247217;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                alt_only_existence,
                long_insert_segmentation,
                te,
                canonical_boundary);

        assert(result.emit_te_call);
        assert(result.emit_unknown_te);
        assert(result.final_qc == "PASS_FINAL_TE_CALL_UNKNOWN");
    }

    {
        EventExistenceEvidence strong_existence;
        strong_existence.alt_struct_reads = 12;
        strong_existence.ref_span_reads = 1;
        strong_existence.score = 3.0;

        EventSegmentationEvidence occupied_segmentation;
        occupied_segmentation.has_insert_seq = true;
        occupied_segmentation.score = 0.20;

        BoundaryEvidence noncanonical_boundary;
        noncanonical_boundary.score = 0.25;

        TEAlignmentEvidence te;
        te.best_family = "PIF-ISL2EU";
        te.best_subfamily = "PIF-ISL2EU-4";
        te.best_identity = 0.543036;
        te.best_query_coverage = 0.987654;
        te.cross_family_margin = 0.0495866;
        te.pass = false;
        te.qc_reason = "TE_ALIGNMENT_LOW_IDENTITY";

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                strong_existence,
                occupied_segmentation,
                te,
                noncanonical_boundary);

        assert(result.emit_te_call);
        assert(result.emit_unknown_te);
        assert(result.final_qc == "PASS_FINAL_TE_CALL_UNKNOWN");
    }

    {
        Pipeline pipeline(PipelineConfig{}, nullptr);

        ComponentCall component;
        component.chrom = "8";
        component.tid = 0;
        component.bin_start = 25088312;
        component.bin_end = 25092851;

        EventReadEvidence event_evidence;
        event_evidence.alt_struct_reads = 12;
        event_evidence.ref_span_reads = 1;

        EventConsensus consensus;
        consensus.consensus_len = 616;
        consensus.qc_reason = "PASS_EVENT_CONSENSUS";

        EventSegmentation segmentation;
        segmentation.left_ref_end = 25090631;
        segmentation.right_ref_start = 25090685;
        segmentation.insert_seq = std::string(486, 'A');
        segmentation.left_flank_align_len = 50;
        segmentation.right_flank_align_len = 80;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "PIF-ISL2EU";
        te.best_subfamily = "PIF-ISL2EU-4";
        te.best_identity = 0.543036;
        te.best_query_coverage = 0.987654;
        te.cross_family_margin = 0.0495866;
        te.qc_reason = "TE_ALIGNMENT_LOW_IDENTITY";

        GenotypeCall genotype;
        genotype.genotype = "0/1";
        genotype.af = 0.923077;
        genotype.gq = 99;

        FinalBoundaryDecision boundary;
        boundary.boundary_type = "NONCANONICAL";
        boundary.boundary_len = 54;
        boundary.qc = "PASS_BOUNDARY_NONCANONICAL_CONSISTENT";

        FinalTeAcceptanceDecision acceptance;
        acceptance.pass = true;
        acceptance.qc = "PASS_FINAL_TE_CALL_UNKNOWN";

        const FinalCall call = pipeline.emit_final_te_call(
            component,
            event_evidence,
            consensus,
            segmentation,
            te,
            genotype,
            boundary,
            acceptance);

        assert(call.te_name == "UNKNOWN");
        assert(call.family == "UNKNOWN");
        assert(call.subfamily == "UNKNOWN");
        assert(call.te_qc == "TE_ALIGNMENT_LOW_IDENTITY");
        assert(call.final_qc == "PASS_FINAL_TE_CALL_UNKNOWN|TE_ALIGNMENT_LOW_IDENTITY|PASS_BOUNDARY_NONCANONICAL_CONSISTENT|PASS_EVENT_CONSENSUS|PASS_EVENT_SEGMENTATION");
    }

    return 0;
}
