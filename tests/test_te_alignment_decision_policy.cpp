#ifdef NDEBUG
#undef NDEBUG
#endif
#include "decision_policy.h"
#define private public
#include "pipeline.h"
#undef private

#include <cassert>
#include <string>

namespace {

void set_precise_bam_breakdown(
    placer::EventExistenceEvidence& evidence,
    int32_t split_reads,
    int32_t indel_reads,
    int32_t left_clip_reads,
    int32_t right_clip_reads) {
    evidence.alt_split_reads = split_reads;
    evidence.alt_indel_reads = indel_reads;
    evidence.alt_left_clip_reads = left_clip_reads;
    evidence.alt_right_clip_reads = right_clip_reads;
}

}  // namespace

int main() {
    using namespace placer;

    EventExistenceEvidence existence;
    existence.score = 2.0;
    existence.alt_struct_reads = 8;
    existence.alt_split_reads = 2;
    existence.alt_indel_reads = 1;
    existence.alt_left_clip_reads = 3;
    existence.alt_right_clip_reads = 4;

    EventSegmentationEvidence segmentation;
    segmentation.has_insert_seq = true;
    segmentation.has_left_flank = true;
    segmentation.has_right_flank = true;
    segmentation.pair_valid = true;
    segmentation.insert_len = 500;
    segmentation.score = 1.5;

    BoundaryEvidence boundary;
    boundary.geometry_defined = true;
    boundary.canonical_pass = true;
    boundary.evidence_consistent = true;
    boundary.score = 0.5;

    {
        TEAlignmentEvidence te;
        assert(te.annotation_confidence == "NA");
        assert(te.annotation_class == "NA");
        assert(te.annotation_order == "NA");
        assert(te.annotation_intervals == "NA");
        assert(te.annotation_residual_fraction == 0.0);
        assert(te.annotation_masked_fraction == 0.0);
        te.best_family = "Gypsy";
        te.best_identity = 0.72;
        te.best_query_coverage = 0.68;
        te.cross_family_margin = 0.08;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;
        te.te_sequence_explanation = explain_te_sequence_structure(
            std::string(35, 'C') +
                "GTCAGTCAGATCGATGCTAGCTAGGACCTAGTCAGATCGATGCTAGCTA" +
                "GATCGATGCTAGCTAGGACCTAGTCAGATCGATGCTAGCTAGTCAGTCAG",
            te.qc_reason,
            te.best_family,
            te.best_subfamily,
            te.best_identity,
            0.10,
            0.90,
            0.0,
            te.cross_family_margin,
            te.second_score,
            "TE_MODEL_OUTLIER",
            -0.80);

        const JointDecisionResult result =
            evaluate_joint_hypotheses(existence, segmentation, te, boundary);

        assert(!result.emit_te_call);
        assert(result.emit_evidence_te_call);
        assert(result.final_qc == "TE_AMBIGUOUS");
        assert(result.lfdr_qc == "TE_LFDR_HIGH");
        assert(result.robust_mechanistic_qc == "TE_LFDR_HIGH");
    }

    {
        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.72;
        te.best_query_coverage = 0.72;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(existence, segmentation, te, boundary);

        assert(result.emit_te_call);
        assert(result.emit_unknown_te);
        assert(result.final_qc == "PASS_TE_CLOSED");
        assert(result.lfdr_qc == "TE_LFDR_HIGH");
    }

    {
        EventExistenceEvidence strong_existence;
        strong_existence.alt_struct_reads = 124;
        set_precise_bam_breakdown(strong_existence, 20, 20, 40, 44);
        strong_existence.ref_span_reads = 89;
        strong_existence.gq = 99;
        strong_existence.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 1355;
        complete_segmentation.score = 1.5;

        BoundaryEvidence tsd_boundary;
        tsd_boundary.geometry_defined = true;
        tsd_boundary.canonical_pass = true;
        tsd_boundary.evidence_consistent = true;
        tsd_boundary.boundary_type = "TSD";
        tsd_boundary.boundary_len = 40;
        tsd_boundary.score = 1.0;

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.653205;
        te.best_query_coverage = 0.99262;
        te.cross_family_margin = 0.122309;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_EDGE";
        te.sequence_model_score = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                strong_existence,
                complete_segmentation,
                te,
                tsd_boundary);

        assert(!result.emit_te_call);
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence strong_existence;
        strong_existence.alt_struct_reads = 28;
        set_precise_bam_breakdown(strong_existence, 4, 4, 10, 10);
        strong_existence.ref_span_reads = 7;
        strong_existence.gq = 99;
        strong_existence.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 164;
        complete_segmentation.score = 1.5;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 36;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.60;
        te.best_query_coverage = 0.841463;
        te.cross_family_margin = 0.504878;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                strong_existence,
                complete_segmentation,
                te,
                small_del_boundary);

        assert(result.emit_te_call);
        assert(result.emit_unknown_te);
        assert(result.final_qc == "PASS_TE_CLOSED");
        assert(result.lfdr_qc == "TE_LFDR_HIGH");
    }

    {
        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.72;
        te.best_query_coverage = 0.72;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_OUTLIER";
        te.sequence_model_score = -0.45;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(existence, segmentation, te, boundary);

        assert(!result.emit_te_call);
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.72;
        te.best_query_coverage = 0.72;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(existence, segmentation, te, boundary);

        assert(result.emit_te_call);
        assert(result.emit_unknown_te);
        assert(result.final_qc == "PASS_TE_CLOSED");
        assert(result.lfdr_qc == "TE_LFDR_HIGH");
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

        assert(!result.emit_te_call);
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

        assert(!result.emit_te_call);
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
        te.annotation_confidence = "LOW";
        te.annotation_class = "DNA";
        te.annotation_order = "TIR";
        te.annotation_intervals = "q=0-480,t=15-495,id=0.543,cov=0.988";
        te.annotation_residual_fraction = 0.012346;
        te.annotation_masked_fraction = 0.0;
        te.qc_reason = "TE_ALIGNMENT_LOW_IDENTITY";

        GenotypeCall genotype;
        genotype.genotype = "0/1";
        genotype.af = 0.923077;
        genotype.gq = 99;

        FinalBoundaryDecision boundary;
        boundary.boundary_type = "NONCANONICAL";
        boundary.boundary_len = 54;
        boundary.qc = "PASS_BOUNDARY_NONCANONICAL_CONSISTENT";

        JointDecisionResult joint;
        joint.emit_unknown_te = true;
        joint.final_qc = "PASS_TE_IMPRECISE";
        joint.best_explanation = "TE";
        joint.explanation_residual = "structural=0;missing=0";
        joint.explanation_path = "NA";

        const FinalCall call = pipeline.emit_final_te_call(
            component,
            event_evidence,
            consensus,
            segmentation,
            te,
            genotype,
            boundary,
            joint);

        assert(call.te_name == "UNKNOWN");
        assert(call.family == "UNKNOWN");
        assert(call.subfamily == "UNKNOWN");
        assert(!call.family_committed);
        assert(call.te_annotation_confidence == "LOW");
        assert(call.te_annotation_class == "DNA");
        assert(call.te_annotation_order == "TIR");
        assert(call.te_annotation_intervals == "q=0-480,t=15-495,id=0.543,cov=0.988");
        assert(call.te_annotation_residual_fraction == 0.012346);
        assert(call.te_annotation_masked_fraction == 0.0);
        assert(call.te_qc == "TE_ALIGNMENT_LOW_IDENTITY");
        assert(call.final_qc == "PASS_TE_IMPRECISE|TE_ALIGNMENT_LOW_IDENTITY|PASS_BOUNDARY_NONCANONICAL_CONSISTENT|PASS_EVENT_CONSENSUS|PASS_EVENT_SEGMENTATION");
        assert(call.best_explanation == "TE");
        assert(call.explanation_residual == "structural=0;missing=0");
        assert(call.explanation_path == "NA");
    }

    {
        Pipeline pipeline(PipelineConfig{}, nullptr);

        ComponentCall component;
        component.chrom = "3";
        component.tid = 0;
        component.anchor_pos = 44440766;
        component.bin_start = 44439000;
        component.bin_end = 44443000;

        EventReadEvidence event_evidence;
        event_evidence.bp_left = 44440762;
        event_evidence.bp_right = 44440771;
        event_evidence.alt_struct_reads = 40;
        event_evidence.ref_span_reads = 1;

        EventConsensus consensus;
        consensus.consensus_len = 431;
        consensus.qc_reason = "PASS_EVENT_CONSENSUS";

        EventSegmentation segmentation;
        segmentation.left_ref_end = 44440762;
        segmentation.right_ref_start = 44440771;
        segmentation.insert_seq = std::string(431, 'A');
        segmentation.left_flank_align_len = 91;
        segmentation.right_flank_align_len = 88;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "hAT-Ac";
        te.best_subfamily = "hAT-Ac-10";
        te.best_identity = 0.655;
        te.best_query_coverage = 0.323;
        te.cross_family_margin = 0.02;
        te.pass = false;
        te.qc_reason = "TE_ALIGNMENT_LOW_COVERAGE";
        te.sequence_model_label = "TE_MODEL_OUTLIER";

        GenotypeCall genotype;
        genotype.genotype = "0/1";
        genotype.af = 0.97561;
        genotype.gq = 99;

        FinalBoundaryDecision boundary;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 10;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        JointDecisionResult joint;
        joint.emit_structural_event_call = true;
        joint.emit_unknown_te = true;
        joint.final_qc = "PASS_STRUCTURAL_INSERTION";
        joint.best_explanation = "INSERTION_NON_TE";
        joint.explanation_residual = "structural=0;missing=0";
        joint.explanation_path = "NA";

        const FinalCall call = pipeline.emit_final_te_call(
            component,
            event_evidence,
            consensus,
            segmentation,
            te,
            genotype,
            boundary,
            joint);

        assert(call.bp_left == 44440762);
        assert(call.bp_right == 44440771);
        assert(call.pos == 44440766);
        assert(call.family == "UNKNOWN");
        assert(call.subfamily == "UNKNOWN");
        assert(call.te_name == "UNKNOWN");
        assert(call.final_qc == "PASS_STRUCTURAL_INSERTION|TE_ALIGNMENT_LOW_COVERAGE|PASS_BOUNDARY_SMALL_DEL|PASS_EVENT_CONSENSUS|PASS_EVENT_SEGMENTATION");
        assert(call.best_explanation == "INSERTION_NON_TE");
    }

    {
        Pipeline pipeline(PipelineConfig{}, nullptr);

        ComponentCall component;
        component.chrom = "7";
        component.tid = 0;
        component.anchor_pos = 67143037;
        component.bin_start = 67143000;
        component.bin_end = 67145000;

        EventReadEvidence event_evidence;
        event_evidence.bp_left = 67144040;
        event_evidence.bp_right = 67144187;
        event_evidence.alt_struct_reads = 13;
        event_evidence.ref_span_reads = 4;

        EventConsensus consensus;
        consensus.consensus_len = 679;
        consensus.qc_reason = "PASS_EVENT_CONSENSUS";

        EventSegmentation segmentation;
        segmentation.left_ref_end = 67144160;
        segmentation.right_ref_start = 67144178;
        segmentation.insert_seq = std::string(546, 'A');
        segmentation.left_flank_align_len = 51;
        segmentation.right_flank_align_len = 82;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.555389;
        te.best_query_coverage = 1.0;
        te.cross_family_margin = 0.527917;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        GenotypeCall genotype;
        genotype.genotype = "0/1";
        genotype.af = 0.764706;
        genotype.gq = 99;

        FinalBoundaryDecision boundary;
        boundary.boundary_type = "SMALL_DEL";
        boundary.boundary_len = 18;
        boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

        JointDecisionResult joint;
        joint.emit_unknown_te = true;
        joint.final_qc = "PASS_TE_IMPRECISE";
        joint.best_explanation = "TE";
        joint.explanation_residual = "structural=0;missing=0";
        joint.explanation_path = "NA";

        const FinalCall call = pipeline.emit_final_te_call(
            component,
            event_evidence,
            consensus,
            segmentation,
            te,
            genotype,
            boundary,
            joint);

        assert(call.bp_left == 67144040);
        assert(call.bp_right == 67144187);
        assert(call.pos == 67144113);
    }

    {
        Pipeline pipeline(PipelineConfig{}, nullptr);

        ComponentCall component;
        component.chrom = "12";
        component.tid = 0;
        component.anchor_pos = 1300;

        EventReadEvidence event_evidence;
        event_evidence.bp_left = 900;
        event_evidence.bp_right = 1600;

        EventConsensus consensus;
        consensus.consensus_len = 500;
        consensus.qc_reason = "PASS_EVENT_CONSENSUS";

        EventSegmentation segmentation;
        segmentation.left_ref_start = 1040;
        segmentation.left_ref_end = 1120;
        segmentation.right_ref_start = 1600;
        segmentation.right_ref_end = 1600;
        segmentation.insert_seq = std::string(420, 'A');
        segmentation.left_flank_align_len = 80;
        segmentation.right_flank_align_len = 0;
        segmentation.pass = true;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION_ONE_SIDED_LEFT";

        TEAlignmentEvidence te;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";

        GenotypeCall genotype;
        genotype.genotype = "0/1";

        FinalBoundaryDecision boundary;
        boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        JointDecisionResult joint;
        joint.emit_evidence_te_call = true;
        joint.emit_unknown_te = true;
        joint.final_qc = "PASS_TE_IMPRECISE";

        const FinalCall call = pipeline.emit_final_te_call(
            component,
            event_evidence,
            consensus,
            segmentation,
            te,
            genotype,
            boundary,
            joint);

        assert(call.bp_left == 1120);
        assert(call.bp_right == 1120);
        assert(call.pos == 1120);
        assert(call.final_qc.find("CONSENSUS_FLANK_ANCHORED_BREAKPOINT") !=
               std::string::npos);
    }

    {
        Pipeline pipeline(PipelineConfig{}, nullptr);

        ComponentCall component;
        component.chrom = "12";
        component.tid = 0;
        component.anchor_pos = 2400;

        EventReadEvidence event_evidence;
        event_evidence.bp_left = 2380;
        event_evidence.bp_right = 2420;

        EventConsensus consensus;
        consensus.consensus_len = 820;
        consensus.qc_reason = "PASS_EVENT_CONSENSUS";

        EventSegmentation segmentation;
        segmentation.left_ref_end = 2395;
        segmentation.right_ref_start = 2405;
        segmentation.insert_seq = std::string(660, 'A');
        segmentation.left_flank_align_len = 80;
        segmentation.right_flank_align_len = 80;
        segmentation.pass = true;
        segmentation.qc_reason = "PASS_EVENT_SEGMENTATION";

        TEAlignmentEvidence te;
        te.best_family = "Gypsy";
        te.best_subfamily = "Gypsy-12";
        te.best_identity = 0.97;
        te.best_query_coverage = 0.94;
        te.cross_family_margin = 0.24;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        te.annotation_confidence = "HIGH";
        te.te_sequence_explanation.status = TeAnnotationStatus::kResolved;
        te.te_sequence_explanation.family = "Gypsy";
        te.te_sequence_explanation.subfamily = "Gypsy-12";

        GenotypeCall genotype;
        genotype.genotype = "0/1";

        FinalBoundaryDecision boundary;
        boundary.qc = "REJECT_BOUNDARY_TSD_RANGE";

        JointDecisionResult joint;
        joint.emit_te_call = true;
        joint.emit_unknown_te = true;
        joint.emit_evidence_te_call = true;
        joint.final_qc = "PASS_TE_IMPRECISE";
        joint.te_structure_path_confidence = 0.97;
        joint.te_core_coverage = 0.94;

        const FinalCall call = pipeline.emit_final_te_call(
            component,
            event_evidence,
            consensus,
            segmentation,
            te,
            genotype,
            boundary,
            joint);

        assert(call.family == "UNKNOWN");
        assert(call.subfamily == "UNKNOWN");
        assert(call.sequence_family_candidate == "Gypsy");
        assert(call.sequence_subfamily_candidate == "Gypsy-12");
        assert(call.sequence_family_commit_eligible);
        assert(call.final_qc.find("SEQUENCE_FAMILY_COMMITTED") ==
               std::string::npos);

        te.cross_family_margin = 0.01;
        const FinalCall ambiguous_call = pipeline.emit_final_te_call(
            component,
            event_evidence,
            consensus,
            segmentation,
            te,
            genotype,
            boundary,
            joint);
        assert(ambiguous_call.family == "UNKNOWN");
        assert(ambiguous_call.subfamily == "UNKNOWN");
        assert(!ambiguous_call.family_committed);
        assert(!ambiguous_call.sequence_family_commit_eligible);
        assert(ambiguous_call.final_qc.find("SEQUENCE_FAMILY_COMMITTED") ==
               std::string::npos);

        te.best_family = "Unknown";
        te.best_subfamily = "Unknown-7";
        te.cross_family_margin = 0.24;
        te.te_sequence_explanation.family = "Unknown";
        te.te_sequence_explanation.subfamily = "Unknown-7";
        joint.emit_unknown_te = false;
        const FinalCall library_unknown_call = pipeline.emit_final_te_call(
            component,
            event_evidence,
            consensus,
            segmentation,
            te,
            genotype,
            boundary,
            joint);
        assert(library_unknown_call.family == "Unknown");
        assert(library_unknown_call.subfamily == "Unknown-7");
        assert(library_unknown_call.family_committed);
    }

    {
        EventExistenceEvidence strong_existence;
        strong_existence.alt_struct_reads = 15;
        set_precise_bam_breakdown(strong_existence, 6, 2, 4, 3);
        strong_existence.ref_span_reads = 2;
        strong_existence.gq = 99;
        strong_existence.score = 3.0;

        EventSegmentationEvidence one_sided_segmentation;
        one_sided_segmentation.has_insert_seq = true;
        one_sided_segmentation.has_left_flank = false;
        one_sided_segmentation.has_right_flank = true;
        one_sided_segmentation.pair_valid = true;
        one_sided_segmentation.insert_len = 263;
        one_sided_segmentation.score = 0.25;

        BoundaryEvidence unclosed_boundary;
        unclosed_boundary.geometry_defined = false;
        unclosed_boundary.canonical_pass = false;
        unclosed_boundary.evidence_consistent = false;

        TEAlignmentEvidence te;
        te.best_family = "5S-Deu-L2";
        te.best_subfamily = "5S-Deu-L2-3";
        te.best_identity = 0.980989;
        te.best_query_coverage = 1.0;
        te.cross_family_margin = 0.12336;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                strong_existence,
                one_sided_segmentation,
                te,
                unclosed_boundary);

        assert(result.emit_te_call);
        assert(result.emit_unknown_te);
        assert(result.emit_evidence_te_call);
        assert(result.final_qc == "PASS_TE_IMPRECISE");
    }

    {
        EventExistenceEvidence strong_existence;
        strong_existence.alt_struct_reads = 7;
        strong_existence.ref_span_reads = 5;
        strong_existence.gq = 99;
        strong_existence.score = 3.0;

        EventSegmentationEvidence one_sided_segmentation;
        one_sided_segmentation.has_insert_seq = true;
        one_sided_segmentation.has_left_flank = true;
        one_sided_segmentation.has_right_flank = false;
        one_sided_segmentation.pair_valid = true;
        one_sided_segmentation.insert_len = 435;
        one_sided_segmentation.score = 0.25;

        BoundaryEvidence unclosed_boundary;
        unclosed_boundary.geometry_defined = false;

        TEAlignmentEvidence te;
        te.best_family = "5S-Deu-L2";
        te.best_subfamily = "5S-Deu-L2-4";
        te.best_identity = 0.896254;
        te.best_query_coverage = 0.793103;
        te.cross_family_margin = 0.158838;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                strong_existence,
                one_sided_segmentation,
                te,
                unclosed_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence strong_existence;
        strong_existence.alt_struct_reads = 5;
        strong_existence.ref_span_reads = 4;
        strong_existence.gq = 96;
        strong_existence.score = 3.0;

        EventSegmentationEvidence one_sided_segmentation;
        one_sided_segmentation.has_insert_seq = true;
        one_sided_segmentation.has_left_flank = true;
        one_sided_segmentation.has_right_flank = false;
        one_sided_segmentation.pair_valid = true;
        one_sided_segmentation.insert_len = 454;
        one_sided_segmentation.score = 0.25;

        BoundaryEvidence unclosed_boundary;
        unclosed_boundary.geometry_defined = false;

        TEAlignmentEvidence te;
        te.best_family = "5S-Deu-L2";
        te.best_subfamily = "NA";
        te.best_identity = 0.827922;
        te.best_query_coverage = 0.665198;
        te.cross_family_margin = 0.489058;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                strong_existence,
                one_sided_segmentation,
                te,
                unclosed_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence alt_only_existence;
        alt_only_existence.alt_struct_reads = 7;
        set_precise_bam_breakdown(alt_only_existence, 1, 0, 3, 3);
        alt_only_existence.ref_span_reads = 0;
        alt_only_existence.gq = 18;
        alt_only_existence.score = 0.0;

        EventSegmentationEvidence one_sided_segmentation;
        one_sided_segmentation.has_insert_seq = true;
        one_sided_segmentation.has_left_flank = true;
        one_sided_segmentation.has_right_flank = false;
        one_sided_segmentation.pair_valid = true;
        one_sided_segmentation.insert_len = 381;
        one_sided_segmentation.score = 0.371321;

        BoundaryEvidence unclosed_boundary;
        unclosed_boundary.geometry_defined = false;

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.55794;
        te.best_query_coverage = 0.994751;
        te.cross_family_margin = 0.0744059;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                alt_only_existence,
                one_sided_segmentation,
                te,
                unclosed_boundary);

        assert(!result.emit_te_call);
        assert(!result.emit_unknown_te);
        assert(result.emit_evidence_te_call);
        assert(result.final_qc == "TE_AMBIGUOUS");
        assert(result.posterior_qc == "TE_POSTERIOR_LOW");
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence weak_one_sided_resolved_annotation;
        weak_one_sided_resolved_annotation.alt_struct_reads = 7;
        weak_one_sided_resolved_annotation.ref_span_reads = 5;
        weak_one_sided_resolved_annotation.gq = 99;
        weak_one_sided_resolved_annotation.score = 3.0;

        EventSegmentationEvidence one_sided_segmentation;
        one_sided_segmentation.has_insert_seq = true;
        one_sided_segmentation.has_left_flank = true;
        one_sided_segmentation.has_right_flank = false;
        one_sided_segmentation.pair_valid = true;
        one_sided_segmentation.insert_len = 435;
        one_sided_segmentation.score = 0.25;

        BoundaryEvidence unclosed_boundary;
        unclosed_boundary.geometry_defined = false;

        TEAlignmentEvidence resolved_annotation;
        resolved_annotation.best_family = "5S-Deu-L2";
        resolved_annotation.best_subfamily = "5S-Deu-L2-4";
        resolved_annotation.best_identity = 0.98;
        resolved_annotation.best_query_coverage = 0.99;
        resolved_annotation.cross_family_margin = 0.40;
        resolved_annotation.pass = true;
        resolved_annotation.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        resolved_annotation.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        resolved_annotation.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                weak_one_sided_resolved_annotation,
                one_sided_segmentation,
                resolved_annotation,
                unclosed_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence softclip_only_existence;
        softclip_only_existence.alt_struct_reads = 38;
        softclip_only_existence.alt_split_reads = 0;
        softclip_only_existence.alt_indel_reads = 0;
        softclip_only_existence.alt_left_clip_reads = 38;
        softclip_only_existence.alt_right_clip_reads = 0;
        softclip_only_existence.ref_span_reads = 10;
        softclip_only_existence.gq = 99;
        softclip_only_existence.score = 3.0;

        EventSegmentationEvidence one_sided_segmentation;
        one_sided_segmentation.has_insert_seq = true;
        one_sided_segmentation.has_left_flank = true;
        one_sided_segmentation.has_right_flank = false;
        one_sided_segmentation.pair_valid = true;
        one_sided_segmentation.insert_len = 240;
        one_sided_segmentation.score = 1.0;

        BoundaryEvidence unclosed_boundary;
        unclosed_boundary.geometry_defined = false;
        unclosed_boundary.canonical_pass = false;
        unclosed_boundary.evidence_consistent = false;

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.556777;
        te.best_query_coverage = 0.979167;
        te.cross_family_margin = 0.02753;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                softclip_only_existence,
                one_sided_segmentation,
                te,
                unclosed_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence weak_clip_existence;
        weak_clip_existence.alt_struct_reads = 3;
        weak_clip_existence.alt_split_reads = 0;
        weak_clip_existence.alt_indel_reads = 0;
        weak_clip_existence.alt_left_clip_reads = 2;
        weak_clip_existence.alt_right_clip_reads = 1;
        weak_clip_existence.ref_span_reads = 10;
        weak_clip_existence.gq = 40;
        weak_clip_existence.score = 1.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 309;
        complete_segmentation.score = 1.0;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 16;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.561983;
        te.best_query_coverage = 0.983819;
        te.cross_family_margin = 0.0334339;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_EDGE";
        te.sequence_model_score = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                weak_clip_existence,
                complete_segmentation,
                te,
                small_del_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence precise_existence;
        precise_existence.alt_struct_reads = 3;
        precise_existence.alt_split_reads = 1;
        precise_existence.alt_indel_reads = 0;
        precise_existence.alt_left_clip_reads = 1;
        precise_existence.alt_right_clip_reads = 1;
        precise_existence.ref_span_reads = 10;
        precise_existence.gq = 40;
        precise_existence.score = 1.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 309;
        complete_segmentation.score = 1.0;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 16;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.561983;
        te.best_query_coverage = 0.983819;
        te.cross_family_margin = 0.0334339;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_EDGE";
        te.sequence_model_score = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                precise_existence,
                complete_segmentation,
                te,
                small_del_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence ref_present_existence;
        ref_present_existence.alt_struct_reads = 7;
        set_precise_bam_breakdown(ref_present_existence, 1, 0, 3, 3);
        ref_present_existence.ref_span_reads = 1;
        ref_present_existence.gq = 18;
        ref_present_existence.score = 0.0;

        EventSegmentationEvidence one_sided_segmentation;
        one_sided_segmentation.has_insert_seq = true;
        one_sided_segmentation.has_left_flank = true;
        one_sided_segmentation.has_right_flank = false;
        one_sided_segmentation.pair_valid = true;
        one_sided_segmentation.insert_len = 381;
        one_sided_segmentation.score = 0.371321;

        BoundaryEvidence unclosed_boundary;
        unclosed_boundary.geometry_defined = false;

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.55794;
        te.best_query_coverage = 0.994751;
        te.cross_family_margin = 0.0744059;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                ref_present_existence,
                one_sided_segmentation,
                te,
                unclosed_boundary);

        assert(!result.emit_te_call);
        assert(!result.emit_unknown_te);
        assert(result.emit_evidence_te_call);
        assert(result.final_qc == "TE_AMBIGUOUS");
        assert(result.posterior_qc == "TE_POSTERIOR_LOW");
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence weak_one_sided_existence;
        weak_one_sided_existence.alt_struct_reads = 8;
        weak_one_sided_existence.alt_split_reads = 1;
        weak_one_sided_existence.alt_indel_reads = 0;
        weak_one_sided_existence.alt_left_clip_reads = 1;
        weak_one_sided_existence.alt_right_clip_reads = 7;
        weak_one_sided_existence.ref_span_reads = 1;
        weak_one_sided_existence.gq = 7;
        weak_one_sided_existence.score = -0.65;

        EventSegmentationEvidence one_sided_segmentation;
        one_sided_segmentation.has_insert_seq = true;
        one_sided_segmentation.has_left_flank = true;
        one_sided_segmentation.has_right_flank = false;
        one_sided_segmentation.pair_valid = true;
        one_sided_segmentation.insert_len = 180;
        one_sided_segmentation.score = 1.34;

        BoundaryEvidence unclosed_boundary;
        unclosed_boundary.geometry_defined = false;
        unclosed_boundary.canonical_pass = false;
        unclosed_boundary.evidence_consistent = false;

        TEAlignmentEvidence te;
        te.best_family = "UNKNOWN";
        te.best_subfamily = "UNKNOWN";
        te.best_identity = 0.586207;
        te.best_query_coverage = 0.955556;
        te.cross_family_margin = 0.00587191;
        te.pass = true;
        te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                weak_one_sided_existence,
                one_sided_segmentation,
                te,
                unclosed_boundary);

        assert(!result.emit_te_call);
        assert(!result.emit_unknown_te);
        assert(result.emit_evidence_te_call);
        assert(result.final_qc == "TE_AMBIGUOUS");
        assert(result.posterior_qc == "TE_POSTERIOR_LOW");
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence small_del_clip_dominated;
        small_del_clip_dominated.alt_struct_reads = 45;
        small_del_clip_dominated.alt_split_reads = 0;
        small_del_clip_dominated.alt_indel_reads = 5;
        small_del_clip_dominated.alt_left_clip_reads = 6;
        small_del_clip_dominated.alt_right_clip_reads = 14;
        small_del_clip_dominated.ref_span_reads = 37;
        small_del_clip_dominated.gq = 99;
        small_del_clip_dominated.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 631;
        complete_segmentation.score = 1.0;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 25;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence resolved_te;
        resolved_te.best_family = "Rex-Babar";
        resolved_te.best_subfamily = "Rex-Babar-11";
        resolved_te.best_identity = 0.961409;
        resolved_te.best_query_coverage = 0.935024;
        resolved_te.cross_family_margin = 0.875169;
        resolved_te.pass = true;
        resolved_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        resolved_te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        resolved_te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                small_del_clip_dominated,
                complete_segmentation,
                resolved_te,
                small_del_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence weak_family_only_small_del;
        weak_family_only_small_del.alt_struct_reads = 43;
        weak_family_only_small_del.alt_split_reads = 0;
        weak_family_only_small_del.alt_indel_reads = 19;
        weak_family_only_small_del.alt_left_clip_reads = 0;
        weak_family_only_small_del.alt_right_clip_reads = 73;
        weak_family_only_small_del.ref_span_reads = 2;
        weak_family_only_small_del.gq = 72;
        weak_family_only_small_del.score = 2.6;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 83;
        complete_segmentation.score = 0.50;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 17;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence family_only_te;
        family_only_te.best_family = "Maverick";
        family_only_te.best_subfamily = "NA";
        family_only_te.best_identity = 0.733945;
        family_only_te.best_query_coverage = 0.987952;
        family_only_te.cross_family_margin = 0.143713;
        family_only_te.pass = true;
        family_only_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY";

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                weak_family_only_small_del,
                complete_segmentation,
                family_only_te,
                small_del_boundary);

        assert(result.emit_te_call);
        assert(!result.emit_unknown_te);
        assert(result.final_qc == "PASS_TE_CLOSED");
        assert(result.lfdr_qc == "TE_LFDR_HIGH");
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence precise_small_del;
        precise_small_del.alt_struct_reads = 80;
        precise_small_del.alt_split_reads = 22;
        precise_small_del.alt_indel_reads = 18;
        precise_small_del.alt_left_clip_reads = 24;
        precise_small_del.alt_right_clip_reads = 26;
        precise_small_del.ref_span_reads = 20;
        precise_small_del.gq = 99;
        precise_small_del.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 900;
        complete_segmentation.score = 1.2;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 12;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence resolved_te;
        resolved_te.best_family = "Rex-Babar";
        resolved_te.best_subfamily = "Rex-Babar-11";
        resolved_te.best_identity = 0.96;
        resolved_te.best_query_coverage = 0.95;
        resolved_te.cross_family_margin = 0.80;
        resolved_te.pass = true;
        resolved_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        resolved_te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        resolved_te.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                precise_small_del,
                complete_segmentation,
                resolved_te,
                small_del_boundary);

        assert(result.emit_te_call);
        assert(!result.emit_unknown_te);
        assert(result.final_qc == "PASS_TE_CLOSED");
        assert(result.robust_mechanistic_qc == "PASS_TE_LFDR");
        assert(result.robust_mechanistic_worst_case_lfdr <= 0.10);
        assert(result.mechanistic_ref_conflict_signal > 0.0);
    }

    {
        EventExistenceEvidence weak_family_only_blunt;
        weak_family_only_blunt.alt_struct_reads = 25;
        weak_family_only_blunt.alt_split_reads = 0;
        weak_family_only_blunt.alt_indel_reads = 4;
        weak_family_only_blunt.alt_left_clip_reads = 6;
        weak_family_only_blunt.alt_right_clip_reads = 15;
        weak_family_only_blunt.ref_span_reads = 38;
        weak_family_only_blunt.gq = 99;
        weak_family_only_blunt.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 1681;
        complete_segmentation.score = 1.6;

        BoundaryEvidence blunt_boundary;
        blunt_boundary.geometry_defined = true;
        blunt_boundary.canonical_pass = true;
        blunt_boundary.evidence_consistent = true;
        blunt_boundary.boundary_type = "BLUNT";
        blunt_boundary.boundary_len = 0;
        blunt_boundary.score = 1.0;

        TEAlignmentEvidence weak_family_only_te;
        weak_family_only_te.best_family = "hAT-Tip100";
        weak_family_only_te.best_subfamily = "NA";
        weak_family_only_te.best_identity = 0.762914;
        weak_family_only_te.best_query_coverage = 0.856633;
        weak_family_only_te.cross_family_margin = 0.139714;
        weak_family_only_te.pass = true;
        weak_family_only_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY";
        weak_family_only_te.sequence_model_label = "TE_MODEL_EDGE";
        weak_family_only_te.sequence_model_score = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                weak_family_only_blunt,
                complete_segmentation,
                weak_family_only_te,
                blunt_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence one_sided_resolved_with_invalid_boundary;
        one_sided_resolved_with_invalid_boundary.alt_struct_reads = 45;
        one_sided_resolved_with_invalid_boundary.alt_split_reads = 0;
        one_sided_resolved_with_invalid_boundary.alt_indel_reads = 3;
        one_sided_resolved_with_invalid_boundary.alt_left_clip_reads = 5;
        one_sided_resolved_with_invalid_boundary.alt_right_clip_reads = 38;
        one_sided_resolved_with_invalid_boundary.ref_span_reads = 5;
        one_sided_resolved_with_invalid_boundary.gq = 23;
        one_sided_resolved_with_invalid_boundary.score = 0.15;

        EventSegmentationEvidence one_sided_segmentation;
        one_sided_segmentation.has_insert_seq = true;
        one_sided_segmentation.has_left_flank = true;
        one_sided_segmentation.has_right_flank = false;
        one_sided_segmentation.pair_valid = true;
        one_sided_segmentation.insert_len = 3075;
        one_sided_segmentation.score = 1.1;

        BoundaryEvidence invalid_boundary;
        invalid_boundary.geometry_defined = false;
        invalid_boundary.canonical_pass = false;
        invalid_boundary.evidence_consistent = false;
        invalid_boundary.boundary_type = "REJECT";
        invalid_boundary.boundary_len = 0;
        invalid_boundary.score = -2.0;
        invalid_boundary.qc = "REJECT_BOUNDARY_INVALID_REF_SEGMENTS";

        TEAlignmentEvidence resolved_te_annotation;
        resolved_te_annotation.best_family = "PIF-Harbinger";
        resolved_te_annotation.best_subfamily = "PIF-Harbinger-1";
        resolved_te_annotation.best_identity = 0.989653;
        resolved_te_annotation.best_query_coverage = 0.973333;
        resolved_te_annotation.cross_family_margin = 0.430067;
        resolved_te_annotation.pass = true;
        resolved_te_annotation.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        resolved_te_annotation.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        resolved_te_annotation.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                one_sided_resolved_with_invalid_boundary,
                one_sided_segmentation,
                resolved_te_annotation,
                invalid_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence structurally_supported_but_model_uninformative;
        structurally_supported_but_model_uninformative.alt_struct_reads = 40;
        structurally_supported_but_model_uninformative.alt_split_reads = 3;
        structurally_supported_but_model_uninformative.alt_indel_reads = 3;
        structurally_supported_but_model_uninformative.alt_left_clip_reads = 20;
        structurally_supported_but_model_uninformative.alt_right_clip_reads = 20;
        structurally_supported_but_model_uninformative.ref_span_reads = 20;
        structurally_supported_but_model_uninformative.gq = 99;
        structurally_supported_but_model_uninformative.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 900;
        complete_segmentation.score = 1.2;

        BoundaryEvidence tsd_boundary;
        tsd_boundary.geometry_defined = true;
        tsd_boundary.canonical_pass = true;
        tsd_boundary.evidence_consistent = true;
        tsd_boundary.boundary_type = "TSD";
        tsd_boundary.boundary_len = 8;
        tsd_boundary.score = 1.0;

        TEAlignmentEvidence resolved_annotation;
        resolved_annotation.best_family = "Rex-Babar";
        resolved_annotation.best_subfamily = "Rex-Babar-11";
        resolved_annotation.best_identity = 0.99;
        resolved_annotation.best_query_coverage = 0.99;
        resolved_annotation.cross_family_margin = 0.80;
        resolved_annotation.pass = true;
        resolved_annotation.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        resolved_annotation.sequence_model_label = "TE_MODEL_UNAVAILABLE";
        resolved_annotation.sequence_model_score = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                structurally_supported_but_model_uninformative,
                complete_segmentation,
                resolved_annotation,
                tsd_boundary);

        assert(!result.emit_te_call);

        resolved_annotation.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        resolved_annotation.sequence_model_score = 0.50;

        const JointDecisionResult model_supported_result =
            evaluate_joint_hypotheses(
                structurally_supported_but_model_uninformative,
                complete_segmentation,
                resolved_annotation,
                tsd_boundary);

        assert(!model_supported_result.emit_te_call);
        assert(model_supported_result.best_explanation == "TE");
        assert(model_supported_result.final_qc == "TE_AMBIGUOUS");
    }

    {
        EventExistenceEvidence existence_with_clip_concordance;
        existence_with_clip_concordance.alt_struct_reads = 3;
        existence_with_clip_concordance.alt_split_reads = 1;
        existence_with_clip_concordance.alt_indel_reads = 1;
        existence_with_clip_concordance.alt_left_clip_reads = 2;
        existence_with_clip_concordance.alt_right_clip_reads = 2;
        existence_with_clip_concordance.ref_span_reads = 1;
        existence_with_clip_concordance.gq = 30;
        existence_with_clip_concordance.score = 0.25;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 240;
        complete_segmentation.score = 0.85;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 8;
        small_del_boundary.score = 0.25;

        TEAlignmentEvidence weak_te_like;
        weak_te_like.best_family = "UNKNOWN";
        weak_te_like.best_subfamily = "UNKNOWN";
        weak_te_like.best_identity = 0.53;
        weak_te_like.best_query_coverage = 0.92;
        weak_te_like.pass = false;
        weak_te_like.qc_reason = "TE_ALIGNMENT_LOW_IDENTITY";
        weak_te_like.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        weak_te_like.sequence_model_score = 0.20;

        const JointDecisionResult without_clip_concordance =
            evaluate_joint_hypotheses(
                existence_with_clip_concordance,
                complete_segmentation,
                weak_te_like,
                small_del_boundary);
        assert(!without_clip_concordance.emit_te_call);

        ClipInsertConcordanceEvidence clip_concordance;
        clip_concordance.pass = true;
        clip_concordance.full_insert_reads = 1;
        clip_concordance.left_clip_reads = 2;
        clip_concordance.right_clip_reads = 2;
        clip_concordance.max_left_identity = 0.94;
        clip_concordance.max_right_identity = 0.95;
        clip_concordance.qc = "PASS_CLIP_INSERT_CONCORDANCE";

        const JointDecisionResult with_clip_concordance =
            evaluate_joint_hypotheses(
                existence_with_clip_concordance,
                complete_segmentation,
                weak_te_like,
                small_del_boundary,
                &clip_concordance);
        assert(!with_clip_concordance.emit_te_call);
        assert(with_clip_concordance.final_qc.find("RESCUE") == std::string::npos);
        assert(
            with_clip_concordance.final_qc.find("PASS_CLIP_INSERT_CONCORDANCE") ==
            std::string::npos);
    }

    {
        EventExistenceEvidence existence_with_clip_concordance;
        existence_with_clip_concordance.alt_struct_reads = 4;
        existence_with_clip_concordance.alt_split_reads = 1;
        existence_with_clip_concordance.alt_indel_reads = 1;
        existence_with_clip_concordance.alt_left_clip_reads = 2;
        existence_with_clip_concordance.alt_right_clip_reads = 2;
        existence_with_clip_concordance.ref_span_reads = 1;
        existence_with_clip_concordance.gq = 40;
        existence_with_clip_concordance.score = 0.50;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 240;
        complete_segmentation.score = 0.90;

        BoundaryEvidence boundary;
        boundary.geometry_defined = true;
        boundary.canonical_pass = true;
        boundary.evidence_consistent = true;
        boundary.score = 0.50;

        TEAlignmentEvidence non_te_like;
        non_te_like.best_identity = 0.20;
        non_te_like.best_query_coverage = 0.15;
        non_te_like.pass = false;
        non_te_like.qc_reason = "NO_TE_ALIGNMENT_MATCH";
        non_te_like.sequence_model_label = "TE_MODEL_OUTLIER";
        non_te_like.sequence_model_score = -0.60;

        ClipInsertConcordanceEvidence clip_concordance;
        clip_concordance.pass = true;
        clip_concordance.full_insert_reads = 1;
        clip_concordance.left_clip_reads = 2;
        clip_concordance.right_clip_reads = 2;
        clip_concordance.max_left_identity = 0.94;
        clip_concordance.max_right_identity = 0.95;
        clip_concordance.qc = "PASS_CLIP_INSERT_CONCORDANCE";

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                existence_with_clip_concordance,
                complete_segmentation,
                non_te_like,
                boundary,
                &clip_concordance);
        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence strict_bam_supported;
        strict_bam_supported.alt_struct_reads = 12;
        strict_bam_supported.alt_split_reads = 2;
        strict_bam_supported.alt_indel_reads = 3;
        strict_bam_supported.alt_left_clip_reads = 4;
        strict_bam_supported.alt_right_clip_reads = 4;
        strict_bam_supported.ref_span_reads = 4;
        strict_bam_supported.gq = 99;
        strict_bam_supported.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 288;
        complete_segmentation.score = 0.85;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 12;
        small_del_boundary.score = 0.50;

        TEAlignmentEvidence model_supported_no_named_hit;
        model_supported_no_named_hit.pass = false;
        model_supported_no_named_hit.qc_reason = "NO_TE_ALIGNMENT";
        model_supported_no_named_hit.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        model_supported_no_named_hit.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                strict_bam_supported,
                complete_segmentation,
                model_supported_no_named_hit,
                small_del_boundary);

        assert(!result.emit_te_call);
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence strict_bam_supported;
        strict_bam_supported.alt_struct_reads = 12;
        strict_bam_supported.alt_split_reads = 2;
        strict_bam_supported.alt_indel_reads = 3;
        strict_bam_supported.alt_left_clip_reads = 4;
        strict_bam_supported.alt_right_clip_reads = 4;
        strict_bam_supported.ref_span_reads = 4;
        strict_bam_supported.gq = 99;
        strict_bam_supported.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 288;
        complete_segmentation.score = 0.85;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 12;
        small_del_boundary.score = 0.50;

        TEAlignmentEvidence model_outlier_no_named_hit;
        model_outlier_no_named_hit.pass = false;
        model_outlier_no_named_hit.qc_reason = "NO_TE_ALIGNMENT";
        model_outlier_no_named_hit.sequence_model_label = "TE_MODEL_OUTLIER";
        model_outlier_no_named_hit.sequence_model_score = -0.50;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                strict_bam_supported,
                complete_segmentation,
                model_outlier_no_named_hit,
                small_del_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence ref_heavy;
        ref_heavy.alt_struct_reads = 12;
        ref_heavy.alt_split_reads = 2;
        ref_heavy.alt_indel_reads = 3;
        ref_heavy.alt_left_clip_reads = 4;
        ref_heavy.alt_right_clip_reads = 4;
        ref_heavy.ref_span_reads = 8;
        ref_heavy.gq = 99;
        ref_heavy.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 288;
        complete_segmentation.score = 0.85;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 12;
        small_del_boundary.score = 0.50;

        TEAlignmentEvidence model_supported_no_named_hit;
        model_supported_no_named_hit.pass = false;
        model_supported_no_named_hit.qc_reason = "NO_TE_ALIGNMENT";
        model_supported_no_named_hit.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        model_supported_no_named_hit.sequence_model_score = 0.40;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                ref_heavy,
                complete_segmentation,
                model_supported_no_named_hit,
                small_del_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence edge_model_supported_unknown;
        edge_model_supported_unknown.alt_struct_reads = 8;
        edge_model_supported_unknown.alt_split_reads = 0;
        edge_model_supported_unknown.alt_indel_reads = 1;
        edge_model_supported_unknown.alt_left_clip_reads = 6;
        edge_model_supported_unknown.alt_right_clip_reads = 1;
        edge_model_supported_unknown.ref_span_reads = 3;
        edge_model_supported_unknown.gq = 99;
        edge_model_supported_unknown.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 1025;
        complete_segmentation.score = 1.24;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 18;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence edge_unknown_te;
        edge_unknown_te.best_family = "UNKNOWN";
        edge_unknown_te.best_subfamily = "UNKNOWN";
        edge_unknown_te.best_identity = 0.546802;
        edge_unknown_te.best_query_coverage = 1.0;
        edge_unknown_te.cross_family_margin = 0.532295;
        edge_unknown_te.pass = true;
        edge_unknown_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        edge_unknown_te.sequence_model_label = "TE_MODEL_EDGE";
        edge_unknown_te.sequence_model_score = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                edge_model_supported_unknown,
                complete_segmentation,
                edge_unknown_te,
                small_del_boundary);

        assert(!result.emit_te_call);
        assert(result.best_explanation == "TE");
        assert(result.final_qc == "TE_AMBIGUOUS");
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence edge_model_supported_resolved;
        edge_model_supported_resolved.alt_struct_reads = 30;
        edge_model_supported_resolved.alt_split_reads = 1;
        edge_model_supported_resolved.alt_indel_reads = 4;
        edge_model_supported_resolved.alt_left_clip_reads = 7;
        edge_model_supported_resolved.alt_right_clip_reads = 19;
        edge_model_supported_resolved.ref_span_reads = 0;
        edge_model_supported_resolved.gq = 99;
        edge_model_supported_resolved.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 1185;
        complete_segmentation.score = 0.875;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 30;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence edge_resolved_te;
        edge_resolved_te.best_family = "Gypsy";
        edge_resolved_te.best_subfamily = "Gypsy-117";
        edge_resolved_te.best_identity = 0.982639;
        edge_resolved_te.best_query_coverage = 0.971308;
        edge_resolved_te.cross_family_margin = 0.942631;
        edge_resolved_te.pass = true;
        edge_resolved_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        edge_resolved_te.sequence_model_label = "TE_MODEL_EDGE";
        edge_resolved_te.sequence_model_score = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                edge_model_supported_resolved,
                complete_segmentation,
                edge_resolved_te,
                small_del_boundary);

        assert(!result.emit_te_call);
        assert(result.best_explanation == "TE");
        assert(result.final_qc == "TE_AMBIGUOUS");
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence edge_model_low_identity;
        edge_model_low_identity.alt_struct_reads = 19;
        edge_model_low_identity.alt_split_reads = 0;
        edge_model_low_identity.alt_indel_reads = 3;
        edge_model_low_identity.alt_left_clip_reads = 9;
        edge_model_low_identity.alt_right_clip_reads = 7;
        edge_model_low_identity.ref_span_reads = 5;
        edge_model_low_identity.gq = 99;
        edge_model_low_identity.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 1078;
        complete_segmentation.score = 1.0;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 28;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence low_identity_edge_te;
        low_identity_edge_te.best_family = "MULE-MuDR";
        low_identity_edge_te.best_subfamily = "MULE-MuDR-3";
        low_identity_edge_te.best_identity = 0.548485;
        low_identity_edge_te.best_query_coverage = 0.995362;
        low_identity_edge_te.cross_family_margin = 0.0105093;
        low_identity_edge_te.pass = false;
        low_identity_edge_te.qc_reason = "TE_ALIGNMENT_LOW_IDENTITY";
        low_identity_edge_te.sequence_model_label = "TE_MODEL_EDGE";
        low_identity_edge_te.sequence_model_score = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                edge_model_low_identity,
                complete_segmentation,
                low_identity_edge_te,
                small_del_boundary);

        assert(!result.emit_te_call);
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        EventExistenceEvidence ref_heavy_edge_unknown;
        ref_heavy_edge_unknown.alt_struct_reads = 8;
        ref_heavy_edge_unknown.alt_split_reads = 0;
        ref_heavy_edge_unknown.alt_indel_reads = 1;
        ref_heavy_edge_unknown.alt_left_clip_reads = 6;
        ref_heavy_edge_unknown.alt_right_clip_reads = 1;
        ref_heavy_edge_unknown.ref_span_reads = 6;
        ref_heavy_edge_unknown.gq = 99;
        ref_heavy_edge_unknown.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 1025;
        complete_segmentation.score = 1.24;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 18;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence edge_unknown_te;
        edge_unknown_te.best_identity = 0.546802;
        edge_unknown_te.best_query_coverage = 1.0;
        edge_unknown_te.cross_family_margin = 0.532295;
        edge_unknown_te.pass = true;
        edge_unknown_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        edge_unknown_te.sequence_model_label = "TE_MODEL_EDGE";
        edge_unknown_te.sequence_model_score = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                ref_heavy_edge_unknown,
                complete_segmentation,
                edge_unknown_te,
                small_del_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence outlier_model_supported_unknown;
        outlier_model_supported_unknown.alt_struct_reads = 8;
        outlier_model_supported_unknown.alt_split_reads = 0;
        outlier_model_supported_unknown.alt_indel_reads = 1;
        outlier_model_supported_unknown.alt_left_clip_reads = 6;
        outlier_model_supported_unknown.alt_right_clip_reads = 1;
        outlier_model_supported_unknown.ref_span_reads = 3;
        outlier_model_supported_unknown.gq = 99;
        outlier_model_supported_unknown.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 1025;
        complete_segmentation.score = 1.24;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 18;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence outlier_unknown_te;
        outlier_unknown_te.best_identity = 0.546802;
        outlier_unknown_te.best_query_coverage = 1.0;
        outlier_unknown_te.cross_family_margin = 0.532295;
        outlier_unknown_te.pass = true;
        outlier_unknown_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        outlier_unknown_te.sequence_model_label = "TE_MODEL_OUTLIER";
        outlier_unknown_te.sequence_model_score = -0.50;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                outlier_model_supported_unknown,
                complete_segmentation,
                outlier_unknown_te,
                small_del_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence edge_model_supported_unknown;
        edge_model_supported_unknown.alt_struct_reads = 8;
        edge_model_supported_unknown.alt_split_reads = 0;
        edge_model_supported_unknown.alt_indel_reads = 1;
        edge_model_supported_unknown.alt_left_clip_reads = 6;
        edge_model_supported_unknown.alt_right_clip_reads = 1;
        edge_model_supported_unknown.ref_span_reads = 3;
        edge_model_supported_unknown.gq = 99;
        edge_model_supported_unknown.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 1025;
        complete_segmentation.score = 1.24;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 18;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence high_risk_low_annotation;
        high_risk_low_annotation.best_family = "UNKNOWN";
        high_risk_low_annotation.best_subfamily = "UNKNOWN";
        high_risk_low_annotation.best_identity = 0.546802;
        high_risk_low_annotation.best_query_coverage = 1.0;
        high_risk_low_annotation.cross_family_margin = 0.532295;
        high_risk_low_annotation.pass = true;
        high_risk_low_annotation.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        high_risk_low_annotation.sequence_model_label = "TE_MODEL_EDGE";
        high_risk_low_annotation.sequence_model_score = 0.0;
        high_risk_low_annotation.annotation_confidence = "LOW";
        high_risk_low_annotation.annotation_residual_fraction = 0.62;
        high_risk_low_annotation.annotation_masked_fraction = 0.71;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                edge_model_supported_unknown,
                complete_segmentation,
                high_risk_low_annotation,
                small_del_boundary);

        assert(!result.emit_te_call);
    }

    {
        EventExistenceEvidence edge_model_supported_unknown;
        edge_model_supported_unknown.alt_struct_reads = 8;
        edge_model_supported_unknown.alt_split_reads = 0;
        edge_model_supported_unknown.alt_indel_reads = 1;
        edge_model_supported_unknown.alt_left_clip_reads = 6;
        edge_model_supported_unknown.alt_right_clip_reads = 1;
        edge_model_supported_unknown.ref_span_reads = 3;
        edge_model_supported_unknown.gq = 99;
        edge_model_supported_unknown.score = 3.0;

        EventSegmentationEvidence complete_segmentation;
        complete_segmentation.has_insert_seq = true;
        complete_segmentation.has_left_flank = true;
        complete_segmentation.has_right_flank = true;
        complete_segmentation.pair_valid = true;
        complete_segmentation.insert_len = 1025;
        complete_segmentation.score = 1.24;

        BoundaryEvidence small_del_boundary;
        small_del_boundary.geometry_defined = true;
        small_del_boundary.canonical_pass = true;
        small_del_boundary.evidence_consistent = true;
        small_del_boundary.boundary_type = "SMALL_DEL";
        small_del_boundary.boundary_len = 18;
        small_del_boundary.score = 1.0;

        TEAlignmentEvidence clean_low_annotation;
        clean_low_annotation.best_family = "UNKNOWN";
        clean_low_annotation.best_subfamily = "UNKNOWN";
        clean_low_annotation.best_identity = 0.546802;
        clean_low_annotation.best_query_coverage = 1.0;
        clean_low_annotation.cross_family_margin = 0.532295;
        clean_low_annotation.pass = true;
        clean_low_annotation.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        clean_low_annotation.sequence_model_label = "TE_MODEL_EDGE";
        clean_low_annotation.sequence_model_score = 0.0;
        clean_low_annotation.annotation_confidence = "LOW";
        clean_low_annotation.annotation_residual_fraction = 0.02;
        clean_low_annotation.annotation_masked_fraction = 0.0;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(
                edge_model_supported_unknown,
                complete_segmentation,
                clean_low_annotation,
                small_del_boundary);

        assert(!result.emit_te_call);
        assert(result.best_explanation == "TE");
        assert(result.final_qc == "TE_AMBIGUOUS");
        assert(result.final_qc.find("RESCUE") == std::string::npos);
    }

    {
        TEAlignmentEvidence resolved_te;
        resolved_te.best_family = "Gypsy";
        resolved_te.best_subfamily = "Gypsy-1";
        resolved_te.pass = true;
        resolved_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        resolved_te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        resolved_te.sequence_model_score = 0.4;

        const ExplanationDecision decision =
            evaluate_event_explanations_for_test(
                existence,
                segmentation,
                resolved_te,
                boundary);

        assert(decision.best.kind == ExplanationKind::kTe);
        assert(decision.final_qc == "PASS_TE_CLOSED");
        assert(decision.emit_te_call);
    }

    {
        EventExistenceEvidence alt_only;
        alt_only.score = 2.0;
        alt_only.alt_struct_reads = 6;
        alt_only.alt_split_reads = 1;
        alt_only.alt_indel_reads = 1;
        alt_only.alt_left_clip_reads = 4;
        alt_only.alt_right_clip_reads = 0;
        alt_only.ref_span_reads = 0;

        EventSegmentationEvidence one_sided;
        one_sided.has_insert_seq = true;
        one_sided.has_left_flank = true;
        one_sided.has_right_flank = false;
        one_sided.pair_valid = false;
        one_sided.insert_len = 420;
        one_sided.score = 0.5;
        one_sided.qc = "PASS_EVENT_SEGMENTATION_ONE_SIDED";

        BoundaryEvidence open_boundary;
        open_boundary.geometry_defined = false;
        open_boundary.canonical_pass = false;
        open_boundary.evidence_consistent = false;
        open_boundary.qc = "REJECT_BOUNDARY_UNSET";

        TEAlignmentEvidence unknown_te;
        unknown_te.best_family = "UNKNOWN";
        unknown_te.best_subfamily = "UNKNOWN";
        unknown_te.pass = true;
        unknown_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        unknown_te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        unknown_te.sequence_model_score = 0.4;

        const ExplanationDecision decision =
            evaluate_event_explanations_for_test(
                alt_only,
                one_sided,
                unknown_te,
                open_boundary);

        assert(decision.best.kind == ExplanationKind::kTe);
        assert(decision.final_qc == "PASS_TE_IMPRECISE");
        assert(decision.emit_te_call);
        assert(decision.emit_unknown_te);
    }

    {
        TEAlignmentEvidence outlier_te;
        outlier_te.best_family = "UNKNOWN";
        outlier_te.best_subfamily = "UNKNOWN";
        outlier_te.pass = true;
        outlier_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_UNKNOWN";
        outlier_te.sequence_model_label = "TE_MODEL_OUTLIER";
        outlier_te.sequence_model_score = -0.8;

        const ExplanationDecision decision =
            evaluate_event_explanations_for_test(
                existence,
                segmentation,
                outlier_te,
                boundary);

        assert(decision.best.kind != ExplanationKind::kTe ||
               decision.final_qc == "TE_AMBIGUOUS");
        assert(!decision.emit_te_call);
    }

    {
        TEAlignmentEvidence resolved_te;
        resolved_te.best_family = "Gypsy";
        resolved_te.best_subfamily = "Gypsy-1";
        resolved_te.pass = true;
        resolved_te.qc_reason = "PASS_INSERT_TE_ALIGNMENT";
        resolved_te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
        resolved_te.sequence_model_score = 0.4;

        const JointDecisionResult result =
            evaluate_joint_hypotheses(existence, segmentation, resolved_te, boundary);
        assert(result.final_qc.find("RESCUE") == std::string::npos);
        assert(result.best_explanation == "TE");
        assert(result.explanation_residual != "NA");
    }

    return 0;
}
