#include "decision_policy.h"
#include "pipeline.h"

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

int main() {
    using namespace placer;

    EventExistenceEvidence existence;
    existence.best_gt = "0/0";
    existence.af = 1.0;
    existence.gq = 0;
    existence.alt_struct_reads = 7;
    existence.alt_split_reads = 0;
    existence.alt_indel_reads = 1;
    existence.alt_left_clip_reads = 3;
    existence.alt_right_clip_reads = 3;
    existence.ref_span_reads = 0;
    existence.depth = 7;
    existence.score = -1.0;

    EventSegmentationEvidence seg;
    seg.has_consensus = true;
    seg.has_left_flank = true;
    seg.has_right_flank = true;
    seg.has_insert_seq = true;
    seg.pair_valid = true;
    seg.left_align_len = 52;
    seg.right_align_len = 80;
    seg.insert_len = 154;
    seg.score = 0.678154;
    seg.qc = "PASS_EVENT_SEGMENTATION";

    TEAlignmentEvidence te;
    te.best_family = "RTE-BovB";
    te.best_identity = 0.75974;
    te.best_query_coverage = 1.0;
    te.cross_family_margin = 0.188312;
    te.pass = true;
    te.qc_reason = "PASS_INSERT_TE_ALIGNMENT_FAMILY_ONLY";
    te.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
    te.sequence_model_score = 0.40;

    BoundaryEvidence boundary;
    boundary.geometry_defined = true;
    boundary.canonical_pass = true;
    boundary.evidence_consistent = true;
    boundary.boundary_type = "SMALL_DEL";
    boundary.boundary_len = 17;
    boundary.score = 1.0;
    boundary.qc = "PASS_BOUNDARY_SMALL_DEL";

    const JointDecisionResult joint =
        evaluate_joint_hypotheses(existence, seg, te, boundary);

    assert(joint.best.kind == FinalHypothesisKind::kTeResolved);
    assert(joint.emit_te_call);
    assert(joint.final_qc == "PASS_TE_CLOSED");
    assert(joint.posterior_qc == "PASS_TE_POSTERIOR");
    assert(joint.lfdr_qc == "PASS_TE_LFDR");
    assert(joint.robust_mechanistic_qc == "PASS_TE_LFDR");
    assert(joint.robust_mechanistic_worst_case_lfdr <= 0.10);
    return 0;
}
