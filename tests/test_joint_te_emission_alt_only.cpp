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
    assert(joint.final_qc == "PASS_FINAL_TE_CALL");
    return 0;
}
