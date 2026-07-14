#ifdef NDEBUG
#undef NDEBUG
#endif
#include "decision_policy.h"
#include "pipeline.h"

#include <cassert>
#include <cstdint>
#include <string>

// Focused unit tests for the two structural-sanity predicates used by
// evaluate_joint_hypotheses to decide whether the TE-unknown (h2) and
// TE-resolved (h3) hypotheses are eligible for ranking. After the move to an
// lFDR-primary emission gate these predicates are necessary conditions only:
// they reject a hypothesis when the event definitionally cannot be a TE call
// (no insert sequence, no TE alignment, annotation quality too low, resolved
// call without a resolved qc_reason). They no longer encode read-count,
// insert-length, segmentation-score, one-sided or independent-support ladders --
// whether an eligible hypothesis is emitted is decided by the calibrated
// worst-case local FDR, tested elsewhere.

namespace {

using placer::BoundaryEvidence;
using placer::compute_te_resolved_hard_veto;
using placer::compute_te_unknown_hard_veto;
using placer::EventExistenceEvidence;
using placer::EventSegmentationEvidence;
using placer::TEAlignmentEvidence;

EventSegmentationEvidence closed_segmentation() {
    EventSegmentationEvidence s;
    s.has_consensus = true;
    s.has_left_flank = true;
    s.has_right_flank = true;
    s.has_insert_seq = true;
    s.pair_valid = true;
    s.insert_len = 300;
    s.score = 0.9;
    s.qc = "PASS_EVENT_SEGMENTATION";
    return s;
}

EventExistenceEvidence some_existence() {
    EventExistenceEvidence e;
    e.alt_struct_reads = 8;
    e.alt_split_reads = 5;
    e.alt_indel_reads = 0;
    e.alt_left_clip_reads = 3;
    e.alt_right_clip_reads = 3;
    e.ref_span_reads = 2;
    e.gq = 60;
    e.score = 1.0;
    return e;
}

TEAlignmentEvidence te_pass(const std::string& qc_reason) {
    TEAlignmentEvidence t;
    t.pass = true;
    t.qc_reason = qc_reason;
    t.sequence_model_label = "TE_MODEL_IN_DISTRIBUTION";
    t.sequence_model_score = 0.9;
    return t;
}

BoundaryEvidence make_boundary(const std::string& type, int32_t len) {
    BoundaryEvidence b;
    b.geometry_defined = true;
    b.canonical_pass = true;
    b.evidence_consistent = true;
    b.boundary_type = type;
    b.boundary_len = len;
    return b;
}

// The resolved (h3) predicate: needs an insert sequence, a passing TE alignment
// with a resolved qc_reason, and adequate annotation quality.
void test_resolved_veto() {
    const auto seg = closed_segmentation();
    const auto existence = some_existence();
    const auto boundary = make_boundary("TSD", 0);
    const auto te = te_pass("PASS_INSERT_TE_ALIGNMENT");

    // Baseline: insert sequence, passing resolved TE, good annotation -> eligible.
    assert(!compute_te_resolved_hard_veto(existence, seg, te, boundary));

    // Missing insert sequence vetoes.
    auto no_insert = seg;
    no_insert.has_insert_seq = false;
    assert(compute_te_resolved_hard_veto(existence, no_insert, te, boundary));

    // A non-passing TE alignment vetoes.
    auto te_fail = te;
    te_fail.pass = false;
    assert(compute_te_resolved_hard_veto(existence, seg, te_fail, boundary));

    // The "unknown" qc_reason is not a resolved call.
    const auto te_unknown = te_pass("PASS_INSERT_TE_ALIGNMENT_UNKNOWN");
    assert(compute_te_resolved_hard_veto(existence, seg, te_unknown, boundary));

    // A SMALL_DEL boundary is no longer a structural-support gate: an otherwise
    // well-formed resolved TE stays eligible regardless of boundary type.
    const auto small_del = make_boundary("SMALL_DEL", 50);
    assert(!compute_te_resolved_hard_veto(existence, seg, te, small_del));
}

// The unknown (h2) predicate: needs an insert sequence, some TE alignment (not
// NO_TE_ALIGNMENT), and adequate annotation quality.
void test_unknown_veto() {
    const auto seg = closed_segmentation();
    const auto existence = some_existence();
    const auto boundary = make_boundary("TSD", 0);
    const auto te = te_pass("PASS_INSERT_TE_ALIGNMENT");

    // Baseline: eligible.
    assert(!compute_te_unknown_hard_veto(existence, seg, te, boundary));

    // Missing insert sequence vetoes.
    auto no_insert = seg;
    no_insert.has_insert_seq = false;
    assert(compute_te_unknown_hard_veto(existence, no_insert, te, boundary));

    // A hard "no TE alignment" reason vetoes.
    const auto te_none = te_pass("NO_TE_ALIGNMENT");
    assert(compute_te_unknown_hard_veto(existence, seg, te_none, boundary));

    // An UNKNOWN-reason call with an insert sequence and good annotation is
    // eligible for ranking; the lFDR gate, not this predicate, decides emission.
    const auto te_unknown = te_pass("PASS_INSERT_TE_ALIGNMENT_UNKNOWN");
    assert(!compute_te_unknown_hard_veto(existence, seg, te_unknown, boundary));

    // Low annotation quality vetoes.
    auto te_low_annotation = te_unknown;
    te_low_annotation.annotation_confidence = "LOW";
    te_low_annotation.annotation_residual_fraction = 0.9;
    te_low_annotation.annotation_masked_fraction = 0.9;
    assert(compute_te_unknown_hard_veto(existence, seg, te_low_annotation, boundary));
}

}  // namespace

int main() {
    test_resolved_veto();
    test_unknown_veto();
    return 0;
}
