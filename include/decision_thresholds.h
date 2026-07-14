#ifndef PLACER_DECISION_THRESHOLDS_H
#define PLACER_DECISION_THRESHOLDS_H

namespace placer {

// Decision-policy constants that survive the move to a calibrated,
// lFDR-primary emission gate.
//
// The read-count / genotype-quality / insert-length / segmentation-score
// *ladders* that used to live here (kAltReads*, kGq*, kInsertLen*, kSegScore*,
// kPrecise*, precise-fraction and alt/ref-ratio floors) have been removed. They
// were hand-drawn decision surfaces fit one site at a time to a benchmark, and
// final TE emission no longer depends on them: it is gated by the robust
// worst-case local-FDR (see evaluate_joint_hypotheses / evaluate_robust_
// mechanistic_lfdr) against a single target risk kTargetQ.
//
// What remains are not feature thresholds but decision-theoretic constants:
//   - kTargetQ           the one policy knob, a target false-call risk (FDR).
//   - the log-evidence cut-offs are Bayes-factor scale points (log e ~ 2 is
//     "decisive" on the Jeffreys/Kass-Raftery scale), not tuned read counts.
//   - the artifact-posterior caps are probabilities, not feature ladders.
//   - kScore* are the linear weights that *rank* competing hypotheses (h0..h3);
//     ranking is diagnostic and does not gate emission.
struct DecisionThresholds {
    // The single emission knob: emit a final TE call iff the robust worst-case
    // local FDR is at or below this target false-call risk. Precision-first
    // default; a strict mode can lower it to 0.05.
    static constexpr double kTargetQ = 0.10;

    // Artifact-posterior guards (probabilities).
    static constexpr double kArtifactPosteriorCap = 0.65;
    static constexpr double kArtifactPosteriorStrictCap = 0.35;

    // Final log-evidence (Bayes-factor) decision cut-offs.
    static constexpr double kStructuralLogEvidence = 2.0;

    // Linear weights that combine the per-signal supports into each competing
    // hypothesis score (h0..h3) in evaluate_joint_hypotheses. This score is the
    // cross-candidate ranking key (joint.best.total selects the winning
    // breakpoint candidate); it is diagnostic and does not gate emission.
    // Negative values are penalties.
    static constexpr double kScoreRefExistenceWeight = -1.0;
    static constexpr double kScoreRefSegmentationWeight = -0.5;
    static constexpr double kScoreNonTeSegmentationWeight = 0.7;
    static constexpr double kScoreNonTeModelWeight = -0.4;
    static constexpr double kScoreTeSegmentationWeight = 0.8;
    static constexpr double kScoreTeBoundaryWeight = 0.2;
};

}  // namespace placer

#endif  // PLACER_DECISION_THRESHOLDS_H
