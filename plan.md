# PLACER TE Calling Improvement Plan (Post-v0.0.1)

## 1. Problem Statement

Real-data behavior shows extremely low recall (about 0.1%) even after threshold relaxation.
Current evidence suggests this is mainly an architecture bottleneck, not only a parameter bottleneck:

- TE label is still coupled too early in the pipeline.
- Components can be filtered before assembly rescue can act.
- Low-quality ONT reads and TE-library mismatch can produce empty/weak TE labels early.

## 2. Evidence from Existing Methods

Common long-read TE callers use a candidate-first architecture:

- TELR: insertion candidate -> local assembly -> TE annotation -> allele/frequency estimate.
- xTea (long-read mode): candidate discovery + local assembly + annotation.
- sTELLeR: split/insertion clustering first, TE matching later.

Design implication:

- TE assignment should be posterior evidence, not a hard early gate.
- Assembly should run for insertion-like candidates even when TE family is uncertain.

## 3. Target Design (Architecture + Statistics)

### A1. TE-late pipeline (core architecture change)

- Stage 1: insertion candidate discovery and component clustering (TE-agnostic).
- Stage 2: local POA assembly using abPOA.
- Stage 3: TE annotation on consensus + fragment evidence.
- Stage 4: genotype / placeability with explicit uncertainty output.

Expected impact:

- Reduces false negatives caused by early TE-label failure.
- Keeps candidates alive for post-assembly rescue and annotation.

### A2. Replace hard TE gates with posterior TE score

For each family `f`, compute:

`log P(f | x) ~ log P(seq | f) + log P(len | f) + log P(hallmark | f) + log pi_f`

Where:

- `P(seq | f)`: k-mer/minimizer support likelihood.
- `P(len | f)`: family-specific inserted-length prior (important for ONT truncation patterns).
- `P(hallmark | f)`: TSD/polyA/orientation consistency likelihood.
- `pi_f`: prior family frequency from reference/population panel.

Output:

- `TE=UNK` allowed when posterior margin is low.
- Keep call with uncertainty instead of hard-drop.

### A3. Length-aware family modeling

Add per-family robust length priors:

- L1: heavy 5' truncation, long-tail length.
- Alu/SVA: shorter and more concentrated.

Use robust distributions (mixture/log-normal or empirical bins with smoothing) to avoid overfitting.

### A4. Robust support modeling (replace rigid support thresholds)

Current hard support filters should be complemented by probabilistic calibration:

- Support count model: beta-binomial or negative-binomial for over-dispersed coverage.
- Breakpoint spread model: robust scale (MAD/Huber) with calibrated cutoff.
- Report posterior confidence or q-value (FDR-controlled) instead of only pass/fail bits.

### A5. Denoising before clustering

For noisy support tracks in repetitive regions:

- Apply total variation / fused-lasso denoising before component split/merge.
- Use SIAM-backed optimization solvers (FISTA or Split Bregman) for stable convergence.

Expected impact:

- Fewer over-split components and fewer false merges in repeat-rich windows.

### A6. Multi-hypothesis TE annotation

Do not collapse to one TE family too early.

- Keep top-K family hypotheses with posterior scores.
- Resolve to top-1 only when margin exceeds threshold.
- Otherwise output uncertain class (`UNK` or family-set), still preserving insertion call.

## 4. Concrete Code Changes in Current Repository

Priority changes in `src/pipeline/pipeline.cpp`:

1. Remove TE-name hard requirement in evidence Stage-A gate.
   - Current behavior ties `pass_te_consistency` to non-empty `te_name`.
2. Remove TE-name hard filtering before assembly read selection.
   - Keep all valid fragments for POA; TE-hit can be used as weight, not hard filter.
3. Allow post-assembly insertion pass when TE is uncertain.
   - Convert `FAIL_NO_TE_LABEL` to `PASS_INSERTION_TE_UNCERTAIN` when non-TE evidence is strong.

Supporting updates:

- Extend `FinalCall` with uncertainty fields (`te_posterior_top1`, `te_posterior_margin`, `te_status`).
- Add per-stage reject counters in pipeline summary:
  - `reject_evidence_min_support`
  - `reject_evidence_breakpoint_mad`
  - `reject_assembly_qc`
  - `reject_te_low_support`
  - `pass_insertion_te_uncertain`

## 5. Milestones

### M3 (Architecture Decoupling)

- TE-agnostic candidate-to-assembly path enabled.
- No candidate dropped only because early `te_name` is empty.

Exit criteria:

- On real data, `assembled_calls` increases substantially from current baseline.
- Stage counter report confirms fewer early TE-driven drops.

### M4 (Posterior TE Scoring + UNK)

- Implement posterior TE scoring and top-K hypothesis output.
- Add `TE=UNK` pathway without dropping insertion call.

Exit criteria:

- Recall increases without catastrophic precision collapse.
- TE uncertain calls are explicitly quantifiable.

### M5 (Robust Statistical Filters)

- Integrate probabilistic support calibration (beta-binomial/NB).
- Add breakpoint robust-denoising option for repeat-rich bins.

Exit criteria:

- Better precision-recall stability across coverage regimes.
- Reduced sensitivity to manual threshold tuning.

### M6 (Parameter Learning)

- Learn/update calibration parameters from truth-labeled sets iteratively.
- Keep per-platform presets (ONT/PacBio) versioned in config docs.

Exit criteria:

- Reproducible gain on held-out datasets.
- Documented training/evaluation protocol and parameter snapshots.

## 6. Evaluation Protocol

Primary metrics:

- Recall, precision, F1.
- Family-level accuracy (when label is not UNK).
- Breakpoint accuracy (50/100/500 bp windows).
- Calibration metrics (Brier score / reliability bins for posterior confidence).

Required outputs:

- Stage-wise funnel table from `components` to `final_calls`.
- Separate accounting of `TE-certain` vs `TE-uncertain` insertion calls.

## 7. Key References

- TELR (NAR 2022): https://doi.org/10.1093/nar/gkac794
- xTea (Nat Commun 2021): https://doi.org/10.1038/s41467-021-24041-8
- sTELLeR (Bioinformatics 2024): https://doi.org/10.1093/bioinformatics/btae686
- Sniffles2 (Nat Biotechnol 2024): https://doi.org/10.1038/s41587-023-02024-y
- Dfam (NAR 2013): https://doi.org/10.1093/nar/gks1265
- FISTA (SIAM J Imaging Sci 2009): https://doi.org/10.1137/080716542
- Split Bregman (SIAM J Imaging Sci 2009): https://doi.org/10.1137/080725891

