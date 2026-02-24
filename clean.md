# PLACER Clean Architecture Plan

## 1. Goal

Keep detection sensitivity under noisy long-read and weak TE-reference conditions, while reducing runtime branch complexity.

Target: make runtime pipeline minimal and stable, move iterative learning/library growth to offline.

## 2. Design Principles

- Separate concerns: insertion existence and TE family identity are different tasks.
- Keep runtime deterministic: avoid deep branch coupling and ad-hoc fallback paths.
- Open-set first: allow `TE_UNCERTAIN`/`NON_TE` without dropping insertion calls.
- Offline iteration: bootstrap library and parameter learning should not complicate hot path.

## 3. Target Pipeline (Runtime-3 + Offline-1)

### Runtime Stage A: TE-agnostic candidate/evidence

- Build components from split/SA/indel signals.
- Apply core evidence hard filters.

### Runtime Stage B: Local assembly + insertion acceptance

- Assemble with abPOA.
- Decide insertion acceptance from placeability quality only.

### Runtime Stage C: Open-set TE classification

- Classify accepted insertions into `TE_CERTAIN` / `TE_UNCERTAIN` / `NON_TE`.
- TE classification must not back-propagate to reject insertion.

### Offline Stage D: Bootstrap refinement

- Export uncertain/non-TE consensus.
- Cluster/rebuild TE library and rerun pass-2 outside runtime.

## 4. Implementation Tasks

- [x] C1. Add explicit insertion-acceptance decision function (placeability-driven).
- [x] C2. Add explicit open-set TE classification function.
- [x] C3. Refactor `process_bin_records` to execute in order:
  - evidence -> assembly -> insertion acceptance -> TE open-set -> genotype/output.
- [x] C4. Keep existing bootstrap export as offline interface.
- [x] C5. Add dedicated insertion QC field in output (separate from TE QC).
- [x] C6. Add focused unit tests for C1/C2 behavior matrix.

## 5. Current Code Mapping

- `src/pipeline/decision_policy.cpp`
  - `evaluate_insertion_acceptance(...)`: insertion pass/high-low confidence.
  - `classify_te_open_set(...)`: open-set TE status decision matrix.
- `src/pipeline/pipeline.cpp`
  - `evaluate_insertion_acceptance_for_pipeline(...)`: adapter from placeability report to policy.
  - `apply_te_open_set_classification(...)`: adapter from policy output to `FinalCall`.
  - `process_bin_records(...)`: decoupled execution order and reduced branch coupling.

## 6. Exit Criteria

- Runtime has no TE-gate hard dependency for insertion acceptance.
- TE classification remains auditable via `te_status`, posterior, and QC tags.
- Existing outputs remain backward compatible for current downstream scripts.
