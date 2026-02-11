# Repository Guidelines

## Project Structure & Module Organization
- `src/` contains phase-based C++ modules: `stream/`, `gate1/`, `component/`, `local_realign/`, `assembly/`, `placeability/`, `genotyping/`, and `te_reverse_index/`. `src/main.cpp` wires the end-to-end pipeline.
- `include/` stores shared public headers (for example, `window_buffer.h`, `placeability.h`). Keep header and implementation names aligned.
- `tests/` contains CTest-registered binaries named `test_*.cpp`.
- `test_data/` holds small local fixtures (for example `te_test.fa`), and `third_party/abPOA/` vendors the abPOA dependency.
- `doc/` and `scripts/` contain planning notes and environment-specific run scripts.

## Build, Test, and Development Commands
- `python3 -m pip install -r requirements.txt` installs Python dependencies (notably `pysam` headers used by CMake).
- `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release` configures the project.
- `cmake --build build -j` compiles `placer` and all module/test targets.
- `ctest --test-dir build --output-on-failure` runs the full CTest suite.
- `./build/placer <input.bam> <ref.fa> [te.fa]` runs the pipeline locally.
- `./build.sh` wraps configure/build/test, but includes a hardcoded demo path; edit before reuse.

## Coding Style & Naming Conventions
- Language baseline: C++17 (and C11 where applicable), defined in `CMakeLists.txt`.
- Follow existing style: 4-space indentation, concise functions, and readable early returns.
- Use `snake_case` for files, functions, and variables; use `PascalCase` for types/classes (for example `PipelineConfig`).
- Add new module code under `src/<module>/` and expose stable interfaces in `include/`.
- No repository-wide formatter/linter is configured; match adjacent code style and include order.

## Testing Guidelines
- Tests are plain C++ executables with `assert` checks, registered in `tests/CMakeLists.txt`.
- Add tests as `tests/test_<feature>.cpp` and register each target with `add_test(...)`.
- Run focused tests during iteration, then run full `ctest` before opening a PR.
- Several current tests use absolute `/mnt/home1/...` paths; switch to local fixtures or guarded paths for portable development.

## Commit & Pull Request Guidelines
- Keep commit subjects short, imperative, and phase/module scoped (history examples: `Add test data`, `Fix phase 1 bugs`).
- Prefer one logical change per commit; avoid mixing refactors and behavior changes.
- PR descriptions should include: purpose, affected modules/phases, and exact validation commands run.
- If output semantics change, include a small `scientific.vcf` excerpt or summary of field-level differences.
- Call out environment assumptions (for example, htslib location, reference indexes) so reviewers can reproduce.
