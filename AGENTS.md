# Repository Guidelines

## Environment Profile
You are working on an Apple Silicon MacBook Pro (M1 Max, 32 GB RAM) with macOS Ventura or later. Prefer ARM64-native binaries, keep Python at ≥3.12, and assume Homebrew, Git, and virtualenv tooling are available. All heavy compute (FFT pipelines, bootstraps) can leverage multiprocessing; monitor memory when loading large CRISPR datasets.

## Project Structure & Module Organization
- `applications/` – CLI entry points (guide designer, scoring, visualizations). Use these for end-to-end demos.
- `experiments/signal_theoretic_crispr/` – research scripts including `breathing_dynamics.py`, `ablation_tests.py`, and benchmark runners.
- `wave_crispr_signal/` & `modules/` – reusable libraries (spectral encoders, geodesic utilities) shared across experiments.
- `proof_pack/` – reproducibility suites with bootstrap/permutation gates; run before publishing numbers.
- `docs/` – mathematical derivations and experiment reports; update when algorithms change.
- `tests/` – pytest-based unit/integration coverage mirroring module layout.

## Build, Test, and Development Commands
- `pip install -r requirements.txt` or `make install` – set up dependencies.
- `python -m pytest -q` – fast unit tests; required before every PR.
- `make test` – runs `run_tests.py` for the broader regression suite.
- `python experiments/test_breathing_dynamics_encoding.py --n-seq 50 --n-runs 5` – reproduce breathing vs. arbitrary benchmarks.
- `python proof_pack/run_validation.py` – full validation (long-running, document runtime and seed).

## Coding Style & Naming Conventions
Follow Ruff defaults (`line-length = 88`). Use type hints and docstrings on public APIs. Keep module names snake_case, classes in PascalCase, constants uppercase (`HELICAL_PERIOD_BP`). Run `ruff check .` and auto-fix before committing. Prefer deterministic RNG seeds and explicit parameter logging in experiments.

## Testing Guidelines
Write pytest tests alongside code (`tests/test_<module>.py`). For stochastic routines, assert statistical bounds and seed RNGs via `numpy.random.default_rng`. Mark any >60 s tests with `pytest.mark.slow`. When adding new scientific gates, replicate patterns in `tests/test_breathing_dynamics_integration.py`.

## Commit & Pull Request Expectations
Use imperative, concise commit messages (“Add CZT harmonic metrics”). PR descriptions must include: summary, validation commands executed, datasets touched, and updated docs/figures. Link related issues, attach result tables or plots when metrics move, and call out long-running proof-pack jobs. Update documentation (`docs/`, `CHARTER.md`) whenever scientific assumptions shift.
