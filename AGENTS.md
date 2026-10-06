# Repository Guidelines

## Project Structure & Module Organization

- `01-scripts/` — primary MATLAB pipeline scripts (run roughly A → I) plus utilities.
  - `01-scripts/HASH/` — HASH inputs/outputs; Fortran sources live in `01-scripts/HASH/code/`.
  - `01-scripts/SKHASH/` — SKHASH Python runner and example control files.
  - `01-scripts/python/` — Jupyter notebooks + helpers for the Y1+Y2 workflow.
- `02-data/` — inputs and intermediate artifacts (mostly generated binaries like `*.mat`, kept out of git).
- `03-output-graphics/` — generated plots (`MyPng/`, `MyPdf/`; gitignored).
- Root files: `FM_buildpath4.m` (MATLAB path bootstrap), `savemyfigureFM4.m` (figure saving), and `CLAUDE.md` (detailed architecture/runbook).

## Build, Test, and Development Commands

- MATLAB (main workflow):
  - `run FM_buildpath4.m`
  - Run stages directly, e.g. `run 01-scripts/G_FM_SKHASH_22OBS.m`
- SKHASH (Python, bundled):
  - `cd 01-scripts/SKHASH && python SKHASH.py examples/hash3/control_file.txt`
- HASH (optional, Fortran build):
  - `cd 01-scripts/HASH/code && make hash_driver3` (requires `gfortran`)

## Coding Style & Naming Conventions

- MATLAB: use 4-space indentation, keep stage prefixes in filenames (`A_...`, `B_...`), and prefer clear structs/tables over `eval` when practical.
- Python: 4-space indentation and `snake_case` modules/functions.
- Paths/config: avoid adding new absolute paths. When touching code with hardcoded paths, centralize them (e.g., `FM_buildpath4.m`, `01-scripts/python/pipeline_config.py`) and prefer repo-root-relative paths (`fullfile`, `pathlib.Path`).

## Testing Guidelines

- No formal unit test suite; validate changes by running the smallest affected stage and confirming expected outputs land in `02-data/<stage>/` or `03-output-graphics/`.
- Keep generated artifacts out of git (`*.mat`, `*.png`, `*.pdf`, `*.mseed`). If a new artifact type is introduced, update `.gitignore` instead of committing it.

## Commit & Pull Request Guidelines

- Existing history uses short, imperative commit messages (e.g., “Add …”). Keep commits focused and scoped to one logical change.
- PRs should include: what stage(s) changed, how you validated (commands + key output filenames), and any reproducibility notes (path changes, external dependencies). Include screenshots only when plot appearance changes.
