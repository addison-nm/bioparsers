# Template Modernization Session

**Date:** 2026-03-27

## Summary

Reviewed and modernized the Python project template. Changes focused on removing legacy tooling, improving shell script robustness, adding linting/formatting, and fleshing out project structure.

## Changes Made

### Tooling
- **Removed `setup.py`** — consolidated everything into `pyproject.toml`
- **Added ruff** (`>=0.9`) as dev dependency with config for pyflakes, pycodestyle, isort, pyupgrade, and flake8-bugbear rules
- **Added `[tool.pytest.ini_options]`** to `pyproject.toml`, removed redundant `pytest_configure` hook from `conftest.py`
- **Added `Makefile`** with `install`, `test`, `lint`, `format`, `check` targets
- **Added `.python-version`** (3.12) for pyenv users

### Bug Fixes
- **Rename script input validation** — rejects names not matching `^[a-z][a-z0-9_]*$`
- **Shell script quoting** — properly quoted `${flist[@]}`, loop variables, and `mv` arguments
- **`TMPDIR` renamed to `TEST_OUTDIR`** in `conftest.py` to avoid shadowing POSIX env var
- **Relative paths in `conftest.py`** — replaced string paths with `Path(__file__).parent / ...` so tests work regardless of working directory

### pyproject.toml
- Set explicit `version = "0.0.1"` (was `dynamic`)
- Bumped `requires-python` to `>= 3.10` to match actual dependency floors
- Added `[tool.setuptools.packages.find]` with `where = ["src"]`

### Project Structure
- Added template directories with READMEs: `data/`, `outputs/`, `scripts/`, `logs/`, `notebooks/`, `docs/`
- Added `docs/troubleshooting.md` with common environment, testing, and CLI issues
- Updated `.gitignore` to exclude contents of `data/`, `outputs/`, `logs/` (keeping READMEs)
- Removed `setup2_make_directories.sh` (directories now ship with template)
- Renamed `setup3_create_environment.sh` to `setup2_create_environment.sh`

### Documentation
- Rewrote `_setup_scripts/README.md` with full usage instructions for each script
- Populated `README.md` with setup, usage, and development sections

## Decisions

- **Kept bash script approach** for project renaming (vs GitHub Actions auto-rename or cookiecutter) — simpler and more transparent for a template repo
- **Did not specify environment name** in `environment.yml` — intentional, left to user
- **LICENSE file** not included in rename script — doesn't contain "myproject"

## Still TODO

- GitHub Actions CI workflow
- Pre-commit hooks (`.pre-commit-config.yaml` with ruff)
