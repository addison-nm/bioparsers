# Setup Scripts

These scripts initialize a new project after creating a repository from this template. They can be deleted once setup is complete.

## Quick Start

From the root of your new repository, run:

```bash
bash _setup_scripts/init.sh
```

This runs the setup scripts below in order.

## Scripts

### 1. Rename Project (`setup1_rename_project.sh`)

Replaces all occurrences of `myproject` with your chosen project name in source files and directories. Should be run first, before installing any dependencies.

```bash
bash _setup_scripts/setup1_rename_project.sh <name>
```

Files modified:

- `pyproject.toml`
- `src/myproject/__main__.py`
- `tests/test_core.py`
- `README.md`
- Renames `src/myproject/` to `src/<name>/`

The name must be a valid Python package name (lowercase, no hyphens, no spaces).

### 2. Create Conda Environment (`setup2_create_environment.sh`)

Creates a local conda environment at `./env` using `environment.yml`.

```bash
bash _setup_scripts/setup2_create_environment.sh
```

After creation, activate with:

```bash
conda activate ./env
```

## Cleanup

Once setup is complete, delete this directory:

```bash
rm -rf _setup_scripts
```
