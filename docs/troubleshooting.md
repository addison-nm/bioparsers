# Troubleshooting

## Environment Issues

### Conda environment fails to create

Make sure conda is installed and available on your `PATH`. Then retry:

```bash
conda env create -p ./env -f environment.yml --yes
```

If dependency resolution fails, try updating conda first:

```bash
conda update -n base conda
```

### Package not found after install

Ensure the environment is activated and the package is installed in editable mode:

```bash
conda activate ./env
pip install -e .[dev]
```

## Testing Issues

### Tests fail with `ModuleNotFoundError`

The package is likely not installed. Run:

```bash
pip install -e .[dev]
```

### Tests pass locally but fail in CI

Check that the Python version and dependencies match between your local environment and CI. The `.python-version` and `environment.yml` files define the expected versions.

## Common Errors

### `command not found: myproject`

The CLI entry point is only available after installing the package:

```bash
pip install -e .[dev]
myproject
```

### `ruff: command not found`

Ruff is included in the dev dependencies. Install with:

```bash
pip install -e .[dev]
```
