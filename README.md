# myproject

## Description

TODO: Describe your project here.

## Setup

Create and activate the conda environment:

```bash
bash _setup_scripts/init.sh
conda activate ./env
```

Or install manually with pip:

```bash
pip install -e .[dev]
```

## Usage

Run the CLI:

```bash
myproject
```

Or as a module:

```bash
python -m myproject
```

## Development

```bash
make test      # run tests
make lint      # check code with ruff
make format    # format code with ruff
make check     # lint + test
```
