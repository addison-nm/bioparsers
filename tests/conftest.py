"""Pytest Configuration File

"""

from pathlib import Path

import pytest
import shutil

DATDIR = Path(__file__).parent / "data"      # data directory for all tests.
TEST_OUTDIR = Path(__file__).parent / "_tmp"  # output directory for all tests.

def remove_dir(dir: str | Path):
    """Helper function to remove a directory recursively."""
    dir = Path(dir).resolve()
    if not dir.is_relative_to(TEST_OUTDIR.resolve()):
        msg = f"Can only use function `remove_dir` from tests.conftest to \
        remove directories in the directory {TEST_OUTDIR}. Got: {dir}"
        raise RuntimeError(msg)
    shutil.rmtree(dir)

#####################
##  Configuration  ##
#####################

def pytest_addoption(parser):
    parser.addoption(
        "--benchmark", action="store_true", default=False, 
        help="run benchmarking tests"
    )

def pytest_collection_modifyitems(config, items):
    benchmark_flag_given = False
    if config.getoption("--benchmark"):
        # --benchmark given in cli: do not skip benchmarking tests
        benchmark_flag_given = True
    skip_benchmark = pytest.mark.skip(reason="need --benchmark option to run")
    for item in items:
        if "benchmark" in item.keywords and not benchmark_flag_given:
            item.add_marker(skip_benchmark)
