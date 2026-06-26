"""Make the test-fixture recipe importable from the builder tests.

The builder tests exercise a representative ``Builder`` from a self-contained
copy under ``tests/builders/_recipes/`` — kept here so the suite does **not**
depend on the live scripts in the top-level ``recipes/`` directory (which are
user-editable and get renamed/reorganized). The fixtures dir is added to
``sys.path`` so the copy imports as a normal module.
"""

import os
import sys

RECIPES = os.path.join(os.path.dirname(__file__), "_recipes")
if RECIPES not in sys.path:
    sys.path.insert(0, RECIPES)
