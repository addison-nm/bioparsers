"""Make the recipe scripts importable from builder tests.

Concrete builders live in the top-level ``recipes/`` scripts (each recipe
demonstrates defining a custom Builder), so the unit tests import them from
there by adding ``recipes/`` to ``sys.path``. (The recipes' shared
Pfam runner now lives in the package — ``bioparsers.builders.uniprot`` —
so no cross-directory import is involved.)
"""

import os
import sys

RECIPES = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "recipes"))
if RECIPES not in sys.path:
    sys.path.insert(0, RECIPES)
