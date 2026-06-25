"""Make the recipe scripts importable from builder tests.

Concrete builders live in the top-level ``recipes/`` scripts (each recipe
demonstrates defining a custom Builder), so the unit tests import them from
there. Adding ``recipes/`` to ``sys.path`` also lets the recipes' own
``from _pfam_runner import ...`` resolve on import.
"""

import os
import sys

RECIPES = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "recipes"))
if RECIPES not in sys.path:
    sys.path.insert(0, RECIPES)
