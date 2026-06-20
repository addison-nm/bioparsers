"""Entry point for ``python -m bioparsers``.

The console-script entry point (``bioparsers`` on the PATH) is wired
directly to :func:`bioparsers.main.main` in ``pyproject.toml``.
"""

import sys

from bioparsers.main import main

if __name__ == "__main__":
    sys.exit(main())
