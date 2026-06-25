"""UniProt-specific builder helpers and filters.

These operate on the UniProt record dict shape emitted by
``bioparsers.parsers.uniprot_dat`` (``comments`` as ``{topic, text}``,
``description.rec_name``, ``status``, ``cross_references``, ...), and on
UniProtKB's text conventions (``{ECO:...}`` evidence, ``(PubMed:...)``
citations). A different source database (e.g. Pfam) would get its own
sibling subpackage with its own helpers/filters; the generic framework
(``Builder``, JSONL I/O) stays at ``bioparsers.builders``.
"""

from bioparsers.builders.uniprot import filters, helpers

__all__ = ["helpers", "filters"]
