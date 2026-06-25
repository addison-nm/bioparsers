"""Abstract base for dataset builders.

A :class:`Builder` turns a stream of parsed UniProt records (plain dicts,
as emitted by the parser's JSONL output) into a stream of curated output
records (also plain dicts). It is the dataset-layer analogue of a parser's
``Record`` subclass: each concrete builder is a single, self-describing
unit.

Every concrete builder must declare two things, enforced at
subclass-definition time (fail-loud, mirroring ``Record``/``SchemaError``):

- ``name``: a stable, versioned identifier (e.g. ``"uniprot_flat_demo"``).
- ``description``: a long-form text description that documents the shape
  of each output record. If a subclass does not set ``description``
  explicitly, its class docstring is used. A builder with neither is a
  bug and raises ``TypeError``.

Builders are streaming-first: :meth:`Builder.build` returns an iterator and
must not materialize the whole input. Configuration (filters, options) is
passed to ``__init__``; the input record stream is passed to ``build``.
"""

from __future__ import annotations

import inspect
from abc import ABC, abstractmethod
from typing import ClassVar, Iterable, Iterator


class Builder(ABC):
    """Base class for a dataset composition.

    Subclass it, set a versioned ``name``, document the output record form
    in a long-form ``description`` (or the class docstring), and implement
    :meth:`build`.
    """

    #: Stable, versioned identifier, e.g. ``"uniprot_flat_demo"``.
    name: ClassVar[str] = ""
    #: Long-form description of the output record form. Defaults to the
    #: class docstring when not set explicitly.
    description: ClassVar[str] = ""

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        # Skip still-abstract intermediates (``build`` not yet implemented);
        # only finished builders must declare their identity and output
        # contract. ``__abstractmethods__`` isn't populated until after
        # ``__init_subclass__`` runs, so test ``build`` directly.
        if getattr(cls.build, "__isabstractmethod__", False):
            return
        if not cls.name:
            raise TypeError(
                f"{cls.__name__} must set a non-empty class-level `name`"
            )
        text = cls.description or (cls.__doc__ or "")
        if not text.strip():
            raise TypeError(
                f"{cls.__name__} must document its output record form via a "
                "`description` or a class docstring"
            )
        cls.description = inspect.cleandoc(text)

    @abstractmethod
    def build(self, records: Iterable[dict]) -> Iterator[dict]:
        """Yield curated output dicts from parsed UniProt record dicts."""

    def __call__(self, records: Iterable[dict]) -> Iterator[dict]:
        return self.build(records)
