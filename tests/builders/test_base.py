"""Tests for the Builder ABC's subclass-definition enforcement."""

import pytest

from bioparsers.builders.base import Builder


class TestBuilderEnforcement:

    def test_missing_name_raises(self):
        with pytest.raises(TypeError, match="name"):
            class NoName(Builder):
                description = "documented"
                def build(self, records):
                    return iter(())

    def test_missing_description_raises(self):
        with pytest.raises(TypeError, match="document"):
            class NoDesc(Builder):
                name = "x"
                def build(self, records):
                    return iter(())

    def test_docstring_used_as_description(self):
        class WithDoc(Builder):
            """This docstring documents the output form."""
            name = "with_doc"
            def build(self, records):
                return iter(())

        assert WithDoc.description == "This docstring documents the output form."

    def test_abstract_intermediate_is_exempt(self):
        # A subclass that does not implement build is still abstract and
        # need not declare name/description.
        class Intermediate(Builder):
            pass

        assert getattr(Intermediate.build, "__isabstractmethod__", False)

    def test_callable_delegates_to_build(self):
        class Echo(Builder):
            """Yields records unchanged."""
            name = "echo"
            def build(self, records):
                yield from records

        recs = [{"a": 1}, {"a": 2}]
        assert list(Echo()(recs)) == recs
