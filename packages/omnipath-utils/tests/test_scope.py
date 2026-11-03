"""Tests for the 007 build-scope resolution (``db/_presets.py``).

Pure Python, no DB — validates the nested tiers and the explicit
taxid/name-list parsing of ``--scope``.
"""

from __future__ import annotations

import pytest

from omnipath_utils.db._presets import (
    SCOPES,
    CORE_MODEL,
    EXTENDED_MODEL,
    MODEL_ORGANISMS,
    resolve_scope,
)


class TestScopeTiers:
    """The named tiers are nested supersets."""

    def test_only_human(self):
        assert resolve_scope('only-human') == [9606]

    def test_core_model(self):
        assert resolve_scope('core-model') == [9606, 10090, 10116]

    def test_nesting_is_superset(self):
        # only-human ⊂ core-model ⊂ extended-model ⊂ model-organisms
        assert set([9606]).issubset(CORE_MODEL)
        assert set(CORE_MODEL).issubset(EXTENDED_MODEL)
        assert set(EXTENDED_MODEL).issubset(MODEL_ORGANISMS)

    def test_model_organisms_membership(self):
        # the 007 clarification set: incl. pig, rhesus, E. coli, dog
        for taxid in (9606, 10090, 10116, 9823, 9544, 8364,
                      7955, 7227, 6239, 4932, 511145, 9615):
            assert taxid in MODEL_ORGANISMS

    def test_complete_is_none(self):
        assert resolve_scope('complete') is None
        assert SCOPES['complete'] is None

    def test_case_insensitive(self):
        assert resolve_scope('Core-Model') == resolve_scope('core-model')


class TestExplicitList:
    """Explicit taxid / name lists."""

    def test_taxid_list(self):
        assert resolve_scope('9606,9823,9598') == [9606, 9823, 9598]

    def test_name_list(self):
        assert resolve_scope('human,pig,chimpanzee') == [9606, 9823, 9598]

    def test_mixed_and_semicolons(self):
        assert resolve_scope('human; 9823 ;chimp') == [9606, 9823, 9598]

    def test_dedup_preserves_order(self):
        assert resolve_scope('human,mouse,human') == [9606, 10090]

    def test_unknown_name_raises(self):
        with pytest.raises(ValueError):
            resolve_scope('human,notanorganism')

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            resolve_scope(',,')

    def test_none_passes_through(self):
        assert resolve_scope(None) is None
