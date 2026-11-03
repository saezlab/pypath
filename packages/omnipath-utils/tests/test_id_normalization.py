"""One canonical form per chemical identifier namespace (spec 011 T024-T026).

WP1's whole join depends on this: the resolver's chemical lookup joins on
the normalized key, and a namespace whose translation-side and build-side
forms disagree never overlaps -- which is exactly the defect that survived
a full cycle unmeasured. These tests are the value_pattern each namespace's
normalized form must match, and the normalize() call that must produce it
from every raw form actually seen in the field.
"""

from __future__ import annotations

import re

import pytest

from omnipath_utils.mapping._id_types import normalize_identifier, value_pattern


def _assert_normalizes(id_type: str, raw: str, expected: str) -> None:
    got = normalize_identifier(id_type, raw)
    assert got == expected, f'{id_type} {raw!r} -> {got!r}, expected {expected!r}'
    pattern = value_pattern(id_type)
    assert pattern is not None, f'{id_type} has no registered value_pattern'
    assert re.fullmatch(pattern, got), (
        f'{id_type} normalized form {got!r} does not match its own '
        f'value_pattern {pattern!r}'
    )


@pytest.mark.parametrize(
    'raw, expected',
    [
        ('15377', 'CHEBI:15377'),
        ('chebi:15377', 'CHEBI:15377'),
        ('CHEBI:15377', 'CHEBI:15377'),
        ('Chebi:15377', 'CHEBI:15377'),
        (' CHEBI:15377 ', 'CHEBI:15377'),
    ],
)
def test_chebi_prefix(raw, expected):
    _assert_normalizes('chebi', raw, expected)


@pytest.mark.parametrize(
    'raw, expected',
    [
        ('HMDB00001', 'HMDB0000001'),  # old 5-digit form
        ('HMDB0000001', 'HMDB0000001'),  # current 7-digit form
        ('hmdb0000001', 'HMDB0000001'),
        ('HMDB1', 'HMDB0000001'),
    ],
)
def test_hmdb_padding(raw, expected):
    _assert_normalizes('hmdb', raw, expected)


@pytest.mark.parametrize(
    'raw, expected',
    [
        ('000000001', 'SLM:000000001'),
        ('SLM:000000001', 'SLM:000000001'),
        ('slm:000000001', 'SLM:000000001'),
    ],
)
def test_swisslipids_prefix(raw, expected):
    _assert_normalizes('swisslipids', raw, expected)


@pytest.mark.parametrize(
    'raw, expected',
    [
        ('lmfa01010001', 'LMFA01010001'),
        ('LMFA01010001', 'LMFA01010001'),
    ],
)
def test_lipidmaps_form(raw, expected):
    _assert_normalizes('lipidmaps', raw, expected)


@pytest.mark.parametrize(
    'raw, expected',
    [
        ('C2', 'C00002'),
        ('c00002', 'C00002'),
        ('C00002', 'C00002'),
    ],
)
def test_kegg_form(raw, expected):
    _assert_normalizes('kegg', raw, expected)


@pytest.mark.parametrize(
    'raw, expected',
    [
        ('0000123', '123'),
        ('123', '123'),
        (' 123 ', '123'),
    ],
)
def test_pubchem_leading_zeros(raw, expected):
    _assert_normalizes('pubchem', raw, expected)


@pytest.mark.parametrize(
    'raw, expected',
    [
        ('chembl123', 'CHEMBL123'),
        ('CHEMBL123', 'CHEMBL123'),
        ('CHEMBL 123', 'CHEMBL123'),
    ],
)
def test_chembl_form(raw, expected):
    _assert_normalizes('chembl', raw, expected)


@pytest.mark.parametrize(
    'raw, expected',
    [
        ('50-00-0', '50-00-0'),
        (' 50-00-0 ', '50-00-0'),
    ],
)
def test_cas_form(raw, expected):
    _assert_normalizes('cas', raw, expected)


@pytest.mark.parametrize(
    'raw, expected',
    [
        ('bqjcrhhnabkakn-kbqpjgbksa-n', 'BQJCRHHNABKAKN-KBQPJGBKSA-N'),
        ('BQJCRHHNABKAKN-KBQPJGBKSA-N', 'BQJCRHHNABKAKN-KBQPJGBKSA-N'),
        ('InChIKey=BQJCRHHNABKAKN-KBQPJGBKSA-N', 'BQJCRHHNABKAKN-KBQPJGBKSA-N'),
    ],
)
def test_structure_key_case_and_prefix(raw, expected):
    _assert_normalizes('inchikey', raw, expected)


def test_unknown_id_type_returns_none():
    assert normalize_identifier('not_a_real_type', 'whatever') is None


def test_empty_value_returns_none():
    assert normalize_identifier('chebi', '') is None
    assert normalize_identifier('chebi', None) is None


def test_normalize_identifier_is_the_one_reusable_implementation():
    """T026: other packages import this, not reimplement it.

    Asserted here as an import-surface contract -- the function omnipath-build
    is meant to call lives at this exact path, with this exact name.
    """
    import omnipath_utils.mapping._id_types as mod

    assert callable(mod.normalize_identifier)
    assert callable(mod.value_pattern)
