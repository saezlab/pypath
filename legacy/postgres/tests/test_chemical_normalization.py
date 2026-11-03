"""normalize_chemical_identifier mirrors omnipath-utils' normalizers, one
per build identifier type (spec 011 T040).
"""

from __future__ import annotations

import pytest

from omnipath_build.resolver.chemical_normalization import (
    is_chemical_identifier_type,
    normalize_chemical_identifier,
)


@pytest.mark.parametrize(
    'identifier_type, raw, expected',
    [
        ('Chebi:MI:0474', 'chebi:15377', 'CHEBI:15377'),
        ('Chebi:MI:0474', '15377', 'CHEBI:15377'),
        ('Hmdb:OM:0004', 'HMDB00001', 'HMDB0000001'),
        ('Kegg Compound:MI:2012', 'c00002', 'C00002'),
        ('Pubchem Compound:OM:0002', '0000123', '123'),
        (
            'Standard Inchi Key:MI:1101',
            'InChIKey=BQJCRHHNABKAKN-KBQPJGBKSA-N',
            'BQJCRHHNABKAKN-KBQPJGBKSA-N',
        ),
    ],
)
def test_normalizes_by_build_identifier_type(identifier_type, raw, expected):
    assert normalize_chemical_identifier(identifier_type, raw) == expected


def test_non_chemical_identifier_type_returns_none():
    assert normalize_chemical_identifier('Uniprot:MI:0486', 'P12345') is None


def test_empty_value_returns_none():
    assert normalize_chemical_identifier('Chebi:MI:0474', '') is None
    assert normalize_chemical_identifier('Chebi:MI:0474', None) is None


def test_malformed_value_returns_none():
    """A value that never matches its namespace's raw form -- the raw
    value must be kept by the caller, not replaced with garbage."""
    assert normalize_chemical_identifier('Chebi:MI:0474', 'not-a-chebi-id') is None


@pytest.mark.parametrize(
    'identifier_type',
    [
        'Chebi:MI:0474', 'Hmdb:OM:0004', 'Kegg Compound:MI:2012',
        'Pubchem Compound:OM:0002', 'Standard Inchi Key:MI:1101',
    ],
)
def test_is_chemical_identifier_type_true_for_chemical_namespaces(identifier_type):
    assert is_chemical_identifier_type(identifier_type)


def test_is_chemical_identifier_type_false_for_non_chemical():
    assert not is_chemical_identifier_type('Uniprot:MI:0486')
    assert not is_chemical_identifier_type(None)
