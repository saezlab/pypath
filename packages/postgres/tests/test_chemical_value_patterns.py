"""The syntax of each chemical namespace, as the vocabulary stores it for the coverage check.

The check is ``canonical_identifier !~ value_pattern``: PostgreSQL's match is unanchored, so a
pattern that is not anchored accepts anything that contains a valid value, such as two
identifiers joined by a pipe.
"""

import re

import pytest

from omnipath_postgres.relational.db import schema


def stored(slug):
    return schema.chemical_value_pattern(slug)


@pytest.mark.parametrize("value", ["CHEBI:15377", "CHEBI:1"])
def test_a_published_chebi_identity_matches(value):
    assert re.search(stored("chebi"), value)


@pytest.mark.parametrize(
    "value",
    ["15377", "chebi:15377", "CHEBI:", "CHEBI:15377 | CHEBI:15378", "74546 | 82922", "xCHEBI:1", "CHEBI:1x"],
)
def test_what_is_not_a_published_chebi_identity_does_not(value):
    assert not re.search(stored("chebi"), value)


def test_two_joined_identifiers_never_match_any_namespace():
    joined = {
        "swisslipids": "SLM:000000239 | SLM:000057063",
        "kegg": "C00001|C00002",
        "hmdb": "HMDB0000001 | HMDB0000002",
        "pubchem": "1 | 2",
        "chembl": "CHEMBL1 | CHEMBL2",
        "lipidmaps": "LMFA01010001 | LMFA01010002",
        "inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-C | AAAAAAAAAAAAAA-BBBBBBBBBB-D",
    }
    for slug, value in joined.items():
        assert not re.search(stored(slug), value), slug


def test_every_pattern_is_anchored_at_both_ends():
    for slug in schema._CHEMICAL_VALUE_PATTERNS:
        pattern = stored(slug)
        assert pattern.startswith("^") and pattern.endswith("$"), slug


def test_the_alternation_of_pubchem_does_not_escape_its_anchors():
    assert re.search(stored("pubchem"), "0") and re.search(stored("pubchem"), "12")
    assert not re.search(stored("pubchem"), "x0") and not re.search(stored("pubchem"), "0x")
