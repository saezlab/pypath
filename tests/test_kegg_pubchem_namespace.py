"""KEGG's conv/pubchem cross-reference is a PubChem substance ID, not a
compound ID (spec 011 T032, spec.md US1 acceptance scenario 3).

KEGG mints its own compound accessions and only cites PubChem through the
`conv/pubchem` endpoint, which returns substance IDs (SID), not compound IDs
(CID). The resource's own ``mints=`` declaration confirms this, as does the
comment directly above the schema this test inspects. Tagging that value as
a PubChem *compound* identifier feeds a wrong structure candidate into the
resolver, since WP1's chemical projection joins compound IDs straight to
InChIKeys. This is a pre-existing, self-documented bug in the schema
declaration, not new field-extraction logic.
"""

from __future__ import annotations

from typing import Any

from pypath.internals.cv_terms import IdentifierNamespaceCv
from pypath.inputs_v2.kegg import reactions_schema


def _term_of(cv) -> Any:
    return cv.term_source.value


def _member_cvs_for(field_name: str):
    """CV specs anywhere in the reactions schema whose value column
    selects ``field_name``."""

    found = []
    for member in reactions_schema.membership.members:
        for cv in member.identifiers.cvs:
            value_source = cv.value_source
            if getattr(value_source, 'selector', None) == field_name:
                found.append(cv)
    return found


def test_reactant_pubchem_is_tagged_as_substance():
    cvs = _member_cvs_for('reactant_pubchem')
    assert cvs, 'no CV reads reactant_pubchem -- has the field been renamed?'
    for cv in cvs:
        assert _term_of(cv) is IdentifierNamespaceCv.PUBCHEM_SUBSTANCE, (
            f'reactant_pubchem tagged as {_term_of(cv)!r}, expected '
            f'PUBCHEM_SUBSTANCE -- KEGG conv/pubchem returns SIDs, not CIDs'
        )


def test_product_pubchem_is_tagged_as_substance():
    cvs = _member_cvs_for('product_pubchem')
    assert cvs, 'no CV reads product_pubchem -- has the field been renamed?'
    for cv in cvs:
        assert _term_of(cv) is IdentifierNamespaceCv.PUBCHEM_SUBSTANCE, (
            f'product_pubchem tagged as {_term_of(cv)!r}, expected '
            f'PUBCHEM_SUBSTANCE -- KEGG conv/pubchem returns SIDs, not CIDs'
        )


def test_kegg_does_not_mint_pubchem_compound_identifiers():
    """The resource's own ``mints=`` declaration already says as much."""

    from pypath.inputs_v2.kegg import config

    assert IdentifierNamespaceCv.PUBCHEM_COMPOUND not in config.mints
