"""Tests for the ``mints`` declaration on ``ResourceConfig``.

A resource description states which identifier namespaces it is the minting
authority for. The authority relation (``identifier_authority`` in
omnipath-build) is populated from this declaration, not inferred from what a
resource's schema happens to tag data as — a resource can carry cross-
references to a namespace without minting it.
"""

from __future__ import annotations

from pypath.inputs_v2.base import ResourceConfig
from pypath.internals.cv_terms import (
    IdentifierNamespaceCv,
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)


def _config(**overrides) -> ResourceConfig:
    defaults = dict(
        id=ResourceCv.CHEBI,
        name='Test Resource',
        url='https://example.org',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        description='A resource used only to exercise the mints field.',
    )
    defaults.update(overrides)
    return ResourceConfig(**defaults)


def test_mints_defaults_to_empty():
    """A resource that declares nothing mints nothing — a pure consumer."""
    config = _config()
    assert config.mints == ()


def test_mints_declares_the_namespaces_a_resource_is_authoritative_for():
    config = _config(mints=(IdentifierNamespaceCv.CHEBI,))
    assert config.mints == (IdentifierNamespaceCv.CHEBI,)


def test_mints_holds_more_than_one_namespace():
    config = _config(
        mints=(
            IdentifierNamespaceCv.KEGG_COMPOUND,
            IdentifierNamespaceCv.KEGG_REACTION,
        ),
    )
    assert IdentifierNamespaceCv.KEGG_COMPOUND in config.mints
    assert IdentifierNamespaceCv.KEGG_REACTION in config.mints


def test_mints_does_not_claim_a_namespace_a_resource_only_cites():
    """The regression this field exists to catch: KEGG cites PubChem

    identifiers it does not mint. Declaring ``mints`` explicitly is what lets
    the build tell "this resource is the authority" apart from "this resource
    also carries a cross-reference".
    """
    kegg = _config(
        id=ResourceCv.KEGG_METABOLIC,
        mints=(IdentifierNamespaceCv.KEGG_COMPOUND,),
    )
    assert IdentifierNamespaceCv.PUBCHEM_COMPOUND not in kegg.mints
