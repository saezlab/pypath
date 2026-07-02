"""Tests for the resource 3-name model (Milestone M).

Names are sourced from the authoritative resources.json (via resource_registry),
not a hand-curated dict.
"""

from __future__ import annotations

import pytest

from pypath.inputs_v2.resource_names import (
    build_filter_index,
    resolve_filter,
    resolve_names,
    resource_registry,
    slugify,
    validate_resource_name,
)


@pytest.mark.parametrize(
    'name, expected',
    [
        ('RaMP', 'ramp'),
        ('CellPhoneDB', 'cellphonedb'),
        ('LIPID MAPS', 'lipidmaps'),
        ('Guide to Pharmacology', 'guidetopharmacology'),
    ],
)
def test_slugify(name, expected):
    assert slugify(name) == expected


def test_validate_accepts_conforming_names():
    assert validate_resource_name(
        'ramp', 'RaMP', 'Relational Database for Metabolomic Pathways'
    ) == []


@pytest.mark.parametrize(
    'slug, short, full',
    [
        ('RaMP', 'RaMP', 'x'),         # slug not lowercase
        ('ra_mp', 'RaMP', 'x'),        # slug underscore
        ('ra mp', 'RaMP', 'x'),        # slug space
        ('ramp', 'Ra_MP', 'x'),        # short underscore
        ('ramp', 'Ra MP', 'x'),        # short space
        ('ramp', 'RaMP', 'a_b'),       # full underscore
    ],
)
def test_validate_rejects_violations(slug, short, full):
    assert validate_resource_name(slug, short, full) != []


def test_registry_loads_from_resources_json():
    """The registry is sourced from the authoritative resources.json."""
    registry = resource_registry()
    assert len(registry) > 50
    # Self-spelled short = the JSON key; full from full_name.
    assert registry['ramp'].short == 'RaMP'
    assert registry['ramp'].full == 'Relational Database for Metabolomic Pathways'
    assert registry['signor'].short == 'SIGNOR'


def test_resolve_names_from_registry():
    names = resolve_names(name='RaMP')
    assert names.slug == 'ramp'
    assert names.short == 'RaMP'
    assert names.full == 'Relational Database for Metabolomic Pathways'


def test_resolve_names_by_canonical_slug_ignores_inconsistent_name():
    """Passing the canonical slug resolves via resources.json, not the long `name`."""
    names = resolve_names(name='Signaling Network Open Resource', slug='signor')
    assert names.slug == 'signor'
    assert names.short == 'SIGNOR'  # the resources.json key, not the long name


def test_resolve_names_derived_for_unknown():
    names = resolve_names(name='FooBarResource')
    assert names.slug == 'foobarresource'
    assert names.short == 'FooBarResource'
    assert names.full == 'FooBarResource'


def test_resolve_names_explicit_override():
    names = resolve_names(name='whatever', slug='myslug', short='MyShort', full='My Full')
    assert (names.slug, names.short, names.full) == ('myslug', 'MyShort', 'My Full')


def test_resolve_filter_slug_short_synonym():
    """SC-019: slug / short / synonym all resolve to the same canonical slug."""
    index = build_filter_index()
    assert resolve_filter('signor', index) == 'signor'
    assert resolve_filter('SIGNOR', index) == 'signor'
    assert resolve_filter('Signor', index) == 'signor'  # resources.json synonym
    assert resolve_filter('nonexistent-xyz', index) is None


def test_resource_config_names_uses_resourcecv_slug():
    """ResourceConfig.names() resolves via the ResourceCv member, not `name`."""
    from pypath.inputs_v2.base import ResourceConfig
    from pypath.internals.cv_terms import (
        LicenseCV,
        ResourceCv,
        UpdateCategoryCV,
    )

    config = ResourceConfig(
        id=ResourceCv.RAMP if hasattr(ResourceCv, 'RAMP') else ResourceCv.CHEMBL,
        name='some inconsistent long name',
        url='x',
        license=LicenseCV.CC_BY_4_0,
        update_category=UpdateCategoryCV.REGULAR,
        description='d',
    )
    names = config.names()
    # slug comes from the ResourceCv member name, not the inconsistent `name`.
    assert names.slug == slugify(config.id.name)
