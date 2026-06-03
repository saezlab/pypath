"""Tests for the resource 3-name model (Milestone M)."""

from __future__ import annotations

import pytest

from pypath.inputs_v2.resource_names import (
    RESOURCE_NAMES,
    build_filter_index,
    resolve_filter,
    resolve_names,
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


def test_every_registry_entry_is_rule_conforming():
    """SC-018/SC-020: the curated audit registry has no rule violations."""
    for names in RESOURCE_NAMES.values():
        assert validate_resource_name(names.slug, names.short, names.full) == [], names


def test_resolve_names_curated():
    names = resolve_names(name='RaMP')
    assert names.slug == 'ramp'
    assert names.short == 'RaMP'
    assert names.full == 'Relational Database for Metabolomic Pathways'
    assert 'RaMP-DB' in names.synonyms


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
    assert resolve_filter('ramp', index) == 'ramp'
    assert resolve_filter('RaMP', index) == 'ramp'
    assert resolve_filter('RaMP-DB', index) == 'ramp'
    assert resolve_filter('rampdb', index) == 'ramp'
    assert resolve_filter('nonexistent-xyz', index) is None
