"""T107-T112 (spec 011 WP7/US7) -- the lipid nomenclature grammar binding.

Pure unit tests, no database: :func:`parse_lipid_name` is a stateless
function over one name. Skipped whole when ``pygoslin`` is not installed
(the ``lipid`` extra) -- T112 covers that absence explicitly via a monkeypatch
rather than an environment without the package.
"""

from __future__ import annotations

import pytest

pygoslin = pytest.importorskip('pygoslin')

from omnipath_utils.mapping._lipid import (  # noqa: E402
    lipid_nomenclature_available,
    parse_lipid_name,
)


def test_t107_normalizes_format_variants_and_records_parse_fields():
    """Two format variants of the same lipid normalize to one name, and the
    parse fields (level, chains_possible, chains_listed) are populated."""
    a = parse_lipid_name('PC(16:0/18:1)')
    b = parse_lipid_name('PC 16:0/18:1')
    assert a is not None and b is not None
    assert a['lipid_name'] == b['lipid_name']
    assert a['lipid_level'] == 'sn_position'
    assert a['chains_possible'] == 2
    assert a['chains_listed'] == 2


@pytest.mark.parametrize(
    ('name', 'chains_possible'),
    [
        ('LPC 16:0', 1),  # one-chain class
        ('PC 16:0/18:1', 2),
        ('TG 18:0/18:1/18:1', 3),
        ('CL 18:0/18:1/18:2/18:3', 4),
    ],
)
def test_t108_a_fully_specified_name_keys_to_itself(name, chains_possible):
    """A fully specified name at 1/2/3/4 chains: chains_listed ==
    chains_possible, and re-parsing its own canonical rendering is a fixed
    point (keys to itself)."""
    parsed = parse_lipid_name(name)
    assert parsed is not None
    assert parsed['chains_possible'] == chains_possible
    assert parsed['chains_listed'] == chains_possible
    again = parse_lipid_name(parsed['lipid_name'])
    assert again['lipid_name'] == parsed['lipid_name']
    assert again['lipid_level'] == parsed['lipid_level']


@pytest.mark.parametrize('name', ['TG 58:12_20:0', 'CL 54:6_18:2'])
def test_t109_a_partially_specified_name_keeps_its_own_identity(name):
    """A name that sums some chains and names the rest does not collapse to
    the species total -- it is its own level, chains_listed strictly between
    0 and chains_possible."""
    parsed = parse_lipid_name(name)
    assert parsed is not None
    assert parsed['lipid_level'] == 'partially_specified'
    assert 0 < parsed['chains_listed'] < parsed['chains_possible']
    # Must not have collapsed to the species-level name (a bare total, one
    # chain token) -- the species rendering of the same total.
    species_name = f"{name.split()[0]} {parsed['total_carbon']}:{parsed['total_db']}"
    species = parse_lipid_name(species_name)
    assert species is not None
    assert species['lipid_level'] == 'species'
    assert parsed['lipid_name'] != species['lipid_name']


def test_t110_both_source_conventions_agree_on_one_canonical_name():
    """The underscore (parts add) and dash-FA (named chain subtracted from
    the total) conventions for the same lipid produce one canonical name,
    and the listed parts add to the species total."""
    underscore = parse_lipid_name('TG 58:12_20:0')
    dash_fa = parse_lipid_name('TG 78:12-FA20:0')
    assert underscore is not None and dash_fa is not None
    assert underscore['lipid_name'] == dash_fa['lipid_name']
    assert underscore['total_carbon'] == 78
    assert underscore['total_db'] == 12


@pytest.mark.parametrize(
    ('name', 'expected_chains'),
    [
        ('MG 18:1', 1),  # one-chain class: naive text-token count would be 3
        ('LPC(16:0/0:0)', 1),  # one-chain class, explicit empty slot: naive would be 2
        ('CE 18:1', 1),  # sterol ester: naive text-token count would be 2
    ],
)
def test_t111_chain_counts_come_from_the_parse_not_rendered_tokens(
    name, expected_chains,
):
    """The classes whose rendering inflates a naive C:D-token count over the
    real chain count -- chains_possible/chains_listed must still be correct,
    because they come from the parse object's fa_list, not string counting."""
    import re

    parsed = parse_lipid_name(name)
    assert parsed is not None
    assert parsed['chains_possible'] == expected_chains
    assert parsed['chains_listed'] == expected_chains
    naive_token_count = len(re.findall(r'\d+:\d+', parsed['lipid_name']))
    assert naive_token_count != expected_chains, (
        f'{name} was supposed to demonstrate the token-count bug, but its '
        f'rendering {parsed["lipid_name"]!r} has {naive_token_count} tokens '
        f'== the real chain count -- pick a different example'
    )


def test_t112_absent_extra_degrades_cleanly(monkeypatch):
    """With the grammar unavailable, parsing returns None (never raises),
    matching _chemistry.py's degradation convention -- callers building
    lipid_name_node/_edge from this get nothing to write, not a crash."""
    lipid_nomenclature_available.cache_clear()
    monkeypatch.setattr(
        'omnipath_utils.mapping._lipid.lipid_nomenclature_available',
        lambda: False,
    )
    assert parse_lipid_name('PC 16:0/18:1') is None
    lipid_nomenclature_available.cache_clear()


def test_non_lipid_name_returns_none():
    assert parse_lipid_name('glucose') is None
    assert parse_lipid_name('') is None
    assert parse_lipid_name(None) is None


def test_capability_reports_available_in_this_environment():
    # pygoslin is installed in this test environment (importorskip above
    # already proved that); the capability function must agree.
    assert lipid_nomenclature_available() is True
