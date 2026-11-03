"""One canonical identity absorbing an unreasonable number of source records
fails the build, by name (spec 011 T031, spec.md US1 acceptance scenario 5).

A hub -- one identity that swallows tens of thousands of unrelated source
records -- repeats the placeholder-merge failure from the lipid hierarchy
work, now in the chemical key space. The check is generic over the
identifier, not specific to any one namespace, since the next hub can
originate from any of them.
"""

from __future__ import annotations

import pytest

from omnipath_build.resolver.preflight import check_absorption_threshold


def test_passes_when_every_identity_is_under_threshold():
    counts = {'CHEBI:1': 5, 'CHEBI:2': 100}
    check_absorption_threshold(counts, threshold=1000)  # must not raise


def test_fails_and_names_the_identifier_over_threshold():
    counts = {'CHEBI:1': 5, 'CHEBI:2': 200000}
    with pytest.raises(RuntimeError, match='CHEBI:2'):
        check_absorption_threshold(counts, threshold=100000)


def test_names_every_offending_identity_at_once():
    counts = {'CHEBI:1': 200000, 'HMDB0000001': 300000}
    with pytest.raises(RuntimeError, match='CHEBI:1') as exc:
        check_absorption_threshold(counts, threshold=100000)
    assert 'HMDB0000001' in str(exc.value)


def test_exactly_at_threshold_is_not_a_failure():
    counts = {'CHEBI:1': 100000}
    check_absorption_threshold(counts, threshold=100000)  # must not raise
