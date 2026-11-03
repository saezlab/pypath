"""A namespace whose evidence and resolver keys never overlap fails the
build loudly, by name (spec 011 T029).

A silent zero-overlap namespace is exactly the shape of the defect this
cycle exists to fix (spec.md US1 acceptance scenario 6). The ChEBI key
never matched between the translation database and the build's own
evidence, and nothing said so. WP1's own normalization landing does not
guard against a *future* namespace drifting the same way, so the check
runs every build, not just this one.
"""

from __future__ import annotations

import pytest

from omnipath_build.resolver.preflight import check_key_space_overlap


def test_passes_when_every_namespace_overlaps():
    evidence_keys = {'chebi': {'CHEBI:1', 'CHEBI:2'}, 'hmdb': {'HMDB0000001'}}
    resolver_keys = {'chebi': {'CHEBI:2', 'CHEBI:3'}, 'hmdb': {'HMDB0000001'}}
    check_key_space_overlap(evidence_keys, resolver_keys)  # must not raise


def test_fails_and_names_the_namespace_with_no_overlap():
    evidence_keys = {'chebi': {'CHEBI:1', 'CHEBI:2'}, 'hmdb': {'HMDB0000001'}}
    resolver_keys = {'chebi': set(), 'hmdb': {'HMDB0000001'}}
    with pytest.raises(RuntimeError, match='chebi'):
        check_key_space_overlap(evidence_keys, resolver_keys)


def test_fails_and_names_a_namespace_missing_from_the_resolver_entirely():
    evidence_keys = {'kegg': {'C00002'}}
    resolver_keys = {}
    with pytest.raises(RuntimeError, match='kegg'):
        check_key_space_overlap(evidence_keys, resolver_keys)


def test_names_every_offending_namespace_at_once():
    evidence_keys = {'chebi': {'CHEBI:1'}, 'kegg': {'C00002'}}
    resolver_keys = {'chebi': set(), 'kegg': set()}
    with pytest.raises(RuntimeError, match='chebi') as exc:
        check_key_space_overlap(evidence_keys, resolver_keys)
    assert 'kegg' in str(exc.value)


def test_a_namespace_with_no_evidence_at_all_is_not_a_failure():
    """Nothing to resolve for a namespace this build never ingested --
    that is not the same defect as a namespace present on both sides that
    fails to match."""
    evidence_keys = {'chebi': set()}
    resolver_keys = {'hmdb': {'HMDB0000001'}}
    check_key_space_overlap(evidence_keys, resolver_keys)  # must not raise
