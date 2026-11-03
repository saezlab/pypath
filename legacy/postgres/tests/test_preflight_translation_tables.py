"""The build's identifier-resolution pre-flight reports which resolver relations
are present and usable, warns about any that are missing or empty, and records
the entity types left unresolvable by a missing relation."""

from __future__ import annotations

import logging

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

from omnipath_build import duckdb_load  # noqa: E402
from omnipath_build.duckdb_load import (  # noqa: E402
    CHEMICAL_ENTITY_TYPE,
    REQUIRED_TRANSLATION_TABLES,
    preflight_translation_tables,
)


def _simulate_utils(monkeypatch, *, absent=(), empty=()):
    """Pretend a resolver database is attached, with the named relations dropped
    or emptied and every other required relation present and non-empty."""
    monkeypatch.setattr(duckdb_load, '_live_utils_attached', lambda con: True)
    monkeypatch.setattr(
        duckdb_load,
        '_attached_utils_relation_exists',
        lambda con, name: name not in absent,
    )
    monkeypatch.setattr(
        duckdb_load,
        '_attached_utils_relation_nonempty',
        lambda con, name: name not in absent and name not in empty,
    )


def test_all_present_reports_no_missing(monkeypatch):
    _simulate_utils(monkeypatch)
    con = duckdb.connect(':memory:')

    report = preflight_translation_tables(con)

    assert {r['relation'] for r in report} == set(REQUIRED_TRANSLATION_TABLES)
    assert all(r['status'] == 'present' for r in report)
    assert con.execute(
        'SELECT count(*) FROM missing_translation_entity_type'
    ).fetchone()[0] == 0


def test_dropped_relation_is_reported_missing_and_warns(monkeypatch, caplog):
    _simulate_utils(monkeypatch, absent=('resolver_chemical',))
    con = duckdb.connect(':memory:')

    with caplog.at_level(logging.WARNING):
        report = preflight_translation_tables(con)

    by_relation = {r['relation']: r['status'] for r in report}
    assert by_relation['resolver_chemical'] == 'absent'
    # A chemical can only be resolved by resolver_chemical, so dropping it leaves
    # chemicals unresolvable.
    blocked = [
        row[0]
        for row in con.execute(
            'SELECT entity_type FROM missing_translation_entity_type'
        ).fetchall()
    ]
    assert blocked == [CHEMICAL_ENTITY_TYPE]
    assert any(
        'resolver_chemical' in record.getMessage()
        and record.levelno == logging.WARNING
        for record in caplog.records
    )


def test_empty_relation_is_reported_empty(monkeypatch):
    _simulate_utils(monkeypatch, empty=('resolver_chemical',))
    con = duckdb.connect(':memory:')

    report = preflight_translation_tables(con)

    by_relation = {r['relation']: r['status'] for r in report}
    assert by_relation['resolver_chemical'] == 'empty'


def test_protein_blocked_only_when_every_covering_relation_gone(monkeypatch):
    # A protein can be resolved by resolver_gene, resolver_gene_protein_global, or
    # resolver_protein — losing just one leaves proteins resolvable.
    _simulate_utils(monkeypatch, absent=('resolver_protein',))
    con = duckdb.connect(':memory:')
    preflight_translation_tables(con)
    assert con.execute(
        'SELECT count(*) FROM missing_translation_entity_type'
    ).fetchone()[0] == 0

    _simulate_utils(
        monkeypatch,
        absent=(
            'resolver_gene',
            'resolver_gene_protein_global',
            'resolver_protein',
        ),
    )
    con2 = duckdb.connect(':memory:')
    preflight_translation_tables(con2)
    blocked = {
        row[0]
        for row in con2.execute(
            'SELECT entity_type FROM missing_translation_entity_type'
        ).fetchall()
    }
    assert 'Protein:MI:0326' in blocked


def test_no_resolver_database_is_a_noop(monkeypatch):
    monkeypatch.setattr(duckdb_load, '_live_utils_attached', lambda con: False)
    con = duckdb.connect(':memory:')

    report = preflight_translation_tables(con)

    assert report == []
    # The table still exists (empty), so canonicalisation can join against it.
    assert con.execute(
        'SELECT count(*) FROM missing_translation_entity_type'
    ).fetchone()[0] == 0
