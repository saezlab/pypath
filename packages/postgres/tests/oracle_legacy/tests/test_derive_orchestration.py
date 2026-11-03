"""The derive orchestration registers the interaction step, and it is not silent.

Three things this covers:

- the interaction projection is a **registered step of its own**, run once per
  build — ``rebuild_derived_tables`` is called with ``interactions=False``, so
  the 20-minute projection is not paid twice;
- a **failing** step is reported at error level in the structured derive shape
  and re-raised, so a broken projection aborts the derive with a non-zero exit
  instead of passing as a successful build (Principle V, no silently skipped
  phase) — unlike the supplementary network views, which stay caught;
- importing the package **configures logging**, so the
  ``--progress`` output actually reaches a sink rather than being swallowed by
  a root logger left at WARNING.

No database is touched: the projection is stubbed.

    uv run --with pytest pytest tests/test_derive_orchestration.py -v
"""

from __future__ import annotations

import inspect
import logging
import os

import pytest

import omnipath_build  # noqa: F401  (imports `_session`, configuring logging)
from omnipath_build import cli
from omnipath_build.db.derived_tables import InteractionDeriveStats

STATS = InteractionDeriveStats(
    interactions=11,
    parties=22,
    records=33,
    rows_by_class={'signaling': 30, 'ligand_receptor': 3},
    sign_conflict={'signed_rows': 10.0, 'both_flags_rows': 1.0},
    seconds=1.25,
    step_seconds={'interaction_header': 0.5, 'interaction_fact_resource': 0.75},
)


class _Conn:
    """The little the derive step asks of a connection."""

    def __init__(self) -> None:
        self.rolled_back = False

    def rollback(self) -> None:
        self.rolled_back = True


def _lines(caplog) -> list[str]:
    return [record.getMessage() for record in caplog.records]


def test_failing_interaction_step_is_reported_and_propagates(monkeypatch, caplog):
    """A deliberately failing step: error line, rollback, and the error re-raised."""

    def _boom(conn, *, schema, progress):
        raise RuntimeError('interaction projection exploded')

    monkeypatch.setattr(cli, 'rebuild_interaction_tables', _boom)
    conn = _Conn()
    caplog.set_level(logging.INFO, logger='omnipath_build.cli')

    with pytest.raises(RuntimeError, match='exploded'):
        cli._derive_interactions(conn, schema='public')

    failures = [
        record
        for record in caplog.records
        if 'event=interactions_failed' in record.getMessage()
    ]
    assert failures, 'a failing derive step was not reported'
    assert failures[0].levelno == logging.ERROR
    assert 'interaction projection exploded' in failures[0].getMessage()
    assert conn.rolled_back, 'the failed step left its transaction open'
    assert not any('event=interactions_done' in line for line in _lines(caplog))


def test_successful_step_logs_the_cost_figures(monkeypatch, caplog):
    """The done line carries the counts the manifest reads out of --progress."""

    monkeypatch.setattr(
        cli,
        'rebuild_interaction_tables',
        lambda conn, *, schema, progress: STATS,
    )
    caplog.set_level(logging.INFO, logger='omnipath_build.cli')

    stats = cli._derive_interactions(_Conn(), schema='public')

    assert stats is STATS
    done = [line for line in _lines(caplog) if 'event=interactions_done' in line]
    assert done, 'the derive step reported no done line'
    line = done[0]
    assert 'records=33' in line
    assert 'interactions=11' in line
    assert 'parties=22' in line
    assert 'seconds=1.250' in line
    assert 'signaling:30' in line
    assert any('event=interactions_start' in l for l in _lines(caplog))


def test_derive_cost_reaches_the_manifest_with_real_numbers():
    """The stats the step returned become the manifest's per-step cost record."""
    cost = cli._interaction_derive_cost(STATS)
    assert cost['interaction_fact_resource'] == {'seconds': 0.75, 'rows': 33}
    assert cost['interaction_party']['rows'] == 22
    assert cost['interaction_header']['rows'] == 11
    # The build leaves one interaction table, so the collapse is not reported
    # as a step that ran unmeasured — it is not reported at all.
    assert 'interaction_fact_combined' not in cost
    # A build that ran no projection records no cost rather than zeros.
    assert cli._interaction_derive_cost(None) is None


def test_projection_is_registered_once_not_twice():
    """`rebuild_derived_tables` is asked not to repeat the registered step."""
    source = inspect.getsource(cli.main)
    assert 'interactions=False' in source, (
        'rebuild_derived_tables would run the interaction projection a second '
        'time; the derive registers it as a step of its own'
    )
    # And it is registered after the classification that fills the map it reads.
    assert source.index('classify_interaction_class(') < source.index(
        '_derive_interactions('
    )


def test_progress_logging_is_configured():
    """Importing the package configures a session, so INFO output is visible."""
    assert logging.getLogger('omnipath_build').level == logging.INFO, (
        'the package logger is not at INFO; pkg_infra leaves root at WARNING, '
        'so the derive --progress output would be silent'
    )
    assert logging.getLogger().handlers, 'pkg_infra installed no handlers'
    assert logging.getLogger(
        'omnipath_build.db.derived_tables'
    ).isEnabledFor(logging.INFO)


# --- The cost really reaches the manifest --------------------------------
# A scratch schema and the small fixture graph, so this runs the registered
# step for real without paying for the full projection.

DATABASE_URL = os.environ.get('DATABASE_URL')
SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_DERIVE_ORCHESTRATION',
    'derive_orchestration_test',
)


@pytest.fixture(scope='module')
def scratch_conn():
    if not DATABASE_URL:
        pytest.skip('DATABASE_URL not set; the manifest test needs a Postgres')
    import psycopg2

    from omnipath_build.db import schema as build_schema
    from tests.fixtures.interaction_graph import build_interaction_fixture

    connection = psycopg2.connect(DATABASE_URL)
    try:
        build_schema.ensure_schema(connection, schema=SCRATCH, drop_existing=True)
        connection.commit()
        build_interaction_fixture(connection, SCRATCH)
        yield connection
    finally:
        with connection.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {SCRATCH} CASCADE')
        connection.commit()
        connection.close()


def test_registered_step_cost_reaches_the_manifest(scratch_conn):
    """The step runs, and its own numbers land in interactions_derive_cost."""
    from omnipath_build.db.resources import emit_build_manifest

    stats = cli._derive_interactions(scratch_conn, schema=SCRATCH)
    assert stats.records > 0, 'the fixture graph projected no records'

    emit_build_manifest(
        scratch_conn,
        schema=SCRATCH,
        derive_cost=cli._interaction_derive_cost(stats),
    )
    with scratch_conn.cursor() as cur:
        cur.execute(f'SELECT interactions_derive_cost FROM {SCRATCH}.build_manifest')
        cost = cur.fetchone()[0]

    assert cost is not None, 'interactions_derive_cost landed NULL'
    by_step = {entry['step']: entry for entry in cost['steps']}
    assert by_step['interaction_fact_resource']['rows'] == stats.records
    assert by_step['interaction_fact_resource']['seconds'] > 0
    assert by_step['interaction_party']['rows'] == stats.parties
    assert by_step['interaction_header']['rows'] == stats.interactions
