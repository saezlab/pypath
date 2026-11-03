"""The catalogue round trip across a deferred load.

The projection loads ``interaction``, ``interaction_party`` and
``interaction_fact_resource`` with their foreign keys and secondary indexes
**dropped**, and puts them back before the step ends. A measurement priced it
at 709.7 s against 1,814.7 s over the four tables that existed then, with the
largest step falling from 674.0 s to 272.2 s — sixty per cent of that step was
constraint and index maintenance rather than the work of building a header.

**The saving is only a saving if the far side of the step is the near side.**
So this file asserts the round trip rather than the seconds: the same
constraints, the same indexes, every definition identical string for string,
and the same rows, class counts and sign-conflict summary as a load that
deferred nothing.

**A restored-but-``NOT VALID`` key is the failure this file exists to catch.**
Postgres will happily take a foreign key back as ``NOT VALID``: it looks like a
constraint in ``\\d``, it is enforced on every row written **after** it, and it
says nothing whatever about the fourteen million rows the load just wrote. That
is not the constraint the schema declares, and a deferral that returns one has
traded a guarantee for time. Every foreign key is therefore asserted
``convalidated`` by the catalogue, not by its presence.

The counts come from the schema rather than from the note that priced the
deferral. That measurement was taken on 2026-08-20 over **four** tables and
recorded 15 constraints and 20 indexes, 13 of them foreign keys and 18
secondary. Removing ``interaction_fact_combined`` — four foreign keys and nine
indexes — left the amended schema holding 11 constraints and 11 indexes, and
the record's index on its header id brought the second number to **12**, of
which **9 foreign keys and 10 secondary indexes** are what the deferral drops
and restores. The two primary keys stay through the load, because the header
insert deduplicates with ``ON CONFLICT (interaction_id) DO NOTHING`` and needs
its unique index while the insert runs.

The counts are asserted rather than derived on purpose: an index added to the
schema and not to the deferral would otherwise be a silent gap, so a schema
change is meant to land here as a failing number that someone recounts.

The fixture graph is a dozen relations, so nothing here measures anything. It
is the round trip that is being asserted, and the round trip is a property of
the mechanism rather than of the row count.

Run::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_interactions_deferral.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import inspect
import os

import psycopg2
import pytest

from tests.fixtures.interaction_graph import build_interaction_fixture

DATABASE_URL = os.environ.get('DATABASE_URL')
SCRATCH_DEFERRED = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_DEFERRAL',
    'interactions_deferral_test',
)
SCRATCH_PLAIN = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_DEFERRAL_BASELINE',
    'interactions_deferral_baseline_test',
)

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the deferral test needs a Postgres',
)

#: The tables the deferral covers. The vocabularies and the license table are
#: seeded rather than loaded, so nothing is deferred on them.
PROJECTION_TABLES = (
    'interaction',
    'interaction_party',
    'interaction_fact_resource',
)

#: What the amended schema holds on those three tables, counted from the
#: catalogue and asserted here so a schema change that adds an object without
#: adding it to the deferral is a failure rather than a silent gap.
#:
#: ``TABLE_CONSTRAINTS`` counts foreign keys and primary keys — the constraint
#: kinds the deferral measurement counted. Postgres 17 and later also list
#: every ``NOT NULL`` in ``pg_constraint``; those are column properties that no
#: load can drop, and counting them would make the number depend on the server
#: version rather than on the schema.
TABLE_CONSTRAINTS = 11
TABLE_INDEXES = 12

#: What the load runs without, and gets back validated.
DEFERRED_FOREIGN_KEYS = 9
DEFERRED_INDEXES = 10

#: The parameter that turns the deferral on. Named once, because the fixture and
#: the guard below must not spell it differently.
DEFER_PARAMETER = 'defer_constraints'


def _require_the_deferral_exists() -> None:
    """Fail with what is missing, rather than erroring out of a TypeError.

    Constitution III puts this test before the implementation, so it has to
    fail for a reason a reader can act on. ``rebuild_interaction_tables``
    without the parameter is the build before the deferral, and saying so is
    more use than an unexpected keyword argument traceback.
    """
    from omnipath_build.db.derived_tables import rebuild_interaction_tables

    parameters = inspect.signature(rebuild_interaction_tables).parameters
    if DEFER_PARAMETER not in parameters:
        pytest.fail(
            f'rebuild_interaction_tables takes no {DEFER_PARAMETER!r}: the '
            f'load still runs through its foreign keys and secondary '
            f'indexes, so there is no deferral for the catalogue to '
            f'round-trip'
        )


def _project(schema: str, *, deferred: bool):
    """Build the fixture graph in ``schema`` and project it, once."""
    from omnipath_build.db import schema as build_schema
    from omnipath_build.db.derived_tables import rebuild_interaction_tables

    connection = psycopg2.connect(DATABASE_URL)
    build_schema.ensure_schema(connection, schema=schema, drop_existing=True)
    connection.commit()
    build_interaction_fixture(connection, schema)
    stats = rebuild_interaction_tables(
        connection,
        schema=schema,
        **{DEFER_PARAMETER: deferred},
    )
    return connection, stats


@pytest.fixture(scope='module')
def arms():
    """Both arms: the same fixture graph loaded with and without the deferral.

    Two schemas rather than two runs in one, because the comparison is between
    catalogues and a second run in one schema would compare a catalogue with
    itself.
    """
    _require_the_deferral_exists()
    deferred = _project(SCRATCH_DEFERRED, deferred=True)
    plain = _project(SCRATCH_PLAIN, deferred=False)
    try:
        yield deferred, plain
    finally:
        for connection, schema in (
            (deferred[0], SCRATCH_DEFERRED),
            (plain[0], SCRATCH_PLAIN),
        ):
            with connection.cursor() as cur:
                cur.execute(f'DROP SCHEMA IF EXISTS {schema} CASCADE')
            connection.commit()
            connection.close()


def _constraints(conn, schema: str) -> dict[str, tuple[str, str, bool]]:
    """Every foreign and primary key on the three tables, by name.

    The definition comes from ``pg_get_constraintdef``, which renders the
    schema it lives in, so the schema name is replaced by a placeholder — the
    two arms are in two schemas and the comparison is about everything else.
    """
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT con.conname,
                   con.contype,
                   pg_get_constraintdef(con.oid),
                   con.convalidated
            FROM pg_constraint con
            JOIN pg_class cls ON cls.oid = con.conrelid
            JOIN pg_namespace ns ON ns.oid = cls.relnamespace
            WHERE ns.nspname = %s
              AND cls.relname = ANY(%s)
              AND con.contype IN ('f', 'p')
            """,
            (schema, list(PROJECTION_TABLES)),
        )
        return {
            name: (contype, definition.replace(schema, '{schema}'), validated)
            for name, contype, definition, validated in cur.fetchall()
        }


def _indexes(conn, schema: str) -> dict[str, str]:
    """Every index on the three tables, by name, with its definition."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT indexname, indexdef
            FROM pg_indexes
            WHERE schemaname = %s AND tablename = ANY(%s)
            """,
            (schema, list(PROJECTION_TABLES)),
        )
        return {
            name: definition.replace(schema, '{schema}')
            for name, definition in cur.fetchall()
        }


def _row_counts(conn, schema: str) -> dict[str, int]:
    counts = {}
    with conn.cursor() as cur:
        for table in PROJECTION_TABLES:
            cur.execute(f'SELECT count(*) FROM {schema}.{table}')
            counts[table] = cur.fetchone()[0]
    return counts


# ---------------------------------------------------------------------------
# The catalogue on the far side is the catalogue on the near side
# ---------------------------------------------------------------------------


def test_the_deferred_load_leaves_the_same_constraints(arms):
    """Name for name, definition for definition — and the count is asserted too.

    The count is what catches an object the deferral drops and forgets to
    restore *and* an object the schema gained without the deferral learning
    about it. The definitions are what catch a restore that rebuilt something
    subtly different — a key without its ``ON DELETE CASCADE``, say, which
    still looks like the same key by name.
    """
    (deferred_conn, _), (plain_conn, _) = arms
    deferred = _constraints(deferred_conn, SCRATCH_DEFERRED)
    plain = _constraints(plain_conn, SCRATCH_PLAIN)

    assert len(plain) == TABLE_CONSTRAINTS, (
        f'the undeferred schema holds {len(plain)} foreign and primary keys on '
        f'{list(PROJECTION_TABLES)}, not {TABLE_CONSTRAINTS}: recount the '
        f'schema and amend the constant, do not weaken the assertion'
    )
    assert deferred == plain


def test_the_deferred_load_leaves_the_same_indexes(arms):
    """Every index comes back, with the definition it went away with."""
    (deferred_conn, _), (plain_conn, _) = arms
    deferred = _indexes(deferred_conn, SCRATCH_DEFERRED)
    plain = _indexes(plain_conn, SCRATCH_PLAIN)

    assert len(plain) == TABLE_INDEXES, (
        f'the undeferred schema holds {len(plain)} indexes on '
        f'{list(PROJECTION_TABLES)}, not {TABLE_INDEXES}'
    )
    assert deferred == plain


def test_every_restored_foreign_key_is_validated(arms):
    """The failure this file exists to catch.

    ``NOT VALID`` is the cheap way to put a key back, and it is the wrong one:
    the constraint then describes rows written after the load and says nothing
    about the rows the load wrote. Revalidating all 229.9 million row checks
    set-based measured 44.6 s against 726.3 s of per-row triggers through the
    load, so the validation is the saving rather than an extra on top of it.
    """
    (deferred_conn, _), _plain = arms
    unvalidated = sorted(
        name
        for name, (contype, _definition, validated) in _constraints(
            deferred_conn, SCRATCH_DEFERRED
        ).items()
        if contype == 'f' and not validated
    )
    assert unvalidated == [], (
        f'foreign keys restored NOT VALID: {unvalidated}. They look like '
        f'constraints and enforce nothing about the rows already loaded'
    )


def test_the_deferral_drops_the_keys_and_indexes_it_says_it_does(arms):
    """The step reports what it deferred, and the numbers are the schema's.

    Reported rather than inferred, because the manifest records it and
    because a deferral that quietly stopped dropping anything would still pass
    every catalogue assertion above — the round trip of a load that deferred
    nothing is trivially unchanged.
    """
    (_conn, stats), _plain = arms
    deferral = getattr(stats, 'deferral', None)
    assert deferral, 'the step reported no deferral record'
    assert deferral['deferred'] is True
    assert deferral['constraints_deferred'] == DEFERRED_FOREIGN_KEYS
    assert deferral['indexes_deferred'] == DEFERRED_INDEXES
    assert deferral['catalogue_unchanged'] is True
    assert deferral['revalidate_seconds'] is not None, (
        'the step recorded no revalidation cost, so nothing says the keys came '
        'back validated rather than NOT VALID'
    )


def test_the_undeferred_arm_defers_nothing(arms):
    """The baseline arm is a baseline: it drops no key and no index.

    Without this the two arms could be the same arm, and every comparison
    above would be a comparison of a run with itself.
    """
    _deferred, (_conn, stats) = arms
    deferral = getattr(stats, 'deferral', None)
    assert deferral, 'the step reported no deferral record'
    assert deferral['deferred'] is False
    assert not deferral['constraints_deferred']
    assert not deferral['indexes_deferred']


# ---------------------------------------------------------------------------
# And the rows on the far side are the rows on the near side
# ---------------------------------------------------------------------------


def test_the_deferred_load_writes_the_same_rows(arms):
    """Digit for digit, on all three tables."""
    (deferred_conn, _), (plain_conn, _) = arms
    assert _row_counts(deferred_conn, SCRATCH_DEFERRED) == _row_counts(
        plain_conn, SCRATCH_PLAIN
    )


def test_the_deferred_load_reports_the_same_class_counts(arms):
    """All eight interaction classes, every one of them, at the same count."""
    (_deferred_conn, deferred_stats), (_plain_conn, plain_stats) = arms
    assert deferred_stats.rows_by_class == plain_stats.rows_by_class
    assert len(plain_stats.rows_by_class) == 8, (
        f'the class vocabulary holds {len(plain_stats.rows_by_class)} classes, '
        f'not the eight the derivation produces, so this comparison covers '
        f'less than it claims'
    )


def test_the_deferred_load_reports_the_same_sign_conflict(arms):
    """The sign-conflict summary measures the rows, so it cannot move.

    It is measured by folding the record, and the fold reads the collapse index
    the deferral drops and restores. A summary that differs between the arms
    would mean the restore changed what the record says, not merely when it
    said it.
    """
    (_deferred_conn, deferred_stats), (_plain_conn, plain_stats) = arms
    assert deferred_stats.sign_conflict == plain_stats.sign_conflict
    assert plain_stats.sign_conflict['signed_rows'] > 0, (
        'the fixture graph produced no signed rows, so the comparison holds '
        'nothing'
    )
