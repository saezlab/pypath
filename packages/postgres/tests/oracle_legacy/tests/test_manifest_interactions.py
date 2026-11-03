"""The manifest records interaction derive cost and the preset inventory.

The build leaves one interaction table, so the cost block names it alone:
``interaction_fact_resource``, the record. Every measured query scope still
reports what materialising it cost — the scopes this build declined to
materialise included — and ``interactions_deferral_cost`` reports what the
deferred load saved and what restoring the catalogue cost.

Both are volatile — per-step seconds and row counts, and the registered preset
list — so they are written next to the identity hash, never inside it. These
tests assert the manifest carries them and that two builds over identical inputs
still agree on ``build_id`` when their timings differ (the cycle-007 lesson).

Run against a build database, e.g. on dev4::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_manifest_interactions.py -v

The manifests are written into a throwaway schema, so the build's own
``public.build_manifest`` row is never rewritten. Skipped without DATABASE_URL.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os

import pytest

from omnipath_build.db.resources import (
    INTERACTION_RECORD_STEP,
    emit_build_manifest,
)
from omnipath_build.network_views import (
    NetworkDefinition,
    ensure_network_registry,
    register_network,
)

DATABASE_URL = os.environ.get('DATABASE_URL')
TEST_SCHEMA = 'manifest_interactions_test'

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; manifest test needs a database',
)

PRESET = NetworkDefinition(
    name='_manifest_preset',
    kind='ligand_receptor',
    included_sources=('connectomedb2025',),
    interaction_class_scope=('ligand_receptor',),
    default_attributes=('endpoints',),
    mandatory_attributes=('label', 'references', 'evidence'),
)

# Two runs of the same build: identical inputs, different timings and row counts
# observed along the way.
COST_FIRST = {
    INTERACTION_RECORD_STEP: {'seconds': 12.5, 'rows': 1000},
    'interaction_assay': {'seconds': 3.25, 'rows': 400},
    'interaction_party': {'seconds': 5.0, 'rows': 2000},
    'reaction_projection': {'seconds': 8.75, 'rows': 700},
    'intercell': {'seconds': 1.5, 'rows': 90},
}
COST_SECOND = {
    step: {'seconds': cost['seconds'] * 2 + 0.5, 'rows': cost['rows']}
    for step, cost in COST_FIRST.items()
}


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture(scope='module')
def build_schema(conn):
    """A throwaway schema holding the inputs a manifest is derived from."""
    with conn.cursor() as cur:
        cur.execute(f'DROP SCHEMA IF EXISTS {TEST_SCHEMA} CASCADE')
        cur.execute(f'CREATE SCHEMA {TEST_SCHEMA}')
        cur.execute(
            f"""
            CREATE TABLE {TEST_SCHEMA}.resources (
              resource_id text PRIMARY KEY,
              entity_count bigint NOT NULL DEFAULT 0,
              interaction_count bigint NOT NULL DEFAULT 0,
              association_count bigint NOT NULL DEFAULT 0,
              identifier_count bigint NOT NULL DEFAULT 0,
              ontology_term_count bigint NOT NULL DEFAULT 0,
              input_module_commit text,
              input_module_dirty boolean NOT NULL DEFAULT false
            )
            """
        )
        cur.execute(
            f"""
            INSERT INTO {TEST_SCHEMA}.resources
              (resource_id, entity_count, interaction_count, input_module_commit)
            VALUES ('signor', 10, 20, 'abc123'), ('connectomedb2025', 5, 7, 'def456')
            """
        )
    conn.commit()
    ensure_network_registry(conn, registry_schema=TEST_SCHEMA)
    register_network(conn, PRESET, registry_schema=TEST_SCHEMA)
    yield TEST_SCHEMA
    with conn.cursor() as cur:
        cur.execute(f'DROP SCHEMA IF EXISTS {TEST_SCHEMA} CASCADE')
    conn.commit()


def _manifest(conn, schema):
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT build_id, package_commits, resources,
                   interactions_derive_cost, network_presets
            FROM {schema}.build_manifest
            """
        )
        return cur.fetchone()


def test_manifest_records_derive_cost_and_preset_inventory(conn, build_schema):
    emit_build_manifest(conn, schema=build_schema, derive_cost=COST_FIRST)
    _build_id, _commits, _resources, cost, presets = _manifest(conn, build_schema)

    assert isinstance(cost, dict), 'no interactions_derive_cost in the manifest'
    by_step = {entry['step']: entry for entry in cost['steps']}
    assert set(by_step) == set(COST_FIRST)
    assert by_step[INTERACTION_RECORD_STEP]['seconds'] == pytest.approx(12.5)
    assert by_step[INTERACTION_RECORD_STEP]['rows'] == 1000

    assert isinstance(presets, list) and presets, 'no network_presets inventory'
    registered = {entry['name']: entry for entry in presets}
    assert PRESET.name in registered
    entry = registered[PRESET.name]
    assert entry['kind'] == PRESET.kind
    assert entry['included_sources'] == list(PRESET.included_sources)
    assert entry['interaction_class_scope'] == list(PRESET.interaction_class_scope)
    assert entry['mandatory_attributes'] == list(PRESET.mandatory_attributes)


def test_build_id_is_stable_across_differing_timings(conn, build_schema):
    """Identical inputs, different timings — the same build identity."""
    first = emit_build_manifest(conn, schema=build_schema, derive_cost=COST_FIRST)
    _id, _c, _r, first_cost, _p = _manifest(conn, build_schema)
    second = emit_build_manifest(conn, schema=build_schema, derive_cost=COST_SECOND)
    _id2, _c2, _r2, second_cost, _p2 = _manifest(conn, build_schema)

    assert first_cost != second_cost, 'the two runs recorded the same timings'
    assert first.build_id == second.build_id


def test_derive_cost_and_presets_stay_out_of_the_hash_payload(conn, build_schema):
    """The identity hash covers the package commits and resources — nothing else."""
    emit_build_manifest(conn, schema=build_schema, derive_cost=COST_SECOND)
    build_id, package_commits, resources, cost, presets = _manifest(
        conn, build_schema
    )
    payload = {'package_commits': package_commits, 'resources': resources}
    expected = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()[:12]
    assert build_id == expected
    # And the volatile fields really were present while that hash was computed.
    assert cost and presets


def test_manifest_omits_cost_when_no_derive_step_ran(conn, build_schema):
    """A build with no interaction derive steps records no cost, and still hashes."""
    stats = emit_build_manifest(conn, schema=build_schema)
    build_id, _commits, _resources, cost, presets = _manifest(conn, build_schema)
    assert cost is None
    assert build_id == stats.build_id
    assert presets is not None


# --- one interaction table, and the deferral over its load ----------------
#
# With the materialisation gone, `interaction_tables` names the record and
# nothing else, and no scope is the all-resources one any more. What the
# collapse used to buy the reader — an attributable share of the projection's
# cost — is now bought by `interactions_deferral_cost`, which says what the
# deferred load saved and what putting the catalogue back cost.

COST_INTERACTION_TABLES = {
    'interaction_party': {'seconds': 679.1, 'rows': 28359137},
    INTERACTION_RECORD_STEP: {'seconds': 138.6, 'rows': 14686404},
    'intercell': {'seconds': 1.5, 'rows': 90},
}

# A build that wrote no interaction table at all: nothing to name apart.
COST_NO_INTERACTION_TABLE = {
    'interaction_party': {'seconds': 5.0, 'rows': 2000},
    'intercell': {'seconds': 1.5, 'rows': 90},
}

SCOPE_COST = {
    'kinase_extra': {'seconds': 4.0, 'rows': 120},
    'academic_license': {
        'table': 'interaction_fact_academic',
        'seconds': 51.5,
        'rows': 900,
        'sources': ('signor',),
    },
    'signor_only': {
        'table': 'interaction_fact_signor',
        'seconds': 7.5,
        'rows': 20,
        'sources': ('signor',),
    },
}

# What the deferral measurement found, in the shape the step hands over.
# Measured on dev4 2026-08-24, after the materialisation came out: the deferral
# drops nine foreign keys and nine secondary indexes, not the 13 and 18 the
# first measurement counted, because `interaction_fact_combined` carried four of
# the keys and nine of the indexes and no longer exists to carry them.
DEFERRAL_COST = {
    'deferred': True,
    'seconds_saved': 737.0,
    'load_seconds': 399.5,
    'drop_seconds': 1.003,
    'restore_seconds': 25.537,
    'revalidate_seconds': 31.910,
    'constraints_deferred': 9,
    'indexes_deferred': 9,
    'catalogue_unchanged': True,
}


def _deferral(conn, schema):
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT interactions_deferral_cost FROM {schema}.build_manifest'
        )
        return cur.fetchone()[0]


def test_interaction_tables_name_the_record_alone(conn, build_schema):
    """One table is built, so one table is named — no collapse, no total."""
    emit_build_manifest(
        conn, schema=build_schema, derive_cost=COST_INTERACTION_TABLES
    )
    _id, _c, _r, cost, _p = _manifest(conn, build_schema)

    tables = cost['interaction_tables']
    assert tables == {
        'record': {
            'table': INTERACTION_RECORD_STEP,
            'seconds': pytest.approx(138.6),
            'rows': 14686404,
        }
    }
    # The removed halves are gone rather than emitted empty: a `collapse: null`
    # would read as a table that ran and was not measured, which is a claim
    # about a table this build does not have.
    assert 'collapse' not in tables
    assert 'total_seconds' not in tables

    # The naming is an extra reading of the same step, not a replacement: the
    # record stays in `steps` too, so the whole-derive total is still readable.
    by_step = {entry['step']: entry for entry in cost['steps']}
    assert INTERACTION_RECORD_STEP in by_step


def test_unmeasured_record_is_null_and_not_zero(conn, build_schema):
    """Absent and free stay distinguishable: nothing reported renders ``null``.

    A record rendered as zero reads as "the table was free", which is a claim
    about the build; one rendered as ``null`` reads as "nobody measured it",
    which is a claim about the measurement.
    """
    emit_build_manifest(
        conn,
        schema=build_schema,
        # The record was written but never timed or counted.
        derive_cost={INTERACTION_RECORD_STEP: {'seconds': None, 'rows': None}},
    )
    _id, _c, _r, cost, _p = _manifest(conn, build_schema)

    record = cost['interaction_tables']['record']
    assert record['table'] == INTERACTION_RECORD_STEP
    assert record['seconds'] is None, 'an unmeasured table was rendered as a number'
    assert record['rows'] is None, 'an uncounted table was rendered as a number'


def test_interaction_tables_absent_when_no_table_reported(conn, build_schema):
    """The block is omitted, not emitted empty, when the record did not run."""
    emit_build_manifest(
        conn, schema=build_schema, derive_cost=COST_NO_INTERACTION_TABLE
    )
    _id, _c, _r, cost, _p = _manifest(conn, build_schema)

    assert 'interaction_tables' not in cost
    assert set(cost) == {'steps'}


def test_scope_cost_orders_materialised_scopes_first(conn, build_schema):
    """Materialised scopes first, then the measured ones, each by name."""
    emit_build_manifest(
        conn,
        schema=build_schema,
        derive_cost=COST_INTERACTION_TABLES,
        scope_cost=SCOPE_COST,
    )
    _id, _c, _r, cost, _p = _manifest(conn, build_schema)

    assert [entry['scope'] for entry in cost['scopes']] == [
        'academic_license',
        'signor_only',
        'kinase_extra',
    ]


def test_scope_measured_without_being_materialised_round_trips(conn, build_schema):
    """The scope this build declined to store still reports its cost.

    Declining to materialise a scope has to stay an available answer to a cost
    overrun, and that argument is made from the number for the scope that was
    declined — so ``table: null`` with ``materialised: false`` must survive the
    round trip through the manifest rather than being dropped as "not a table".
    """
    emit_build_manifest(
        conn,
        schema=build_schema,
        derive_cost=COST_INTERACTION_TABLES,
        scope_cost=SCOPE_COST,
    )
    _id, _c, _r, cost, _p = _manifest(conn, build_schema)
    by_scope = {entry['scope']: entry for entry in cost['scopes']}

    declined = by_scope['kinase_extra']
    assert declined['table'] is None
    assert declined['materialised'] is False
    assert declined['seconds'] == pytest.approx(4.0)
    assert declined['rows'] == 120
    assert declined['sources'] is None

    # A scope that named a table is materialised without having to say so.
    stored = by_scope['academic_license']
    assert stored['table'] == 'interaction_fact_academic'
    assert stored['materialised'] is True
    assert stored['sources'] == ['signor']


def test_scope_record_without_a_name_is_warned_about_and_dropped(
    conn,
    build_schema,
    caplog,
):
    """A malformed scope record costs a warning, never the build."""
    scope_cost = [
        {'seconds': 4.0, 'rows': 120},  # no `scope`: unattributable
        {
            'scope': 'signor_only',
            'table': 'interaction_fact_signor',
            'seconds': 7.5,
            'rows': 20,
        },
    ]
    with caplog.at_level(logging.WARNING, logger='omnipath_build.db.resources'):
        stats = emit_build_manifest(
            conn,
            schema=build_schema,
            derive_cost=COST_INTERACTION_TABLES,
            scope_cost=scope_cost,
        )
    _id, _c, _r, cost, _p = _manifest(conn, build_schema)

    assert [entry['scope'] for entry in cost['scopes']] == ['signor_only']
    assert stats.build_id
    assert any(
        'scope-cost record without a scope' in record.getMessage()
        for record in caplog.records
    )


def test_deferral_cost_is_recorded(conn, build_schema):
    """What the deferral saved and what the restore cost, per build."""
    emit_build_manifest(
        conn,
        schema=build_schema,
        derive_cost=COST_INTERACTION_TABLES,
        deferral_cost=DEFERRAL_COST,
    )
    deferral = _deferral(conn, build_schema)

    assert deferral == {
        'deferred': True,
        'seconds_saved': pytest.approx(737.0),
        'load_seconds': pytest.approx(399.5),
        'drop_seconds': pytest.approx(1.003),
        'restore_seconds': pytest.approx(25.537),
        'revalidate_seconds': pytest.approx(31.910),
        'constraints_deferred': 9,
        'indexes_deferred': 9,
        'catalogue_unchanged': True,
    }

    # It is its own column, not a nested reading of the derive cost: the
    # deferral is a property of how the load ran, not of what a step cost.
    _id, _c, _r, cost, _p = _manifest(conn, build_schema)
    assert 'deferral' not in cost


def test_deferral_cost_is_null_when_the_deferral_did_not_run(conn, build_schema):
    """A build without the deferral reports ``null``, never zero.

    Zero seconds saved is a measurement — a deferral that bought nothing. No
    deferral at all is the absence of one, and the two must not read alike,
    because the first is a regression to chase and the second is not.
    """
    emit_build_manifest(
        conn, schema=build_schema, derive_cost=COST_INTERACTION_TABLES
    )
    assert _deferral(conn, build_schema) is None

    # And a deferral that measured nothing at all says nothing, rather than
    # claiming a run with every number missing.
    emit_build_manifest(
        conn,
        schema=build_schema,
        derive_cost=COST_INTERACTION_TABLES,
        deferral_cost={},
    )
    assert _deferral(conn, build_schema) is None


def test_deferral_cost_tolerates_a_partial_producer(conn, build_schema):
    """The step reports what it measured; the rest stays ``null``, not zero."""
    emit_build_manifest(
        conn,
        schema=build_schema,
        derive_cost=COST_INTERACTION_TABLES,
        deferral_cost={'seconds_saved': 1105.0, 'revalidate_seconds': 44.6},
    )
    deferral = _deferral(conn, build_schema)

    assert deferral['seconds_saved'] == pytest.approx(1105.0)
    assert deferral['revalidate_seconds'] == pytest.approx(44.6)
    assert deferral['drop_seconds'] is None
    assert deferral['restore_seconds'] is None
    assert deferral['constraints_deferred'] is None
    assert deferral['indexes_deferred'] is None
    assert deferral['catalogue_unchanged'] is None
    # A producer that measured a saving deferred something, whether or not it
    # said so — the same inference `materialised` makes from `table`.
    assert deferral['deferred'] is True


def test_malformed_deferral_cost_is_warned_about_and_dropped(
    conn,
    build_schema,
    caplog,
):
    """An unreadable number costs a warning, never the build."""
    with caplog.at_level(logging.WARNING, logger='omnipath_build.db.resources'):
        stats = emit_build_manifest(
            conn,
            schema=build_schema,
            derive_cost=COST_INTERACTION_TABLES,
            deferral_cost={
                'seconds_saved': 'eighteen minutes',
                'revalidate_seconds': 44.6,
            },
        )
    deferral = _deferral(conn, build_schema)

    assert stats.build_id
    assert deferral['seconds_saved'] is None
    assert deferral['revalidate_seconds'] == pytest.approx(44.6)
    assert any(
        'deferral-cost field' in record.getMessage() for record in caplog.records
    )

    # A record with nothing readable in it at all is dropped, not stored empty.
    with caplog.at_level(logging.WARNING, logger='omnipath_build.db.resources'):
        emit_build_manifest(
            conn,
            schema=build_schema,
            derive_cost=COST_INTERACTION_TABLES,
            deferral_cost={'how_long': 'a while'},
        )
    assert _deferral(conn, build_schema) is None


def test_table_scope_and_deferral_cost_do_not_move_the_build_id(
    conn, build_schema
):
    """Two builds differing only in these costs share one identity.

    The amendment added fields to the manifest; the cycle-007 lesson says a
    field describing a run may never move the identity of what the run built.
    Asserted here so a later edit cannot fold them into the hashed payload
    unnoticed.
    """
    first = emit_build_manifest(
        conn,
        schema=build_schema,
        derive_cost=COST_INTERACTION_TABLES,
        scope_cost=SCOPE_COST,
        deferral_cost=DEFERRAL_COST,
    )
    _id, commits, resources, first_cost, _p = _manifest(conn, build_schema)
    first_deferral = _deferral(conn, build_schema)

    slower = {
        step: {'seconds': cost['seconds'] * 3 + 1.0, 'rows': cost['rows']}
        for step, cost in COST_INTERACTION_TABLES.items()
    }
    fewer_scopes = {
        'signor_only': {
            'table': None,
            'materialised': False,
            'seconds': 9999.0,
            'rows': 1,
        },
    }
    second = emit_build_manifest(
        conn,
        schema=build_schema,
        derive_cost=slower,
        scope_cost=fewer_scopes,
        # This run took the deferral and bought nothing by it.
        deferral_cost={'deferred': True, 'seconds_saved': 0.0},
    )
    _id2, _c2, _r2, second_cost, _p2 = _manifest(conn, build_schema)

    assert first_cost['interaction_tables'] != second_cost['interaction_tables']
    assert first_cost['scopes'] != second_cost['scopes']
    assert first_deferral != _deferral(conn, build_schema)
    assert first.build_id == second.build_id

    # And the identity is still the hash of the two names it always covered.
    payload = {'package_commits': commits, 'resources': resources}
    assert first.build_id == hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()[:12]
