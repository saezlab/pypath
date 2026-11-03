"""A dataset that names a resource the build has not loaded.

This happens constantly and legitimately. A dataset is registered once and the
resources behind it arrive over several cycles; a capped development build
loads a fraction of them; a source is renamed and the old name is briefly
stale. So the question is not how to prevent it but what a build should do when
it meets it, and there are exactly three answers:

1. **Fail the dataset.** Wrong: one unloaded resource takes down eleven that
   are present, and the dataset that fails is usually the one being worked on.
2. **Contribute rows anyway.** Wrong in a way that does not announce itself:
   there are no rows to contribute, so the only way to produce any is to invent
   them.
3. **Contribute nothing, and say which resource and why.** The one this build
   takes.

The third answer is only different from the second if somebody hears it. A
silent empty contribution and a resource that genuinely has no interactions are
the same observation, and a developer who cannot tell them apart debugs the
wrong thing. So the log line is asserted as hard as the row count, and it must
name the resource.

Run against a build database, e.g. on dev4::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_preset_missing_resource.py -v
"""

from __future__ import annotations

import logging
import os

import pytest

from omnipath_build.network_views import (
    NetworkDefinition,
    ensure_network_registry,
    register_network,
)

DATABASE_URL = os.environ.get('DATABASE_URL')
TEST_SCHEMA = 'preset_missing_resource_test'

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the missing-resource rule needs a database',
)

# A name no loader will ever produce, so the test cannot pass because the
# resource quietly turned up.
ABSENT = '_no_such_resource_loaded'

# Present in the fixture, so the dataset has something to keep serving.
PRESENT = 'signor'


@pytest.fixture(scope='module')
def conn():
    """An open connection to the built database."""
    psycopg2 = pytest.importorskip('psycopg2')
    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture(scope='module')
def sandbox(conn):
    """A throwaway schema holding one loaded resource and one record row."""
    with conn.cursor() as cur:
        cur.execute(f'DROP SCHEMA IF EXISTS {TEST_SCHEMA} CASCADE')
        cur.execute(f'CREATE SCHEMA {TEST_SCHEMA}')
    conn.commit()
    ensure_network_registry(conn, registry_schema=TEST_SCHEMA)
    with conn.cursor() as cur:
        cur.execute(
            f"""
            CREATE TABLE {TEST_SCHEMA}.data_source (
              source_id bigint PRIMARY KEY, name text NOT NULL
            );
            INSERT INTO {TEST_SCHEMA}.data_source VALUES (1, '{PRESENT}');
            CREATE TABLE {TEST_SCHEMA}.interaction_fact_resource (
              subject_entity_id bigint NOT NULL,
              object_entity_id bigint NOT NULL,
              interaction_class_id int NOT NULL,
              source_id bigint NOT NULL
            );
            INSERT INTO {TEST_SCHEMA}.interaction_fact_resource VALUES (10, 20, 1, 1);
            """
        )
    conn.commit()
    try:
        yield TEST_SCHEMA
    finally:
        with conn.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {TEST_SCHEMA} CASCADE')
        conn.commit()


@pytest.fixture
def mixed_preset():
    """A dataset naming one loaded resource and one that is not."""
    return NetworkDefinition(
        name='_missing_resource_preset',
        kind='ligand_receptor',
        included_sources=(PRESENT, ABSENT),
        interaction_class_scope=('ligand_receptor',),
    )


def _contribution(conn, schema, name):
    """Record rows a resource contributes under a schema.

    Args:
        conn: An open connection.
        schema: The schema holding the record.
        name: The resource name.

    Returns:
        The row count, which is zero for a resource that is not loaded.
    """
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT count(*) FROM {schema}.interaction_fact_resource f
            JOIN {schema}.data_source d USING (source_id)
            WHERE d.name = %s
            """,
            [name],
        )
        return cur.fetchone()[0]


def test_the_dataset_registers(conn, sandbox, mixed_preset):
    """One unloaded resource does not take the dataset down with it."""
    register_network(conn, mixed_preset, registry_schema=sandbox)
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT included_sources FROM {sandbox}.network_registry WHERE name = %s',
            [mixed_preset.name],
        )
        row = cur.fetchone()
    assert row is not None, (
        'the dataset failed to register because one of its resources is not '
        'loaded; a dataset assembled over several cycles would never register '
        'at all'
    )
    assert ABSENT in row[0], (
        'the unloaded resource was dropped from the stored scope. It must stay '
        'there: it is what the dataset is waiting for, and removing it turns a '
        'known gap into a forgotten one'
    )


def test_the_missing_resource_contributes_nothing(conn, sandbox, mixed_preset):
    """No row is invented for a resource that has none."""
    register_network(conn, mixed_preset, registry_schema=sandbox)
    assert _contribution(conn, sandbox, ABSENT) == 0, (
        f'{ABSENT} contributes rows to the record and is not loaded; those '
        f'rows came from nowhere'
    )


def test_the_loaded_resource_still_contributes(conn, sandbox, mixed_preset):
    """The rest of the dataset keeps working."""
    register_network(conn, mixed_preset, registry_schema=sandbox)
    assert _contribution(conn, sandbox, PRESENT), (
        f'{PRESENT} contributes nothing either; the unloaded resource took the '
        f'whole dataset with it'
    )


def test_the_warning_names_the_resource(conn, sandbox, mixed_preset, caplog):
    """The empty contribution is announced, and announced specifically.

    Without the name the message is worth nothing: a developer reading "some
    resource is missing" has to work out which one by hand, which is exactly
    the work the message exists to save.
    """
    with caplog.at_level(logging.WARNING, logger='omnipath_build'):
        register_network(conn, mixed_preset, registry_schema=sandbox)
    warnings = [
        record.getMessage() for record in caplog.records
        if record.levelno >= logging.WARNING
    ]
    naming = [message for message in warnings if ABSENT in message]
    assert naming, (
        f'registering a dataset over an unloaded resource logged {warnings}, '
        f'none of which names {ABSENT}. An unannounced empty contribution is '
        f'indistinguishable from a resource with nothing to say'
    )


def test_no_warning_when_every_resource_is_loaded(conn, sandbox, caplog):
    """The warning means something, so it is not emitted for every dataset."""
    whole = NetworkDefinition(
        name='_fully_loaded_preset',
        kind='ligand_receptor',
        included_sources=(PRESENT,),
    )
    with caplog.at_level(logging.WARNING, logger='omnipath_build'):
        register_network(conn, whole, registry_schema=sandbox)
    noisy = [
        record.getMessage() for record in caplog.records
        if record.levelno >= logging.WARNING and PRESENT in record.getMessage()
    ]
    assert not noisy, (
        f'a dataset whose resources are all loaded warned anyway: {noisy}. A '
        f'warning every build emits is a warning nobody reads'
    )
