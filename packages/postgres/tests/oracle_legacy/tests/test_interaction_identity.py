"""The interaction header id is content-addressed.

``interaction.interaction_id`` is a deterministic hash of the **sorted
participant multiset** and the interaction class, so it is endpoint-independent:
the same participants handed over in a different order name the same
interaction, and a rebuild over unchanged inputs mints no new headers. That is
what makes ``ON CONFLICT DO NOTHING`` a safe, additive rebuild and what lets a
caller ask for an interaction from either endpoint.

The identity is defined twice — as a DuckDB macro next to ``content_uuid``
(``duckdb_load.py``, the load side) and as a SQL expression in the derive step
(``db/derived_tables.py``, the Postgres side). Both must produce the same uuid
for the same content, so the cross-engine agreement is asserted here rather than
assumed.

Run::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_interaction_identity.py -v

The macro tests need no database; the projection tests are skipped without one.
"""

from __future__ import annotations

import os

import pytest

from tests.fixtures.interaction_graph import build_interaction_fixture

DATABASE_URL = os.environ.get('DATABASE_URL')
SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_INTERACTION_IDENTITY',
    'interaction_identity_test',
)

ENTITY_ONE = 'a0000000-0000-4000-8000-000000000001'
ENTITY_TWO = 'b0000000-0000-4000-8000-000000000002'


@pytest.fixture(scope='module')
def duck():
    """A DuckDB connection carrying the load-side identity macros."""
    import duckdb

    from omnipath_build.duckdb_load import (
        _create_duckdb_content_uuid_macro,
    )

    connection = duckdb.connect()
    _create_duckdb_content_uuid_macro(connection)
    try:
        yield connection
    finally:
        connection.close()


def _duck_uuid(duck, participants: list[str], interaction_class: str) -> str:
    literal = ', '.join(f"'{value}'::UUID" for value in participants)
    return str(
        duck.sql(
            f"SELECT interaction_content_uuid([{literal}], "
            f"'{interaction_class}')"
        ).fetchone()[0]
    )


def test_participant_order_does_not_change_the_id(duck):
    """Sorting is part of the identity: input order is not."""
    forward = _duck_uuid(duck, [ENTITY_ONE, ENTITY_TWO], 'signaling')
    reverse = _duck_uuid(duck, [ENTITY_TWO, ENTITY_ONE], 'signaling')
    assert forward == reverse


def test_the_class_is_part_of_the_identity(duck):
    """Same participants, different class, different interaction."""
    signaling = _duck_uuid(duck, [ENTITY_ONE, ENTITY_TWO], 'signaling')
    transport = _duck_uuid(duck, [ENTITY_ONE, ENTITY_TWO], 'transport')
    assert signaling != transport


def test_the_participants_are_a_multiset_not_a_set(duck):
    """A homodimer keeps both of its participants, so it is its own interaction."""
    single = _duck_uuid(duck, [ENTITY_ONE], 'signaling')
    doubled = _duck_uuid(duck, [ENTITY_ONE, ENTITY_ONE], 'signaling')
    assert single != doubled


@pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the cross-engine check needs a Postgres',
)
def test_duckdb_and_postgres_mint_the_same_id(duck):
    """One identity, two engines: the load side and the derive side agree."""
    import psycopg2

    from omnipath_build.db.derived_tables import interaction_content_uuid_sql

    expression = interaction_content_uuid_sql(
        participants='ARRAY[%s::uuid, %s::uuid]',
        interaction_class='%s',
    )
    # The class comes first in the payload, then the participants — and the
    # participants are handed over in the opposite order to the macro call
    # below, which is the point of the assertion.
    connection = psycopg2.connect(DATABASE_URL)
    connection.autocommit = True
    try:
        with connection.cursor() as cur:
            cur.execute(
                f'SELECT {expression}',
                ['signaling', ENTITY_TWO, ENTITY_ONE],
            )
            from_postgres = str(cur.fetchone()[0])
    finally:
        connection.close()
    assert from_postgres == _duck_uuid(
        duck,
        [ENTITY_ONE, ENTITY_TWO],
        'signaling',
    )


@pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the projection test needs a Postgres',
)
class TestHeadersOverTheFixtureGraph:
    """The derive step's headers, over the hand-built graph."""

    @pytest.fixture(scope='class')
    def built(self):
        import psycopg2

        from omnipath_build.db import schema as build_schema
        from omnipath_build.db.derived_tables import rebuild_interaction_tables

        connection = psycopg2.connect(DATABASE_URL)
        try:
            build_schema.ensure_schema(
                connection,
                schema=SCRATCH,
                drop_existing=True,
            )
            connection.commit()
            build_interaction_fixture(connection, SCRATCH)
            first = rebuild_interaction_tables(connection, schema=SCRATCH)
            second = rebuild_interaction_tables(connection, schema=SCRATCH)
            yield connection, first, second
        finally:
            with connection.cursor() as cur:
                cur.execute(f'DROP SCHEMA IF EXISTS {SCRATCH} CASCADE')
            connection.commit()
            connection.close()

    def test_the_rebuild_is_idempotent(self, built):
        """A second run over unchanged inputs adds no header rows."""
        connection, first, second = built
        assert first.interactions == second.interactions
        with connection.cursor() as cur:
            cur.execute(f'SELECT count(*) FROM {SCRATCH}.interaction')
            assert cur.fetchone()[0] == first.interactions

    def test_both_directions_share_one_header(self, built):
        """A→B and B→A are two facts of one endpoint-independent interaction."""
        connection, _first, _second = built
        # Counted over the record's collapse key rather than over a collapsed
        # table: the materialisation went, and the claim — two facts,
        # one header — is a property of the key and the header link, both of
        # which the record carries per contributing resource.
        with connection.cursor() as cur:
            cur.execute(
                f"""
                SELECT count(DISTINCT r.interaction_id),
                       count(DISTINCT (r.subject_entity_id,
                                       r.object_entity_id,
                                       r.interaction_class_id))
                FROM {SCRATCH}.interaction_fact_resource r
                JOIN {SCRATCH}.entity subject
                  ON subject.entity_id = r.subject_entity_id
                JOIN {SCRATCH}.entity object
                  ON object.entity_id = r.object_entity_id
                WHERE subject.canonical_identifier IN ('FIXTURE_g', 'FIXTURE_h')
                  AND object.canonical_identifier IN ('FIXTURE_g', 'FIXTURE_h')
                """
            )
            headers, facts = cur.fetchone()
        assert facts == 2, 'the opposite-direction pair lost a row'
        assert headers == 1, 'the header is not endpoint-independent'

    def test_every_header_id_matches_its_participants(self, built):
        """The stored id is the hash of the participants actually recorded."""
        connection, _first, _second = built
        from omnipath_build.db.derived_tables import (
            interaction_content_uuid_sql,
        )

        expression = interaction_content_uuid_sql(
            participants='party.participants',
            interaction_class='vic.name',
        )
        with connection.cursor() as cur:
            cur.execute(
                f"""
                SELECT count(*)
                FROM {SCRATCH}.interaction i
                JOIN {SCRATCH}.vocab_interaction_class vic
                  ON vic.interaction_class_id = i.interaction_class_id
                JOIN (
                  SELECT interaction_id, array_agg(entity_id) AS participants
                  FROM {SCRATCH}.interaction_party
                  GROUP BY interaction_id
                ) party ON party.interaction_id = i.interaction_id
                WHERE i.interaction_id <> {expression}
                """
            )
            assert cur.fetchone()[0] == 0
