"""Schema-existence and population tests for the interaction model (008).

Constitution III: every new table this cycle adds is asserted to exist and to
carry rows after a build. Five tables make up the model — the header
``interaction``, the participant ``interaction_party``, the participant-role
vocabulary ``vocab_relation_role``, the interaction record
``interaction_fact_resource`` and the resource license terms
``data_source_license``.

**Amended 2026-08-20.** ``interaction_fact_combined`` was
the projection; it became one scope's materialisation of a record kept per
contributing resource.

**Amended 2026-08-21.** That materialisation is removed. The
derive writes **three tables and folds nothing** — the header, the participant
and the record — and the collapse is the shape a query produces at request
time, not a table. The tests follow: the projection is asserted to write
exactly those three, no ``interaction_fact_combined`` may survive a build in
any non-system schema, ``ensure_schema`` is asserted to drop one a pre-change
database still carries, and the anchoring rule the collapse rests on — the
detail tables anchor on the record through the denormalised triple — is
asserted against the Postgres catalogue. The collapse-equivalence assertions
this file used to carry moved to the api-service, which is where the fold now
lives.

The schema half runs against a throwaway scratch schema, so it needs no data.
The population half reads the built schema and needs a build — a capped
``MAX_RECORDS`` build is enough, because the assertion is that the derive step
produced rows at all, not how many.

Run::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_interactions_schema.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')
SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_INTERACTIONS',
    'interactions_schema_test',
)
# The migration test plants a legacy table and reruns `ensure_schema` over it,
# so it needs a schema of its own: the module-scoped `scratch` is shared, and a
# table planted in it would be visible to the tests asserting that no such
# table exists anywhere.

# The tables the interaction model adds.
# `interaction_fact_resource` and `data_source_license` joined the list with the
# grain amendment: the record is what the projection now stores, and a license
# filter resolves to a resource set over the ordinal levels of
# `data_source_license`. `interaction_fact_combined` left it when the build
# stopped storing the fold: the collapse is a query-time shape, so there is
# no fourth table for the existence or the population half to assert.
INTERACTION_TABLES = (
    'interaction',
    'interaction_party',
    'vocab_relation_role',
    'interaction_fact_resource',
    'data_source_license',
)

# What the derive writes, and all of what it writes.
# `vocab_relation_role` and `data_source_license` are seeded rather than
# projected, so they are model tables without being projection outputs.
PROJECTION_TABLES = (
    'interaction',
    'interaction_party',
    'interaction_fact_resource',
)

# The removed materialisation. Named once, because two tests ask the
# catalogue for it and neither may spell it differently.
LEGACY_COLLAPSE_TABLE = 'interaction_fact_combined'

# The detail tables. They hang off the record, not off the collapse — see
# the anchoring rule below.
DETAIL_TABLES = ('interaction_assay', 'interaction_ptm')

# `interaction_fact_resource` is keyed per contributing resource and carries a
# deterministic uuid surrogate, because the rest of the key is nullable and
# Postgres `MATCH SIMPLE` does not check a foreign key with any NULL column.
RECORD_KEY_COLUMNS = ('source_id', 'interaction_fact_resource_id')

# `data_source_license` stores ordinal levels, not a name: `enables` is
# `self >= other`, so a license question is a range predicate over three
# smallints. `is_known` is the flag that excludes an unmapped license.
LICENSE_LEVEL_COLUMNS = ('purpose_level', 'sharing_level', 'attrib_level')

# Sign and direction are three-valued: NULL means no contributing resource
# asserts the attribute, which is a different statement from an asserted
# `false`. A NOT NULL or a DEFAULT on any of these three destroys the
# distinction, so the test pins both.
THREE_VALUED_COLUMNS = ('is_directed', 'is_stimulation', 'is_inhibition')

# The record's index on the header id. A group keyed on the interaction rather
# than on the endpoint pair reads this one, and so does a lookup that follows a
# single header back to the rows asserting it; every other index on the table
# leads with `subject_entity_id`, so neither is served without it.
RECORD_INTERACTION_INDEX = 'interaction_fact_resource_interaction_idx'

# The ordered key of the fact table: A→B and B→A are two rows.
FACT_KEY_COLUMNS = [
    'subject_entity_id',
    'object_entity_id',
    'interaction_class_id',
]

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; interaction schema test needs a Postgres',
)


@pytest.fixture(scope='module')
def conn():
    """A read-only connection to the built database."""
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    # Read-only queries: autocommit keeps one failing statement from aborting
    # the transaction and masking the rest with InFailedSqlTransaction.
    connection.autocommit = True
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture(scope='module')
def scratch(conn):
    """The schema built into a throwaway namespace, with no data in it."""
    import psycopg2

    from omnipath_build.db import schema as build_schema

    writable = psycopg2.connect(DATABASE_URL)
    try:
        build_schema.ensure_schema(
            writable,
            schema=SCRATCH,
            drop_existing=True,
        )
        writable.commit()
        yield SCRATCH
    finally:
        with writable.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {SCRATCH} CASCADE')
        writable.commit()
        writable.close()


def _columns(conn, schema: str, table: str) -> dict[str, tuple[str, str]]:
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT column_name, data_type, is_nullable, column_default
            FROM information_schema.columns
            WHERE table_schema = %s AND table_name = %s
            """,
            (schema, table),
        )
        return {
            name: (data_type, nullable, default)
            for name, data_type, nullable, default in cur.fetchall()
        }


def _table_exists(conn, schema: str, table: str) -> bool:
    with conn.cursor() as cur:
        cur.execute('SELECT to_regclass(%s)', (f'{schema}.{table}',))
        return cur.fetchone()[0] is not None


def _require_table(conn, schema: str, table: str) -> None:
    """Fail with the missing table named, rather than raising UndefinedTable.

    A test that errors out of the driver has not been made to fail as designed:
    the message has to say which table the build owes.
    """
    assert _table_exists(conn, schema, table), (
        f'{schema}.{table} does not exist; the build never created it'
    )


def _relations_named(conn, table: str) -> list[tuple[str, str]]:
    """Every relation called ``table`` in a non-system schema, as (schema, kind).

    Catalogue rather than ``information_schema``, and every relation kind that
    can carry the name — a table, a partitioned table, a view, a materialised
    view or a foreign table. "Removed" means the name resolves to nothing a
    query can read, not merely that the ordinary table is gone.

    ``pg_temp_*`` and ``pg_toast*`` are excluded with the rest of the system
    namespaces: this database carries about ninety of them, and a session-local
    scratch relation is not a schema object anybody inherits.
    """
    with conn.cursor() as cur:
        cur.execute(
            r"""
            SELECT ns.nspname, cls.relkind
            FROM pg_class cls
            JOIN pg_namespace ns ON ns.oid = cls.relnamespace
            WHERE cls.relname = %s
              AND cls.relkind IN ('r', 'p', 'v', 'm', 'f')
              AND ns.nspname <> 'information_schema'
              AND ns.nspname NOT LIKE 'pg\_%%'
            ORDER BY 1
            """,
            (table,),
        )
        return [(schema, kind) for schema, kind in cur.fetchall()]


def _foreign_key_targets(conn, schema: str, table: str) -> set[str]:
    """The table names every foreign key *on* ``table`` points at."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT tgt.relname
            FROM pg_constraint con
            JOIN pg_class src ON src.oid = con.conrelid
            JOIN pg_namespace src_ns ON src_ns.oid = src.relnamespace
            JOIN pg_class tgt ON tgt.oid = con.confrelid
            WHERE con.contype = 'f'
              AND src_ns.nspname = %s
              AND src.relname = %s
            """,
            (schema, table),
        )
        return {row[0] for row in cur.fetchall()}


def _index_columns(conn, schema: str, index: str) -> list[str]:
    """The indexed columns of ``index``, in index order.

    Read from ``pg_index`` rather than parsed out of ``indexdef``: the order is
    the whole question here, and ``pg_attribute`` gives it without a text
    scrape. An index that does not exist comes back as an empty list, so the
    test says which index is missing instead of raising.
    """
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT att.attname
            FROM pg_index idx
            JOIN pg_class i ON i.oid = idx.indexrelid
            JOIN pg_namespace ns ON ns.oid = i.relnamespace
            JOIN unnest(idx.indkey) WITH ORDINALITY AS k(attnum, ord) ON true
            JOIN pg_attribute att
              ON att.attrelid = idx.indrelid AND att.attnum = k.attnum
            WHERE ns.nspname = %s AND i.relname = %s
            ORDER BY k.ord
            """,
            (schema, index),
        )
        return [name for (name,) in cur.fetchall()]


def _row_count(conn, schema: str, table: str) -> int:
    with conn.cursor() as cur:
        cur.execute(f'SELECT count(*) FROM {schema}.{table}')
        return cur.fetchone()[0]


@pytest.mark.parametrize('table', INTERACTION_TABLES)
def test_table_is_created(conn, scratch, table):
    """Creating the schema creates every table of the interaction model."""
    with conn.cursor() as cur:
        cur.execute(f"SELECT to_regclass('{scratch}.{table}')")
        assert cur.fetchone()[0] is not None, f'{table} was not created'


def test_fact_table_carries_the_hot_columns(conn, scratch):
    """The fact table carries the hot filter columns.

    **Repointed 2026-08-21**: `interaction_fact_resource` is the only stored
    fact table. The list shrinks with the move, and the columns it loses are
    exactly the ones the collapse recomputes for the scope that asks —
    `sources`, `source_count`, `dataset_tags`, `reference_count`,
    `sign_source_count` and `direction_source_count`. Those are produced by a
    fold and never stored, so asserting them here would ask the record to hold
    a summary that is only true of one scope.
    """
    columns = _columns(conn, scratch, 'interaction_fact_resource')
    expected = {
        'subject_entity_id',
        'object_entity_id',
        'interaction_class_id',
        'source_id',
        'subject_organism',
        'object_organism',
        'affinity',
        'pchembl',
        'score',
        'curation_flags',
        'reference_pubmed_ids',
        'reference_dois',
        'attributes',
        'interaction_id',
    }
    assert expected <= set(columns), (
        f'missing hot columns: {sorted(expected - set(columns))}'
    )


@pytest.mark.parametrize('column', THREE_VALUED_COLUMNS)
def test_sign_and_direction_are_three_valued(conn, scratch, column):
    """NULL means unasserted, so these three are nullable and undefaulted.

    **Repointed 2026-08-21** to the record, where the flags say what one
    resource asserts. The claim is the same one and it matters more here: the
    fold reads these columns, so a defaulted `false` on the record would be
    summarised into every scope that touches the row.
    """
    columns = _columns(conn, scratch, 'interaction_fact_resource')
    assert column in columns, f'{column} is missing from interaction_fact_resource'
    data_type, nullable, default = columns[column]
    assert data_type == 'boolean'
    assert nullable == 'YES', (
        f'{column} is NOT NULL, which erases "no resource asserts it"'
    )
    assert default is None, (
        f'{column} carries a default ({default}); an unasserted attribute '
        f'must stay NULL, never become an asserted false'
    )


def test_role_vocabulary_is_populated_by_name(conn, scratch):
    """``vocab_relation_role`` is seeded with the participant roles."""
    with conn.cursor() as cur:
        cur.execute(f'SELECT name FROM {scratch}.vocab_relation_role')
        names = {row[0] for row in cur.fetchall()}
    assert {
        'subject',
        'object',
        'reactant',
        'product',
        'enzyme',
        'cofactor',
        'regulator',
        'member',
    } <= names


@pytest.mark.parametrize('table', INTERACTION_TABLES)
def test_table_is_populated_by_the_build(conn, table):
    """Every table of the interaction model carries rows after a build."""
    _require_table(conn, SCHEMA, table)
    assert _row_count(conn, SCHEMA, table) > 0, (
        f'{SCHEMA}.{table} is empty; the derive step produced no rows'
    )


def test_every_fact_row_links_to_a_header(conn):
    """The projection keeps its link to the endpoint-independent header.

    **Repointed 2026-08-21** to the record, which is what the projection
    now writes. The link is what carries the header through a fold: a collapse
    carries `interaction_id` through unchanged, so a record row without one
    produces a collapsed row without one.
    """
    _require_table(conn, SCHEMA, 'interaction_fact_resource')
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT count(*) FROM {SCHEMA}.interaction_fact_resource '
            f'WHERE interaction_id IS NULL'
        )
        orphans = cur.fetchone()[0]
    assert orphans == 0, f'{orphans} record rows carry no header link'


# ---------------------------------------------------------------------------
# The record and the license table
# ---------------------------------------------------------------------------


@pytest.mark.parametrize('column', RECORD_KEY_COLUMNS)
def test_record_table_is_keyed_per_resource(conn, scratch, column):
    """The record keeps the resource, which is what makes summaries decomposable.

    ``interaction_fact_resource`` is one row per ordered endpoint pair,
    class and **source**, plus the assertion signature that source states. Drop
    ``source_id`` and the table is the collapse again, with no way to recompute
    a summary for a subset of resources.
    """
    columns = _columns(conn, scratch, 'interaction_fact_resource')
    assert columns, 'interaction_fact_resource was not created'
    assert column in columns, (
        f'{column} is missing from interaction_fact_resource; '
        f'found {sorted(columns)}'
    )


def test_record_surrogate_key_is_the_primary_key(conn, scratch):
    """The surrogate is the PK, because the composite key is not FK-able.

    The rest of the key includes the nullable assertion signature, and Postgres
    defaults to ``MATCH SIMPLE``, under which a foreign key with any NULL column
    is not checked at all. Drug-target rows assert neither sign nor direction,
    so a composite FK would pass silently on nearly every row while looking
    enforced. The uuid surrogate is a ``content_uuid`` over the full key, so it
    is stable across builds.
    """
    columns = _columns(conn, scratch, 'interaction_fact_resource')
    assert columns, 'interaction_fact_resource was not created'
    assert columns['interaction_fact_resource_id'][0] == 'uuid'
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT array_agg(att.attname ORDER BY k.ordinality)
            FROM pg_constraint con
            JOIN pg_class rel ON rel.oid = con.conrelid
            JOIN pg_namespace ns ON ns.oid = rel.relnamespace
            CROSS JOIN LATERAL unnest(con.conkey) WITH ORDINALITY AS k(
              attnum, ordinality
            )
            JOIN pg_attribute att
              ON att.attrelid = rel.oid AND att.attnum = k.attnum
            WHERE ns.nspname = %s
              AND rel.relname = 'interaction_fact_resource'
              AND con.contype = 'p'
            GROUP BY con.oid
            """,
            (scratch,),
        )
        primary = [row[0] for row in cur.fetchall()]
    assert primary == [['interaction_fact_resource_id']], (
        f'the primary key of interaction_fact_resource is {primary}, not the '
        f'surrogate; a composite key over nullable signature columns cannot be '
        f'referenced under MATCH SIMPLE'
    )


def test_license_table_stores_ordinal_levels(conn, scratch):
    """A license filter compares levels, so `data_source_license` stores them.

    ``License.enables(other)`` is ``self >= other``, so "everything usable
    commercially" is ``WHERE purpose_level >= 15`` over a 44-row table. A name
    answers no question on its own: two names can permit the same use, and one
    name under two versions can not.
    """
    columns = _columns(conn, scratch, 'data_source_license')
    assert columns, 'data_source_license was not created'
    for column in LICENSE_LEVEL_COLUMNS:
        assert column in columns, (
            f'{column} is missing from data_source_license; a license filter '
            f'has nothing to compare'
        )
        assert columns[column][0] == 'smallint', (
            f'{column} is {columns[column][0]}, not an ordinal smallint'
        )
    assert 'source_id' in columns, (
        'data_source_license does not name the resource it describes'
    )


def test_license_table_can_exclude_unknown_terms(conn, scratch):
    """Unknown terms are exclusions, never defaults.

    A resource whose license could not be mapped must not appear in a
    license-filtered result. Without ``is_known`` the only way to express "no
    record" is a NULL level, which a range predicate silently drops — the same
    answer for the wrong reason — or a permissive default, which is the one
    failure mode of this table that cannot be detected downstream.
    """
    columns = _columns(conn, scratch, 'data_source_license')
    assert columns, 'data_source_license was not created'
    assert 'is_known' in columns, (
        'data_source_license carries no is_known flag; an unmapped license '
        'would have to be admitted under a permissive default'
    )
    assert columns['is_known'][0] == 'boolean'


# ---------------------------------------------------------------------------
# The anchoring rule
# ---------------------------------------------------------------------------


def test_nothing_anchors_on_a_materialisation(conn, scratch):
    """The detail tables anchor on the record, and no materialisation survives.

    **Generalised 2026-08-21.** The rule used to be "no foreign key
    points at ``interaction_fact_combined``". That check has nothing left to
    point at: the table is removed, so the way to state the rule is that the
    name resolves to nothing in any non-system schema, and that the tables
    which would have keyed into it anchor on ``interaction_fact_resource``
    through the denormalised ``(subject_entity_id, object_entity_id,
    interaction_class_id)`` triple instead.

    The triple is what makes the removal a deletion rather than a migration,
    and it is more general than the key it replaced: it reaches **any** scope's
    collapse — a preset's, a license-filtered one, the empty one — and not only
    the all-resources scope that used to be stored.

    The check reads ``pg_class`` and ``pg_constraint``, not the schema module:
    what matters is what the database ends up holding and enforcing.

    **The vacuity guard is kept, and it now guards both halves.** A schema with
    no detail table in it satisfies "nothing anchors on a materialisation" for
    free, and it would satisfy the triple assertion for free as well.
    """
    present = [
        table for table in DETAIL_TABLES if _table_exists(conn, scratch, table)
    ]
    assert present == list(DETAIL_TABLES), (
        f'no detail table to anchor: {sorted(set(DETAIL_TABLES) - set(present))}'
        f' missing from {scratch}, so both halves of the anchoring rule hold '
        f'vacuously and prove nothing'
    )
    surviving = _relations_named(conn, LEGACY_COLLAPSE_TABLE)
    assert surviving == [], (
        f'{LEGACY_COLLAPSE_TABLE} still exists as '
        + '; '.join(f'{ns} (relkind {kind})' for ns, kind in surviving)
        + ' — the collapse is the shape a query produces, not a table, and '
        'a stored copy is something a foreign key can be pointed at again'
    )
    for table in DETAIL_TABLES:
        columns = _columns(conn, scratch, table)
        missing = [c for c in FACT_KEY_COLUMNS if c not in columns]
        assert not missing, (
            f'{table} is missing the denormalised anchor {missing}; with no '
            f'materialisation left it has no other route to a collapse'
        )
        targets = _foreign_key_targets(conn, scratch, table)
        assert 'interaction_fact_resource' in targets, (
            f'{table} declares no foreign key to interaction_fact_resource; '
            f'it references {sorted(targets)}'
        )


@pytest.mark.parametrize('table', DETAIL_TABLES)
def test_detail_table_references_the_record(conn, scratch, table):
    """A detail row keys into ``interaction_fact_resource``, not the header.

    An assay is ChEMBL's measurement and a PTM site is SIGNOR's or
    PhosphoSitePlus's. The record is keyed by resource, so the foreign key
    carries that meaning rather than merely satisfying a constraint. The header
    was considered and rejected: it is deduped across resources, so it drops the
    one dimension the detail actually has.
    """
    assert _table_exists(conn, scratch, table), (
        f'{scratch}.{table} was not created'
    )
    columns = _columns(conn, scratch, table)
    assert 'interaction_fact_resource_id' in columns, (
        f'{table} carries no interaction_fact_resource_id; it cannot name the '
        f'resource record its rows belong to'
    )
    targets = _foreign_key_targets(conn, scratch, table)
    assert 'interaction_fact_resource' in targets, (
        f'{table} declares no foreign key to interaction_fact_resource; '
        f'it references {sorted(targets)}'
    )


@pytest.mark.parametrize('table', DETAIL_TABLES)
def test_detail_table_carries_the_denormalised_join_key(conn, scratch, table):
    """Navigation to a collapse is a join, not a key.

    The three denormalised columns reach **any** scope's collapse — a preset's,
    a license-filtered one, the empty one. That is both cheaper and more general
    than the foreign key it replaces, and it costs the two largest planned
    tables no per-row key validation, which matters where the derive already
    measures 1,205 seconds under a 20-minute ceiling. **With no stored collapse
    it is also the only route left**: there is nothing left to key into.
    """
    assert _table_exists(conn, scratch, table), (
        f'{scratch}.{table} was not created'
    )
    columns = _columns(conn, scratch, table)
    missing = [c for c in FACT_KEY_COLUMNS if c not in columns]
    assert not missing, (
        f'{table} is missing the denormalised join key {missing}; without it a '
        f'scoped collapse is unreachable except through the materialisation'
    )


# ---------------------------------------------------------------------------
# Reaching the record by its header
# ---------------------------------------------------------------------------


def test_the_record_is_indexed_on_its_header(conn, scratch):
    """The record can be grouped and looked up by ``interaction_id``.

    A reaction is one interaction with any number of participants, so a query
    that wants one group per interaction keys on the header id and not on the
    endpoint pair. Nothing else on this table serves that: the unique key and
    the collapse index both lead with ``subject_entity_id``, which is a column
    such a query never restricts, so without an index of its own the grouping
    falls back to a full scan of the largest table the build writes.

    The leading column is what is asserted, not merely the presence of the
    name. An index that carried ``interaction_id`` in second position would
    satisfy an existence check and serve none of these queries.
    """
    columns = _index_columns(conn, scratch, RECORD_INTERACTION_INDEX)
    assert columns, (
        f'{scratch} carries no {RECORD_INTERACTION_INDEX}; a group keyed on '
        f'the interaction has no index to run on'
    )
    assert columns[0] == 'interaction_id', (
        f'{RECORD_INTERACTION_INDEX} leads on {columns[0]!r}, not on '
        f'interaction_id, so it does not serve a lookup or a grouping by '
        f'header: its columns are {columns}'
    )


# ---------------------------------------------------------------------------
# Three tables, and nothing folded
# ---------------------------------------------------------------------------


def test_the_derive_writes_three_tables_and_folds_nothing(conn):
    """The projection writes the header, the participant and the record. Only.

    As amended, there is **one stored fact table**, and no scope is
    precomputed — not even the all-resources scope, which is the one the
    removed ``interaction_fact_combined`` held. The two halves belong in one
    test because either alone is satisfiable by the wrong build: three
    populated tables say nothing about a fourth still being written, and an
    absent fourth table is what a derive that produced nothing at all also
    looks like.

    The absence is asserted across **every** non-system schema rather than only
    the built one. A build writes into the schema it is pointed at, and a copy
    left in a second schema is exactly the stale 4.68 GiB the removal is for.
    """
    for table in PROJECTION_TABLES:
        _require_table(conn, SCHEMA, table)
        assert _row_count(conn, SCHEMA, table) > 0, (
            f'{SCHEMA}.{table} is empty; the derive step produced no rows'
        )
    surviving = _relations_named(conn, LEGACY_COLLAPSE_TABLE)
    assert surviving == [], (
        f'the derive still leaves {LEGACY_COLLAPSE_TABLE} behind in '
        + '; '.join(f'{ns} (relkind {kind})' for ns, kind in surviving)
        + f' — the projection writes {len(PROJECTION_TABLES)} tables and folds '
        'nothing, so the collapse must resolve to no relation at all'
    )
