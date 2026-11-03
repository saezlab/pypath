"""Sign and direction, per resource and per scope.

**Amended 2026-08-20.** Sign and direction are asserted **per resource**, on
``interaction_fact_resource`` — the record. The collapse of that record over a
resource set is one scope's fold rather than the stored grain. So the file
holds three layers of assertion: what a single resource says on its own row,
what the all-resources collapse says, and what a collapse restricted to a
resource subset says. A resource that is silent leaves NULL on its own row and
never inherits what a neighbour on the same endpoint pair asserted. A scoped
collapse reports the kept resources' count, references and sign flags alone:
selecting a wider collapse's row and filtering it with ``sources &&
ARRAY[...]`` returns the right interaction carrying the wrong numbers, and the
scope tests here are built on a three-resource pair precisely so that
implementation fails them.

**Amended 2026-08-21.** The collapse is no longer a table and no longer a build
routine: ``interaction_fact_combined`` is dropped and
``collapse_interaction_scope`` is deleted, because the derive has no caller for
a fold and the query path has one, in the api-service. The collapse assertions
here read ``tests.fixtures.collapse`` instead — the fold written once on the
test side, as a view for the all-resources scope and as a statement for a
scoped one. Every assertion below is the one it was: they were about what a
collapse must report, never about where it was stored.

``is_directed``, ``is_stimulation`` and ``is_inhibition`` are **three-valued**.
NULL means *no contributing resource asserts the attribute*, which is a
different statement from an asserted ``false`` — 97.4 per cent of the graph's
edges carry no sign at all, so a design that defaulted them to ``false`` would
mislabel 13.9 million edges as "known not stimulatory". Where resources
disagree, **both** sign flags are true: the disagreement is represented, never
resolved by a silent winner. And ``sources``/``references`` list **every**
contributor, including the ones asserting neither sign nor direction, so
``sign_source_count <= cardinality(sources)`` is a real inequality.

The second half covers the sign-conflict summary the derive step measures on
every build — signed rows, both-flags-true rows and their percentage, and the
single-resource versus cross-resource split — is produced and reaches the build
manifest **without** entering the ``build_id`` hash, because it is a
measurement of the build, not part of its identity.

Run::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_interactions_sign_direction.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import hashlib
import json
import os

import pytest

from tests.fixtures.collapse import (
    COLLAPSE_VIEW,
    collapse_sql,
    create_collapse_view,
)
from tests.fixtures.interaction_graph import build_interaction_fixture

DATABASE_URL = os.environ.get('DATABASE_URL')
SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_INTERACTIONS_SIGN',
    'interactions_sign_test',
)

#: The record table: one row per endpoint pair, class, resource and assertion
#: signature. What a resource asserts lives here.
RECORD_TABLE = 'interaction_fact_resource'

#: The build routine that used to fold the record for a resource scope.
#: Dropping the materialisation it filled left the derive with no caller for it,
#: so it went too: one body with one caller, in the package that calls it. A
#: build that still exports it has not finished the removal.
REMOVED_COLLAPSE_ROUTINE = 'collapse_interaction_scope'

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the projection test needs a Postgres',
)


@pytest.fixture(scope='module')
def built():
    """The fixture graph, projected by the derive step, in a scratch schema."""
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
        stats = rebuild_interaction_tables(connection, schema=SCRATCH)
        # The all-resources collapse, as a view over the record. The derive
        # no longer writes such a table, and the assertions below are about
        # what the fold reports rather than about what the build stores.
        create_collapse_view(connection, SCRATCH)
        yield connection, stats
    finally:
        with connection.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {SCRATCH} CASCADE')
        connection.commit()
        connection.close()


#: The collapsed shape, in the order both helpers below select it.
_FACT_KEYS = (
    'is_directed',
    'is_stimulation',
    'is_inhibition',
    'sign_source_count',
    'direction_source_count',
    'sources',
    'source_count',
    'reference_pubmed_ids',
    'reference_count',
)


def _row(conn, subject: str, object_: str) -> dict[str, object]:
    """One row of the all-resources collapse, by endpoint pair."""
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT f.is_directed, f.is_stimulation, f.is_inhibition,
                   f.sign_source_count, f.direction_source_count,
                   f.sources, f.source_count,
                   f.reference_pubmed_ids, f.reference_count
            FROM {SCRATCH}.{COLLAPSE_VIEW} f
            JOIN {SCRATCH}.entity subject
              ON subject.entity_id = f.subject_entity_id
            JOIN {SCRATCH}.entity object
              ON object.entity_id = f.object_entity_id
            WHERE subject.canonical_identifier = %s
              AND object.canonical_identifier = %s
            """,
            (f'FIXTURE_{subject}', f'FIXTURE_{object_}'),
        )
        rows = cur.fetchall()
    assert len(rows) == 1, f'{subject}->{object_} produced {len(rows)} rows'
    return dict(zip(_FACT_KEYS, rows[0]))


def test_an_unasserted_sign_stays_null(built):
    """No resource asserts a sign, so the columns are NULL — never `false`."""
    conn, _stats = built
    row = _row(conn, 'e', 'f')
    assert row['is_stimulation'] is None, (
        'an unasserted sign became an asserted value'
    )
    assert row['is_inhibition'] is None
    assert row['sign_source_count'] == 0


def test_disagreeing_resources_set_both_flags(built):
    """One resource says stimulation, another inhibition: both surface."""
    conn, _stats = built
    row = _row(conn, 'c', 'd')
    assert row['is_stimulation'] is True
    assert row['is_inhibition'] is True


def test_a_contributor_asserting_no_sign_still_counts_as_provenance(built):
    """Attribute-poor evidence stays in `sources`."""
    conn, _stats = built
    row = _row(conn, 'c', 'd')
    assert 'fixture_res_c' in row['sources'], (
        'the resource asserting neither sign nor direction was dropped'
    )
    assert row['sign_source_count'] == 2
    assert row['sign_source_count'] <= len(row['sources'])
    assert row['source_count'] == len(row['sources'])


def test_the_sign_source_count_never_exceeds_the_sources(built):
    """The inequality holds on every row, not only the constructed one."""
    conn, _stats = built
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT count(*) FROM {SCRATCH}.{COLLAPSE_VIEW}
            WHERE sign_source_count > cardinality(sources)
               OR direction_source_count > cardinality(sources)
            """
        )
        assert cur.fetchone()[0] == 0


def test_an_opposite_direction_pair_is_two_rows(built):
    """A→B and B→A are never merged into one bidirectional record."""
    conn, _stats = built
    forward = _row(conn, 'g', 'h')
    reverse = _row(conn, 'h', 'g')
    assert forward['is_directed'] is True
    assert reverse['is_directed'] is True


def test_direction_is_null_where_no_predicate_speaks_to_it(built):
    """`has_member` says nothing about direction, so the column stays NULL."""
    conn, _stats = built
    row = _row(conn, 'q', 'r')
    assert row['is_directed'] is None
    assert row['direction_source_count'] == 0



def test_a_symmetric_predicate_does_not_assert_undirectedness(built):
    """`interacts_with` leaves `is_directed` NULL, it does not write false.

    Decided 2026-08-18, reverting a first pass that read a symmetric predicate as
    an assertion of undirectedness and wrote false on 8.2M rows. The predicate
    vocabulary is a coarse ontology layer the resources did not choose per
    interaction, so a symmetric verb is not a resource saying "this interaction
    has no direction". An unasserted attribute never becomes an asserted false.

    The pair is `e`-`f`, the symmetric verb over a class that says nothing about
    order. The ligand-receptor pair used to stand here and no longer can: its
    class asserts the order on its own, which the test below is about.
    """
    conn, _stats = built
    row = _row(conn, 'e', 'f')
    assert row['is_directed'] is None
    assert row['direction_source_count'] == 0


def test_a_ligand_receptor_row_is_directed_under_a_symmetric_predicate(built):
    """The class names its two endpoints asymmetrically, so the order is fixed.

    `a`-`b` carries `interacts_with`, the same verb as the pair above, and the
    same silence from its resource about direction. What differs is the class:
    the participant annotations make it ligand-receptor, and the projection
    stores the ligand as the subject. A row whose class is an ordered pair of
    roles is directed however coarse the verb the ingest layer gave it.
    """
    conn, _stats = built
    row = _row(conn, 'a', 'b')
    assert row['is_directed'] is True


def test_no_collapsed_row_ever_asserts_a_false_sign_or_direction(built):
    """The whole collapse, not one case: every flag is true or NULL.

    Widened 2026-08-20 from direction alone to both signs: never defaulting an
    unasserted attribute is one rule over three columns, and a defaulted
    `false` is as wrong on a sign as it is on a direction.
    """
    conn, _stats = built
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT count(*) FROM {SCRATCH}.{COLLAPSE_VIEW}
            WHERE is_directed IS FALSE
               OR is_stimulation IS FALSE
               OR is_inhibition IS FALSE
            """
        )
        assert cur.fetchone()[0] == 0

# --- Reopened 2026-08-20: the assertion is per resource --------------------


def _require_record_table(conn) -> None:
    """Say why the record grain is missing, rather than raising UndefinedTable."""
    with conn.cursor() as cur:
        cur.execute('SELECT to_regclass(%s)', (f'{SCRATCH}.{RECORD_TABLE}',))
        present = cur.fetchone()[0]
    if present is None:
        pytest.fail(
            f'{SCRATCH}.{RECORD_TABLE} does not exist: the derive still folds '
            f'every resource into one row, so what a single resource asserts '
            f'has nowhere to live'
        )


def _record(conn, subject: str, object_: str) -> dict[str, dict[str, object]]:
    """The record rows for one ordered endpoint pair, keyed by resource name."""
    _require_record_table(conn)
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT source.name, r.is_directed, r.is_stimulation,
                   r.is_inhibition, r.reference_pubmed_ids
            FROM {SCRATCH}.{RECORD_TABLE} r
            JOIN {SCRATCH}.entity subject
              ON subject.entity_id = r.subject_entity_id
            JOIN {SCRATCH}.entity object
              ON object.entity_id = r.object_entity_id
            JOIN {SCRATCH}.data_source source
              ON source.source_id = r.source_id
            WHERE subject.canonical_identifier = %s
              AND object.canonical_identifier = %s
            """,
            (f'FIXTURE_{subject}', f'FIXTURE_{object_}'),
        )
        rows = cur.fetchall()
    keys = ('is_directed', 'is_stimulation', 'is_inhibition', 'references')
    return {row[0]: dict(zip(keys, row[1:])) for row in rows}


def test_the_record_carries_one_row_per_contributing_resource(built):
    """Three resources report c->d, so the record holds three rows for it."""
    conn, _stats = built
    rows = _record(conn, 'c', 'd')
    assert set(rows) == {'fixture_res_a', 'fixture_res_b', 'fixture_res_c'}, (
        'the record grain is not per resource'
    )


def test_a_silent_resource_leaves_null_on_its_own_row(built):
    """The core of the per-resource amendment: silence is not inherited.

    `fixture_res_c` reports c->d and asserts no sign. Its neighbours on the
    same endpoint pair assert opposite signs. A row that carried the group's
    summary would show `fixture_res_c` as both stimulating and inhibiting,
    which is a claim that resource never made.
    """
    conn, _stats = built
    rows = _record(conn, 'c', 'd')
    silent = rows['fixture_res_c']
    assert silent['is_stimulation'] is None, (
        'a silent resource inherited a neighbour sign assertion'
    )
    assert silent['is_inhibition'] is None, (
        'a silent resource inherited a neighbour sign assertion'
    )
    # And the two that do assert keep their own sign, not the group summary.
    assert rows['fixture_res_a']['is_stimulation'] is True
    assert rows['fixture_res_a']['is_inhibition'] is None
    assert rows['fixture_res_b']['is_inhibition'] is True
    assert rows['fixture_res_b']['is_stimulation'] is None


def test_a_resource_silent_on_direction_leaves_null_on_its_own_row(built):
    """Direction is per resource too, and NULL where none is asserted."""
    conn, _stats = built
    unsigned = _record(conn, 'e', 'f')['fixture_res_c']
    assert unsigned['is_directed'] is None
    directed = _record(conn, 'g', 'h')['fixture_res_a']
    assert directed['is_directed'] is True


def test_the_record_holds_only_that_resources_references(built):
    """A reference belongs to the resource that supplied it, not to the pair."""
    conn, _stats = built
    rows = _record(conn, 'c', 'd')
    assert rows['fixture_res_a']['references'] == ['11111111']
    assert rows['fixture_res_c']['references'] == ['33333333']
    assert not rows['fixture_res_b']['references']


def test_no_record_row_ever_asserts_a_false_sign_or_direction(built):
    """At the record grain too: unasserted stays NULL, never `false`."""
    conn, _stats = built
    _require_record_table(conn)
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT count(*) FROM {SCRATCH}.{RECORD_TABLE}
            WHERE is_directed IS FALSE
               OR is_stimulation IS FALSE
               OR is_inhibition IS FALSE
            """
        )
        assert cur.fetchone()[0] == 0, (
            'a record row defaulted an unasserted attribute to false'
        )


def test_an_opposite_direction_pair_is_two_records(built):
    """A→B and B→A stay two rows at the record grain, as at the collapse."""
    conn, _stats = built
    forward = _record(conn, 'g', 'h')
    reverse = _record(conn, 'h', 'g')
    assert set(forward) == set(reverse) == {'fixture_res_a'}
    assert forward['fixture_res_a']['is_directed'] is True
    assert reverse['fixture_res_a']['is_directed'] is True


# --- The scope case: a collapse restricted to a subset of resources ---------


def _collapse_query(sources) -> tuple[str, tuple[object, ...]]:
    """The scoped collapse, as a statement the test embeds as a subquery."""
    statement, parameters = collapse_sql(schema=SCRATCH, sources=list(sources))
    return statement, tuple(parameters or ())


def test_the_build_exports_no_collapse_routine():
    """The fold left the build, and it left no stub behind.

    Without the materialisation the derive has no caller for a fold, and the
    query path has one, in the api-service. A routine still exported from
    ``omnipath_build.db.derived_tables`` would be the second body the removal
    exists to prevent — and one the api-service cannot import in any case: its
    ``pyproject.toml`` declares no build package, and the import only ever
    resolved because both checkouts sit on one path (Principle I/II).
    """
    from omnipath_build.db import derived_tables

    assert not hasattr(derived_tables, REMOVED_COLLAPSE_ROUTINE), (
        f'omnipath_build.db.derived_tables.{REMOVED_COLLAPSE_ROUTINE} still '
        f'exists: the fold has one caller, and that caller is not in '
        f'this repository'
    )


def _scoped(conn, subject: str, object_: str, sources) -> dict[str, object]:
    """One collapsed row for an endpoint pair, over `sources` alone."""
    statement, parameters = _collapse_query(sources)
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT f.is_directed, f.is_stimulation, f.is_inhibition,
                   f.sign_source_count, f.direction_source_count,
                   f.sources, f.source_count,
                   f.reference_pubmed_ids, f.reference_count
            FROM ({statement}) f
            JOIN {SCRATCH}.entity subject
              ON subject.entity_id = f.subject_entity_id
            JOIN {SCRATCH}.entity object
              ON object.entity_id = f.object_entity_id
            WHERE subject.canonical_identifier = %s
              AND object.canonical_identifier = %s
            """,
            parameters + (f'FIXTURE_{subject}', f'FIXTURE_{object_}'),
        )
        rows = cur.fetchall()
    assert len(rows) == 1, (
        f'{subject}->{object_} scoped to {list(sources)} produced '
        f'{len(rows)} rows'
    )
    return dict(zip(_FACT_KEYS, rows[0]))


def test_the_scope_fixture_is_not_vacuous(built):
    """The scope guard: c->d must be a row several resources report.

    97.43 per cent of the built rows carry a single resource and pass the
    scope tests whether or not the rule is implemented. This pair carries
    three, and the all-resources collapse disagrees with the single-resource
    scope on the count, the references and both sign flags — so reading the
    collapsed row and filtering it with `sources && ARRAY[...]` cannot pass
    the three tests below.
    """
    conn, _stats = built
    combined = _row(conn, 'c', 'd')
    assert combined['source_count'] == 3
    assert combined['is_stimulation'] is True
    assert combined['is_inhibition'] is True
    assert sorted(combined['reference_pubmed_ids']) == ['11111111', '33333333']


def test_a_scoped_collapse_reports_only_that_resources_assertion(built):
    """Restricted to one resource, the collapse reports its assertion alone."""
    conn, _stats = built
    row = _scoped(conn, 'c', 'd', ['fixture_res_b'])
    assert row['is_inhibition'] is True
    assert row['is_stimulation'] is None, (
        'the scoped row reports a sign asserted by a resource outside the scope'
    )
    assert row['sign_source_count'] == 1
    assert row['sources'] == ['fixture_res_b']


def test_a_resource_scoped_query_counts_only_that_resource(built):
    """`source_count` is one, not the three of the wider fold."""
    conn, _stats = built
    row = _scoped(conn, 'c', 'd', ['fixture_res_a'])
    assert row['source_count'] == 1
    assert row['sources'] == ['fixture_res_a']


def test_a_resource_scoped_query_returns_only_that_resources_references(built):
    """The references of the kept resource, and no other."""
    conn, _stats = built
    row = _scoped(conn, 'c', 'd', ['fixture_res_a'])
    assert row['reference_pubmed_ids'] == ['11111111']
    assert row['reference_count'] == 1


def test_a_resource_scoped_query_reports_only_that_resources_sign(built):
    """The flags say what `fixture_res_a` asserts, nothing more."""
    conn, _stats = built
    row = _scoped(conn, 'c', 'd', ['fixture_res_a'])
    assert row['is_stimulation'] is True
    assert row['is_inhibition'] is None, (
        'the scoped row carries the inhibition another resource asserted'
    )
    assert row['sign_source_count'] == 1


def test_a_scope_on_the_silent_resource_reports_no_sign_at_all(built):
    """The other half of the scope rule: scoping to the contributor that asserts
    nothing yields NULL sign flags, not the pair's summary."""
    conn, _stats = built
    row = _scoped(conn, 'c', 'd', ['fixture_res_c'])
    assert row['is_stimulation'] is None
    assert row['is_inhibition'] is None
    assert row['sign_source_count'] == 0
    assert row['source_count'] == 1
    assert row['reference_pubmed_ids'] == ['33333333']


# --- The sign-conflict summary reaches the manifest -------------------------


def test_the_step_measures_the_sign_conflict_rate(built):
    """Every build reports how often both sign flags land on one row."""
    _conn, stats = built
    conflict = stats.sign_conflict
    # Five rows carry a sign: c->d and i->j, plus the three orthosteric pairs,
    # whose `Agonist` annotation is both a class and a positive sign. There is
    # no orthosteric fixture row that is not also signed, and that is the data
    # rather than the fixture: every term that names the class — agonist,
    # antagonist, inhibitor, activator — is a sign term as well.
    assert conflict['signed_rows'] == 5
    assert conflict['both_flags_rows'] == 2
    assert conflict['both_flags_percent'] == pytest.approx(100.0 * 2 / 5, abs=0.01)
    # c->d is two resources disagreeing; i->j is one resource asserting both
    # under two predicates that share the `signaling` class.
    assert conflict['single_resource_rows'] == 1
    assert conflict['cross_resource_rows'] == 1


def test_the_summary_is_recorded_next_to_the_derive_cost(built):
    """The measurement reaches the manifest, outside the identity hash."""
    from omnipath_build.db.resources import emit_build_manifest

    conn, _stats = built
    with conn.cursor() as cur:
        cur.execute(
            f"""
            CREATE TABLE IF NOT EXISTS {SCRATCH}.resources (
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
            f"INSERT INTO {SCRATCH}.resources (resource_id) VALUES ('fixture') "
            f'ON CONFLICT DO NOTHING'
        )
    conn.commit()

    stats = emit_build_manifest(
        conn,
        schema=SCRATCH,
        derive_cost={'interaction_fact_resource': {'seconds': 1.0, 'rows': 12}},
    )
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT build_id, package_commits, resources, '
            f'interactions_derive_cost FROM {SCRATCH}.build_manifest'
        )
        build_id, commits, resources, cost = cur.fetchone()

    conflict = cost['sign_conflict']
    assert set(conflict) >= {
        'signed_rows',
        'both_flags_rows',
        'both_flags_percent',
        'single_resource_rows',
        'cross_resource_rows',
    }
    assert conflict['signed_rows'] == 5

    # And it stayed out of the identity: the hash covers the package commits
    # and the resource inventory, nothing else.
    payload = {'package_commits': commits, 'resources': resources}
    expected = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()[:12]
    assert build_id == expected == stats.build_id
