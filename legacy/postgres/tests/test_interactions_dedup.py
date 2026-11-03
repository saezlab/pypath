"""The interaction grain and its collapse.

**Amended 2026-08-20.** The dedup claim moved. One row per
ordered ``(subject_entity_id, object_entity_id, interaction_class_id)`` is the
contract of the **collapsed output for a stated scope**, not the stored grain.
So this file asserts two things that used to be one:

* ``interaction_fact_resource`` is the **record**: one row per
  ``(subject_entity_id, object_entity_id, interaction_class_id, source_id)``
  **plus the assertion signature** ``(is_directed, is_stimulation,
  is_inhibition)`` that resource states. A resource asserting two contradicting
  signs for the same endpoints under different predicates therefore keeps
  **two** rows, because the signature is part of the key.
* **the collapse** holds one row per ordered triple over the scope it was
  folded for. **Amended 2026-08-21**: it is no longer a table.
  ``interaction_fact_combined`` materialised the all-resources scope and
  is removed, so the collapse assertions here read
  ``tests.fixtures.collapse`` — the fold written once on the test side — through
  a view over the record. What they assert is unchanged, which is the evidence
  that they were always about the collapse and never about its storage. The
  pair of tests that compared the fold against the materialisation went with
  the table: there is nothing left for the fold to disagree with, and the
  equivalence the api-service owes its own fold is asserted there.

Endpoints stay ordered on both grains, so A→B and B→A remain distinct rows.

The class is what makes the key discriminate, and the derive step reads it off
the resource annotations: participant-role evidence first, then
interaction-level annotation, then the predicate, then ``other``. Resolving it from the predicate
alone put 93.5 per cent of the graph in ``other`` and left five classes empty,
so the precedence is asserted here class by class.

**And the grain is asserted class by class too.** The fold is written once and
the class is only a column in it, so a rule shown to hold for the class the
first slice exercised has not been shown to hold for the rest. The fixture
carries a coverage pair per class — two resources, two references, one ordered
key — and the dedup and provenance assertions run over all of them. A class the
graph evidences but the projection never names is the failure these catch, and
it is a real one: a whole resource can arrive under a verb the curated map has
no entry for, and be served as ``other`` with nothing in the per-class counts
to say so.

Run::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_interactions_dedup.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os

import psycopg2
import pytest

from tests.fixtures.collapse import COLLAPSE_VIEW, create_collapse_view
from tests.fixtures.interaction_graph import (
    COVERED_CLASSES,
    REFERENCES,
    SOURCE_A,
    SOURCE_B,
    SOURCE_NAMES,
    build_interaction_fixture,
)

DATABASE_URL = os.environ.get('DATABASE_URL')
SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_INTERACTIONS_DEDUP',
    'interactions_dedup_test',
)

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the projection test needs a Postgres',
)

@pytest.fixture(scope='module')
def built():
    """The fixture graph, projected by the derive step, in a scratch schema."""
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
        # The collapse is a query-time shape, so the tests below read it
        # through a view over the record rather than through a table the derive
        # writes — there is no such table any more.
        create_collapse_view(connection, SCRATCH)
        yield connection, stats
    finally:
        with connection.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {SCRATCH} CASCADE')
        connection.commit()
        connection.close()


def _query(conn, statement: str, params=None) -> list[tuple]:
    """Run a read query, keeping a failure from poisoning the next test.

    The fixture is module-scoped and hands every test the same connection, so
    one statement erroring — over a table this cycle has not built yet, say —
    would leave the transaction aborted and make every later test fail for a
    reason that is not its own. Roll back and re-raise instead.
    """
    try:
        with conn.cursor() as cur:
            cur.execute(statement, params)
            return cur.fetchall()
    except psycopg2.Error:
        conn.rollback()
        raise


def _facts(conn, subject: str, object_: str) -> list[dict[str, object]]:
    """Every collapsed row for a fixture endpoint pair, by name.

    An endpoint pair reaches one collapsed row **per class**, so this returns a
    list: the class is part of the key, and a pair two resources report under
    two different classes is two rows rather than one row with two labels.
    """
    rows = _query(
        conn,
        f"""
        SELECT vic.name, f.sources, f.source_count, f.reference_pubmed_ids,
               f.reference_count, f.interaction_id, f.is_stimulation,
               f.is_inhibition, f.sign_source_count
        FROM {SCRATCH}.{COLLAPSE_VIEW} f
        JOIN {SCRATCH}.entity subject
          ON subject.entity_id = f.subject_entity_id
        JOIN {SCRATCH}.entity object
          ON object.entity_id = f.object_entity_id
        LEFT JOIN {SCRATCH}.vocab_interaction_class vic
          ON vic.interaction_class_id = f.interaction_class_id
        WHERE subject.canonical_identifier = %s
          AND object.canonical_identifier = %s
        """,
        (f'FIXTURE_{subject}', f'FIXTURE_{object_}'),
    )
    return [
        {
            'class': name,
            'sources': sorted(sources or []),
            'source_count': source_count,
            'pubmed': sorted(pubmed or []),
            'reference_count': reference_count,
            'interaction_id': interaction_id,
            'is_stimulation': is_stimulation,
            'is_inhibition': is_inhibition,
            'sign_source_count': sign_source_count,
        }
        for (
            name,
            sources,
            source_count,
            pubmed,
            reference_count,
            interaction_id,
            is_stimulation,
            is_inhibition,
            sign_source_count,
        ) in rows
    ]


def _fact(conn, subject: str, object_: str) -> dict[str, object] | None:
    """The single collapsed row for a fixture endpoint pair, by name.

    Most fixture rows carry one class, and this is the shape they read in.
    """
    facts = _facts(conn, subject, object_)
    assert len(facts) <= 1, f'{subject}->{object_} produced {len(facts)} rows'
    return facts[0] if facts else None


def _records(conn, subject: str, object_: str) -> list[dict[str, object]]:
    """Every record row for a fixture endpoint pair.

    One per contributing resource **and** assertion signature.
    """
    rows = _query(
        conn,
        f"""
        SELECT ds.name, vic.name, r.is_directed, r.is_stimulation,
               r.is_inhibition, r.reference_pubmed_ids,
               r.interaction_fact_resource_id, r.interaction_id
        FROM {SCRATCH}.interaction_fact_resource r
        JOIN {SCRATCH}.entity subject
          ON subject.entity_id = r.subject_entity_id
        JOIN {SCRATCH}.entity object
          ON object.entity_id = r.object_entity_id
        JOIN {SCRATCH}.data_source ds ON ds.source_id = r.source_id
        LEFT JOIN {SCRATCH}.vocab_interaction_class vic
          ON vic.interaction_class_id = r.interaction_class_id
        WHERE subject.canonical_identifier = %s
          AND object.canonical_identifier = %s
        ORDER BY ds.name, r.is_stimulation, r.is_inhibition
        """,
        (f'FIXTURE_{subject}', f'FIXTURE_{object_}'),
    )
    return [
        {
            'source': source,
            'class': class_name,
            'is_directed': is_directed,
            'is_stimulation': is_stimulation,
            'is_inhibition': is_inhibition,
            'pubmed': sorted(pubmed or []),
            'record_id': record_id,
            'interaction_id': interaction_id,
        }
        for (
            source,
            class_name,
            is_directed,
            is_stimulation,
            is_inhibition,
            pubmed,
            record_id,
            interaction_id,
        ) in rows
    ]


# ---------------------------------------------------------------------------
# The record grain — `interaction_fact_resource`
# ---------------------------------------------------------------------------


def test_the_record_key_is_unique(built):
    """One row per endpoint pair, class, resource **and** assertion signature.

    Postgres `GROUP BY` treats NULLs as equal, which is the comparison this key
    needs: the signature columns are nullable, and two silent rows from one
    resource are a duplicate rather than two distinct assertions.
    """
    conn, _stats = built
    rows = _query(
        conn,
        f"""
        SELECT count(*) FROM (
          SELECT 1 FROM {SCRATCH}.interaction_fact_resource
          GROUP BY subject_entity_id, object_entity_id, interaction_class_id,
                   source_id, is_directed, is_stimulation, is_inhibition
          HAVING count(*) > 1
        ) duplicated
        """,
    )
    assert rows[0][0] == 0


def test_the_record_keeps_one_row_per_contributing_resource(built):
    """Three resources on one endpoint pair are three record rows, not one.

    Each carries what *that* resource asserts, so the silent third contributor
    keeps NULL sign columns instead of inheriting a neighbour's assertion.
    """
    conn, _stats = built
    records = _records(conn, 'c', 'd')
    assert [record['source'] for record in records] == [
        'fixture_res_a',
        'fixture_res_b',
        'fixture_res_c',
    ], 'the signalling trio did not keep one record row per resource'
    assert {record['class'] for record in records} == {'signaling'}
    by_source = {record['source']: record for record in records}
    assert by_source['fixture_res_a']['is_stimulation'] is True
    assert by_source['fixture_res_a']['is_inhibition'] is None
    assert by_source['fixture_res_b']['is_inhibition'] is True
    assert by_source['fixture_res_b']['is_stimulation'] is None
    assert by_source['fixture_res_c']['is_stimulation'] is None
    assert by_source['fixture_res_c']['is_inhibition'] is None


def test_contradicting_signs_from_one_resource_keep_two_rows(built):
    """The assertion signature is part of the key, so one resource stating both
    signs under two predicates keeps **two** record rows."""
    conn, _stats = built
    records = _records(conn, 'i', 'j')
    assert len(records) == 2, (
        'one resource asserting both signs must keep two record rows; got '
        f'{len(records)}'
    )
    assert {record['source'] for record in records} == {'fixture_res_a'}
    signatures = {
        (record['is_stimulation'], record['is_inhibition'])
        for record in records
    }
    assert signatures == {(True, None), (None, True)}
    assert len({record['record_id'] for record in records}) == 2, (
        'the two rows must carry distinct surrogate keys'
    )


def test_the_contradiction_collapses_to_one_row_with_both_flags(built):
    """Collapsed over the whole resource set the same pair is one row, with
    both sign flags true and a `sign_source_count` of one."""
    conn, _stats = built
    fact = _fact(conn, 'i', 'j')
    assert fact is not None, 'the single-resource conflict produced no row'
    assert fact['is_stimulation'] is True
    assert fact['is_inhibition'] is True
    assert fact['source_count'] == 1
    assert fact['sign_source_count'] == 1


def test_the_record_keeps_each_resources_references_to_itself(built):
    """A record row carries **this** resource's references, not the pool.

    Pooling them at record grain is the defect the amendment names: the value
    would then describe resources a scoped query excluded.
    """
    conn, _stats = built
    by_source = {
        record['source']: record for record in _records(conn, 'c', 'd')
    }
    assert by_source['fixture_res_a']['pubmed'] == ['11111111']
    assert by_source['fixture_res_b']['pubmed'] == []
    assert by_source['fixture_res_c']['pubmed'] == ['33333333']


def test_the_record_has_a_unique_surrogate_key(built):
    """`interaction_fact_resource_id` is present on every row and unique.

    The composite key includes nullable signature columns and Postgres defaults
    to `MATCH SIMPLE`, under which a foreign key with any NULL column is not
    checked at all — so the detail tables can only anchor on the surrogate.
    """
    conn, _stats = built
    rows = _query(
        conn,
        f"""
        SELECT count(*),
               count(interaction_fact_resource_id),
               count(DISTINCT interaction_fact_resource_id)
        FROM {SCRATCH}.interaction_fact_resource
        """,
    )
    total, present, distinct = rows[0]
    assert total > 0, 'the record table is empty'
    assert present == total, 'some record rows carry no surrogate key'
    assert distinct == total, 'the surrogate key is not unique'


def test_the_surrogate_key_is_deterministic_across_rebuilds(built):
    """A re-run reproduces the same ids, so a child row's key survives a build.

    The id is a `content_uuid` over the full key rather than a serial, which is
    what makes it stable; a rebuild that reshuffled it would turn every detail
    table into a migration.
    """
    conn, _stats = built
    from omnipath_build.db.derived_tables import rebuild_interaction_tables

    def keyed_ids() -> dict[tuple, str]:
        rows = _query(
            conn,
            f"""
            SELECT subject_entity_id, object_entity_id, interaction_class_id,
                   source_id, is_directed, is_stimulation, is_inhibition,
                   interaction_fact_resource_id
            FROM {SCRATCH}.interaction_fact_resource
            """,
        )
        return {row[:-1]: row[-1] for row in rows}

    before = keyed_ids()
    rebuild_interaction_tables(conn, schema=SCRATCH)
    assert keyed_ids() == before, (
        'the record surrogate key changed on a rebuild of the same content'
    )


def test_opposite_directions_stay_two_record_rows(built):
    """The record key is ordered too, so A→B and B→A never merge."""
    conn, _stats = built
    forward = _records(conn, 'g', 'h')
    reverse = _records(conn, 'h', 'g')
    assert len(forward) == 1 and len(reverse) == 1
    assert forward[0]['record_id'] != reverse[0]['record_id']


def test_every_record_row_links_to_its_header(built):
    """The record carries the header FK; the materialisation inherits it."""
    conn, _stats = built
    rows = _query(
        conn,
        f'SELECT count(*) FROM {SCRATCH}.interaction_fact_resource '
        f'WHERE interaction_id IS NULL',
    )
    assert rows[0][0] == 0
    rows = _query(
        conn,
        f"""
        SELECT count(*)
        FROM {SCRATCH}.interaction_fact_resource r
        LEFT JOIN {SCRATCH}.interaction i
          ON i.interaction_id = r.interaction_id
        WHERE i.interaction_id IS NULL
        """,
    )
    assert rows[0][0] == 0


def test_no_record_row_is_left_without_a_class(built):
    """Every record row reaches at least `other`, so the key never goes NULL."""
    conn, _stats = built
    rows = _query(
        conn,
        f'SELECT count(*) FROM {SCRATCH}.interaction_fact_resource '
        f'WHERE interaction_class_id IS NULL',
    )
    assert rows[0][0] == 0


# ---------------------------------------------------------------------------
# The all-resources collapse
# ---------------------------------------------------------------------------


def test_the_carried_columns_come_through_the_collapse_unchanged(built):
    """`subject_organism`, `object_organism` and `interaction_id` are carried
    rather than recomputed, so a record row must state each."""
    conn, _stats = built
    rows = _query(
        conn,
        f"""
        SELECT count(*)
        FROM {SCRATCH}.{COLLAPSE_VIEW} f
        WHERE NOT EXISTS (
          SELECT 1 FROM {SCRATCH}.interaction_fact_resource r
          WHERE r.subject_entity_id = f.subject_entity_id
            AND r.object_entity_id = f.object_entity_id
            AND r.interaction_class_id = f.interaction_class_id
            AND r.subject_organism IS NOT DISTINCT FROM f.subject_organism
            AND r.object_organism IS NOT DISTINCT FROM f.object_organism
            AND r.interaction_id IS NOT DISTINCT FROM f.interaction_id
        )
        """,
    )
    assert rows[0][0] == 0, (
        'a collapsed row carries a carried-through value no record row states'
    )


def test_the_collapse_invents_no_summary_value(built):
    """Values the collapse summarises without a fixed rule — affinity, pchembl,
    score, curation flags — must still come from a contributing record row."""
    conn, _stats = built
    for column in ('affinity', 'pchembl', 'score'):
        rows = _query(
            conn,
            f"""
            SELECT count(*)
            FROM {SCRATCH}.{COLLAPSE_VIEW} f
            WHERE f.{column} IS NOT NULL
              AND NOT EXISTS (
                SELECT 1 FROM {SCRATCH}.interaction_fact_resource r
                WHERE r.subject_entity_id = f.subject_entity_id
                  AND r.object_entity_id = f.object_entity_id
                  AND r.interaction_class_id = f.interaction_class_id
                  AND r.{column} = f.{column}
              )
            """,
        )
        assert rows[0][0] == 0, f'{column} on a collapsed row has no record row'
    rows = _query(
        conn,
        f"""
        SELECT count(*)
        FROM {SCRATCH}.{COLLAPSE_VIEW} f
        CROSS JOIN LATERAL unnest(coalesce(f.curation_flags, '{{}}'::text[]))
          AS c(flag)
        WHERE NOT EXISTS (
          SELECT 1 FROM {SCRATCH}.interaction_fact_resource r
          WHERE r.subject_entity_id = f.subject_entity_id
            AND r.object_entity_id = f.object_entity_id
            AND r.interaction_class_id = f.interaction_class_id
            AND c.flag = ANY (r.curation_flags)
        )
        """,
    )
    assert rows[0][0] == 0, 'a curation flag on a collapsed row has no record'


def test_the_ordered_key_is_unique(built):
    """The collapsed output holds one row per ordered endpoint pair and class."""
    conn, _stats = built
    rows = _query(
        conn,
        f"""
        SELECT count(*) FROM (
          SELECT 1 FROM {SCRATCH}.{COLLAPSE_VIEW}
          GROUP BY subject_entity_id, object_entity_id,
                   interaction_class_id
          HAVING count(*) > 1
        ) duplicated
        """,
    )
    assert rows[0][0] == 0


def test_three_resources_fold_into_one_row(built):
    """One collapsed row per endpoint pair and class, with all provenance."""
    conn, _stats = built
    fact = _fact(conn, 'c', 'd')
    assert fact is not None, 'the signalling pair produced no fact row'
    assert fact['class'] == 'signaling'
    assert fact['sources'] == [
        'fixture_res_a',
        'fixture_res_b',
        'fixture_res_c',
    ]
    assert fact['source_count'] == 3


def test_references_aggregate_across_the_contributors(built):
    """Every contributing resource's references reach the collapsed row."""
    conn, _stats = built
    fact = _fact(conn, 'c', 'd')
    assert fact['pubmed'] == ['11111111', '33333333']
    assert fact['reference_count'] == 2


def test_opposite_directions_stay_two_rows(built):
    """The key is ordered, so A→B and B→A are never merged."""
    conn, _stats = built
    forward = _fact(conn, 'g', 'h')
    reverse = _fact(conn, 'h', 'g')
    assert forward is not None and reverse is not None
    assert forward['class'] == reverse['class'] == 'signaling'


def test_every_fact_row_links_to_its_header(built):
    """The projection keeps the header link the general API reads."""
    conn, _stats = built
    rows = _query(
        conn,
        f'SELECT count(*) FROM {SCRATCH}.{COLLAPSE_VIEW} '
        f'WHERE interaction_id IS NULL',
    )
    assert rows[0][0] == 0
    rows = _query(
        conn,
        f"""
        SELECT count(*)
        FROM {SCRATCH}.{COLLAPSE_VIEW} f
        LEFT JOIN {SCRATCH}.interaction i
          ON i.interaction_id = f.interaction_id
        WHERE i.interaction_id IS NULL
        """,
    )
    assert rows[0][0] == 0


def test_every_header_has_its_participants(built):
    """Each header's ``interaction_party`` rows match its stated arity."""
    conn, _stats = built
    rows = _query(
        conn,
        f"""
        SELECT count(*)
        FROM {SCRATCH}.interaction i
        LEFT JOIN (
          SELECT interaction_id, count(*) AS parties
          FROM {SCRATCH}.interaction_party
          GROUP BY interaction_id
        ) p ON p.interaction_id = i.interaction_id
        WHERE p.parties IS DISTINCT FROM i.arity
        """,
    )
    assert rows[0][0] == 0


@pytest.mark.parametrize(
    ('subject', 'object_', 'expected'),
    [
        # Tier 1, participant-role evidence: a ligand on one side and a
        # receptor on the other, under the generic `interacts_with` verb the
        # predicate vocabulary cannot classify.
        ('a', 'b', 'ligand_receptor'),
        # Tier 2, interaction-level annotation, again under `interacts_with`.
        ('k', 'l', 'orthosteric'),
        ('m', 'n', 'allosteric'),
        # Tier 3, the predicate.
        ('c', 'd', 'signaling'),
        ('o', 'p', 'transport'),
        ('maturation_s', 'maturation_o', 'maturation'),
        # The fallback, which stays a real class rather than a dumping ground.
        ('q', 'r', 'other'),
    ],
)
def test_the_class_follows_the_r18_precedence(built, subject, object_, expected):
    """The class comes from what the resource asserts, not from the verb alone."""
    conn, _stats = built
    fact = _fact(conn, subject, object_)
    assert fact is not None, f'{subject}->{object_} produced no fact row'
    assert fact['class'] == expected


def test_no_fact_row_is_left_without_a_class(built):
    """Every row reaches at least `other`, so the class key never goes NULL."""
    conn, _stats = built
    rows = _query(
        conn,
        f'SELECT count(*) FROM {SCRATCH}.{COLLAPSE_VIEW} '
        f'WHERE interaction_class_id IS NULL',
    )
    assert rows[0][0] == 0


def test_the_step_reports_its_per_class_counts(built):
    """Deriving the class from annotations is fragile, so every run counts it.

    Every class in the vocabulary is named, including the ones at zero: a count
    that is absent reads as "not measured", and the whole point of the report is
    that a class collapsing back to nothing is visible in the build output.
    """
    conn, _stats = built
    _connection, stats = built
    vocabulary = {
        name
        for (name,) in _query(
            conn,
            f'SELECT name FROM {SCRATCH}.vocab_interaction_class',
        )
    }
    assert set(stats.rows_by_class) == vocabulary
    assert stats.rows_by_class['signaling'] >= 3


@pytest.mark.parametrize('class_name', COVERED_CLASSES)
def test_every_class_the_graph_evidences_reaches_the_projection(
    built,
    class_name,
):
    """A class the graph evidences is a class the projection produces rows for.

    The projection is not the ligand-receptor slice with a fallback beside it.
    Every class some resource states, by a participant role, by an
    interaction-level annotation or by the verb it publishes under, has to
    arrive in the fact table under its own name — otherwise it is served as
    `other`, which is the one label that cannot be filtered on usefully.
    """
    _conn, stats = built
    assert stats.rows_by_class[class_name] > 0, (
        f'the graph evidences {class_name} and the projection produced no row '
        f'for it; the rows are in the fallback class instead'
    )


def test_a_class_no_resource_evidences_is_reported_at_zero(built):
    """A class with no evidence stays empty and stays visible.

    Transcription-factor targets are a missing-resource gap rather than a
    derivation to force: nothing in the graph asserts one. The report names the
    class at zero, which is how a gap is surfaced instead of being filled with
    guesses.
    """
    _conn, stats = built
    assert stats.rows_by_class['tf_target'] == 0


# ---------------------------------------------------------------------------
# The dedup and aggregation rule, class by class
# ---------------------------------------------------------------------------


def _coverage_pair(class_name: str) -> tuple[str, str]:
    return f'{class_name}_s', f'{class_name}_o'


def _coverage_references(class_name: str) -> list[str]:
    key = f'pair_{class_name}'
    return sorted(
        reference
        for (relation, _source), reference in REFERENCES.items()
        if relation == key
    )


@pytest.mark.parametrize('class_name', COVERED_CLASSES)
def test_the_ordered_triple_folds_to_one_row_in_every_class(built, class_name):
    """Two resources on one ordered pair and class fold to a single row.

    Asserted per class rather than on `signaling` alone: the fold is written
    once and the classes are only a column in it, so a rule that held for the
    class the first slice exercised has to be shown holding for the rest.
    """
    conn, _stats = built
    subject, object_ = _coverage_pair(class_name)
    facts = _facts(conn, subject, object_)
    assert len(facts) == 1, (
        f'{class_name} folded to {len(facts)} rows, not one'
    )
    assert facts[0]['class'] == class_name


@pytest.mark.parametrize('class_name', COVERED_CLASSES)
def test_provenance_aggregates_across_resources_in_every_class(
    built,
    class_name,
):
    """The folded row names both contributors and both their references."""
    conn, _stats = built
    subject, object_ = _coverage_pair(class_name)
    fact = _fact(conn, subject, object_)
    assert fact is not None, f'{class_name} produced no folded row'
    assert fact['sources'] == sorted(
        (SOURCE_NAMES[SOURCE_A], SOURCE_NAMES[SOURCE_B])
    )
    assert fact['source_count'] == 2
    assert fact['pubmed'] == _coverage_references(class_name)
    assert fact['reference_count'] == 2


@pytest.mark.parametrize('class_name', COVERED_CLASSES)
def test_the_record_keeps_one_row_per_resource_in_every_class(
    built,
    class_name,
):
    """Below the fold the grain is per resource, whatever the class.

    The stored table keeps the contributors apart so a query scoped to one of
    them can recompute the summary from its own rows. That is the property the
    fold depends on, and it has to hold in every class rather than in the ones
    the first slice happened to cover.
    """
    conn, _stats = built
    subject, object_ = _coverage_pair(class_name)
    records = _records(conn, subject, object_)
    assert [record['source'] for record in records] == [
        SOURCE_NAMES[SOURCE_A],
        SOURCE_NAMES[SOURCE_B],
    ]
    assert {record['class'] for record in records} == {class_name}
    by_source = {record['source']: record for record in records}
    assert by_source[SOURCE_NAMES[SOURCE_A]]['pubmed'] == [
        REFERENCES[(f'pair_{class_name}', SOURCE_A)]
    ]
    assert by_source[SOURCE_NAMES[SOURCE_B]]['pubmed'] == [
        REFERENCES[(f'pair_{class_name}', SOURCE_B)]
    ]


def test_the_class_is_part_of_the_collapse_key(built):
    """One endpoint pair reported under two classes stays two folded rows.

    Merging them would answer "what kind of interaction is this?" with a list,
    and would pool the provenance of a signalling claim with that of a
    pharmacological one. Each row carries the resource that made its own claim.
    """
    conn, _stats = built
    facts = _facts(conn, 'dual_s', 'dual_o')
    assert len(facts) == 2, (
        f'one pair under two classes folded to {len(facts)} rows'
    )
    by_class = {fact['class']: fact for fact in facts}
    assert set(by_class) == {'signaling', 'orthosteric'}
    assert by_class['signaling']['sources'] == [SOURCE_NAMES[SOURCE_A]]
    assert by_class['orthosteric']['sources'] == [SOURCE_NAMES[SOURCE_B]]
    assert by_class['signaling']['source_count'] == 1
    assert by_class['orthosteric']['source_count'] == 1


def test_the_ordered_key_holds_outside_the_signalling_class(built):
    """A→B and B→A stay two rows in a class the first slice never exercised."""
    conn, _stats = built
    forward = _fact(conn, 'maturation_s', 'maturation_o')
    reverse = _fact(conn, 'maturation_o', 'maturation_s')
    assert forward is not None and reverse is not None
    assert forward['class'] == reverse['class'] == 'maturation'
    assert forward['sources'] != reverse['sources']


def test_the_step_reports_what_the_fallback_class_is_made_of(built):
    """The fallback is broken down by the verb its rows arrived under.

    `other` is a real class and also the place every unrecognised verb lands,
    so a per-class count alone cannot tell the two apart: a resource publishing
    a whole class under a predicate no rule maps shows up as more of the same
    large number. Reporting the fallback per predicate is what makes that
    visible, and it is the only report that would have caught it.
    """
    _conn, stats = built
    fallback = stats.fallback_predicates
    assert fallback, 'the step reported no breakdown of the fallback class'
    assert 'has_member' in fallback, (
        'the verb the fallback rows arrived under is not named'
    )
    assert fallback['has_member'] > 0
    assert 'transports' not in fallback, (
        'a verb that reaches a class of its own is not part of the fallback'
    )


# ---------------------------------------------------------------------------
# The `source_count` histogram
# ---------------------------------------------------------------------------


def test_the_histogram_counts_every_collapse_key_once(built):
    """`interaction_source_count_histogram` partitions the collapse keys.

    One row per **observed** level, the count of keys at each, and nothing
    outside: the levels sum to the number of collapse keys, so a key is counted
    once and no key is missed. Written to the database rather than only
    returned, because the consumer is the api-service's guardrail and
    it reads the build's output, not the build's process.

    Nine levels on dev4 at present cardinality, which is a property of the data
    and not of the table — this fixture holds three resources, so it shows the
    levels it has.
    """
    conn, _stats = built
    histogram = _query(
        conn,
        f'SELECT source_count, keys FROM '
        f'{SCRATCH}.interaction_source_count_histogram ORDER BY source_count',
    )
    assert histogram, 'the derive recorded no source_count histogram'
    keys = _query(
        conn,
        f"""
        SELECT count(*) FROM (
          SELECT 1 FROM {SCRATCH}.interaction_fact_resource
          GROUP BY subject_entity_id, object_entity_id, interaction_class_id
        ) folded
        """,
    )[0][0]
    assert sum(count for _level, count in histogram) == keys
    assert all(level >= 1 for level, _count in histogram), (
        'a collapse key with no contributing resource is not a collapse key'
    )


def test_the_histogram_counts_resources_and_not_record_rows(built):
    """The level is `count(DISTINCT source_id)`, never `count(*)`.

    One resource asserting two contradicting signs keeps **two** record rows
    under the per-resource grain, and it is still one resource contributing. A
    histogram built on `count(*)` would price `source_count >= 2` from a
    population that includes single-resource keys, which is the guardrail
    estimating a filter from the wrong distribution.
    """
    conn, _stats = built
    # i->j is the single-resource conflict: two record rows, one resource.
    rows = _query(
        conn,
        f"""
        SELECT count(*), count(DISTINCT r.source_id)
        FROM {SCRATCH}.interaction_fact_resource r
        JOIN {SCRATCH}.entity subject
          ON subject.entity_id = r.subject_entity_id
        JOIN {SCRATCH}.entity object
          ON object.entity_id = r.object_entity_id
        WHERE subject.canonical_identifier = 'FIXTURE_i'
          AND object.canonical_identifier = 'FIXTURE_j'
        """,
    )
    record_rows, resources = rows[0]
    assert record_rows == 2 and resources == 1, (
        'the single-resource conflict fixture no longer has two record rows '
        'from one resource, so this test proves nothing'
    )
    level = _query(
        conn,
        f"""
        SELECT count(DISTINCT r.source_id)
        FROM {SCRATCH}.interaction_fact_resource r
        JOIN {SCRATCH}.entity subject
          ON subject.entity_id = r.subject_entity_id
        JOIN {SCRATCH}.entity object
          ON object.entity_id = r.object_entity_id
        WHERE subject.canonical_identifier = 'FIXTURE_i'
          AND object.canonical_identifier = 'FIXTURE_j'
        GROUP BY r.subject_entity_id, r.object_entity_id,
                 r.interaction_class_id
        """,
    )[0][0]
    assert level == 1
    at_one = _query(
        conn,
        f'SELECT keys FROM {SCRATCH}.interaction_source_count_histogram '
        f'WHERE source_count = 1',
    )
    assert at_one and at_one[0][0] >= 1, (
        'the one-resource level is absent from the histogram, so the two-row '
        'key was counted as two resources'
    )
