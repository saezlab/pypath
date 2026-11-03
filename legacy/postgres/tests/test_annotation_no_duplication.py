"""An entity's attributes are stored once, whatever its degree.

That a ligand is secreted is a fact about the ligand. It is the same fact in
every interaction the ligand takes part in, and a hub entity takes part in
thousands. Storing it per interaction multiplies one fact by the degree of the
node, and the multiplier is not small: this build's most connected entity has
five figures of interactions behind it.

The rule matters more once reactions are native. A reaction is a hyperedge with
many participants, so a per-interaction annotation store would multiply each
participant's attributes by the size of the reaction *and* by the number of
reactions it appears in. That is the arrangement this test forbids, and it
forbids it structurally rather than by counting rows: the annotation store's
key must not contain an interaction.

The row-count invariant is asserted alongside it, because a table can carry the
right key and still be filled per interaction if something upstream inserts one
row per pair.

Run against a build database, e.g. on dev4::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_annotation_no_duplication.py -v
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the no-duplication rule needs a built database',
)

# The per-entity annotation store the node projection reads.
STORE = 'entity_evidence_annotation'

# Column names that would tie an entity's attributes to one interaction.
INTERACTION_KEYS = (
    'interaction_id',
    'interaction_fact_resource_id',
    'subject_entity_id',
    'object_entity_id',
)

# How many of the most connected entities the degree check reads. Enough for
# the degree to vary by orders of magnitude, small enough to stay cheap.
SAMPLE = 25


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
def busy_entities(conn):
    """The most connected annotated entities, with degree and annotation count.

    Args:
        conn: An open connection to the built database.

    Returns:
        `[(entity_id, degree, annotations), …]`, busiest first.
    """
    with conn.cursor() as cur:
        cur.execute(
            f"""
            WITH degree AS (
              SELECT entity_id, count(*) AS interactions FROM (
                SELECT subject_entity_id AS entity_id
                FROM {SCHEMA}.interaction_fact_resource
                UNION ALL
                SELECT object_entity_id FROM {SCHEMA}.interaction_fact_resource
              ) ends
              GROUP BY entity_id
              ORDER BY 2 DESC
              LIMIT %s
            )
            SELECT d.entity_id, d.interactions,
                   count(link.annotation_key) AS annotations
            FROM degree d
            LEFT JOIN {SCHEMA}.entity_evidence_resolution res
              ON res.entity_id = d.entity_id
            LEFT JOIN {SCHEMA}.{STORE} link
              ON link.entity_evidence_id = res.entity_evidence_id
             AND link.source_id = res.source_id
            GROUP BY 1, 2
            ORDER BY 2 DESC
            """,
            [SAMPLE],
        )
        return cur.fetchall()


def test_the_annotation_store_has_no_interaction_in_its_key(conn):
    """Structurally per entity: no column of the store names an interaction."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT column_name FROM information_schema.columns
            WHERE table_schema = %s AND table_name = %s
            """,
            [SCHEMA, STORE],
        )
        columns = {name for (name,) in cur.fetchall()}
    assert columns, f'{SCHEMA}.{STORE} does not exist'
    offending = sorted(columns & set(INTERACTION_KEYS))
    assert not offending, (
        f'{STORE} carries {offending}; an entity attribute keyed by an '
        f'interaction is stored once per interaction, which for a hyperedge '
        f'means once per participant per reaction'
    )


def test_the_store_holds_one_row_per_entity_evidence_and_term(conn):
    """The grain is the key, and the key is not the pair."""
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT count(*) AS rows,
                   count(DISTINCT (source_id, entity_evidence_id, annotation_key))
                     AS distinct_keys
            FROM {SCHEMA}.{STORE}
            """
        )
        rows, distinct_keys = cur.fetchone()
    assert rows, f'{STORE} is empty; the invariant would hold vacuously'
    assert rows == distinct_keys, (
        f'{STORE} holds {rows} rows over {distinct_keys} distinct keys. The '
        f'duplicates are the multiplication this rule forbids, arriving through '
        f'an insert rather than through a column'
    )


def test_annotation_count_does_not_track_degree(busy_entities):
    """A hub carries no more attribute rows than a leaf with the same terms.

    The busiest entities of the build differ in degree by orders of magnitude.
    If their attributes were stored per interaction, the annotation count would
    follow that spread. It does not, because a per-entity store cannot.
    """
    annotated = [row for row in busy_entities if row[2]]
    assert annotated, (
        'none of the most connected entities carries an annotation, so the '
        'comparison has nothing to compare'
    )
    degrees = [row[1] for row in annotated]
    counts = [row[2] for row in annotated]
    assert max(degrees) > min(degrees), (
        'the sampled entities all have the same degree; widen the sample '
        'before trusting this test'
    )
    for entity_id, degree, annotations in annotated:
        assert annotations <= degree or degree < SAMPLE, (
            f'entity {entity_id} has degree {degree} and {annotations} '
            f'annotation rows'
        )
    assert max(counts) < max(degrees), (
        f'the annotation counts {sorted(counts)[-3:]} reach the degrees '
        f'{sorted(degrees)[-3:]}; the attributes are being stored per '
        f'interaction after all'
    )
