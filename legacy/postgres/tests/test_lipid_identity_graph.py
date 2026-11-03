"""spec 011 WP7/US7 T118/T119/T124 -- the lipid generalization graph.

Runs against a built main DB with the lipid identity graph populated
(``_populate_lipid_identity_graph``, called from ``rebuild_derived_tables``).
Skipped when DATABASE_URL/OMNIPATH_UTILS_PG_URL (utils DB, for the
population step) is unset.
"""

from __future__ import annotations

import os

import pytest

DB_URL = os.environ.get('OMNIPATH_DB_URL') or os.environ.get('DATABASE_URL')
UTILS_URL = os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DB_URL or not UTILS_URL,
    reason='need DATABASE_URL and OMNIPATH_BUILD_UTILS_PG_URL to populate + verify',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DB_URL)
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture(scope='module')
def graph(conn):
    from omnipath_build.db.derived_tables import _populate_lipid_identity_graph

    with conn.cursor() as cur:
        nodes, edges = _populate_lipid_identity_graph(cur, SCHEMA, UTILS_URL)
    conn.commit()
    return nodes, edges


def _scalar(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchone()[0]


def test_the_graph_is_populated(conn, graph):
    nodes, edges = graph
    assert nodes > 0
    assert edges > 0
    assert (
        _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.lipid_name_node') == nodes
    )
    assert (
        _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.lipid_name_edge') == edges
    )


def test_no_self_edges(conn, graph):
    """T119: no self-edge -- also a DB CHECK constraint, this confirms the
    population logic itself never even attempts one."""
    n = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.lipid_name_edge
        WHERE child_name = parent_name AND child_level = parent_level
          AND child_chains_listed = parent_chains_listed
          AND child_chains_possible = parent_chains_possible
        """,
    )
    assert n == 0


def test_no_parent_exceeds_the_hub_threshold(conn, graph):
    """T119: the check that would have caught the SwissLipids placeholder
    merge (one node absorbing 184,509 identifiers) from the hierarchy side."""
    from omnipath_build.db.derived_tables import LIPID_HUB_THRESHOLD

    worst = _scalar(
        conn,
        f"""
        SELECT max(children) FROM (
          SELECT count(*) AS children
          FROM {SCHEMA}.lipid_name_edge
          GROUP BY parent_name, parent_level, parent_chains_listed,
                   parent_chains_possible
        ) x
        """,
    )
    assert worst is not None
    assert worst <= LIPID_HUB_THRESHOLD


def test_every_non_species_node_reaches_a_species_edge(conn, graph):
    """T118: every node above species level has a generated is_a edge to its
    species node."""
    orphaned = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.lipid_name_node n
        WHERE n.lipid_level <> 'species'
          AND n.chains_listed > 0
          AND n.lipid_class IS NOT NULL
          AND n.total_carbon IS NOT NULL
          AND n.total_db IS NOT NULL
          AND NOT EXISTS (
            SELECT 1 FROM {SCHEMA}.lipid_name_edge e
            WHERE e.child_name = n.lipid_name AND e.child_level = n.lipid_level
              AND e.child_chains_listed = n.chains_listed
              AND e.child_chains_possible = n.chains_possible
          )
        """,
    )
    assert orphaned == 0


def test_partially_specified_names_keep_their_own_identity(conn, graph):
    """T124: the partial specifications the old species-downgrade guard used
    to discard (research R9: 3,468 distinct) still exist as their own
    nodes, distinct from any species node."""
    n = _scalar(
        conn,
        f"SELECT count(*) FROM {SCHEMA}.lipid_name_node "
        f"WHERE lipid_level = 'partially_specified'",
    )
    # Not an exact match to R9's 3,468 -- this build's live corpus, not the
    # exact one research.md measured -- but must be the same order of
    # magnitude, not zero (the collapse-to-species bug this replaces would
    # produce zero here).
    assert n > 3000


def test_no_two_chemically_different_chain_assignments_share_a_node(conn, graph):
    """T124: the independent test spec.md's own acceptance criterion names
    -- no standardized name identifies two lipids differing in chain
    assignment. The full (lipid_name, lipid_level, chains_listed,
    chains_possible) key is the primary key already (enforced by Postgres);
    this confirms lipid_name ALONE is genuinely not safe as an identity (the
    reason the composite key exists at all), so a caller must never
    silently drop the other three columns."""
    ambiguous = _scalar(
        conn,
        f"""
        SELECT count(*) FROM (
          SELECT lipid_name FROM {SCHEMA}.lipid_name_node
          GROUP BY lipid_name
          HAVING count(DISTINCT (lipid_level, chains_listed, chains_possible)) > 1
        ) x
        """,
    )
    assert ambiguous > 0, (
        'expected real lipid_name collisions across level/chains -- if this '
        'is ever 0, the composite-key design may no longer be necessary, '
        'worth re-checking rather than assuming'
    )
