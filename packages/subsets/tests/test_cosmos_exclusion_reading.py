"""COSMOS reads the exclusion table; it is not changed by it."""

from omnipath_subsets.cosmos import build, translate
from omnipath_subsets.exclusion_metabolites.build import _sql_text as exclusion_sql
import pytest
from test_cosmos_gene_reference import product_input  # noqa: F401  (fixture)

COLUMNS = (
    "source_label, target_label, source_entity_id::text, target_entity_id::text, "
    "source_type, target_type, source_compartment, target_compartment, mor, "
    "interaction_type, reaction_entity_id::text, interaction_id::text, reaction_index, "
    "orphan, reverse, direction, source_id_type, target_id_type"
)


def project(cur):
    """Rebuild the edge table whole, as the product does, and return its rows."""
    # The fixture predates this table; label staging reads it and finds nothing.
    cur.execute(
        "CREATE TABLE IF NOT EXISTS entity_reference_context("
        "entity_id uuid, reference_entity_key text)"
    )
    cur.execute(build._sql_text("cosmos_edge_table.sql"))
    cur.execute("TRUNCATE cosmos_edge RESTART IDENTITY")
    translate.stage_label_translation(cur)
    cur.execute(
        build._sql_text("project_edges.sql"),
        {
            "build_id": "fixture",
            "reaction_entity_types": list(build.REACTION_ENTITY_TYPES),
            "reaction_sources": list(build.REACTION_SOURCES),
            "direction_term": build.DIRECTION_TERM,
            "reaction_type_term": build.REACTION_TYPE_TERM,
            "transport_types": build.TRANSPORT_TYPES,
        },
    )
    cur.execute(build._sql_text("project_connectors.sql"), {"build_id": "fixture"})
    cur.execute(f"SELECT {COLUMNS} FROM cosmos_edge ORDER BY cosmos_edge_id")
    rows = cur.fetchall()
    translate.drop_label_translation(cur)
    return rows


def publish_exclusion_rows(cur, ids, classes):
    """What a finished exclusion product leaves in the release schema."""
    cur.execute("ALTER TABLE entity ADD COLUMN IF NOT EXISTS label text")
    cur.execute("ALTER TABLE identifier_evidence ADD COLUMN IF NOT EXISTS value_normalized text")
    cur.execute(exclusion_sql("exclusion_tables.sql"))
    cur.executemany(
        "INSERT INTO exclusion_metabolite(entity_id, exclusion_class, assignment_source, origin,"
        " match_kind, evidence, list_version, build_id)"
        " VALUES (%s, %s, 'project_seed', 'curated', 'entity', 'fixture', 'fixture', 'fixture')",
        [(ids[name], klass) for name, klass in classes],
    )


@pytest.mark.integration
def test_cosmos_publishes_the_same_edges_with_or_without_the_exclusion_table(product_input):  # noqa: F811
    cur, ids = product_input
    before = project(cur)
    assert before, "the fixture must produce edges for the comparison to mean anything"

    publish_exclusion_rows(cur, ids, [("input", "small_molecules"), ("output", "cofactors")])
    after = project(cur)

    assert after == before


@pytest.mark.integration
def test_the_product_does_not_read_the_exclusion_table():
    # The cosmos SQL is a closed set of files; none of them names the table.
    names = ("cosmos_edge_table.sql", "project_edges.sql", "project_connectors.sql")
    for name in names:
        assert "exclusion_metabolite" not in build._sql_text(name), name


@pytest.mark.integration
def test_the_contract_join_counts_participants_per_class(product_input):  # noqa: F811
    cur, ids = product_input
    project(cur)
    publish_exclusion_rows(
        cur, ids, [("input", "small_molecules"), ("input", "cofactors"), ("output", "ions")]
    )
    cur.execute(
        """
        SELECT x.exclusion_class, count(*) AS edge_rows
        FROM cosmos_edge p
        JOIN exclusion_metabolite x
          ON x.entity_id IN (p.source_entity_id, p.target_entity_id)
        GROUP BY 1 ORDER BY 2 DESC, 1
        """
    )
    counts = dict(cur.fetchall())
    # One reactant and one product, each in two directions, over three catalysts.
    assert counts == {"cofactors": 6, "small_molecules": 6, "ions": 6}
