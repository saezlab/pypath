"""Build the MetSigDB membership table as Parquet from a published release.

The step reads the pinned resources' Parquet through DuckDB: it parses no upstream
resource and downloads nothing. Each resource runs the same phases. An extraction query
stages one row per published pair, a shared projection reads the metabolite and set sides
off the release's entities, and a shared insert writes the rows. Only the extraction
differs per resource.
"""

from __future__ import annotations

import logging
from pathlib import Path

import duckdb

from omnipath_core.subset_layout import BuildContext
from omnipath_subsets.duckdb_script import run_script
from omnipath_subsets.metsigdb.mapping import KEGG_OVERVIEW_MAP_NUMBERS, RESOURCES, ResourceRule
from omnipath_subsets.release_source import (
    CHEMICAL_ENTITY_TYPES,
    open_release,
    pinned_resources,
    resource_code_revision,
)

_SQL_DIR = Path(__file__).with_name("duckdb_sql")

logger = logging.getLogger(__name__)

TABLE = "metsigdb_membership"

# The extraction file of each rule. Reactome and WikiPathways share the one-hop shape;
# MACdb is the same shape plus the trait type of its set.
EXTRACTIONS = {
    "extract_onehop.sql": "extract_onehop",
    "extract_macdb.sql": "extract_onehop",
    "extract_kegg.sql": "extract_kegg",
    "extract_classyfire.sql": "extract_classyfire",
}


def _sql(name: str) -> str:
    return (_SQL_DIR / f"{name}.sql").read_text(encoding="utf-8")


def _provenance_source(data_root, resources: dict[str, str], source: str) -> str:
    """The input module and commit that parsed the resource, as its build recorded them."""
    module = f"omnipath_resources.inputs_v2.{source}"
    revision = (
        resource_code_revision(data_root, source, resources[source])
        if source in resources
        else None
    )
    return f"{module}@{revision[:12]}" if revision else module


def load_resource(
    con: duckdb.DuckDBPyConnection,
    rule: ResourceRule,
    *,
    build_id: str,
    provenance_source: str,
) -> int:
    """Extract, project and publish one resource into the membership table."""
    parameters = {
        "source": rule.source_name,
        "hierarchy_source": rule.hierarchy_source_name,
        "set_entity_type": rule.set_entity_type,
        "chemical_entity_types": list(CHEMICAL_ENTITY_TYPES),
        "is_pathway": rule.set_type == "pathway",
        "trait_types": rule.name == "MACdb",
        "overview_map_numbers": list(KEGG_OVERVIEW_MAP_NUMBERS),
        "organism_free_pattern": rule.organism_free_pattern,
        "set_id_prefix": rule.set_id_prefix,
        "resource": rule.name,
        "set_type": rule.set_type,
        "provenance_source": provenance_source,
        "build_id": build_id,
    }
    run_script(con, _sql(EXTRACTIONS[rule.extraction]), parameters)
    run_script(con, _sql("project"), parameters)
    run_script(con, _sql("publish"), parameters)
    for name in ("metsigdb_stage", "metsigdb_projection", "metsigdb_set_attrs"):
        con.execute(f"DROP TABLE IF EXISTS {name}")
    return con.execute(f"SELECT count(*) FROM {TABLE} WHERE resource = ?", [rule.name]).fetchone()[
        0
    ]


def _exclusion_members(con: duckdb.DuckDBPyConnection, exclusion_table: Path | None) -> None:
    """The members of the exclusion product, as the view the rows read their classes from."""
    if exclusion_table is None:
        con.execute(
            "CREATE OR REPLACE TEMP VIEW exclusion_member AS "
            "SELECT NULL::VARCHAR AS entity_key, NULL::VARCHAR AS exclusion_class WHERE false"
        )
        return
    literal = str(exclusion_table).replace("'", "''")
    con.execute(
        "CREATE OR REPLACE TEMP VIEW exclusion_member AS "
        f"SELECT entity_key, exclusion_class FROM read_parquet('{literal}')"
    )


def build_table(
    con: duckdb.DuckDBPyConnection,
    *,
    data_root: Path,
    release: str,
    build_id: str,
    exclusion_table: Path | None = None,
    resources: tuple[ResourceRule, ...] = RESOURCES,
) -> dict:
    """Fill ``metsigdb_membership`` in ``con``; returns the rows each resource published."""
    pinned = pinned_resources(data_root, release)
    run_script(con, _sql("tables"))
    _exclusion_members(con, exclusion_table)
    counts = {}
    for rule in resources:
        if rule.source_name not in pinned:
            logger.warning("metsigdb: %s is not pinned by release %s", rule.source_name, release)
            counts[rule.name] = 0
            continue
        counts[rule.name] = load_resource(
            con,
            rule,
            build_id=build_id,
            provenance_source=_provenance_source(data_root, pinned, rule.source_name),
        )
    return counts


COVERAGE_IDENTIFIERS = ("hmdb", "chebi", "kegg", "pubchem", "inchikey")


def coverage(con: duckdb.DuckDBPyConnection) -> dict[str, dict]:
    """What each resource published, and the share of its metabolites with each identifier.

    Measured on the table just built, so a consumer reads the coverage of the release it is
    serving and not a figure from an earlier build: an identifier query returns nothing for
    a resource whose metabolites lack that identifier.
    """
    shares = ", ".join(
        f"count(DISTINCT metabolite_entity_key) FILTER (WHERE {name} IS NOT NULL) AS {name}"
        for name in COVERAGE_IDENTIFIERS
    )
    rows = con.execute(
        f"SELECT resource, count(*) AS n_rows, count(DISTINCT set_source_id) AS sets, "
        f"count(DISTINCT metabolite_entity_key) AS metabolites, {shares} "
        f"FROM {TABLE} GROUP BY resource ORDER BY resource"
    ).fetchall()
    found = {}
    for resource, n_rows, sets, metabolites, *counts in rows:
        found[resource] = {
            "rows": n_rows,
            "sets": sets,
            "metabolites": metabolites,
            "identifiers": {
                name: round(count / metabolites, 4) if metabolites else 0.0
                for name, count in zip(COVERAGE_IDENTIFIERS, counts)
            },
        }
    return found


def write_metsigdb(directory: Path, context: BuildContext, *, seed_directory=None) -> dict:
    """Write the membership table of ``context.release`` into ``directory``."""
    con = open_release(context.data_root, context.release)
    try:
        counts = build_table(
            con,
            data_root=context.data_root,
            release=context.release,
            build_id=context.release_manifest_sha256,
            exclusion_table=context.dependencies["exclusion_metabolites"]
            / "exclusion_metabolite.parquet",
        )
        target = str(directory / f"{TABLE}.parquet").replace("'", "''")
        measured = coverage(con)
        con.execute(
            f"COPY (SELECT * FROM {TABLE} ORDER BY resource, set_source_id, "
            f"metabolite_entity_key) TO '{target}' (FORMAT parquet, COMPRESSION zstd)"
        )
    finally:
        con.close()
    logger.info("metsigdb: %s", counts)
    return {"rows": sum(counts.values()), "resources": counts, "coverage": measured}
