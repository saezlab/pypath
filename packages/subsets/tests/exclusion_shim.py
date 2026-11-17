"""Build the exclusion product the way the pipeline does, and read it as PostgreSQL tables.

The product tests describe a small release as rows of PostgreSQL tables and assert on the
rows the product publishes. The product is built from Parquet: this writes the release as
Parquet, builds the subset, and copies its tables into the same schema, so the tests read
what a loaded release holds.
"""

from pathlib import Path
import tempfile

import duckdb
from psycopg2.extras import Json, execute_values

from omnipath_postgres.subset_loading import create_tables
from omnipath_subsets.exclusion_metabolites import parquet_build
from omnipath_subsets.parquet_products import build_release_subsets
from omnipath_subsets.release_source import open_release
from parquet_release import write_release_from_schema


def build_exclusion_metabolites(connection, *, schema, seed_directory=None):
    """The product of the release held in ``schema``, loaded back into it; returns its counts."""
    with connection.cursor() as cur:
        cur.execute(f'SET search_path = "{schema}", public')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            version = write_release_from_schema(cur, root)
            manifests = build_release_subsets(
                root, version, ["exclusion_metabolites"], seed_directory=seed_directory
            )
            product = root / "subsets" / version / "exclusion_metabolites"
            local = duckdb.connect()
            members = local.execute(
                "SELECT entity_key, exclusion_class, assignment_source, origin, match_kind, "
                "seed_namespace, seed_identifier, evidence, via_entity_key, list_version, "
                f"build_id FROM '{product}/exclusion_metabolite.parquet'"
            ).fetchall()
            report = local.execute(
                "SELECT kind, exclusion_class, namespace, resolved, fallback_only, unresolved, "
                f"detail, build_id FROM '{product}/exclusion_metabolite_report.parquet'"
            ).fetchall()
        cur.execute(
            "DROP TABLE IF EXISTS exclusion_metabolite, exclusion_metabolite_report CASCADE"
        )
        # The tables of a loaded release: the PostgreSQL package's own definition.
        create_tables(connection, schema, "exclusion_metabolites")
        cur.execute(f'SET search_path = "{schema}", public')
        execute_values(
            cur,
            "INSERT INTO exclusion_metabolite VALUES %s",
            members,
            template="(%s::uuid, %s, %s, %s, %s, %s, %s, %s, %s::uuid, %s, %s)",
        )
        execute_values(
            cur,
            "INSERT INTO exclusion_metabolite_report VALUES %s",
            [
                row[:6]
                + (Json(None if row[6] is None else __import__("json").loads(row[6])), row[7])
                for row in report
            ],
        )
    connection.commit()
    return manifests["exclusion_metabolites"]["stats"]


def resolve_seeds(cur, wanted):
    """The seeds resolved against the release held in the cursor's schema."""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        version = write_release_from_schema(cur, root)
        con = open_release(root, version)
        try:
            return parquet_build.resolve_seeds(con, wanted)
        finally:
            con.close()
