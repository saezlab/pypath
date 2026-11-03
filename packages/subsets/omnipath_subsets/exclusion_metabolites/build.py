"""Build the exclusion metabolite tables from a loaded release.

The product reads what the release already published. It resolves no
identifier, calls no service and downloads nothing. Both tables are replaced
whole inside the product transaction, so the same release and the same seed
files always give the same rows.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

from psycopg2 import sql
from psycopg2.extras import Json, execute_values

from omnipath_subsets.exclusion_metabolites import resolve, seeds, version

_SQL_DIR = Path(__file__).with_name("sql")

logger = logging.getLogger(__name__)

# The ChEBI term each derived class is the descendants of. The roots are named
# here, so the size of a derived class is bounded by a term a person can read.
ROOTS = {"ions": "CHEBI:24867", "atoms": "CHEBI:33250"}

# The ChEBI role whose carriers are the derived extended cofactors. The class
# `cofactors` is the reviewed core and its variants; the role adds to it only in
# `cofactors_extended`, so a consumer can leave the role-derived layer out.
ROLE_ROOTS = {"cofactors_extended": "CHEBI:23357"}

ORIGIN = {
    "metaproviz_seed": "curated",
    "project_seed": "curated",
    "override_add": "curated",
    "chebi_hierarchy": "derived",
    "chebi_role": "derived",
}


def _sql_text(name: str) -> str:
    return (_SQL_DIR / name).read_text(encoding="utf-8")


def ensure_tables(conn, *, schema: str = "public") -> None:
    """Create both tables and the review view if they are absent.

    The DDL is idempotent. Applying it to populated tables leaves the shape and
    the rows untouched.
    """
    with conn.cursor() as cur:
        # Unqualified names in the DDL land in the target schema. public.*
        # still resolves for anything the view reads from the canonical layer.
        cur.execute(sql.SQL("SET search_path = {}, public").format(sql.Identifier(schema)))
        cur.execute(_sql_text("exclusion_tables.sql"))
    conn.commit()


def _build_id(cur) -> str:
    cur.execute("SELECT manifest_sha256 FROM parquet_release WHERE singleton")
    row = cur.fetchone()
    if row is None:
        raise LookupError("exclusion_metabolites: parquet_release is empty; load a release first")
    return row[0]


def _coverage(resolutions: list[resolve.Resolution]) -> list[dict]:
    """Seed counts per class and namespace, with the seeds behind each count."""
    groups: dict[tuple[str, str], list[resolve.Resolution]] = {}
    for resolution in resolutions:
        key = (resolution.seed.exclusion_class, resolution.seed.namespace)
        groups.setdefault(key, []).append(resolution)
    rows = []
    for (exclusion_class, namespace), members in sorted(groups.items()):
        unresolved = [r.seed.identifier for r in members if not r.resolved]
        fallback_only = [r.seed.identifier for r in members if r.fallback_only]
        rows.append(
            {
                "exclusion_class": exclusion_class,
                "namespace": namespace,
                "resolved": sum(1 for r in members if r.resolved),
                "fallback_only": len(fallback_only),
                "unresolved": len(unresolved),
                "detail": {"unresolved": unresolved, "fallback_only": fallback_only},
            }
        )
    return rows


def _resolve_roots(cur, terms: dict[str, str]) -> list[tuple[str, str, str]]:
    """The entity behind each ChEBI term, as (class, entity, term).

    A term the release lacks is left out.
    """
    wanted = [
        seeds.Seed(name, "chebi", curie, "metaproviz_seed", "roots")
        for name, curie in terms.items()
    ]
    return [
        (r.seed.exclusion_class, entity_id, r.seed.identifier)
        for r in resolve.resolve_seeds(cur, wanted)
        for entity_id in r.entity_ids
    ]


def _curated_rows(resolutions, *, list_version: str, build_id: str) -> list[tuple]:
    """One row per entity, class and source; the first seed to name it keeps the row."""
    rows: dict[tuple, tuple] = {}
    for resolution in resolutions:
        seed = resolution.seed
        for entity_id in resolution.entity_ids:
            key = (entity_id, seed.exclusion_class, seed.assignment_source)
            rows.setdefault(
                key,
                (
                    entity_id,
                    seed.exclusion_class,
                    seed.assignment_source,
                    ORIGIN[seed.assignment_source],
                    "entity",
                    seed.namespace,
                    seed.identifier,
                    seed.line,
                    None,
                    list_version,
                    build_id,
                ),
            )
    return [rows[key] for key in sorted(rows)]


def build_exclusion_metabolites(
    conn,
    *,
    schema: str = "public",
    seed_directory: Path | str | None = None,
) -> dict:
    """Publish the exclusion tables for the release loaded in ``schema``."""
    ensure_tables(conn, schema=schema)
    seed_set = seeds.read_seeds(seed_directory or seeds.SEED_DIR)
    listed = version.list_version(seed_set.paths)
    with conn.cursor() as cur:
        cur.execute(sql.SQL("SET search_path = {}, public").format(sql.Identifier(schema)))
        build_id = _build_id(cur)
        resolutions = resolve.resolve_seeds(cur, seed_set.seeds)
        rows = _curated_rows(resolutions, list_version=listed.value, build_id=build_id)
        coverage = _coverage(resolutions)

        roots = _resolve_roots(cur, ROOTS)
        role_roots = _resolve_roots(cur, ROLE_ROOTS)

        cur.execute("TRUNCATE exclusion_metabolite, exclusion_metabolite_report")
        execute_values(
            cur,
            "INSERT INTO exclusion_metabolite (entity_id, exclusion_class, assignment_source, "
            "origin, match_kind, seed_namespace, seed_identifier, evidence, via_entity_id, "
            "list_version, build_id) VALUES %s",
            rows,
        )
        cur.execute(_sql_text("stage_edges.sql"))
        cur.execute(
            _sql_text("derive_closure.sql"),
            {
                "root_classes": [r[0] for r in roots],
                "root_entities": [r[1] for r in roots],
                "root_curies": [r[2] for r in roots],
                "entity_types": list(resolve.CHEMICAL_ENTITY_TYPES),
                "list_version": listed.value,
                "build_id": build_id,
            },
        )
        cur.execute(
            _sql_text("derive_roles.sql"),
            {
                "role_classes": [r[0] for r in role_roots],
                "role_entities": [r[1] for r in role_roots],
                "role_curies": [r[2] for r in role_roots],
                "entity_types": list(resolve.CHEMICAL_ENTITY_TYPES),
                "list_version": listed.value,
                "build_id": build_id,
            },
        )
        cur.execute("DROP TABLE IF EXISTS exclusion_subclass_edge")

        vetoes = resolve.resolve_seeds(cur, seed_set.vetoes)
        cur.execute(_sql_text("stage_vetoes.sql"))
        execute_values(
            cur,
            "INSERT INTO exclusion_veto (exclusion_class, namespace, identifier, entity_id) "
            "VALUES %s",
            [
                (v.seed.exclusion_class, v.seed.namespace, v.seed.identifier, entity_id)
                for v in vetoes
                for entity_id in v.entity_ids
            ],
        )
        cur.execute(_sql_text("apply_vetoes.sql"))
        cur.execute(_sql_text("expand_structure_key.sql"))
        cur.execute(_sql_text("apply_vetoes.sql"))
        cur.execute(
            _sql_text("report.sql"),
            {
                "build_id": build_id,
                "coverage": json.dumps(coverage),
                "derived_classes": sorted({r[0] for r in roots}),
                "vetoes": json.dumps(
                    [
                        {
                            "exclusion_class": v.seed.exclusion_class,
                            "namespace": v.seed.namespace,
                            "identifier": v.seed.identifier,
                            "note": v.seed.note,
                            "resolved": v.resolved,
                            "entities": list(v.entity_ids),
                        }
                        for v in vetoes
                    ]
                ),
            },
        )
        if listed.uncommitted:
            cur.execute(
                "INSERT INTO exclusion_metabolite_report (kind, detail, build_id) "
                "VALUES ('uncommitted_inputs', %s, %s)",
                (
                    Json(
                        {
                            "list_version": listed.value,
                            "files": [path.name for path in seed_set.paths],
                        }
                    ),
                    build_id,
                ),
            )
        # Reaction counts need the participant table, which a bare release may lack.
        cur.execute("SELECT to_regclass('interaction_party') IS NOT NULL")
        if cur.fetchone()[0]:
            cur.execute(_sql_text("cofactor_reactions.sql"), {"build_id": build_id})
        cur.execute("ANALYZE exclusion_metabolite")
        cur.execute("ANALYZE exclusion_metabolite_report")
        cur.execute(
            "SELECT exclusion_class, count(DISTINCT entity_id), count(*) "
            "FROM exclusion_metabolite GROUP BY 1"
        )
        sizes = {name: (entities, total) for name, entities, total in cur.fetchall()}

    stats = {
        "build_id": build_id,
        "list_version": listed.value,
        "uncommitted_inputs": listed.uncommitted,
        "rows": sum(total for _, total in sizes.values()),
        "classes": {name: sizes.get(name, (0, 0))[0] for name in seeds.CLASSES},
        "unresolved_roots": sorted(
            (set(ROOTS) | set(ROLE_ROOTS)) - {r[0] for r in roots} - {r[0] for r in role_roots}
        ),
        "resolved_seeds": sum(c["resolved"] for c in coverage),
        "fallback_only_seeds": sum(c["fallback_only"] for c in coverage),
        "unresolved_seeds": sum(c["unresolved"] for c in coverage),
    }
    logger.info("exclusion_metabolites: %s", stats)
    return stats
