"""Build the exclusion metabolite tables as Parquet from a published release.

The product reads the pinned resources' Parquet through DuckDB. It resolves no
identifier, calls no service and downloads nothing, so the same release and the
same seed files always give the same rows.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import duckdb
import pyarrow as pa

from omnipath_subsets.exclusion_metabolites import resolve, seeds, version
from omnipath_core.subset_layout import BuildContext
from omnipath_subsets.release_source import open_release

_SQL_DIR = Path(__file__).with_name("duckdb_sql")

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


logger = logging.getLogger(__name__)

MEMBER_TABLE = "exclusion_metabolite"
REPORT_TABLE = "exclusion_metabolite_report"


def _sql(name: str) -> str:
    return (_SQL_DIR / f"{name}.sql").read_text(encoding="utf-8")


def resolve_seeds(
    con: duckdb.DuckDBPyConnection, wanted: list[seeds.Seed]
) -> list[resolve.Resolution]:
    """Resolve every seed in one query; the result follows the order of ``wanted``."""
    spellings = [
        (index, seed.namespace, spelling)
        for index, seed in enumerate(wanted)
        for spelling in resolve.spellings(seed.namespace, seed.identifier)
    ]
    con.register(
        "seed_spellings",
        pa.table(
            {
                "seed_index": pa.array([s[0] for s in spellings], pa.int32()),
                "ns": pa.array([s[1] for s in spellings], pa.string()),
                "spelling": pa.array([s[2] for s in spellings], pa.string()),
            }
        ),
    )
    try:
        rows = con.execute(
            """
            SELECT DISTINCT w.seed_index, e.entity_key,
                   coalesce(list_contains($standard, e.namespace), false) AS standard_identity
            FROM seed_spellings w
            JOIN release_identifier i ON i.ns = w.ns AND i.id = w.spelling
            JOIN release_entity e ON e.entity_key = i.entity_key
            WHERE list_contains($chemical, e.entity_type)
            """,
            {
                "standard": list(resolve.STANDARD_NAMESPACES),
                "chemical": list(resolve.CHEMICAL_ENTITY_TYPES),
            },
        ).fetchall()
    finally:
        con.unregister("seed_spellings")
    found: dict[int, list[tuple[str, bool]]] = {}
    for index, entity_key, standard in rows:
        found.setdefault(index, []).append((entity_key, standard))
    resolutions = []
    for index, seed in enumerate(wanted):
        matches = sorted(found.get(index, []))
        resolutions.append(
            resolve.Resolution(
                seed,
                tuple(entity_key for entity_key, _ in matches),
                bool(matches) and not any(standard for _, standard in matches),
            )
        )
    return resolutions


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


def _insert(con, statement: str, rows: list) -> None:
    if rows:
        con.executemany(statement, rows)


def _resolve_roots(con, terms: dict[str, str]) -> list[tuple[str, str, str]]:
    """The entity behind each ChEBI term, as (class, entity, term); an absent term is left out."""
    wanted = [
        seeds.Seed(name, "chebi", curie, "metaproviz_seed", "roots")
        for name, curie in terms.items()
    ]
    return [
        (r.seed.exclusion_class, entity_key, r.seed.identifier)
        for r in resolve_seeds(con, wanted)
        for entity_key in r.entity_ids
    ]


def build_tables(
    con: duckdb.DuckDBPyConnection,
    *,
    build_id: str,
    seed_directory: Path | str | None = None,
) -> dict:
    """Fill ``exclusion_metabolite`` and its report in ``con``; returns the counts."""
    seed_set = seeds.read_seeds(seed_directory or seeds.SEED_DIR)
    listed = version.list_version(seed_set.paths)
    con.execute(_sql("tables"))
    chemical = ", ".join(f"'{name}'" for name in resolve.CHEMICAL_ENTITY_TYPES)
    con.execute(
        "CREATE OR REPLACE TEMP VIEW structure_group AS "
        "SELECT entity_key, group_connectivity AS group_key FROM release_entity "
        f"WHERE group_connectivity IS NOT NULL AND entity_type IN ({chemical})"
    )

    resolutions = resolve_seeds(con, seed_set.seeds)
    rows = _curated_rows(resolutions, list_version=listed.value, build_id=build_id)
    coverage = _coverage(resolutions)
    _insert(con, f"INSERT INTO {MEMBER_TABLE} VALUES (?,?,?,?,?,?,?,?,?,?,?)", rows)
    _insert(
        con,
        f"INSERT INTO {REPORT_TABLE} VALUES ('seed_coverage',?,?,?,?,?,?,?)",
        [
            (
                c["exclusion_class"],
                c["namespace"],
                c["resolved"],
                c["fallback_only"],
                c["unresolved"],
                json.dumps(c["detail"]),
                build_id,
            )
            for c in coverage
        ],
    )

    roots = _resolve_roots(con, ROOTS)
    role_roots = _resolve_roots(con, ROLE_ROOTS)
    _insert(con, "INSERT INTO exclusion_root VALUES (?,?,?)", roots)
    _insert(con, "INSERT INTO exclusion_role_root VALUES (?,?,?)", role_roots)
    parameters = {
        "entity_types": list(resolve.CHEMICAL_ENTITY_TYPES),
        "list_version": listed.value,
        "build_id": build_id,
    }
    con.execute(_sql("stage_edges"))
    con.execute(_sql("derive_closure"), parameters)
    con.execute(_sql("derive_roles"), parameters)

    vetoes = resolve_seeds(con, seed_set.vetoes)
    veto_rows = [
        (v.seed.exclusion_class, v.seed.namespace, v.seed.identifier, entity_key)
        for v in vetoes
        for entity_key in v.entity_ids
    ]
    _insert(con, "INSERT INTO exclusion_veto VALUES (?,?,?,?)", veto_rows)
    _insert(
        con,
        "INSERT INTO exclusion_veto_listing VALUES (?,?,?,?,?,?)",
        [
            (
                v.seed.exclusion_class,
                v.seed.namespace,
                v.seed.identifier,
                v.seed.note,
                v.resolved,
                json.dumps(list(v.entity_ids)),
            )
            for v in vetoes
        ],
    )
    for step in ("record_vetoes", "remove_vetoed"):
        con.execute(_sql(step))
    con.execute(_sql("expand_structure_key"))
    for step in ("record_vetoes", "remove_vetoed"):
        con.execute(_sql(step))

    derived = sorted({r[0] for r in roots})
    for report in (
        "report_class_size",
        "report_seed_not_derived",
        "report_derived_not_seeded",
        "report_shared_key",
        "report_no_key",
        "report_veto_applied",
    ):
        con.execute(
            _sql(report),
            {"build_id": build_id, "derived_classes": derived}
            if report in ("report_seed_not_derived", "report_derived_not_seeded")
            else {"build_id": build_id},
        )
    if listed.uncommitted:
        con.execute(
            f"INSERT INTO {REPORT_TABLE} (kind, detail, build_id) VALUES ('uncommitted_inputs', ?, ?)",
            [
                json.dumps(
                    {"list_version": listed.value, "files": [path.name for path in seed_set.paths]}
                ),
                build_id,
            ],
        )
    con.execute(_sql("cofactor_reactions"), {"build_id": build_id})

    duplicated = con.execute(
        "SELECT count(*) FROM (SELECT entity_key, exclusion_class, assignment_source "
        f"FROM {MEMBER_TABLE} GROUP BY ALL HAVING count(*) > 1)"
    ).fetchone()[0]
    if duplicated:
        raise ValueError(f"{duplicated} entities are listed twice in one class by one source")

    sizes = {
        name: (entities, total)
        for name, entities, total in con.execute(
            f"SELECT exclusion_class, count(DISTINCT entity_key), count(*) FROM {MEMBER_TABLE} "
            "GROUP BY 1"
        ).fetchall()
    }
    return {
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


def write_exclusion_metabolites(
    directory: Path,
    context: BuildContext,
    *,
    seed_directory: Path | str | None = None,
) -> dict:
    """Write both tables of the product for ``context.release`` into ``directory``."""
    con = open_release(context.data_root, context.release)
    try:
        stats = build_tables(
            con, build_id=context.release_manifest_sha256, seed_directory=seed_directory
        )
        for table, order in (
            (MEMBER_TABLE, "exclusion_class, entity_key, assignment_source"),
            (REPORT_TABLE, "kind, exclusion_class NULLS FIRST, namespace NULLS FIRST, detail"),
        ):
            target = str(directory / f"{table}.parquet").replace("'", "''")
            con.execute(
                f"COPY (SELECT * FROM {table} ORDER BY {order}) TO '{target}' "
                "(FORMAT parquet, COMPRESSION zstd)"
            )
    finally:
        con.close()
    logger.info("exclusion_metabolites: %s", stats)
    return stats
