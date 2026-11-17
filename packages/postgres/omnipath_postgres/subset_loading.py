"""Load the Parquet subsets of a release into its PostgreSQL schema.

The exclusion metabolites and the MetSigDB memberships are built once, from the published
release, as Parquet (``omnipath-subsets build``). This loads them through the loader's
staging and COPY path: DDL, COPY from the subset's Parquet, indexes, and the
``subset_build_metadata`` row and phase of the product. No product logic runs here.

A subset is loaded only for the release it was built from: its manifest names the digest of
the published release manifest, and the schema's release must be that release. Its entities
must be entities of the loaded release.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from time import perf_counter

import duckdb
from psycopg2 import sql

from omnipath_core.subset_layout import SubsetError, read_subset
from omnipath_subsets.runner import LOADED_PRODUCTS, _metadata, _publish

from .projection import AP_UUID_MACRO, CopyQuery, entity_id_sql

logger = logging.getLogger(__name__)

__all__ = ["LOADED_PRODUCTS", "load_subsets"]

_DDL_DIR = Path(__file__).with_name("relational") / "subsets"
_INDEXES = "-- @indexes"

# The tables of each product, in load order, with the files that hold them. Each is the
# Parquet file of the same name in the subset's directory.
_TABLES = {
    "exclusion_metabolites": ("exclusion_metabolite", "exclusion_metabolite_report"),
    "metsigdb": ("metsigdb_membership",),
}
_DDL = {
    "exclusion_metabolites": "exclusion_tables.sql",
    "metsigdb": "membership_table.sql",
}
# What a table's entity columns refer to, which the loaded release must have.
_ENTITY_COLUMNS = {
    "exclusion_metabolite": ("entity_id", "via_entity_id"),
    "metsigdb_membership": ("metabolite_entity_id", "set_entity_id"),
}


def _literal(value) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _copy_queries(product: str, directory: Path) -> list[CopyQuery]:
    """What to stage from each Parquet table: entity keys become the ids the release gave them."""
    source = lambda table: f"read_parquet({_literal(directory / (table + '.parquet'))})"  # noqa: E731
    entity = lambda column: entity_id_sql(column)  # noqa: E731
    if product == "exclusion_metabolites":
        return [
            CopyQuery(
                "exclusion_metabolite",
                (
                    "entity_id",
                    "exclusion_class",
                    "assignment_source",
                    "origin",
                    "match_kind",
                    "seed_namespace",
                    "seed_identifier",
                    "evidence",
                    "via_entity_id",
                    "list_version",
                    "build_id",
                ),
                f"SELECT {entity('entity_key')}, exclusion_class, assignment_source, origin, "
                f"match_kind, seed_namespace, seed_identifier, evidence, "
                f"CASE WHEN via_entity_key IS NULL THEN NULL ELSE {entity('via_entity_key')} END, "
                f"list_version, build_id FROM {source('exclusion_metabolite')}",
            ),
            CopyQuery(
                "exclusion_metabolite_report",
                (
                    "kind",
                    "exclusion_class",
                    "namespace",
                    "resolved",
                    "fallback_only",
                    "unresolved",
                    "detail",
                    "build_id",
                ),
                "SELECT kind, exclusion_class, namespace, resolved, fallback_only, unresolved, "
                f"detail, build_id FROM {source('exclusion_metabolite_report')}",
            ),
        ]
    return [
        CopyQuery(
            "metsigdb_membership",
            (
                "resource",
                "set_source_id",
                "metabolite_entity_id",
                "metabolite_label",
                "metabolite_entity_type",
                "inchikey",
                "metabolite_structure_key",
                "smiles",
                "hmdb",
                "pubchem",
                "chebi",
                "kegg",
                "has_structure",
                "exclusion_classes",
                "set_entity_id",
                "set_label",
                "set_type",
                "set_sub_type",
                "organism",
                "set_size",
                "set_context",
                "provenance_source",
                "provenance_record",
                "build_id",
            ),
            f"SELECT resource, set_source_id, {entity('metabolite_entity_key')}, "
            "metabolite_label, metabolite_entity_type, inchikey, metabolite_structure_key, "
            "smiles, hmdb, pubchem, chebi, kegg, has_structure, "
            # A PostgreSQL array literal.
            "'{' || array_to_string(list_transform(exclusion_classes, "
            "lambda c: '\"' || c || '\"'), ',') || '}', "
            f"{entity('set_entity_key')}, set_label, set_type, set_sub_type, organism, "
            "set_size, set_context, provenance_source, provenance_record, build_id "
            f"FROM {source('metsigdb_membership')}",
        )
    ]


def _check_release(data_root: Path, identity: tuple[str, str]) -> None:
    """The release loaded in the schema is the published release the subsets were built for."""
    version, digest = identity
    path = data_root / "releases" / f"{version}.json"
    if not path.is_file():
        raise ValueError(
            f"Release {version} is not published under {data_root / 'releases'}, so its "
            "subsets cannot be loaded"
        )
    canonical = json.dumps(
        json.loads(path.read_text()),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )
    if hashlib.sha256(canonical.encode("utf-8")).hexdigest() != digest:
        raise ValueError(
            f"The loaded release is not the published release {version} ({path}): subsets "
            "are loaded for the release they were built from"
        )


def _split(ddl: str) -> tuple[str, str]:
    head, _, tail = ddl.partition(_INDEXES)
    return head, tail.partition("\n")[2]


def ddl(product: str) -> str:
    """The SQL that creates the tables, indexes and views of a product, unqualified."""
    head, tail = _split((_DDL_DIR / _DDL[product]).read_text())
    return head + tail


def create_tables(conn, schema: str, product: str) -> None:
    """The tables, indexes and views of a product in ``schema``, empty (or as they are)."""
    head, tail = _split((_DDL_DIR / _DDL[product]).read_text())
    with conn.cursor() as cur:
        cur.execute(sql.SQL("SET search_path = {}, public").format(sql.Identifier(schema)))
        cur.execute(head)
        cur.execute(tail)
    conn.commit()


def _entity_check(conn, schema: str, table: str) -> None:
    """Every entity the table names is an entity of the loaded release."""
    for column in _ENTITY_COLUMNS[table]:
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL(
                    "SELECT count(*) FROM {t} x WHERE x.{c} IS NOT NULL AND NOT EXISTS "
                    "(SELECT 1 FROM {s}.entity e WHERE e.entity_id = x.{c})"
                ).format(
                    t=sql.Identifier(schema, table),
                    s=sql.Identifier(schema),
                    c=sql.Identifier(column),
                )
            )
            missing = cur.fetchone()[0]
        if missing:
            raise ValueError(
                f"{table}.{column} names {missing} entities the loaded release does not have; "
                "the subset was not built from this release"
            )


def _load_product(
    conn, schema, data_root, identity, product, *, spool, copy, memory_limit, threads
):
    try:
        subset = read_subset(data_root, identity[0], product)
    except SubsetError as error:
        raise ValueError(f"The {product} subset cannot be loaded: {error}") from error
    if subset is None:
        raise ValueError(
            f"Release {identity[0]} has no {product} subset under {data_root / 'subsets'}; "
            f"build it with `omnipath-subsets build --release {identity[0]} --data-root "
            f"{data_root}` (make build-subsets) or leave the product out"
        )
    manifest = subset.manifest
    started = perf_counter()
    head, tail = _split((_DDL_DIR / _DDL[product]).read_text())
    tables = _TABLES[product]
    with conn.cursor() as cur:
        cur.execute(sql.SQL("SET LOCAL search_path = {}, public").format(sql.Identifier(schema)))
        # A product is loaded whole: what an earlier load left is replaced.
        for table in reversed(tables):
            cur.execute(
                sql.SQL("DROP TABLE IF EXISTS {} CASCADE").format(sql.Identifier(schema, table))
            )
        cur.execute(head)
    loaded = {}
    with duckdb.connect() as con:
        con.execute(AP_UUID_MACRO)
        con.execute("SET memory_limit = ?", [memory_limit])
        con.execute("SET threads = ?", [threads])
        for item in _copy_queries(product, subset.directory):
            copied = copy(conn, schema, item, con, spool)
            expected = manifest["tables"][item.table + ".parquet"]["rows"]
            if copied.rows != expected:
                raise ValueError(
                    f"{item.table}: {copied.rows} rows loaded, the subset's manifest says {expected}"
                )
            loaded[item.table] = copied.rows
    for table in tables:
        if table in _ENTITY_COLUMNS:
            _entity_check(conn, schema, table)
    with conn.cursor() as cur:
        cur.execute(tail)
        for table in tables:
            cur.execute(sql.SQL("ANALYZE {}").format(sql.Identifier(schema, table)))
    stats = {
        "release": identity[0],
        "subset_manifest_sha256": subset.manifest_sha256,
        "row_contract": manifest["row_contract"],
        "partial": manifest["partial"],
        "built_at": manifest["built_at"],
        "code_version": manifest["code_version"],
        "rows": loaded,
        "stats": manifest.get("stats", {}),
    }
    return stats, round(perf_counter() - started, 6)


def load_subsets(
    conn,
    schema: str,
    data_root: str | Path,
    identity: tuple[str, str],
    products,
    *,
    spool: Path,
    observer=None,
    memory_limit: str = "1GB",
    duckdb_threads: int = 1,
) -> dict[str, dict]:
    """Load ``products`` for the release in ``schema``, each in one transaction.

    The caller holds the schema lock. A failure rolls back the product in progress and
    leaves the ones already committed.
    """
    from .loader import _copy, _emit

    root = Path(data_root).resolve()
    results = {}
    _check_release(root, identity)
    for product in products:
        if product not in LOADED_PRODUCTS:
            raise ValueError(f"{product} is built in PostgreSQL, not loaded from Parquet")
        _emit(observer, "phase_start", phase=product)
        try:
            _metadata(conn, schema)
            stats, seconds = _load_product(
                conn,
                schema,
                root,
                identity,
                product,
                spool=spool,
                copy=_copy,
                memory_limit=memory_limit,
                threads=duckdb_threads,
            )
            _publish(
                conn,
                schema,
                product,
                identity,
                {"phase_seconds": seconds, "result": stats},
                seconds,
            )
            conn.commit()
        except BaseException:
            conn.rollback()
            raise
        results[product] = stats
        _emit(observer, "phase_committed", phase=product, seconds=seconds, result=stats)
    return results
