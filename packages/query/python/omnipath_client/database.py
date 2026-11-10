"""The DuckDB database of a release: native tables over every resource, written anew.

One database format serves the hosted API (one database per published release, written
on the server by ``build_release_database``) and a local installation (written by
``init()`` into its data directory). Each published and lookup table of the resources
becomes one native table (``entity``, ``relation``, ...) with the provenance columns
``_resource`` and ``_resource_version`` and a global id (``_rid``, ``_eid``: the
resource's index times 2^32 plus its own id). The same tables per resource are views in
a schema named after it (``signor.relation``); ``resources`` lists what the database
holds and ``database_info`` describes it.

Rows are read by id through the global ids, written into the SQL as literal lists: then
DuckDB uses the ART indexes the contract names (``indexes``), and filters a sorted
column by its zone maps. The database's derived tables (``derived_tables``) hold what is
read by a key the tables are not sorted by. A database is always written completely and
put in place when done; nothing is updated in place.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from datetime import datetime, timezone
import hashlib
from importlib.resources import files
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
from typing import Any

import duckdb

from .contracts import CONTRACT, SERVING, TABLES, Installation, ResourceRecord

DATABASE = "omnipath.duckdb"
DATABASE_FORMAT: int = CONTRACT["database_format"]
_PAYLOADS = "evidence_payloads"
_SQL = files("omnipath_client").joinpath("sql")
_GLOBAL_IDS: dict[str, list[str]] = CONTRACT["global_ids"]
# Changes with the contract or the query SQL: a database written by another client
# version gets its query layer written again on connect.
QUERY_LAYER = hashlib.sha256(
    b"".join(_SQL.joinpath(name).read_bytes() for name in ["contract.json", *CONTRACT["macro_files"]])
).hexdigest()[:16]


def _literal(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _name(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


class Source:
    """One resource version to write: its record and the directory of its files."""

    def __init__(self, record: ResourceRecord, directory: Path) -> None:
        self.record = record
        self.directory = Path(directory)
        self.resource_id = record["resource_id"]
        self.version = record["version"]

    def tables(self, payloads: bool = True) -> list[str]:
        """The tables it has: raw evidence only when it was installed and wanted."""
        names = {f["name"] for f in self.record["files"]}
        return [
            t
            for t in TABLES + SERVING
            if t != _PAYLOADS or (payloads and f"{_PAYLOADS}.parquet" in names)
        ]

    def select(self, table: str, index: int | None = None) -> str:
        """Its table read from its file, with the provenance columns; with the resource's
        ``index``, also with the table's global id."""
        path = self.directory / f"{table}.parquet"
        extra = ""
        if index is not None and table in _GLOBAL_IDS:
            column, global_id = _GLOBAL_IDS[table]
            extra = f", ({int(index)}::BIGINT << 32) + {column} AS {global_id}"
        return (
            f"SELECT *, {_literal(self.resource_id)} AS _resource, "
            f"{_literal(self.version)} AS _resource_version{extra} "
            f"FROM read_parquet({_literal(str(path))})"
        )


def _state(sources: Iterable[Source], payloads: bool) -> list[tuple[str, str, bool]]:
    return sorted(
        (s.resource_id, s.version, _PAYLOADS in s.tables(payloads)) for s in sources
    )


def _describe(db, sources: list[Source], payloads: bool, info: dict[str, Any]) -> None:
    db.execute(
        "CREATE OR REPLACE TABLE database_info AS SELECT ?::VARCHAR AS data_dir, "
        "true AS materialized, ?::VARCHAR AS query_layer, ?::VARCHAR AS release, "
        "?::INTEGER AS database_format, ?::TIMESTAMPTZ AS built_at, ?::VARCHAR AS built_by, "
        "?::VARCHAR AS taxonomy",
        [
            info.get("data_dir"),
            QUERY_LAYER,
            info.get("release"),
            DATABASE_FORMAT,
            info.get("built_at"),
            info.get("built_by"),
            info.get("taxonomy"),
        ],
    )
    db.execute(
        "CREATE OR REPLACE TABLE resources "
        "(resource_id VARCHAR, version VARCHAR, has_evidence BOOLEAN)"
    )
    for row in _state(sources, payloads):
        db.execute("INSERT INTO resources VALUES (?, ?, ?)", row)


def _resource_schema(db, source: Source, payloads: bool) -> None:
    schema = _name(source.resource_id)
    db.execute(f"CREATE SCHEMA IF NOT EXISTS {schema}")
    for table in source.tables(payloads):
        db.execute(
            f"CREATE OR REPLACE VIEW {schema}.{table} AS "
            f"SELECT * FROM main.{table} WHERE _resource = {_literal(source.resource_id)}"
        )


def native_lookups() -> list[str]:
    """Lookups over the native tables: one filter on one table, by global id or by key."""
    statements = []
    for table, column, kind in CONTRACT["lookups"]:
        global_id = _GLOBAL_IDS[table][1]
        if kind == "at":
            parameter = "pairs"
            body = (
                f"SELECT * EXCLUDE ({global_id}) FROM {table} WHERE {global_id} IN ("
                "SELECT (i.idx << 32) + w.p.v FROM (SELECT unnest(pairs) AS p) w "
                "JOIN op_resource_index i ON i._resource = w.p.r)"
            )
        else:
            parameter = "vals"
            body = f"SELECT * EXCLUDE ({global_id}) FROM {table} WHERE {column} IN (SELECT unnest(vals))"
        statements.append(f"CREATE OR REPLACE MACRO op_{table}_{kind}({parameter}) AS TABLE {body}")
    return statements


def _query_layer(db) -> None:
    """Biolink names, the lookups and the macros; tables before what reads them."""
    db.execute("CREATE OR REPLACE TABLE op_names (kind VARCHAR, alias VARCHAR, name VARCHAR)")
    db.execute(
        "INSERT INTO op_names VALUES "
        + ", ".join(
            f"({_literal(k)}, {_literal(a)}, {_literal(n)})" for k, a, n in CONTRACT["biolink_names"]
        )
    )
    present = {
        row[0]
        for row in db.execute(
            "SELECT table_name FROM duckdb_tables() WHERE schema_name = 'main' "
            "AND database_name = current_database()"
        ).fetchall()
    }
    for statement in native_lookups():
        table = re.match(r"CREATE OR REPLACE MACRO op_(\w+)_(at|with)\(", statement).group(1)
        if table in present:
            db.execute(statement)
    # The rule tables the macros read (display names, entity groups).
    for table, spec in CONTRACT.get("query_tables", {}).items():
        db.execute(f"CREATE OR REPLACE TABLE {table} ({', '.join(spec['columns'])})")
        if spec["rows"]:
            db.execute(
                f"INSERT INTO {table} VALUES "
                + ", ".join("(" + ", ".join(_value_literal(v) for v in row) + ")" for row in spec["rows"])
            )
    for name in CONTRACT["macro_files"]:
        db.execute(_SQL.joinpath(name).read_text())
    db.execute("UPDATE database_info SET query_layer = ?", [QUERY_LAYER])


def _value_literal(value) -> str:
    return str(int(value)) if isinstance(value, int) and not isinstance(value, bool) else _literal(str(value))


# Derived tables that replace a written table of the same name: {name: its select}.
REPLACED = {
    name: query["select"] for name, query in CONTRACT["derived_tables"].items() if query.get("replaces")
}


def _written_to(table: str) -> str | None:
    """Where a resource table is written: the database, the scratch database ``work`` for
    a table a derived table replaces (its space is never part of the database), or
    nowhere (None) for one replaced without being read."""
    if table not in REPLACED:
        return ""
    return "work." if re.search(rf"\b{table}\b", REPLACED[table]) else None


def _derived(db, present: set[str]) -> None:
    """The derived tables (each sorted) and the indexes, once the tables are written. A
    table that replaces a written one reads it from ``work`` (on the search path)."""
    for name, query in CONTRACT["derived_tables"].items():
        db.execute(
            f"CREATE TABLE omnipath.main.{name} AS SELECT * FROM ({query['select']}) ORDER BY {query['order']}"
        )
    for name, table, column in CONTRACT["indexes"]:
        if table in present:
            db.execute(f"CREATE INDEX {name} ON {table}({column})")


def _taxonomy(db, path: Path | None) -> None:
    """The release's taxon names (common name, else scientific name), empty without its
    taxonomy reference."""
    if path is not None and Path(path).is_file():
        db.execute(
            f"""CREATE OR REPLACE TABLE op_taxonomy AS SELECT DISTINCT taxon_id::VARCHAR AS taxon,
                coalesce(nullif(common_name, ''), nullif(scientific_name, '')) AS name
            FROM read_parquet({_literal(str(path))})
            WHERE coalesce(nullif(common_name, ''), nullif(scientific_name, '')) IS NOT NULL
            ORDER BY taxon"""
        )
    else:
        db.execute("CREATE OR REPLACE TABLE op_taxonomy (taxon VARCHAR, name VARCHAR)")


def taxonomy_file(root: Path, version: str) -> Path:
    """Where a taxonomy reference (named by its SHA-256) is kept, on a server and in a
    local installation alike."""
    return Path(root) / "references" / "taxonomy" / version / "taxonomy.parquet"


def taxonomy_reference(root: Path, manifest: dict[str, Any]) -> Path | None:
    """A release's taxonomy reference file under ``root``, if the release names one."""
    reference = (manifest.get("references") or {}).get("taxonomy")
    if not reference:
        return None
    path = taxonomy_file(root, reference["version"])
    return path if path.is_file() else None


def write_database(
    target: str | Path,
    sources: list[Source],
    *,
    info: dict[str, Any] | None = None,
    taxonomy: str | Path | None = None,
    payloads: bool = True,
    memory_limit: str = "2GB",
    threads: int | None = None,
    temp_directory: str | Path | None = None,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Write the database of ``sources`` to ``target`` and put it in place when complete.

    Resources are numbered in name order for the global ids, and each table's rows
    follow them. ``payloads`` includes raw evidence where a resource has it. Returns the
    target, its size and the seconds taken.
    """
    target = Path(target)
    say = progress or (lambda _message: None)
    started = time.perf_counter()
    info = dict(info or {})
    info.setdefault("built_at", datetime.now(timezone.utc).isoformat())
    staging_dir = Path(tempfile.mkdtemp(dir=target.parent, prefix=f".{target.name}-"))
    staging = staging_dir / target.name
    temporary = Path(temp_directory) if temp_directory else staging_dir / "tmp"
    try:
        db = duckdb.connect()
        try:
            db.execute(
                f"ATTACH {_literal(str(staging))} AS omnipath "
                f"(STORAGE_VERSION {_literal(CONTRACT['storage_version'])})"
            )
            db.execute(f"ATTACH {_literal(str(staging_dir / 'work.duckdb'))} AS work")
            db.execute("USE omnipath")
            db.execute("SET search_path = 'omnipath.main,work.main'")
            db.execute("SET memory_limit = ?", [memory_limit])
            db.execute("SET temp_directory = ?", [str(temporary)])
            if threads is not None:
                db.execute("SET threads = ?", [int(threads)])
            ordered = sorted(sources, key=lambda s: s.resource_id)
            index = {s.resource_id: i for i, s in enumerate(ordered)}
            db.execute("CREATE TABLE op_resource_index (_resource VARCHAR, idx BIGINT)")
            for name, i in index.items():
                db.execute("INSERT INTO op_resource_index VALUES (?, ?)", [name, i])
            present = [t for t in TABLES + SERVING if any(t in s.tables(payloads) for s in ordered)]
            for table in present:
                where = _written_to(table)
                if where is None:
                    continue
                having = [s for s in ordered if table in s.tables(payloads)]
                say(f"Writing {table} ({len(having)} resources)")
                # The columns of every resource (union by name), then each resource's
                # rows in its file's order, so rows follow their global ids.
                union = " UNION ALL BY NAME ".join(s.select(table, index[s.resource_id]) for s in having)
                db.execute(f"CREATE TABLE {where}{table} AS SELECT * FROM ({union}) LIMIT 0")
                for source in having:
                    db.execute(f"INSERT INTO {where}{table} BY NAME " + source.select(table, index[source.resource_id]))
            say("Writing the derived tables and indexes")
            _derived(db, set(present))
            for source in ordered:
                _resource_schema(db, source, payloads)
            _taxonomy(db, taxonomy)
            _describe(db, ordered, payloads, info)
            _query_layer(db)
            db.execute("USE memory")
            db.execute("DETACH work")
            db.execute("DETACH omnipath")
        finally:
            db.close()
        os.replace(staging, target)
    finally:
        shutil.rmtree(staging_dir, ignore_errors=True)
    return {
        "path": str(target),
        "size_bytes": target.stat().st_size,
        "seconds": round(time.perf_counter() - started, 1),
        "resources": len(sources),
    }


def release_sources(data_root: str | Path, resources: dict[str, str]) -> list[Source]:
    """The resource versions of a release manifest, as published under ``data_root``."""
    root = Path(data_root)
    sources = []
    for resource_id, version in sorted(resources.items()):
        directory = root / "resources" / resource_id / version
        names = sorted(p.name for p in directory.glob("*.parquet"))
        missing = [f"{t}.parquet" for t in TABLES + SERVING if t != _PAYLOADS and f"{t}.parquet" not in names]
        if missing:
            raise FileNotFoundError(f"{resource_id}/{version} lacks {', '.join(missing)}")
        record = ResourceRecord(
            resource_id=resource_id,
            version=version,
            files=[
                {"name": n, "size_bytes": (directory / n).stat().st_size}
                for n in names
                if Path(n).stem in TABLES
            ],
        )
        sources.append(Source(record, directory))
    return sources


def build_release_database(
    data_root: str | Path,
    release: str,
    target: str | Path,
    *,
    memory_limit: str = "4GB",
    threads: int | None = None,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Write the database of a published release (``releases/<release>.json``) from its
    resources' files under ``data_root``. Raw evidence payloads stay in their files."""
    root = Path(data_root)
    manifest = json.loads((root / "releases" / f"{release}.json").read_text())
    if manifest.get("version") != release:
        raise ValueError(f"Release file and version disagree: {release}")
    return write_database(
        target,
        release_sources(root, manifest["resources"]),
        info={
            "release": release,
            "built_by": "build_release_database",
            "taxonomy": ((manifest.get("references") or {}).get("taxonomy") or {}).get("version"),
        },
        taxonomy=taxonomy_reference(root, manifest),
        payloads=False,
        memory_limit=memory_limit,
        threads=threads,
        progress=progress,
    )


# The database of a local installation ----------------------------------------------


def _installation_sources(root: Path, records: list[ResourceRecord]) -> list[Source]:
    return [
        Source(r, root / "resources" / r["resource_id"] / r["version"]) for r in records
    ]


def _info(root: Path) -> tuple | None:
    """(query layer, database format, resources, taxonomy version) of the database, or
    None."""
    path = root / DATABASE
    if not path.is_file():
        return None
    try:
        db = duckdb.connect(str(path), read_only=True)
    except duckdb.Error:
        return None
    try:
        columns = {row[0] for row in db.execute("DESCRIBE database_info").fetchall()}
        if "database_format" not in columns:
            return None
        taxonomy = "taxonomy" if "taxonomy" in columns else "NULL"
        layer, database_format, version = db.execute(
            f"SELECT query_layer, database_format, {taxonomy} FROM database_info"
        ).fetchone()
        state = db.execute("SELECT * FROM resources ORDER BY ALL").fetchall()
    except duckdb.Error:
        return None
    finally:
        db.close()
    return layer, database_format, state, version


def _tables_in_sync(info: tuple | None, installation: Installation) -> bool:
    if info is None:
        return False
    _, database_format, state, _ = info
    sources = _installation_sources(Path("."), installation["resources"])
    return database_format == DATABASE_FORMAT and state == _state(sources, True)


def in_sync(root: Path, installation: Installation) -> bool:
    """Whether the database holds exactly the installed resources and taxon names, in
    this client's database format and with its query layer."""
    info = _info(root)
    return (
        _tables_in_sync(info, installation)
        and info[0] == QUERY_LAYER
        and info[3] == installation.get("taxonomy")
    )


def synchronize(root: Path, installation: Installation, memory_limit: str = "2GB") -> None:
    """Bring the database in line with the installation; nothing when it already is.

    Changed resources mean a new database, written completely. When only the query
    layer or the taxon names differ, they alone are rewritten.
    """
    info = _info(root)
    taxonomy = installation.get("taxonomy")
    if _tables_in_sync(info, installation):
        if info[0] != QUERY_LAYER:
            _layer_only(root)
        if info[3] != taxonomy:
            _taxonomy_only(root, taxonomy)
        return
    temporary = root / ".duckdb-tmp"
    try:
        write_database(
            root / DATABASE,
            _installation_sources(root, installation["resources"]),
            info={
                "data_dir": str(root),
                "release": installation["release"],
                "built_by": "init",
                "taxonomy": taxonomy,
            },
            taxonomy=taxonomy_file(root, taxonomy) if taxonomy else None,
            memory_limit=memory_limit,
            temp_directory=temporary,
        )
    finally:
        if temporary.is_dir() and not any(temporary.iterdir()):
            temporary.rmdir()


def _taxonomy_only(root: Path, version: str | None) -> None:
    """Write the taxon names of an otherwise current database."""
    db = duckdb.connect(str(root / DATABASE))
    try:
        db.execute("BEGIN")
        _taxonomy(db, taxonomy_file(root, version) if version else None)
        db.execute("ALTER TABLE database_info ADD COLUMN IF NOT EXISTS taxonomy VARCHAR")
        db.execute("UPDATE database_info SET taxonomy = ?", [version])
        db.execute("COMMIT")
    finally:
        db.close()


def _layer_only(root: Path) -> None:
    """Write the query layer of an otherwise current database."""
    db = duckdb.connect(str(root / DATABASE))
    try:
        db.execute("BEGIN")
        _query_layer(db)
        db.execute("COMMIT")
    finally:
        db.close()
