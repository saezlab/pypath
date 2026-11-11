"""The DuckDB database of a release: every entity and relation merged over its resources.

One database format serves the hosted API (one database per published release, written
on the server by ``build_release_database``) and a local installation (written by
``init()`` into its data directory). The resources' tables are written to the build's
scratch database; the database's own tables are derived from them (the contract's
``derived_tables``, sql/build/): one row per entity (``entity``, numbered ``entity_n``
in key order) and per relation (``relation``, ``relation_n``) with the resources that
state it, their merged identifiers and annotations, every evidence record with its
resource, and lookup tables sorted for the reads they serve. ``resources`` lists what
the database holds and ``database_info`` describes it.

Rows are read by number or key, written into the SQL as literal lists: then DuckDB uses
the ART indexes the contract names (``indexes``), and filters a sorted column by its
zone maps. A database is always written completely and put in place when done; nothing
is updated in place.
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

from . import steps
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


def _names(db) -> None:
    """The Biolink names (op_names), which the shared macros read."""
    db.execute("CREATE OR REPLACE TABLE op_names (kind VARCHAR, alias VARCHAR, name VARCHAR)")
    db.execute(
        "INSERT INTO op_names VALUES "
        + ", ".join(
            f"({_literal(k)}, {_literal(a)}, {_literal(n)})" for k, a, n in CONTRACT["biolink_names"]
        )
    )


def _query_tables(db) -> None:
    """The rule tables the macros read (display names, entity groups)."""
    for table, spec in CONTRACT.get("query_tables", {}).items():
        db.execute(f"CREATE OR REPLACE TABLE {table} ({', '.join(spec['columns'])})")
        if spec["rows"]:
            db.execute(
                f"INSERT INTO {table} VALUES "
                + ", ".join("(" + ", ".join(_value_literal(v) for v in row) + ")" for row in spec["rows"])
            )


def _query_layer(db) -> None:
    """Biolink names, the rule tables and the macros; tables before what reads them."""
    _names(db)
    _query_tables(db)
    for name in CONTRACT["macro_files"]:
        db.execute(_SQL.joinpath(name).read_text())
    db.execute("UPDATE database_info SET query_layer = ?", [QUERY_LAYER])


def _value_literal(value) -> str:
    return str(int(value)) if isinstance(value, int) and not isinstance(value, bool) else _literal(str(value))


# Where each resource table is written: the scratch database (``work``) for those the
# build reads, the database itself for raw evidence (installations only).
_WORK = set(CONTRACT["work_tables"])
_VIEWS = set(CONTRACT.get("work_views", []))
_PARTITION = re.compile(r"\{\{partition ([\w.]+)\}\}")
# Keys are SHA-256 hex: their first characters, then any other.
_SLICES = list("0123456789abcdef")


def _slice(select: str, first: str | None) -> str:
    """A build statement for the keys that start with ``first`` (``""``: with another
    character; None: none)."""
    def condition(m):
        if first is None:
            return "false"
        if first == "":
            return f"left(coalesce({m[1]}, ''), 1) NOT IN ({', '.join(_literal(c) for c in _SLICES)})"
        return f"left({m[1]}, 1) = {_literal(first)}"

    return _PARTITION.sub(condition, select)


def _derived(db, say: Callable[[str], None]) -> None:
    """The derived tables (each sorted), in build order, and the indexes. The stages of
    the build ("scratch") go to the scratch database."""
    for name, query in CONTRACT["derived_tables"].items():
        say(f"Deriving {name}")
        target = f"{'work' if query.get('scratch') else 'omnipath'}.main.{name}"
        select, order = query["select"], query.get("order")
        if _PARTITION.search(select) is None:
            db.execute(f"CREATE TABLE {target} AS SELECT * FROM ({select})" + (f" ORDER BY {order}" if order else ""))
            continue
        # In slices of keys by their first character, in key order: the table is sorted
        # as a whole (its order follows the keys), and no slice needs the memory of all.
        db.execute(f"CREATE TABLE {target} AS SELECT * FROM ({_slice(select, None)}) LIMIT 0")
        for first in [*_SLICES, ""]:
            db.execute(f"INSERT INTO {target} SELECT * FROM ({_slice(select, first)}) ORDER BY {order}")
    for name, table, column in CONTRACT["indexes"]:
        db.execute(f"CREATE INDEX {name} ON omnipath.main.{table}({column})")


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
            present = [t for t in TABLES + SERVING if any(t in s.tables(payloads) for s in ordered)]
            for table in present:
                if table in _WORK:
                    where = "work."
                elif table == _PAYLOADS:
                    where = ""
                else:
                    continue
                having = [s for s in ordered if table in s.tables(payloads)]
                # The columns of every resource (union by name), then each resource's
                # rows in its file's order, so rows follow their global ids.
                union = " UNION ALL BY NAME ".join(s.select(table, index[s.resource_id]) for s in having)
                if table in _VIEWS:
                    db.execute(f"CREATE VIEW work.main.{table} AS {union}")
                    continue
                say(f"Writing {table} ({len(having)} resources)")
                db.execute(f"CREATE TABLE {where}{table} AS SELECT * FROM ({union}) LIMIT 0")
                for source in having:
                    db.execute(f"INSERT INTO {where}{table} BY NAME " + source.select(table, index[source.resource_id]))
            # The merged entities' summaries name taxa and use the shared macros and rules.
            _taxonomy(db, taxonomy)
            _names(db)
            _query_tables(db)
            steps.prepare(db)
            _derived(db, say)
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

    Changed resources or taxon names (the merged entities' summaries name taxa) mean a
    new database, written completely. When only the query layer differs, it alone is
    rewritten.
    """
    info = _info(root)
    taxonomy = installation.get("taxonomy")
    if _tables_in_sync(info, installation) and info[3] == taxonomy:
        if info[0] != QUERY_LAYER:
            _layer_only(root)
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


def _layer_only(root: Path) -> None:
    """Write the query layer of an otherwise current database."""
    db = duckdb.connect(str(root / DATABASE))
    try:
        db.execute("BEGIN")
        _query_layer(db)
        db.execute("COMMIT")
    finally:
        db.close()
