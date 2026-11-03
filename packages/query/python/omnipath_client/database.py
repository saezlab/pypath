"""The DuckDB database of an installation: its tables and its query layer.

The tables are views over the installation's Parquet files, or native tables copied
from them: one table per published and lookup table over every installed resource
(``entity``, ``relation``, ...), the same tables per resource in a schema named after
it (``signor.relation``), and ``resources`` listing what the database holds.

The query layer is what ``search``, ``relations`` and ``entity`` run, in SQL any DuckDB
client can call: the macros of ``sql/*.sql``, Biolink names (``op_names``) and the
installation's lookup macros, which read rows of one table by id or key. Over views a
lookup has one branch per resource; native tables carry global ids (``_rid``, ``_eid``:
the resource's index times 2^32 plus its own id), so a lookup is one filter.
"""

from __future__ import annotations

import hashlib
from importlib.resources import files
import os
from pathlib import Path
import re
import tempfile

import duckdb

from .contracts import CONTRACT, SERVING, TABLES, Installation, ResourceRecord

DATABASE = "omnipath.duckdb"
_PAYLOADS = "evidence_payloads"
_SQL = files("omnipath_client").joinpath("sql")
_GLOBAL_IDS: dict[str, list[str]] = CONTRACT["global_ids"]
# Changes with the contract or a macro: a database written by another client version
# gets its query layer written again on connect.
QUERY_LAYER = hashlib.sha256(
    b"".join(_SQL.joinpath(name).read_bytes() for name in ["contract.json", *CONTRACT["macro_files"]])
).hexdigest()[:16]


def _literal(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _name(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _tables(record: ResourceRecord) -> list[str]:
    """The tables a resource has: raw evidence only when it was installed."""
    names = {f["name"] for f in record["files"]}
    return [t for t in TABLES + SERVING if t != _PAYLOADS or f"{_PAYLOADS}.parquet" in names]


def _source(root: Path, record: ResourceRecord, table: str, index: int | None = None) -> str:
    """One resource's table read from its file, with its provenance columns; with the
    resource's ``index``, also with the table's global id (native tables)."""
    path = root / "resources" / record["resource_id"] / record["version"] / f"{table}.parquet"
    extra = ""
    if index is not None and table in _GLOBAL_IDS:
        column, global_id = _GLOBAL_IDS[table]
        extra = f", ({int(index)}::BIGINT << 32) + {column} AS {global_id}"
    return (
        f"SELECT *, {_literal(record['resource_id'])} AS _resource, "
        f"{_literal(record['version'])} AS _resource_version{extra} "
        f"FROM read_parquet({_literal(str(path))})"
    )


def _state(records: list[ResourceRecord]) -> list[tuple[str, str, bool]]:
    return sorted((r["resource_id"], r["version"], _PAYLOADS in _tables(r)) for r in records)


def _describe(db, root: Path, materialized: bool, records: list[ResourceRecord]) -> None:
    db.execute(
        "CREATE OR REPLACE TABLE database_info AS SELECT ?::VARCHAR AS data_dir, "
        "?::BOOLEAN AS materialized, ?::VARCHAR AS query_layer",
        [str(root), materialized, QUERY_LAYER],
    )
    db.execute(
        "CREATE OR REPLACE TABLE resources "
        "(resource_id VARCHAR, version VARCHAR, has_evidence BOOLEAN)"
    )
    for row in _state(records):
        db.execute("INSERT INTO resources VALUES (?, ?, ?)", row)


def _resource_schema(db, root: Path, record: ResourceRecord, materialized: bool) -> None:
    schema = _name(record["resource_id"])
    db.execute(f"CREATE SCHEMA IF NOT EXISTS {schema}")
    for table in _tables(record):
        source = (
            f"SELECT * FROM main.{table} WHERE _resource = {_literal(record['resource_id'])}"
            if materialized
            else _source(root, record, table)
        )
        db.execute(f"CREATE OR REPLACE VIEW {schema}.{table} AS {source}")


def _configure(db, root: Path, memory_limit: str) -> None:
    db.execute("SET memory_limit = ?", [memory_limit])
    db.execute("SET temp_directory = ?", [str(root / ".duckdb-tmp")])


def view_lookups(resources: list[str]) -> list[str]:
    """Lookups over the per-resource views: a UNION ALL with one branch per resource,
    each reading only the row groups that hold the wanted values."""
    statements = []
    for table, column, kind in CONTRACT["lookups"]:
        parameter, wanted, condition = _lookup_parts(kind, column)
        parts = [
            f"SELECT * FROM {_name(r)}.{table} WHERE " + condition.replace("{r}", r)
            for r in sorted(resources)
        ]
        statements.append(
            f"CREATE OR REPLACE MACRO op_{table}_{kind}({parameter}) AS TABLE "
            f"WITH wanted AS MATERIALIZED ({wanted}) " + " UNION ALL BY NAME ".join(parts)
        )
    return statements


def _lookup_parts(kind: str, column: str) -> tuple[str, str, str]:
    # The argument is evaluated once, into ``wanted``: a macro argument appears wherever
    # its parameter does, and a list read from the views over every resource, repeated
    # in every resource's branch, takes seconds to plan.
    if kind == "at":
        return (
            "pairs",
            "SELECT p.r AS r, p.v AS v FROM (SELECT unnest(pairs) AS p)",
            f"{column} IN (SELECT v FROM wanted WHERE r = '{{r}}')",
        )
    return "vals", "SELECT unnest(vals) AS v", f"{column} IN (SELECT v FROM wanted)"


def native_lookups() -> list[str]:
    """Lookups over native tables: one filter on one table, by global id or by key."""
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


def _query_layer(db, records: list[ResourceRecord], materialized: bool) -> None:
    """Biolink names, the lookups and the macros; tables before what reads them."""
    db.execute("CREATE OR REPLACE TABLE op_names (kind VARCHAR, alias VARCHAR, name VARCHAR)")
    db.execute(
        "INSERT INTO op_names VALUES "
        + ", ".join(
            f"({_literal(k)}, {_literal(a)}, {_literal(n)})" for k, a, n in CONTRACT["biolink_names"]
        )
    )
    present = {
        row[0] for row in db.execute("SELECT table_name FROM duckdb_tables() WHERE schema_name = 'main'").fetchall()
    } | {row[0] for row in db.execute("SELECT view_name FROM duckdb_views() WHERE schema_name = 'main'").fetchall()}
    lookups = (
        native_lookups() if materialized else view_lookups([r["resource_id"] for r in records])
    )
    for statement in lookups:
        table = re.match(r"CREATE OR REPLACE MACRO op_(\w+)_(at|with)\(", statement).group(1)
        if table in present:
            db.execute(statement)
    for name in CONTRACT["macro_files"]:
        db.execute(_SQL.joinpath(name).read_text())
    db.execute("UPDATE database_info SET query_layer = ?", [QUERY_LAYER])


def _index(records: list[ResourceRecord], existing: dict[str, int]) -> dict[str, int]:
    """Each resource's index for global ids: kept for installed ones, new ones after."""
    index = {name: i for name, i in existing.items() if name in {r["resource_id"] for r in records}}
    following = max(existing.values(), default=-1) + 1
    for record in sorted(records, key=lambda r: r["resource_id"]):
        if record["resource_id"] not in index:
            index[record["resource_id"]] = following
            following += 1
    return index


def _write_index(db, index: dict[str, int]) -> None:
    db.execute("CREATE OR REPLACE TABLE op_resource_index (_resource VARCHAR, idx BIGINT)")
    for name, i in sorted(index.items()):
        db.execute("INSERT INTO op_resource_index VALUES (?, ?)", [name, i])


def _rebuild(root: Path, records: list[ResourceRecord], materialized: bool, memory_limit: str):
    """Write the database anew and put it in place when complete."""
    staging = Path(tempfile.mkdtemp(dir=root, prefix=".database-")) / DATABASE
    try:
        db = duckdb.connect()
        try:
            db.execute(
                f"ATTACH {_literal(str(staging))} AS omnipath "
                f"(STORAGE_VERSION {_literal(CONTRACT['storage_version'])})"
            )
            db.execute("USE omnipath")
            _configure(db, root, memory_limit)
            ordered = sorted(records, key=lambda r: r["resource_id"])
            present = [t for t in TABLES + SERVING if any(t in _tables(r) for r in ordered)]
            if materialized:
                index = _index(ordered, {})
                _write_index(db, index)
                for table in present:
                    having = [r for r in ordered if table in _tables(r)]
                    # The columns of every resource (union by name), then each resource's
                    # rows in its file's order, so rows follow their global ids.
                    union = " UNION ALL BY NAME ".join(
                        _source(root, r, table, index[r["resource_id"]]) for r in having
                    )
                    db.execute(f"CREATE TABLE {table} AS SELECT * FROM ({union}) LIMIT 0")
                    for record in having:
                        db.execute(
                            f"INSERT INTO {table} BY NAME "
                            + _source(root, record, table, index[record["resource_id"]])
                        )
            for record in ordered:
                _resource_schema(db, root, record, materialized)
            if not materialized:
                for table in present:
                    union = " UNION ALL BY NAME ".join(
                        f"SELECT * FROM {_name(r['resource_id'])}.{table}"
                        for r in ordered
                        if table in _tables(r)
                    )
                    db.execute(f"CREATE VIEW {table} AS {union}")
            _describe(db, root, materialized, records)
            _query_layer(db, records, materialized)
            db.execute("USE memory")
            db.execute("DETACH omnipath")
        finally:
            db.close()
        os.replace(staging, root / DATABASE)
    finally:
        for leftover in staging.parent.iterdir():
            leftover.unlink()
        staging.parent.rmdir()


def _update(root: Path, records: list[ResourceRecord], memory_limit: str) -> bool:
    """Replace the rows of added, changed and removed resources in native tables, in
    one transaction. False when the database has to be written anew instead."""
    try:
        db = duckdb.connect(str(root / DATABASE))
    except duckdb.Error:
        return False
    try:
        _configure(db, root, memory_limit)
        if db.execute("SELECT materialized FROM database_info").fetchone() != (True,):
            return False
        current = set(db.execute("SELECT * FROM resources").fetchall())
        stale = sorted({row[0] for row in current ^ set(_state(records))})
        by_name = {r["resource_id"]: r for r in records}
        existing = {row[0] for row in db.execute("SHOW TABLES").fetchall()}
        index = _index(records, dict(db.execute("SELECT * FROM op_resource_index").fetchall()))
        db.execute("BEGIN")
        _write_index(db, index)
        for table in TABLES + SERVING:
            fresh = [by_name[n] for n in stale if n in by_name and table in _tables(by_name[n])]
            if table not in existing:
                if fresh:
                    return False
                continue
            db.execute(f"DELETE FROM {table} WHERE _resource IN (SELECT unnest(?))", [stale])
            for record in fresh:
                db.execute(
                    f"INSERT INTO {table} BY NAME "
                    + _source(root, record, table, index[record["resource_id"]])
                )
        for name in stale:
            if name in by_name:
                _resource_schema(db, root, by_name[name], True)
            else:
                db.execute(f"DROP SCHEMA IF EXISTS {_name(name)} CASCADE")
        _describe(db, root, True, records)
        _query_layer(db, records, True)
        db.execute("COMMIT")
        return True
    except duckdb.Error:
        # A resource with a new column or type, or a database without global ids: the
        # tables take their schema anew.
        return False
    finally:
        db.close()


def _info(root: Path) -> tuple | None:
    """(data_dir, materialized, query layer, resources) of the database, or None."""
    path = root / DATABASE
    if not path.is_file():
        return None
    try:
        db = duckdb.connect(str(path), read_only=True)
    except duckdb.Error:
        return None
    try:
        columns = {row[0] for row in db.execute("DESCRIBE database_info").fetchall()}
        layer = "query_layer" if "query_layer" in columns else "NULL"
        data_dir, materialized, query_layer = db.execute(
            f"SELECT data_dir, materialized, {layer} FROM database_info"
        ).fetchone()
        state = db.execute("SELECT * FROM resources ORDER BY ALL").fetchall()
    except duckdb.Error:
        return None
    finally:
        db.close()
    return data_dir, materialized, query_layer, state


def _tables_in_sync(info: tuple | None, root: Path, installation: Installation) -> bool:
    if info is None:
        return False
    data_dir, materialized, _, state = info
    return (
        materialized == installation["materialize"]
        # Views name their files; native tables do not depend on the location.
        and (materialized or data_dir == str(root))
        and state == _state(installation["resources"])
    )


def in_sync(root: Path, installation: Installation) -> bool:
    """Whether the database describes exactly the installed resources, in the
    installation's mode, at this location, with this client's query layer."""
    info = _info(root)
    return _tables_in_sync(info, root, installation) and info[2] == QUERY_LAYER


def synchronize(root: Path, installation: Installation, memory_limit: str = "2GB") -> None:
    """Bring the database in line with the installation; nothing when it already is.

    When only the query layer differs (another client version wrote the database), the
    query layer alone is written again.
    """
    info = _info(root)
    if _tables_in_sync(info, root, installation) and info[2] == QUERY_LAYER:
        return
    records = installation["resources"]
    try:
        if _tables_in_sync(info, root, installation) and _layer_only(root, installation):
            return
        if not (installation["materialize"] and _update(root, records, memory_limit)):
            _rebuild(root, records, installation["materialize"], memory_limit)
    finally:
        temporary = root / ".duckdb-tmp"
        if temporary.is_dir() and not any(temporary.iterdir()):
            temporary.rmdir()


def _layer_only(root: Path, installation: Installation) -> bool:
    """Write the query layer of an otherwise current database; False when it cannot be
    (native tables without global ids, from a client before them)."""
    db = duckdb.connect(str(root / DATABASE))
    try:
        if installation["materialize"]:
            tables = {row[0] for row in db.execute("SHOW TABLES").fetchall()}
            if "op_resource_index" not in tables:
                return False
        if "query_layer" not in {row[0] for row in db.execute("DESCRIBE database_info").fetchall()}:
            db.execute("ALTER TABLE database_info ADD COLUMN query_layer VARCHAR")
        db.execute("BEGIN")
        _query_layer(db, installation["resources"], installation["materialize"])
        db.execute("COMMIT")
        return True
    except duckdb.Error:
        return False
    finally:
        db.close()
