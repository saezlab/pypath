"""The API's queries over the release databases of a serving root, with DuckDB.

Every published release served by the API has one DuckDB database,
``<data_root>/serving/<version>.duckdb``, written by ``omnipath_client.database`` from
the release's resource files (see ``omnipath_api.release_database``). A request reads
one release: the one it names, or by default (``latest``) the newest release with a
database. Each query connection attaches every served database read-only.

The databases hold one table per published and lookup table over all of a release's
resources, with ``_resource``/``_resource_version`` and a global id (``_rid``,
``_eid``: the resource's index times 2^32 plus its own id). Queries see the tables as
before, with ``resource`` = ``source/version``; rows are read by id through the global
ids, written into the SQL as literal lists, so that DuckDB uses the databases' indexes.
Raw evidence payloads are not in the databases and stay in the resources' files.
"""

from __future__ import annotations

from collections import defaultdict
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
import json
import logging
import os
from pathlib import Path
import re
import tempfile
import threading
from types import MappingProxyType
from typing import Any
from urllib.parse import urlsplit

import duckdb

from omnipath_api.queries.store.inventory import ReleaseStore
from omnipath_api.queries.store.query_cache import QueryCache, StepCache
from omnipath_api.queries.store.query_pool import QueryPool
from omnipath_api.queries.store.connection import get_connection, sql_literal
from omnipath_core.schema import PUBLISHED_TABLES, SERVING_TABLES
from omnipath_core.versioning import RESOURCE_FILES
from omnipath_client import steps
from omnipath_client.steps import apply_settings

from omnipath_api.queries.services.resources import ResourcesQueries
from omnipath_api.queries.services.relations import RelationsQueries
from omnipath_api.queries.services.entities import EntitiesQueries
from omnipath_api.queries.services.facets import FacetsQueries
from omnipath_api.queries.services.export import ExportQueries
from omnipath_api.queries.services.evidence import EvidenceQueries
from omnipath_api.queries.services.selection import SelectionQueries
from omnipath_api.queries.services.ontology import OntologyQueries
from omnipath_api.queries.services.stats import StatsQueries

logger = logging.getLogger(__name__)

SERVING_DIR = "serving"
PAYLOAD_FILE = "evidence_payloads.parquet"
_NUMERIC = re.compile(r"[0-9]+(?:\.[0-9]+)*")


def _version_key(version: str) -> tuple[int, ...]:
    return tuple(map(int, version.split(".")))


_SIZE = re.compile(r"\s*(\d+(?:\.\d+)?)\s*([KMGT]?)(i?)B?\s*", re.I)


def _scaled(size: str, factor: int) -> str:
    """A DuckDB memory size times ``factor`` (``3500MB`` times 3 is ``10500MB``)."""
    match = _SIZE.fullmatch(size)
    if match is None:
        raise ValueError(f"Not a memory size: {size!r}")
    number, unit, binary = match.groups()
    unit = unit.upper()
    return f"{float(number) * factor:g}{unit}{'i' if binary else ''}B"


class ReleaseNotServed(LookupError):
    """A published release whose database is not (or no longer) served."""

    status_code = 410


class NothingServed(LookupError):
    """No release database is served at all."""

    status_code = 503


@dataclass(frozen=True)
class ServedRelease:
    """One served release database: what queries need to know about it."""

    version: str
    path: Path
    alias: str
    identity: tuple
    built_at: str
    # {"source/version": inventory record} of the release's resources, with the paths
    # of their files (downloads and raw payloads read those).
    resources: Any
    # {source: index of its global ids}
    index: Any
    # {table: column names} of the database's tables
    columns: Any
    by_index: Any = field(default=None)


@dataclass(frozen=True)
class ServedState:
    """The served release databases, replaced as a whole when one appears or goes."""

    fingerprint: tuple
    releases: Any
    latest: str | None

    @property
    def resources(self):
        """The resources of the newest served release (read-only)."""
        return self.releases[self.latest].resources if self.latest else MappingProxyType({})


def _resource_record(data_root: Path, resource: str, version: str) -> dict[str, Any]:
    directory = data_root / "resources" / resource / version
    tables = {
        name: directory / f"{name}.parquet"
        for name in (*PUBLISHED_TABLES, *SERVING_TABLES)
        if f"{name}.parquet" != PAYLOAD_FILE or (directory / PAYLOAD_FILE).is_file()
    }
    return {
        "resource": resource,
        "version": version,
        "key": f"{resource}/{version}",
        "tables": tables,
        "mtime": 0,
    }


class ReleaseQueries(
    ResourcesQueries,
    RelationsQueries,
    EntitiesQueries,
    FacetsQueries,
    ExportQueries,
    EvidenceQueries,
    SelectionQueries,
    OntologyQueries,
    StatsQueries,
):
    """Queries over the served release databases of a serving root."""

    ENGINE = "duckdb-release"

    def __init__(
        self,
        data_root: str | Path = "data",
        *,
        public_data_url: str | None = None,
        memory_limit: str = "1GB",
        threads: int = 2,
    ) -> None:
        """``memory_limit`` and ``threads`` are per query (OMNIPATH_DUCKDB_MEMORY_LIMIT and
        OMNIPATH_DUCKDB_THREADS override them): the one DuckDB database every query runs in
        gets them times the number of concurrent queries (OMNIPATH_QUERY_CONCURRENCY).
        OMNIPATH_DUCKDB_TEMP_DIRECTORY names a writable directory for spills."""
        self.data_root = Path(data_root).resolve()
        self._memory_limit, self._threads = memory_limit, threads
        self._temp_directory = os.getenv("OMNIPATH_DUCKDB_TEMP_DIRECTORY") or str(
            Path(tempfile.gettempdir()) / "omnipath-duckdb"
        )
        self.public_data_url = (
            public_data_url
            if public_data_url is not None
            else os.getenv("OMNIPATH_PUBLIC_DATA_URL", "")
        ).rstrip("/")
        if self.public_data_url:
            url = urlsplit(self.public_data_url)
            if (
                url.scheme not in {"http", "https"}
                or not url.netloc
                or url.query
                or url.fragment
                or url.username
            ):
                raise ValueError(
                    "OMNIPATH_PUBLIC_DATA_URL must be an HTTP(S) base URL without credentials, query or fragment"
                )
        self.read_only = False
        self.releases = ReleaseStore(self.data_root)
        self._release_scope: ContextVar[dict[str, Any] | None] = ContextVar(
            "release_scope", default=None
        )
        self._served_scope: ContextVar[ServedRelease | None] = ContextVar(
            "served_scope", default=None
        )
        self._db_local = threading.local()
        # The DuckDB database the queries run in (one cache and thread pool for all of
        # them), with its generation; each query gets a cursor on it.
        self._database: tuple[int, duckdb.DuckDBPyConnection] | None = None
        self._database_lock = threading.Lock()
        self._connections = []
        self._connections_lock = threading.Lock()
        self._db_generation = 0
        self._inventory_lock = threading.Lock()
        self._inventory_fingerprint: tuple | None = None
        self._state = ServedState((), MappingProxyType({}), None)
        self._example_cache = QueryCache(max_entries=16, max_bytes=1024**2)
        self._catalog_cache = QueryCache(max_entries=8, max_bytes=64 * 1024**2)
        # Answers of the step-file queries, by release, query and parameters.
        self._query_cache = QueryCache(max_entries=512, max_bytes=128 * 1024**2)
        # Results of the steps marked "cache" (a release's ontology terms, say).
        self._step_cache = StepCache()
        self._query_pool = QueryPool(
            self._open_query_connection,
            lambda: self._db_generation,
            size=int(os.getenv("OMNIPATH_QUERY_CONCURRENCY", "2")),
            timeout=float(os.getenv("OMNIPATH_QUERY_WAIT_SECONDS", "30")),
        )
        self._refresh_inventory_if_stale(force=True)

    # Served databases ----------------------------------------------------------------

    @property
    def serving_dir(self) -> Path:
        return self.data_root / SERVING_DIR

    def _scan_served(self) -> tuple:
        """(version, size, mtime) of every served release database, newest first."""
        if not self.serving_dir.is_dir():
            return ()
        items = []
        for path in self.serving_dir.glob("*.duckdb"):
            version = path.stem
            if not _NUMERIC.fullmatch(version) or not (self.data_root / "releases" / f"{version}.json").is_file():
                continue
            try:
                stat = path.stat()
            except OSError:
                continue
            items.append((version, stat.st_size, stat.st_mtime_ns))
        return tuple(sorted(items, key=lambda item: _version_key(item[0]), reverse=True))

    def _describe_served(self, version: str, identity: tuple) -> ServedRelease:
        path = self.serving_dir / f"{version}.duckdb"
        with duckdb.connect(str(path), read_only=True) as db:
            built_at = db.execute("SELECT built_at::VARCHAR FROM database_info").fetchone()[0]
            index = dict(db.execute("SELECT _resource, idx FROM op_resource_index").fetchall())
            rows = db.execute("SELECT resource_id, version FROM resources ORDER BY ALL").fetchall()
            columns = defaultdict(list)
            for table, column in db.execute(
                "SELECT table_name, column_name FROM duckdb_columns() "
                "WHERE schema_name = 'main' AND database_name = current_database() "
                "ORDER BY table_name, column_index"
            ).fetchall():
                columns[table].append(column)
        resources = {}
        for resource, resource_version in rows:
            record = _resource_record(self.data_root, resource, resource_version)
            resources[record["key"]] = MappingProxyType(record)
        return ServedRelease(
            version=version,
            path=path,
            alias="release_" + version.replace(".", "_"),
            identity=identity,
            built_at=str(built_at),
            resources=MappingProxyType(resources),
            index=MappingProxyType(index),
            columns=MappingProxyType({t: tuple(c) for t, c in columns.items()}),
            by_index=MappingProxyType(
                {i: f"{r}/{v}" for r, v in rows for name, i in index.items() if name == r}
            ),
        )

    def _refresh_inventory_if_stale(self, *, force: bool = False) -> bool:
        """Pick up release databases that appeared, went or were rewritten. Returns True
        when the served set changed."""
        fingerprint = self._scan_served()
        if not force and fingerprint == self._inventory_fingerprint:
            return False
        with self._inventory_lock:
            fingerprint = self._scan_served()
            if not force and fingerprint == self._inventory_fingerprint:
                return False
            previous = self._state.releases
            releases = {}
            for item in fingerprint:
                version = item[0]
                known = previous.get(version)
                try:
                    releases[version] = (
                        known if known is not None and known.identity == item
                        else self._describe_served(version, item)
                    )
                except duckdb.Error:
                    logger.exception("Skipping unreadable release database %s", version)
            latest = next(iter(releases), None)
            self._state = ServedState(fingerprint, MappingProxyType(releases), latest)
            if self._inventory_fingerprint is not None:
                self._reset_query_connection()
            self._inventory_fingerprint = fingerprint
            return True

    def reload_resources(self) -> None:
        """Look for new or removed release databases now."""
        self._refresh_inventory_if_stale(force=True)

    def inventory_snapshot(self) -> ServedState:
        self._refresh_inventory_if_stale()
        return self._state

    def _served(self) -> ServedRelease:
        """The release of the current scope, else the newest served one."""
        scoped = self._served_scope.get()
        if scoped is not None:
            return scoped
        self._refresh_inventory_if_stale()
        state = self._state
        if state.latest is None:
            raise NothingServed(
                f"No release database is served from {self.serving_dir}; publish a release "
                "and build its database (omnipath-api build-release-db)"
            )
        return state.releases[state.latest]

    @property
    def resources(self):
        return self._served().resources

    @property
    def _latest_by_resource(self):
        return {info["resource"]: info for info in self.resources.values()}

    @contextmanager
    def release_scope(self, version: str | None = None, *, snapshot: ServedState | None = None):
        """Pin a request to one served release database."""
        state = snapshot if isinstance(snapshot, ServedState) else self.inventory_snapshot()
        version = version or self.releases.default()
        if version in ("latest", "working"):
            if state.latest is None:
                raise NothingServed(
                    f"No release database is served from {self.serving_dir}"
                )
            version = state.latest
        served = state.releases.get(version)
        if served is None:
            self.releases.get(version)  # FileNotFoundError for an unknown release
            raise ReleaseNotServed(
                f"OmniPath {version} is no longer served by the API; its files remain "
                f"available for download (GET /resources?release={version})"
            )
        manifest = self.releases.get(version)
        served_token = self._served_scope.set(served)
        token = self._release_scope.set(manifest)
        try:
            yield version
        finally:
            self._release_scope.reset(token)
            self._served_scope.reset(served_token)

    # Step-file queries ------------------------------------------------------------

    def run_query(
        self,
        name: str,
        params: dict[str, Any] | None = None,
        resources: list[str] | None = None,
        *,
        output: str = "json",
    ):
        """A step-file query (``omnipath_client/sql/queries/<name>.sql``, the same the
        clients run) on the served release: JSON text (None for no answer) for
        ``output="json"``, an Arrow table for ``"table"``. The selected resources, less
        those the query policy disables, are its ``resources`` parameter."""
        served = self._served()
        query = steps.load(name)
        params = self._query_params(query, params, resources)

        def compute():
            connection = self._db
            connection.execute(f"USE {served.alias}")
            try:
                return steps.run(
                    connection, query, params, output=output,
                    cache=self._step_cache.get, scope=[served.identity, served.built_at],
                )
            except steps.QueryError as exc:
                raise ValueError(str(exc)) from None

        if output != "json":
            return compute()
        return self._query_cache.get(self._answer_key(name, params), compute)

    def _query_params(self, query, params, resources) -> dict[str, Any]:
        """The parameters with the selected resources, less those the query policy
        disables, as ``resources``."""
        params = dict(params or {})
        if "resources" in query.params:
            infos = self._selected_resource_infos(resources)
            selected = sorted({info["resource"] for info in infos})
            # NULL selects every resource; a name no resource has selects none.
            params["resources"] = None if len(infos) == len(self.resources) else selected or [""]
        return params

    def _answer_key(self, name: str, params: dict[str, Any]) -> str:
        return json.dumps([self._served().identity, name, params], sort_keys=True, default=str)

    def remember_answer(
        self, name: str, params: dict[str, Any], resources: list[str] | None, text: str
    ) -> None:
        """Cache a query's JSON answer that another query's answer holds (the first
        relationships page of an entity page, say), so asking for it reads nothing."""
        params = self._query_params(steps.load(name), params, resources)
        self._query_cache.get(self._answer_key(name, params), lambda: text)

    # Connections -----------------------------------------------------------------

    def _open_database(self) -> duckdb.DuckDBPyConnection:
        """The database the queries share, with the served releases attached: one buffer
        cache and one thread pool, sized for the concurrent queries. (A database per
        query kept three copies of the same blocks, each in a third of the memory.)"""
        size = self._query_pool.size
        memory = _scaled(os.getenv("OMNIPATH_DUCKDB_MEMORY_LIMIT", self._memory_limit), size)
        threads = int(os.getenv("OMNIPATH_DUCKDB_THREADS", str(self._threads))) * size
        database = get_connection(memory_limit=memory, threads=threads)
        try:
            Path(self._temp_directory).mkdir(parents=True, exist_ok=True)
            database.execute("SET temp_directory = ?", [self._temp_directory])
            for served in self._state.releases.values():
                database.execute(
                    f"ATTACH {sql_literal(str(served.path))} AS {served.alias} (READ_ONLY)"
                )
        except Exception:
            database.close()
            raise
        return database

    def _open_query_connection(self) -> duckdb.DuckDBPyConnection:
        """A connection (cursor) on the shared database, opened anew when the served
        releases change; cursors of the previous one keep it until they are discarded."""
        with self._database_lock:
            generation = self._db_generation
            if self._database is None or self._database[0] != generation:
                previous = self._database
                self._database = (generation, self._open_database())
                if previous is not None:
                    previous[1].close()
            connection = self._database[1].cursor()
        try:
            # The settings every query connection uses (the index reads, no
            # materializing of whole tables): the same as the clients'.
            connection.execute("SET arrow_large_buffer_size=true")
            apply_settings(connection)
        except Exception:
            connection.close()
            raise
        return connection

    @contextmanager
    def query_scope(self):
        """Lease one bounded database for a complete query, including exports."""
        if getattr(self._db_local, "query_connection", None) is not None:
            yield
            return
        self._refresh_inventory_if_stale()
        with self._query_pool.acquire() as connection:
            self._db_local.query_connection = connection
            try:
                yield
            finally:
                self._db_local.query_connection = None

    def _close_direct_connection(self, connection):
        try:
            connection.close()
        finally:
            with self._connections_lock:
                if connection in self._connections:
                    self._connections.remove(connection)

    @property
    def _db(self) -> duckdb.DuckDBPyConnection:
        """One DuckDB connection per thread so relation search and facets can overlap."""
        scoped = getattr(self._db_local, "query_connection", None)
        if scoped is not None:
            return scoped
        gen = self._db_generation
        conn: duckdb.DuckDBPyConnection | None = getattr(self._db_local, "conn", None)
        conn_gen = getattr(self._db_local, "gen", None)
        if conn is None or conn_gen != gen:
            if conn is not None:
                try:
                    self._close_direct_connection(conn)
                except duckdb.Error:
                    logger.exception("Failed to close stale DuckDB connection")
            conn = self._open_query_connection()
            with self._connections_lock:
                self._connections.append(conn)
            self._db_local.conn = conn
            self._db_local.gen = gen
        return conn

    def _reset_query_connection(self) -> None:
        self._db_generation += 1
        self._query_cache.clear()
        self._step_cache.clear()
        conn = getattr(self._db_local, "conn", None)
        if conn is not None:
            try:
                self._close_direct_connection(conn)
            except duckdb.Error:
                logger.exception("Failed to reset DuckDB connection")
            self._db_local.conn = None
            self._db_local.gen = None

    def close(self):
        """Release all connections after the HTTP worker has drained requests."""
        self._query_pool.close()
        for connection in tuple(self._connections):
            self._close_direct_connection(connection)
        with self._database_lock:
            if self._database is not None:
                self._database[1].close()
                self._database = None

    # Resources of the release ----------------------------------------------------

    def _visible_resources(self) -> dict[str, dict[str, Any]]:
        return self.resources

    def _visible_latest(self) -> dict[str, dict[str, Any]]:
        return self._latest_by_resource

    def _disabled_query_resources(self) -> set[str]:
        """Live serving policy, independent of immutable releases and downloads."""
        try:
            policy = json.loads((self.data_root / "query_policy.json").read_text())
        except FileNotFoundError:
            return set()
        disabled = policy.get("disabled_resources", [])
        if not isinstance(disabled, list) or any(not isinstance(item, str) for item in disabled):
            raise ValueError(
                "query_policy.json disabled_resources must be a list of resource slugs"
            )
        return set(disabled)

    def _selected_resource_infos(
        self, selected_resources: list[str] | None = None, *, include_disabled: bool = False
    ) -> list[dict[str, Any]]:
        visible = self._visible_resources()
        latest = self._visible_latest()
        disabled = set() if include_disabled else self._disabled_query_resources()
        if not visible:
            return []
        if not selected_resources:
            return [info for info in latest.values() if info["resource"] not in disabled]
        infos: list[dict[str, Any]] = []
        seen: set[str] = set()
        for res in selected_resources:
            info = visible.get(res) or latest.get(res)
            if info is None:
                raise ValueError(
                    f"Resource {res!r} is unavailable in the selected OmniPath release"
                )
            if info["resource"] in disabled:
                raise ValueError(f"Resource {res!r} is disabled for queries")
            marker = info["key"]
            if marker not in seen:
                infos.append(info)
                seen.add(marker)
        return infos

    def _table_paths(self, table: str, selected_resources: list[str] | None = None) -> dict:
        """{resource key (source/version): path} of one table's files for the selected
        resources that have it (downloads and raw payloads)."""
        return {
            info["key"]: str(info["tables"][table])
            for info in self._selected_resource_infos(selected_resources)
            if table in info["tables"]
        }

    # Tables ----------------------------------------------------------

    def _in_database(self, table: str) -> bool:
        return table in self._served().columns

    def _projection(self, table: str) -> str:
        """The table's columns as the resource files have them, then ``resource``."""
        columns = self._served().columns[table]
        # A derived table's stored ``resource`` is computed here as for every table; an
        # entity row's key rank (in the sorted copy, standing for its key) is no field.
        helpers = ("_resource", "_resource_version", "_rid", "_eid", "resource", "key_rank")
        hidden = [c for c in helpers if c in columns]
        return (
            f"* EXCLUDE ({', '.join(hidden)}), _resource || '/' || _resource_version AS resource"
        )

    def _source(self, table: str) -> str:
        return f"{self._served().alias}.main.{table}"

    def _file_table(self, table: str, selected_resources: list[str] | None = None) -> str:
        """A table that is only in the resources' files (raw payloads), read from them."""
        paths = self._table_paths(table, selected_resources)
        if not paths:
            return self._empty_table(table)
        listed = ", ".join(sql_literal(path) for path in paths.values())
        return (
            f"(SELECT * EXCLUDE (filename), regexp_extract(filename, '([^/]+/[^/]+)/[^/]+$', 1) "
            f"AS resource FROM read_parquet([{listed}], union_by_name=true, filename=true))"
        )

    @staticmethod
    def _empty_table(table: str) -> str:
        schema = {**PUBLISHED_TABLES, **SERVING_TABLES}[table]
        columns = ", ".join(f"NULL AS {name}" for name in schema.names)
        return f"(SELECT {columns}, NULL::VARCHAR AS resource WHERE FALSE)"

    def _table(self, table: str, selected_resources: list[str] | None = None) -> str:
        """SQL source of one table over the selected resources, with a ``resource`` column."""
        if not self._in_database(table):
            return self._file_table(table, selected_resources)
        infos = self._selected_resource_infos(selected_resources)
        if not infos:
            return self._empty_table(table)
        where = ""
        if len(infos) < len(self.resources):
            names = ", ".join(sql_literal(info["resource"]) for info in infos)
            where = f" WHERE _resource IN ({names})"
        return f"(SELECT {self._projection(table)} FROM {self._source(table)}{where})"

    @staticmethod
    def _jsonify(value: Any) -> Any:
        if value is None or isinstance(value, (str, int, float, bool)):
            return value
        if isinstance(value, dict):
            return {str(k): ReleaseQueries._jsonify(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [ReleaseQueries._jsonify(v) for v in value]
        try:
            import numpy as np

            if isinstance(value, np.ndarray):
                return ReleaseQueries._jsonify(value.tolist())
            if isinstance(value, np.integer):
                return int(value)
            if isinstance(value, np.floating):
                return float(value)
            if isinstance(value, np.bool_):
                return bool(value)
        except ModuleNotFoundError as exc:
            if exc.name != "numpy":
                raise
        return str(value)

    def _fetch_dicts(self, sql: str, params: list[Any] | None = None) -> list[dict[str, Any]]:
        cursor = self._db.execute(sql, params or [])
        cols = [col[0] for col in cursor.description]
        return [self._jsonify(dict(zip(cols, row))) for row in cursor.fetchall()]

    # Downloads are the published tables; serving tables are lookup indexes. Raw source
    # payloads are downloaded on request.
    RESOURCE_PARQUET_FILES = RESOURCE_FILES
    PUBLIC_DOWNLOAD_FILES = tuple(f for f in RESOURCE_FILES if f != PAYLOAD_FILE)
