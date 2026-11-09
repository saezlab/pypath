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
from omnipath_api.queries.store.query_cache import QueryCache
from omnipath_api.queries.store.query_pool import QueryPool
from omnipath_api.queries.store.connection import get_connection, sql_literal
from omnipath_core.measurements import QUANTITY_STRUCT
from omnipath_core.schema import PUBLISHED_TABLES, RELATION_QUALIFIERS, SERVING_TABLES
from omnipath_core.versioning import RESOURCE_FILES
from omnipath_client.contracts import CONTRACT

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
_QUANTITY = [field.name for field in QUANTITY_STRUCT]
# {table: (its own id column, the global id column)} of the release databases.
GLOBAL_IDS: dict[str, tuple[str, str]] = {t: tuple(c) for t, c in CONTRACT["global_ids"].items()}
# Up to this many ids go into the SQL as literals (and through the indexes); more are a
# join, which scans.
LITERAL_LIMIT = 100_000
_HELPERS = ("_resource", "_resource_version", "_rid", "_eid")
# Tables the database also holds sorted across every resource:
# {table: (its sorted copy, the columns it is sorted by)}. Reads by those columns (and
# scans) use the copy; reads by global id the table itself, which is in id order.
SORTED_COPIES = {"entity": ("entity_by_key", ("entity_key", "key_rank"))}
# The entity columns a relation's endpoint shows (shape_entity_summary).
ENDPOINT_FIELDS = (
    "entity_id, entity_key, namespace, identifier, label, entity_type, taxon, has_hierarchy, "
    "parent_count, child_count, reference_entity_key, gene_reference_keys, relation_count"
)
_HEX_KEY = re.compile(r"[0-9a-f]{64}")
_NUMERIC = re.compile(r"[0-9]+(?:\.[0-9]+)*")


def _nest_quantity(row):
    """A table row's ``quantity_*`` columns as one ``quantity`` struct (None when empty)."""
    if "quantity_" + _QUANTITY[0] not in row:
        return row
    quantity = {name: row.pop("quantity_" + name) for name in _QUANTITY}
    row["quantity"] = quantity if any(v is not None for v in quantity.values()) else None
    return row


def _version_key(version: str) -> tuple[int, ...]:
    return tuple(map(int, version.split(".")))


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
    _ENTITY_KEY_RE = re.compile(r"^[0-9a-fA-F]{64}$")
    _ENTITY_SCALAR_COLS = (
        "entity_key, entity_type, namespace, identifier, taxon, label, "
        "has_hierarchy, parent_count, child_count, reference_entity_key, gene_reference_keys"
    )

    def __init__(
        self,
        data_root: str | Path = "data",
        *,
        public_data_url: str | None = None,
        memory_limit: str = "1GB",
        threads: int = 2,
    ) -> None:
        """``memory_limit`` and ``threads`` apply to each query connection unless the
        OMNIPATH_DUCKDB_MEMORY_LIMIT and OMNIPATH_DUCKDB_THREADS variables set them;
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
        self._taxonomy_scope: ContextVar[dict[str, str] | None] = ContextVar(
            "taxonomy_scope", default=None
        )
        self._db_local = threading.local()
        self._connections = []
        self._connections_lock = threading.Lock()
        self._db_generation = 0
        self._inventory_lock = threading.Lock()
        self._inventory_fingerprint: tuple | None = None
        self._state = ServedState((), MappingProxyType({}), None)
        self._facet_cache = QueryCache()
        # Page order of endpoint selections: a gene's relations are sorted once for all pages.
        self._order_cache = QueryCache(max_entries=256, max_bytes=64 * 1024**2)
        self._detail_cache = QueryCache(max_entries=256, max_bytes=32 * 1024**2)
        self._relationship_cache = QueryCache(max_entries=128, max_bytes=32 * 1024**2)
        self._example_cache = QueryCache(max_entries=16, max_bytes=1024**2)
        self._catalog_cache = QueryCache(max_entries=8, max_bytes=64 * 1024**2)
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
        taxonomy_token = None
        try:
            taxonomy_token = self._taxonomy_scope.set(self._resolve_taxonomy_names())
            yield version
        finally:
            if taxonomy_token is not None:
                self._taxonomy_scope.reset(taxonomy_token)
            self._release_scope.reset(token)
            self._served_scope.reset(served_token)

    # Connections -----------------------------------------------------------------

    def _open_query_connection(self) -> duckdb.DuckDBPyConnection:
        connection = get_connection(memory_limit=self._memory_limit, threads=self._threads)
        try:
            Path(self._temp_directory).mkdir(parents=True, exist_ok=True)
            connection.execute("SET temp_directory = ?", [self._temp_directory])
            # Lists of up to LITERAL_LIMIT ids read through the indexes.
            connection.execute(f"SET index_scan_max_count = {LITERAL_LIMIT}")
            connection.execute("SET index_scan_percentage = 0.5")
            # A query reading one table twice (two lookups on entity_term, say) would
            # have its table materialized whole, once, before any filter: 13 s instead
            # of 0.3 s for a gene's entity keys in a full release.
            connection.execute("SET disabled_optimizers = 'common_subplan'")
            for name, value in CONTRACT["settings"].items():
                connection.execute(f"SET {name} = {int(value)}")
            for served in self._state.releases.values():
                connection.execute(
                    f"ATTACH {sql_literal(str(served.path))} AS {served.alias} (READ_ONLY)"
                )
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
        self._facet_cache.clear()
        self._order_cache.clear()
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

    # Resources of the release ----------------------------------------------------

    @staticmethod
    def _is_entity_key(token: str) -> bool:
        return bool(ReleaseQueries._ENTITY_KEY_RE.fullmatch(token))

    def _split_entity_tokens(self, tokens: list[str]) -> tuple[list[str], list[str]]:
        keys: list[str] = []
        other: list[str] = []
        for token in tokens:
            (keys if self._is_entity_key(token) else other).append(token)
        return keys, other

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

    # Tables and lookups ----------------------------------------------------------

    def _in_database(self, table: str) -> bool:
        return table in self._served().columns

    def _derived(self, table: str, resources: list[str] | None = None) -> tuple[str, str] | None:
        """(source, condition on its ``resource``) of a derived table over the selected
        resources, or None when the database has no such table. Derived tables hold
        every resource's rows, so a selection or the query policy is a filter."""
        if not self._in_database(table):
            return None
        infos = self._selected_resource_infos(resources)
        if len(infos) == len(self.resources):
            condition = "TRUE"
        elif not infos:
            condition = "FALSE"
        else:
            condition = f"resource IN ({', '.join(sql_literal(info['key']) for info in infos)})"
        return self._source(table), condition

    def _physical(self, table: str, column: str | None = None) -> str:
        """The database table a read of ``table`` (by ``column``, or a scan) uses."""
        copy = SORTED_COPIES.get(table)
        if copy is not None and copy[0] in self._served().columns:
            if column is None or column in copy[1]:
                return copy[0]
        return table

    def _projection(self, table: str, column: str | None = None) -> str:
        """The table's columns as the resource files have them, then ``resource``."""
        columns = self._served().columns[self._physical(table, column)]
        # A derived table's stored ``resource`` is computed here as for every table; an
        # entity row's key rank (in the sorted copy, standing for its key) is no field.
        helpers = (*_HELPERS, "resource", *(("key_rank",) if table == "entity" else ()))
        hidden = [c for c in helpers if c in columns]
        return (
            f"* EXCLUDE ({', '.join(hidden)}), _resource || '/' || _resource_version AS resource"
        )

    def _source(self, table: str, column: str | None = None) -> str:
        served = self._served()
        return f"{served.alias}.main.{self._physical(table, column)}"

    def _sorted_by(self, table: str, column: str) -> bool:
        """Whether a read of ``table`` by ``column`` reads rows sorted by it over every
        resource: a sorted copy, or the table itself rewritten so (entity_term with its
        entities' key ranks)."""
        copy = SORTED_COPIES.get(table)
        if copy is not None and column in copy[1] and self._physical(table, column) == copy[0]:
            return True
        return (table, column) == ("entity_term", "term") and self._has_term_keys()

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
    def _in(column, values):
        """``column IN values`` that keeps row-group skipping on a sorted column: a
        literal list when short, else the values' range around an unnested list."""
        values = list(dict.fromkeys(values))
        if len(values) <= 1000:
            return f"{column} IN ({','.join('?' for _ in values)})", values
        return (
            f"{column} BETWEEN ? AND ? AND {column} IN (SELECT unnest(?))",
            [min(values), max(values), values],
        )

    def _by_resource(self, pairs):
        """{resource: values} of (resource, value) pairs, for resources of the release."""
        values: dict[str, list] = {}
        resources = self.resources
        for resource, value in pairs:
            if resource in resources:
                values.setdefault(resource, []).append(value)
        return dict(sorted(values.items()))

    def _global_ids(self, table: str, by_resource: dict[str, list]) -> list[int]:
        index = self._served().index
        return sorted(
            {
                (index[resource.split("/", 1)[0]] << 32) + int(value)
                for resource, values in by_resource.items()
                for value in values
            }
        )

    def _matching(self, table: str, column: str, by_resource: dict[str, list]) -> tuple[str, list]:
        """The condition selecting the rows of (resource, value) pairs, with the global
        ids or keys written into it as literals (the indexes apply only to those)."""
        ids = GLOBAL_IDS.get(table)
        if ids is not None and column == ids[0]:
            global_ids = self._global_ids(table, by_resource)
            if len(global_ids) <= LITERAL_LIMIT:
                return f"{ids[1]} IN ({','.join(map(str, global_ids))})", []
            return f"{ids[1]} IN (SELECT unnest(?::BIGINT[]))", [global_ids]
        literal = lambda v: (  # noqa: E731
            str(int(v)) if isinstance(v, int) and not isinstance(v, bool) else sql_literal(str(v))
        )
        if self._sorted_by(table, column):
            # A copy sorted by the column: its values are a few row groups; the pairs
            # keep each resource's own rows.
            values = sorted({str(v) for vs in by_resource.values() for v in vs})
            pairs = ", ".join(
                f"({sql_literal(resource.split('/', 1)[0])}, {literal(v)})"
                for resource, vs in by_resource.items()
                for v in dict.fromkeys(vs)
            )
            return (
                f"{column} IN ({', '.join(literal(v) for v in values)}) "
                f"AND (_resource, {column}) IN (SELECT * FROM (VALUES {pairs}))",
                [],
            )
        if ids is not None and len(by_resource) > 3:
            # By key from more than a few resources, in two steps: the global ids of the wanted (resource, key) pairs
            # from a narrow read of the key column, then the rows by those ids. One
            # condition per resource, OR-ed, would scan the whole table.
            wanted = {
                (resource.split("/", 1)[0], str(value))
                for resource, values in by_resource.items()
                for value in values
            }
            values = sorted({value for _, value in wanted})
            found = self._db.execute(
                f"SELECT _resource, {column}, {ids[1]} FROM {self._source(table)} "
                f"WHERE {column} IN ({', '.join(literal(v) for v in values)})"
            ).fetchall()
            global_ids = sorted({g for r, v, g in found if (r, str(v)) in wanted})
            if not global_ids:
                return "FALSE", []
            if len(global_ids) <= LITERAL_LIMIT:
                return f"{ids[1]} IN ({','.join(map(str, global_ids))})", []
            return f"{ids[1]} IN (SELECT unnest(?::BIGINT[]))", [global_ids]
        parts = []
        for resource, values in by_resource.items():
            literals = ", ".join(literal(v) for v in dict.fromkeys(values))
            parts.append(
                f"(_resource = {sql_literal(resource.split('/', 1)[0])} AND {column} IN ({literals}))"
            )
        return " OR ".join(parts), []

    def _rows_sql(
        self, table: str, column: str, by_resource: dict[str, list], fields: str = "*"
    ) -> tuple[str, list]:
        """A MATERIALIZED read of the rows of (resource, value) pairs: further
        conditions in the same WHERE would make DuckDB scan instead of using an index.
        With plain column names as ``fields``, only those columns (and the row order's)
        are read."""
        condition, params = self._matching(table, column, by_resource)
        projection = self._projection(table, column)
        names = [name.strip() for name in fields.split(",")]
        if fields.strip() != "*" and all(re.fullmatch(r"[a-z_][a-z0-9_]*", n) for n in names):
            available = self._served().columns[self._physical(table, column)]
            ids = GLOBAL_IDS.get(table)
            ordered = [c for c in (ids[0] if ids else None, "ordinal") if c and c in available]
            wanted = list(dict.fromkeys([*names, *ordered]))
            projection = ", ".join(wanted) + ", _resource || '/' || _resource_version AS resource"
        return (
            f"WITH __rows AS MATERIALIZED (SELECT {projection} FROM "
            f"{self._source(table, column)} WHERE {condition})",
            params,
        )

    def _lookup_sql(self, table, column, pairs, where="TRUE", params=()):
        """SQL and parameters of the ``table`` rows (with ``resource``) for (resource,
        value) pairs."""
        by_resource = self._by_resource(pairs)
        if not by_resource:
            return f"SELECT * FROM {self._table(table, [])} WHERE FALSE", []
        if not self._in_database(table):
            return self._file_lookup_sql(table, column, by_resource, where, params)
        cte, values = self._rows_sql(table, column, by_resource)
        return f"{cte} SELECT * FROM __rows WHERE {where}", [*values, *params]

    def _file_lookup_sql(self, table, column, by_resource, where, params):
        parts, values = [], []
        for resource, items in by_resource.items():
            condition, items = self._in(column, items)
            path = sql_literal(str(self.resources[resource]["tables"][table]))
            parts.append(
                f"SELECT *, ?::VARCHAR AS resource FROM read_parquet({path}) "
                f"WHERE {condition} AND {where}"
            )
            values += [resource, *items, *params]
        return " UNION ALL ".join(parts), values

    def _lookup(self, table, column, pairs, fields="*", extra="", params=(), *, narrow=False):
        """Rows of ``table`` for (resource, value) pairs, in resource order:
        ``extra`` (conditions, order, limit) applies to each resource's rows.

        Without a LIMIT in ``extra`` this is one query over every resource; with one,
        one query per resource. ``narrow`` reads only the ``fields`` (plain column
        names) and the row order's columns, for an ``extra`` that refers to no others.
        """
        by_resource = self._by_resource(pairs)
        if by_resource and self._in_database(table) and not re.search(r"\bLIMIT\b", extra, re.I):
            cte, values = self._rows_sql(table, column, by_resource, fields if narrow else "*")
            selected = "*" if fields.strip() == "*" else f"{fields}, resource"
            order = re.search(r"\bORDER\s+BY\b", extra, re.I)
            # Each resource's rows in its file's order, as one read of its file gave them.
            ids = GLOBAL_IDS.get(table)
            native = ", ".join(
                c for c in (ids[0] if ids else None, "ordinal")
                if c and c in self._served().columns[table]
            )
            extra = (
                extra[: order.start()] + "ORDER BY resource, " + extra[order.end():]
                if order
                else f"{extra} ORDER BY resource" + (f", {native}" if native else "")
            )
            return self._fetch_dicts(
                f"{cte} SELECT {selected} FROM __rows WHERE TRUE {extra}", [*values, *params]
            )
        rows = []
        for resource, items in by_resource.items():
            if not self._in_database(table):
                condition, items = self._in(column, items)
                path = sql_literal(str(self.resources[resource]["tables"][table]))
                rows += self._fetch_dicts(
                    f"SELECT {fields}, ? AS resource FROM read_parquet({path}) WHERE {condition} {extra}",
                    [resource, *items, *params],
                )
                continue
            cte, values = self._rows_sql(table, column, {resource: items}, fields if narrow else "*")
            selected = "*" if fields.strip() == "*" else f"{fields}, resource"
            rows += self._fetch_dicts(
                f"{cte} SELECT {selected} FROM __rows WHERE TRUE {extra}", [*values, *params]
            )
        return rows

    def _relations_by_keys(self, keys, resources=None) -> list[dict[str, Any]]:
        """(resource, relation_id, relation_key) of the selected resources' relations
        with these keys, ordered by resource, through the database's key table."""
        keys = list(dict.fromkeys(str(k) for k in keys))
        if not keys:
            return []
        served = self._served()
        selected = {info["key"] for info in self._selected_resource_infos(resources)}
        literals = ", ".join(sql_literal(k) for k in keys)
        found = []
        if "key_prefix" in served.columns["relation_key_lookup"]:
            # By the keys' first 60 bits, then the full key on the relations found; a
            # key that is no sha256 hex is looked for in the relation table itself.
            prefixes = sorted({int(k[:15], 16) for k in keys if _HEX_KEY.fullmatch(k)})
            ids = [
                r[0] for r in self._db.execute(
                    f"SELECT _rid FROM {served.alias}.main.relation_key_lookup "
                    f"WHERE key_prefix IN ({','.join(map(str, prefixes)) or 'NULL'})"
                ).fetchall()
            ]
            rows = self._db.execute(
                f"WITH r AS MATERIALIZED (SELECT relation_key, _rid FROM {served.alias}.main.relation "
                f"WHERE _rid IN ({','.join(map(str, ids)) or 'NULL'})) "
                f"SELECT relation_key, _rid FROM r WHERE relation_key IN ({literals})"
            ).fetchall()
            other = [k for k in keys if not _HEX_KEY.fullmatch(k)]
            if other:
                rows += self._db.execute(
                    f"SELECT relation_key, _rid FROM {served.alias}.main.relation "
                    f"WHERE relation_key IN ({', '.join(sql_literal(k) for k in other)})"
                ).fetchall()
        else:
            rows = self._db.execute(
                f"SELECT relation_key, _rid FROM {served.alias}.main.relation_key_lookup "
                f"WHERE relation_key IN ({literals})"
            ).fetchall()
        for relation_key, global_id in rows:
            resource = served.by_index.get(global_id >> 32)
            if resource in selected:
                found.append(
                    {"resource": resource, "relation_id": global_id & 0xFFFFFFFF, "relation_key": relation_key}
                )
        return sorted(found, key=lambda r: (r["resource"], r["relation_id"]))

    def _children(self, table, parents, field, extra=""):
        """Attach each parent's rows of a child table, in ordinal order, as ``field``.

        Parents are entity or relation rows with ``resource`` and their id. Child rows
        drop the parent id and ordinal; a flattened quantity becomes one struct again.
        """
        column = "relation_id" if table.startswith("relation") else "entity_id"
        found = defaultdict(list)
        if parents:
            for row in self._lookup(
                table,
                column,
                [(p["resource"], p[column]) for p in parents],
                extra=f"{extra} ORDER BY {column}, ordinal",
            ):
                found[(row.pop("resource"), row.pop(column))].append(_nest_quantity(row))
                row.pop("ordinal")
        for parent in parents:
            parent[field] = found.get((parent["resource"], parent[column]), [])
        return parents

    def _with_entity_children(self, rows, *, evidence=False):
        """Entity rows with their identifiers and annotations (and evidence) nested."""
        self._children("entity_identifier", rows, "identifiers")
        self._children("entity_annotation", rows, "annotations")
        if evidence:
            self._children("entity_evidence", rows, "evidence")
        return rows

    def _entity_rows(self, where, params=(), resources=None, *, nested=True, evidence=False):
        """Every resource copy of the matching entities, ordered by resource."""
        rows = self._fetch_dicts(
            f"SELECT * FROM {self._table('entity', resources)} WHERE {where} ORDER BY resource, entity_id",
            list(params),
        )
        return self._with_entity_children(rows, evidence=evidence) if nested else rows

    def _endpoint_rows(self, relations, nested=False):
        """{(resource, entity key): entity row} for both endpoints of relation rows.

        A relation's endpoints are entities of its own resource, read there by key.
        """
        pairs = [
            (r["resource"], r[side])
            for r in relations
            for side in ("subject_entity_key", "object_entity_key")
            if r.get(side)
        ]
        # What an endpoint shows (shape_entity_summary), not the whole row.
        rows = self._lookup("entity", "entity_key", pairs, ENDPOINT_FIELDS, narrow=True) if pairs else []
        if nested:
            self._with_entity_children(rows)
        return {(r["resource"], r["entity_key"]): r for r in rows}

    @staticmethod
    def _with_qualifiers(row):
        """Relation-scope qualifier annotations from the relation's qualifier columns."""
        row["annotations"] = [
            dict(term=term, value=value, scope="relation")
            for term in RELATION_QUALIFIERS
            for value in row.get(term) or []
        ]
        return row

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

    @staticmethod
    def _as_list(value: Any) -> list[Any]:
        if value is None:
            return []
        if isinstance(value, (list, tuple, set)):
            return [item for item in value if item is not None and str(item).strip() != ""]
        return [value]

    def _fetch_dicts(self, sql: str, params: list[Any] | None = None) -> list[dict[str, Any]]:
        cursor = self._db.execute(sql, params or [])
        cols = [col[0] for col in cursor.description]
        return [self._jsonify(dict(zip(cols, row))) for row in cursor.fetchall()]

    # Downloads are the published tables; serving tables are lookup indexes. Raw source
    # payloads are downloaded on request.
    RESOURCE_PARQUET_FILES = RESOURCE_FILES
    PUBLIC_DOWNLOAD_FILES = tuple(f for f in RESOURCE_FILES if f != PAYLOAD_FILE)
