"""The API's analytical queries over the resource Parquet files, with DuckDB.
Supports multi-resource union, predicate pushdown, faceted counts, and binary exports.
"""

from __future__ import annotations

import logging
import json
import os
import re
import threading
from collections import defaultdict
from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit


import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from omnipath_api.queries.store.inventory import ReleaseStore, InventorySnapshot
from omnipath_api.queries.store.query_cache import QueryCache
from omnipath_api.queries.store.query_pool import QueryPool
from omnipath_api.queries.store.connection import get_connection, sql_literal
from omnipath_core.measurements import QUANTITY_STRUCT
from omnipath_core.schema import PUBLISHED_TABLES, RELATION_QUALIFIERS, SERVING_TABLES
from omnipath_core.versioning import RESOURCE_FILES, SERVING_FILES


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

PAYLOAD_FILE = "evidence_payloads.parquet"
_QUANTITY = [field.name for field in QUANTITY_STRUCT]


def _nest_quantity(row):
    """A table row's ``quantity_*`` columns as one ``quantity`` struct (None when empty)."""
    if "quantity_" + _QUANTITY[0] not in row:
        return row
    quantity = {name: row.pop("quantity_" + name) for name in _QUANTITY}
    row["quantity"] = quantity if any(v is not None for v in quantity.values()) else None
    return row


class ParquetQueries(
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
    """Out-of-core streaming queries over a serving root, powered by DuckDB."""

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
        OMNIPATH_DUCKDB_MEMORY_LIMIT and OMNIPATH_DUCKDB_THREADS variables set them."""
        self.data_root = Path(data_root).resolve()
        self._memory_limit, self._threads = memory_limit, threads
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
        self._taxonomy_scope: ContextVar[dict[str, str] | None] = ContextVar(
            "taxonomy_scope", default=None
        )
        self._db_local = threading.local()
        self._connections = []
        self._connections_lock = threading.Lock()
        self._db_generation = 0
        self._inventory_lock = threading.Lock()
        self._inventory_fingerprint: tuple[tuple[Any, ...], ...] | None = None
        self._facet_cache = QueryCache()
        self._detail_cache = QueryCache(max_entries=256, max_bytes=32 * 1024**2)
        self._relationship_cache = QueryCache(max_entries=128, max_bytes=32 * 1024**2)
        self._example_cache = QueryCache(max_entries=16, max_bytes=1024**2)
        self._catalog_cache = QueryCache(max_entries=8, max_bytes=64 * 1024**2)
        self._inventory = InventorySnapshot({}, {}, ())
        self._inventory_scope: ContextVar[InventorySnapshot | None] = ContextVar(
            "inventory_scope", default=None
        )
        self._query_pool = QueryPool(
            self._open_query_connection,
            lambda: self._db_generation,
            size=int(os.getenv("OMNIPATH_QUERY_CONCURRENCY", "2")),
            timeout=float(os.getenv("OMNIPATH_QUERY_WAIT_SECONDS", "30")),
        )
        with self.query_scope():
            self._refresh_inventory_if_stale(force=True)

    @staticmethod
    def _is_entity_key(token: str) -> bool:
        return bool(ParquetQueries._ENTITY_KEY_RE.fullmatch(token))

    def _split_entity_tokens(self, tokens: list[str]) -> tuple[list[str], list[str]]:
        keys: list[str] = []
        other: list[str] = []
        for token in tokens:
            (keys if self._is_entity_key(token) else other).append(token)
        return keys, other

    @classmethod
    def _version_is_ready(cls, ver_dir: Path) -> bool:
        """A published version has every resource and serving table; raw source
        payloads are optional (a client downloads them only on request)."""
        return all(
            (ver_dir / name).is_file()
            for name in RESOURCE_FILES + SERVING_FILES
            if name != PAYLOAD_FILE
        )

    def _scan_inventory_fingerprint(self) -> tuple[tuple[Any, ...], ...]:
        """Cheap disk snapshot of complete resource versions (path, mtime, size)."""
        resources_dir = self.data_root / "resources"
        if not resources_dir.exists():
            return ()
        items: list[tuple[Any, ...]] = []
        for res_dir in resources_dir.iterdir():
            if not res_dir.is_dir():
                continue
            for ver_dir in res_dir.iterdir():
                if not ver_dir.is_dir() or not self._version_is_ready(ver_dir):
                    continue
                ent_file = ver_dir / "entity.parquet"
                rel_file = ver_dir / "relation.parquet"
                stats_file = ver_dir / "resolution_stats.json"
                try:
                    ent_stat = ent_file.stat()
                    rel_stat = rel_file.stat()
                    stats_mtime = stats_file.stat().st_mtime_ns if stats_file.is_file() else 0
                except OSError:
                    continue
                items.append(
                    (
                        f"{res_dir.name}/{ver_dir.name}",
                        ent_stat.st_mtime_ns,
                        ent_stat.st_size,
                        rel_stat.st_mtime_ns,
                        rel_stat.st_size,
                        stats_mtime,
                    )
                )
        return tuple(sorted(items))

    def _open_query_connection(self) -> duckdb.DuckDBPyConnection:
        return get_connection(memory_limit=self._memory_limit, threads=self._threads)

    @contextmanager
    def query_scope(self):
        """Lease one bounded database for a complete query, including exports."""
        if getattr(self._db_local, "query_connection", None) is not None:
            yield
            return
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
        conn = getattr(self._db_local, "conn", None)
        if conn is not None:
            try:
                self._close_direct_connection(conn)
            except duckdb.Error:
                logger.exception("Failed to reset DuckDB connection")
            self._db_local.conn = None
            self._db_local.gen = None

    def _refresh_inventory_if_stale(self, *, force: bool = False) -> bool:
        """Re-index parquet resources when the data directory changes. Returns True on reload."""
        if self._inventory_scope.get() is not None and not force:
            return False
        fingerprint = self._scan_inventory_fingerprint()
        if not force and fingerprint == self._inventory_fingerprint:
            return False
        with self._inventory_lock:
            fingerprint = self._scan_inventory_fingerprint()
            if not force and fingerprint == self._inventory_fingerprint:
                return False
            if self._inventory_fingerprint is not None:
                self._reset_query_connection()
            self._discover_inventory()
            self._inventory_fingerprint = fingerprint
            return True

    def _discover_inventory(self) -> None:
        """Find all available resource parquet files and release manifests."""
        resources = {}
        resources_dir = self.data_root / "resources"
        if resources_dir.exists():
            for res_dir in resources_dir.iterdir():
                if not res_dir.is_dir():
                    continue
                res_name = res_dir.name
                for ver_dir in res_dir.iterdir():
                    if not ver_dir.is_dir() or not self._version_is_ready(ver_dir):
                        continue
                    tables = {
                        name: ver_dir / f"{name}.parquet"
                        for name in (*PUBLISHED_TABLES, *SERVING_TABLES)
                        if f"{name}.parquet" != PAYLOAD_FILE or (ver_dir / PAYLOAD_FILE).is_file()
                    }
                    ent_file, rel_file = tables["entity"], tables["relation"]
                    key = f"{res_name}/{ver_dir.name}"
                    try:
                        rel_rows = pq.read_metadata(rel_file).num_rows
                    except (OSError, pa.ArrowException):
                        logger.exception("Skipping unreadable resource %s", ver_dir)
                        continue
                    try:
                        label_quality = float(
                            self._db.execute(
                                f"SELECT coalesce(avg((label IS DISTINCT FROM identifier)::int), 0) FROM read_parquet({sql_literal(str(ent_file))})"
                            ).fetchone()[0]
                            or 0
                        )
                    except duckdb.Error:
                        logger.warning("Cannot score entity labels in %s", ent_file, exc_info=True)
                        label_quality = 0.0
                    mtime = max(ent_file.stat().st_mtime, rel_file.stat().st_mtime)
                    resources[key] = {
                        "resource": res_name,
                        "version": ver_dir.name,
                        "key": key,
                        "tables": tables,
                        "relations_count": rel_rows,
                        "label_quality": label_quality,
                        "mtime": mtime,
                    }
        self._inventory = InventorySnapshot(
            resources, self._pick_latest_resources(resources), self._scan_inventory_fingerprint()
        )

    def reload_resources(self) -> None:
        """Refresh inventory of parquet resources from disk."""
        with self.query_scope():
            self._refresh_inventory_if_stale(force=True)

    @property
    def resources(self):
        return (self._inventory_scope.get() or self._inventory).resources

    @property
    def _latest_by_resource(self):
        return (self._inventory_scope.get() or self._inventory).latest

    @contextmanager
    def release_scope(self, version: str | None = None, *, snapshot=None):
        """Pin all storage lookups to one atomically published inventory generation."""
        if snapshot is None:
            self._refresh_inventory_if_stale()
            snapshot = self._inventory_scope.get() or self._inventory
        inventory_token = self._inventory_scope.set(snapshot)
        token = taxonomy_token = None
        try:
            version = version or self.releases.default()
            if version == "working":
                version = "latest"
            manifest = None if version == "latest" else self.releases.get(version)
            if manifest is not None:
                for source, resource_version in manifest["resources"].items():
                    key = f"{source}/{resource_version}"
                    if key not in self.resources:
                        raise ValueError(f"OmniPath {version} requires unavailable resource {key}")
            token = self._release_scope.set(manifest)
            taxonomy_token = self._taxonomy_scope.set(self._resolve_taxonomy_names())
            yield version
        finally:
            if taxonomy_token is not None:
                self._taxonomy_scope.reset(taxonomy_token)
            if token is not None:
                self._release_scope.reset(token)
            self._inventory_scope.reset(inventory_token)

    def inventory_snapshot(self):
        self._refresh_inventory_if_stale()
        return self._inventory

    def close(self):
        """Release all connections after the HTTP worker has drained requests."""
        self._query_pool.close()
        for connection in tuple(self._connections):
            self._close_direct_connection(connection)

    def _visible_resources(self) -> dict[str, dict[str, Any]]:
        manifest = self._release_scope.get()
        if manifest is None:
            return self.resources
        return {
            f"{source}/{version}": self.resources[f"{source}/{version}"]
            for source, version in manifest["resources"].items()
        }

    def _visible_latest(self) -> dict[str, dict[str, Any]]:
        if self._release_scope.get() is None:
            return self._latest_by_resource
        return {info["resource"]: info for info in self._visible_resources().values()}

    def _pick_latest_resources(self, resources=None) -> dict[str, dict[str, Any]]:
        """Prefer the highest numeric version; legacy-only resources use file time."""

        def rank(info: dict[str, Any]) -> tuple:
            version = info["version"]
            if re.fullmatch(r"[0-9]+(?:\.[0-9]+)*", version):
                return (1, tuple(map(int, version.split("."))), version)
            return (0, info.get("mtime", 0), version)

        latest: dict[str, dict[str, Any]] = {}
        for info in (resources if resources is not None else self.resources).values():
            current = latest.get(info["resource"])
            if current is None or rank(info) > rank(current):
                latest[info["resource"]] = info
        return latest

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
        self._refresh_inventory_if_stale()
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
        """{resource key (source/version): path} of one table for the selected resources
        that have it."""
        return {
            info["key"]: str(info["tables"][table])
            for info in self._selected_resource_infos(selected_resources)
            if table in info["tables"]
        }

    def _table(self, table: str, selected_resources: list[str] | None = None) -> str:
        """SQL source of one table over the selected resources, with a ``resource`` column."""
        paths = self._table_paths(table, selected_resources)
        if not paths:
            schema = {**PUBLISHED_TABLES, **SERVING_TABLES}[table]
            columns = ", ".join(f"NULL AS {name}" for name in schema.names)
            return f"(SELECT {columns}, NULL::VARCHAR AS resource WHERE FALSE)"
        listed = ", ".join(sql_literal(path) for path in paths.values())
        return (
            f"(SELECT * EXCLUDE (filename), regexp_extract(filename, '([^/]+/[^/]+)/[^/]+$', 1) "
            f"AS resource FROM read_parquet([{listed}], union_by_name=true, filename=true))"
        )

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
        """{resource: values} of (resource, value) pairs, for resources in the inventory."""
        values: dict[str, list] = {}
        for resource, value in pairs:
            if resource in self.resources:
                values.setdefault(resource, []).append(value)
        return dict(sorted(values.items()))

    def _lookup_sql(self, table, column, pairs, where="TRUE", params=()):
        """SQL and parameters of the ``table`` rows (with ``resource``) for (resource,
        value) pairs, reading each resource's file by its sorted ``column``."""
        parts, values = [], []
        for resource, items in self._by_resource(pairs).items():
            condition, items = self._in(column, items)
            path = sql_literal(str(self.resources[resource]["tables"][table]))
            parts.append(
                f"SELECT *, ?::VARCHAR AS resource FROM read_parquet({path}) "
                f"WHERE {condition} AND {where}"
            )
            values += [resource, *items, *params]
        if not parts:
            return f"SELECT * FROM {self._table(table, [])} WHERE FALSE", []
        return " UNION ALL ".join(parts), values

    def _lookup(self, table, column, pairs, fields="*", extra="", params=()):
        """Rows of ``table`` for (resource, value) pairs, one query per resource's file.

        Keys are hashes spread over every file's range, so a multi-key filter across all
        resources reads most files; per resource it reads only the matching row groups.
        """
        rows = []
        for resource, items in self._by_resource(pairs).items():
            condition, items = self._in(column, items)
            path = sql_literal(str(self.resources[resource]["tables"][table]))
            rows += self._fetch_dicts(
                f"SELECT {fields}, ? AS resource FROM read_parquet({path}) WHERE {condition} {extra}",
                [resource, *items, *params],
            )
        return rows

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
        rows = self._lookup("entity", "entity_key", pairs) if pairs else []
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
            return {str(k): ParquetQueries._jsonify(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [ParquetQueries._jsonify(v) for v in value]
        try:
            import numpy as np

            if isinstance(value, np.ndarray):
                return ParquetQueries._jsonify(value.tolist())
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
