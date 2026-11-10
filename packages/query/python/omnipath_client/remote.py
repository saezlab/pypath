"""The hosted API as a client: ``search``, ``relations`` and ``entity`` without a download.

``connect()`` without a data directory returns an ``ApiClient``. It answers like the
local ``Client`` (the same columns, rows and JSON, from the API serving the same
release database) but runs no SQL of its own: ``sql()`` needs a local installation
(``init()``, then ``connect(data_dir)``).
"""

from __future__ import annotations

from collections.abc import Sequence
from copy import deepcopy
from importlib.resources import files
import os
import re
from typing import Any
from urllib.parse import quote

import itertools
from pathlib import Path
import tempfile

import duckdb
import httpx

from . import steps
from .contracts import CONTRACT, validate_catalog, validate_url
from .errors import ClientError

DEFAULT_API_URL = "https://omnipath-metabo-dev.schaul.click/api"
RELATION_EXPORT_LIMIT: int = CONTRACT["relation_export_limit"]
# The local search's columns and their types.
SEARCH_COLUMNS = (
    ("entity_key", "VARCHAR"),
    ("label", "VARCHAR"),
    ("entity_type", "VARCHAR"),
    ("namespace", "VARCHAR"),
    ("identifier", "VARCHAR"),
    ("taxon", "VARCHAR"),
    ("relation_count", "BIGINT"),
    ("reference_entity_key", "VARCHAR"),
    ("gene_reference_keys", "VARCHAR[]"),
)
_answers = itertools.count()
_SCALARS = {"gene_mode", "molecular_endpoint_mode", "is_directed", "sign", "min_evidence_count"}


def _parameters(query: str) -> list[str]:
    """The parameters of one of the shared queries (sql/queries/<name>.sql): the filters
    it takes, in both modes."""
    try:
        return sorted(set(steps.load(query).params) - {"resources"})
    except KeyError:
        raise ClientError(f"No shared query {query}") from None


def _value(name: str, value: Any) -> Any:
    if name in _SCALARS or value is None:
        return value
    if isinstance(value, (str, int)) and not isinstance(value, bool):
        return [value]
    return list(value)


class ApiClient:
    """The hosted API's answers for one release; use as a context manager.

    ``search``, ``relations`` and ``entity`` return what the local client returns, as
    DuckDB relations (search, relations) and a dict (entity). Nothing is downloaded but
    the answers; for SQL, install the resources (``init()``) and ``connect(data_dir)``.
    """

    def __init__(
        self,
        api_url: str | None = None,
        *,
        release: str = "latest",
        timeout: float = 120,
        memory_limit: str = "2GB",
        threads: int | None = None,
        http_client: httpx.Client | None = None,
    ) -> None:
        if not isinstance(release, str) or not release.strip():
            raise ValueError("release must not be empty")
        self.api_url = validate_url(
            api_url or os.environ.get("OMNIPATH_API_URL") or DEFAULT_API_URL
        ).rstrip("/")
        self._release = release
        self._http = http_client or httpx.Client(timeout=timeout, follow_redirects=True)
        self._own_http = http_client is None
        # Answers become relations of an in-memory database, as local ones are.
        self._db = duckdb.connect()
        self._db.execute("SET memory_limit = ?", [memory_limit])
        if threads is not None:
            self._db.execute("SET threads = ?", [threads])
        self._closed = False
        self._filters_of: dict[str, list[str]] = {}

    @property
    def release(self) -> str:
        """The release asked for (``latest``: the newest release the API serves)."""
        return self._release

    @property
    def mode(self) -> str:
        return "api"

    def _ensure_open(self) -> None:
        if self._closed:
            raise ClientError("Client is closed")

    def _request(self, method: str, path: str, **kwargs: Any) -> httpx.Response:
        self._ensure_open()
        params = {**kwargs.pop("params", {}), "release": self._release}
        try:
            response = self._http.request(method, f"{self.api_url}{path}", params=params, **kwargs)
        except httpx.HTTPError as exc:
            raise ClientError(f"The API at {self.api_url} is unreachable: {exc}") from exc
        if response.status_code in (400, 422):
            raise ValueError(_detail(response))
        if response.status_code == 404 and path.startswith("/entities/"):
            return response
        if response.status_code >= 400:
            raise ClientError(f"{method} {path}: HTTP {response.status_code}: {_detail(response)}")
        return response

    def _filters(self, macro: str, internal: set[str], filters: dict[str, Any]) -> dict[str, Any]:
        if macro not in self._filters_of:
            self._filters_of[macro] = sorted(set(_parameters(macro)) - internal)
        available = self._filters_of[macro]
        unknown = set(filters) - set(available)
        if unknown:
            raise TypeError(f"Unknown filters: {sorted(unknown)}; available: {available}")
        return {
            name: _value(name, value) for name, value in filters.items() if value is not None
        }

    def resources(self) -> list[dict[str, Any]]:
        """The release's resources with their catalog metadata."""
        response = self._request("GET", "/resources", params={"shape": "svelte"})
        return deepcopy(list(validate_catalog(response.json().get("resources")).values()))

    @property
    def connection(self) -> duckdb.DuckDBPyConnection:
        raise ClientError(
            "API mode has no database: install resources with init() and "
            "connect(data_dir) for SQL"
        )

    def sql(self, query: str, parameters: Sequence[Any] | None = None) -> duckdb.DuckDBPyRelation:
        """Not available in API mode: SQL runs on a local installation."""
        raise ClientError(
            "sql() needs a local installation: run omnipath_client.init(), then "
            "connect(data_dir) (or connect(mode='local'))"
        )

    def search(self, query: str, *, limit: int = 20, **filters: Any) -> duckdb.DuckDBPyRelation:
        """Entities matching a name or identifier, merged across resources and ranked
        as in the explorer's search; the filters of the local ``search``."""
        named = self._filters("entity_search", {"query", "limit", "cursor_phase", "cursor_offset", "cursor_key"}, filters)
        response = self._request(
            "POST", "/entities/search", json={"query": query, "limit": int(limit), "filters": named}
        )
        rows = [
            (
                e["entityPk"],
                e.get("label"),
                e.get("entityType"),
                e.get("canonicalIdentifierType"),
                e.get("canonicalIdentifier"),
                e.get("taxonomyId"),
                e.get("relationCount"),
                e.get("referenceEntityKey"),
                e.get("geneReferenceKeys") or [],
            )
            for e in response.json().get("entities") or []
        ]
        table = f"answer_{next(_answers)}"
        columns = ", ".join(f"{name} {kind}" for name, kind in SEARCH_COLUMNS)
        self._db.execute(f"CREATE TABLE {table} ({columns})")
        if rows:
            self._db.executemany(
                f"INSERT INTO {table} VALUES ({', '.join('?' for _ in SEARCH_COLUMNS)})", rows
            )
        return self._db.table(table)

    def relations(
        self,
        entities: str | Sequence[str] | None = None,
        *,
        limit: int = RELATION_EXPORT_LIMIT,
        **filters: Any,
    ) -> duckdb.DuckDBPyRelation:
        """Relations as the explorer exports them, with their annotations, evidence and
        the product records their evidence names; the filters of the local
        ``relations``. At most 100,000 relations."""
        if not 1 <= limit <= RELATION_EXPORT_LIMIT:
            raise ValueError(f"limit must be between 1 and {RELATION_EXPORT_LIMIT}")
        if "entity_ids" in filters:
            if entities is not None:
                raise TypeError("Give entities once: positionally or as entity_ids")
            entities = filters.pop("entity_ids")
        named = self._filters("relations", {"entities", "max_relations"}, filters)
        if entities is not None:
            named["entity_ids"] = _value("entities", entities)
        response = self._request("POST", "/export", json={"filters": named, "limit": int(limit)})
        table = f"answer_{next(_answers)}"
        with tempfile.TemporaryDirectory(prefix="omnipath-") as directory:
            path = Path(directory) / "relations.parquet"
            path.write_bytes(response.content)
            self._db.execute(
                f"CREATE TABLE {table} AS SELECT * FROM read_parquet('{str(path).replace(chr(39), chr(39) * 2)}')"
            )
        return self._db.table(table)

    def entity(self, entity_key: str) -> dict[str, Any] | None:
        """One entity by its exact key, with its identifiers, annotations and
        relationship summary as the explorer's detail page shows them; None when the
        release has no such entity."""
        response = self._request(
            "GET", f"/entities/{quote(str(entity_key), safe='')}", params={"paged": "false"}
        )
        if response.status_code == 404:
            return None
        return response.json()

    def close(self) -> None:
        if not self._closed:
            self._db.close()
            if self._own_http:
                self._http.close()
            self._closed = True

    def __enter__(self) -> ApiClient:
        self._ensure_open()
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()


def _detail(response: httpx.Response) -> str:
    try:
        detail = response.json().get("detail")
    except ValueError:
        detail = None
    return str(detail or response.text[:300])
