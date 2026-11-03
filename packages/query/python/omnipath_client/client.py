"""The installed resources as a DuckDB database, with the explorer's queries as SQL."""

from __future__ import annotations

from collections.abc import Sequence
from copy import deepcopy
import json
from pathlib import Path
from typing import Any

import duckdb

from .contracts import CONTRACT, SERVING, ResourceRecord
from .database import DATABASE, in_sync, synchronize
from .errors import ClientError
from .installation import data_directory, read_installation

RELATION_EXPORT_LIMIT: int = CONTRACT["relation_export_limit"]
# Macro parameters that take one value; every other filter takes a list, and one value
# stands for a list of one.
_SCALARS = {"gene_mode", "molecular_endpoint_mode", "is_directed", "sign", "min_evidence_count"}
# Filters of Biolink names, by the kind op_biolink checks them as.
_BIOLINK = {
    "predicates": "predicate",
    "interaction_types": "class",
    "entity_types": "class",
    "partner_types": "class",
}


def _literal(value: Any) -> str:
    """A SQL literal of a filter value: NULL, a boolean, a number, a string or a list.

    The prebuilt queries pass values as literals, not bound parameters: a relation
    from ``sql()`` with parameters is planned without their values, which made a
    query ten times slower.
    """
    if value is None:
        return "NULL"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(int(value))
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_literal(item) for item in value) + "]"
    raise TypeError(f"Unsupported filter value: {value!r}")


def _value(name: str, value: Any) -> Any:
    if name in _SCALARS or value is None:
        return value
    if isinstance(value, (str, int)) and not isinstance(value, bool):
        return [value]
    return list(value)


class Client:
    """The installed resources of one data directory; use as a context manager.

    Everything is read from local files: run ``init()`` first. ``sql`` queries the
    directory's DuckDB database; ``search``, ``relations`` and ``entity`` answer like
    the explorer, through the database's macros ``op_search``, ``op_relations`` and
    ``op_entity``, which any DuckDB client can call too. Returned relations remain
    usable until this client is closed. Instances are intended for one thread. Connect
    again after ``init()`` changed the directory.
    """

    def __init__(
        self,
        data_dir: str | Path | None = None,
        *,
        memory_limit: str = "2GB",
        threads: int | None = None,
    ) -> None:
        self.data_dir = data_directory(data_dir)
        installation = read_installation(self.data_dir)
        if installation is None or not installation["resources"]:
            raise ClientError(
                f"No OmniPath data in {self.data_dir}; run omnipath_client.init() first"
            )
        self._release = installation["release"]
        self._catalog: dict[str, ResourceRecord] = {
            record["resource_id"]: record for record in installation["resources"]
        }
        for record in self._catalog.values():
            for name in (*(f["name"] for f in record["files"]), *(f"{t}.parquet" for t in SERVING)):
                self._check(record, name)
        if not in_sync(self.data_dir, installation):
            # A moved directory, an interrupted init() or another client version: the
            # database (or only its query layer) is written again at once.
            try:
                synchronize(self.data_dir, installation, memory_limit)
            except (OSError, duckdb.Error) as exc:
                raise ClientError(
                    f"The database in {self.data_dir} is out of date; run init() again"
                ) from exc
        self.path = self.data_dir / DATABASE
        self._closed = False
        self._parameters: dict[str, list[str]] = {}
        connection = duckdb.connect(str(self.path), read_only=True)
        try:
            connection.execute("SET memory_limit = ?", [memory_limit])
            if threads is not None:
                connection.execute("SET threads = ?", [threads])
            for name, value in CONTRACT["settings"].items():
                connection.execute(f"SET {name} = {int(value)}")
        except duckdb.Error:
            connection.close()
            raise
        self._connection = connection

    @property
    def release(self) -> str:
        """The release of the last ``init()``."""
        return self._release

    def _ensure_open(self) -> None:
        if self._closed:
            raise ClientError("Client is closed")

    def _check(self, record: dict[str, Any], name: str) -> None:
        """One file of an installed resource; missing or truncated files are errors."""
        path = self.data_dir / "resources" / record["resource_id"] / record["version"] / name
        if not path.resolve().is_relative_to(self.data_dir.resolve()):
            raise ClientError("Installed file escapes the data directory")
        size = next((f["size_bytes"] for f in record["files"] if f["name"] == name), None)
        if not path.is_file() or (size is not None and path.stat().st_size != size):
            raise ClientError(f"Missing or truncated file: {path}; run init() again")

    def resources(self) -> list[dict[str, Any]]:
        """The installed resources with their catalog metadata; an independent copy."""
        self._ensure_open()
        return deepcopy(list(self._catalog.values()))

    @property
    def connection(self) -> duckdb.DuckDBPyConnection:
        """The read-only DuckDB connection to the database, for DuckDB's own interface."""
        self._ensure_open()
        return self._connection

    def sql(self, query: str, parameters: Sequence[Any] | None = None) -> duckdb.DuckDBPyRelation:
        """Run SQL on the database: ``entity``, ``relation`` and the other tables over
        every installed resource (with ``_resource`` and ``_resource_version``), or one
        resource's tables as ``signor.relation``; the macros ``op_search``,
        ``op_relations`` and ``op_entity`` too. Caller-supplied local SQL; not a
        sandbox."""
        return self.connection.sql(
            query, params=list(parameters) if parameters is not None else None
        )

    def _filters(self, macro: str, internal: set[str], filters: dict[str, Any]) -> dict[str, Any]:
        """The filters by name, checked against the macro's parameters."""
        if macro not in self._parameters:
            row = self.connection.execute(
                "SELECT parameters FROM duckdb_functions() WHERE function_name = ? LIMIT 1",
                [macro],
            ).fetchone()
            self._parameters[macro] = sorted(set(row[0]) - internal) if row else []
        available = self._parameters[macro]
        unknown = set(filters) - set(available)
        if unknown:
            raise TypeError(f"Unknown filters: {sorted(unknown)}; available: {available}")
        sources = _value("sources", filters.get("sources")) or []
        missing = sorted(set(sources) - set(self._catalog))
        if missing:
            raise ValueError(
                f"Resources unavailable here: {missing}; installed: {sorted(self._catalog)}"
            )
        values = {name: _value(name, value) for name, value in filters.items()}
        # Biolink names now, not when the rows are read: the macros check them with
        # the same op_biolink, which raises an error naming an unknown one.
        for name, kind in _BIOLINK.items():
            if values.get(name):
                try:
                    self.connection.execute(
                        f"SELECT list(name) FROM op_biolink({_literal(kind)}, {_literal(values[name])})"
                    ).fetchall()
                except duckdb.InvalidInputException as exc:
                    raise ValueError(str(exc).removeprefix("Invalid Input Error: ")) from None
        return values

    def _call(self, macro: str, positional: list[Any], named: dict[str, Any]) -> duckdb.DuckDBPyRelation:
        arguments = [_literal(value) for value in positional] + [
            f"{name} := {_literal(value)}" for name, value in named.items()
        ]
        return self.connection.sql(f"SELECT * FROM {macro}({', '.join(arguments)})")

    def search(self, query: str, *, limit: int = 20, **filters: Any) -> duckdb.DuckDBPyRelation:
        """Entities matching a name or identifier, merged across resources and ranked
        as in the explorer's search. ``filters`` are the explorer's: ``sources``,
        ``entity_types``, ``taxonomy_ids``, ``reference_entity_keys`` and
        ``entity_ids``. The macro ``op_search``."""
        named = self._filters("op_search", {"query", "max_results"}, filters)
        return self._call("op_search", [query], {**named, "max_results": int(limit)})

    def relations(
        self,
        entities: str | Sequence[str] | None = None,
        *,
        limit: int = RELATION_EXPORT_LIMIT,
        **filters: Any,
    ) -> duckdb.DuckDBPyRelation:
        """Relations as the explorer exports them: one row per relation with its
        annotations, evidence and the product records its evidence names.

        ``entities`` are exact names, identifiers, entity keys or group keys such as
        ``gene:entrez:7157``, at either endpoint. ``filters`` are the explorer's, such
        as ``sources``, ``predicates``, ``taxonomy_ids``, ``signs``,
        ``reference_entity_keys``, ``protein_entity_keys`` and ``isoform_identifiers``.
        At most 100,000 relations; larger selections are SQL on ``relation``. The
        macro ``op_relations``."""
        if not 1 <= limit <= RELATION_EXPORT_LIMIT:
            raise ValueError(f"limit must be between 1 and {RELATION_EXPORT_LIMIT}")
        if "entity_ids" in filters:
            if entities is not None:
                raise TypeError("Give entities once: positionally or as entity_ids")
            entities = filters.pop("entity_ids")
        named = self._filters("op_relations", {"entities", "max_relations"}, filters)
        if entities is not None:
            named["entities"] = _value("entities", entities)
        return self._call("op_relations", [], {**named, "max_relations": int(limit)})

    def entity(self, entity_key: str) -> dict[str, Any] | None:
        """One entity by its exact key (as ``search`` returns it), with its identifiers,
        annotations and relationship summary as the explorer's detail page shows them;
        None when no installed resource has it. The macro ``op_entity``."""
        row = self.connection.execute(
            "SELECT details FROM op_entity(?)", [str(entity_key)]
        ).fetchone()
        return json.loads(row[0]) if row else None

    def close(self) -> None:
        if not self._closed:
            self._connection.close()
            self._closed = True

    def __enter__(self) -> Client:
        self._ensure_open()
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()


def connect(
    data_dir: str | Path | None = None,
    *,
    memory_limit: str = "2GB",
    threads: int | None = None,
) -> Client:
    """Open the data directory's database: the directory given, else OMNIPATH_DATA_DIR,
    else the per-user default."""
    return Client(data_dir, memory_limit=memory_limit, threads=threads)
