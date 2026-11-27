"""The Parquet tables of a published release, as views in one DuckDB connection.

A release pins one version of each resource. The subset products read the
pinned resources' tables as they are published; nothing is loaded anywhere
first. Entities are identified across resources by ``entity_key``, so the views
merge a key's rows from every resource that holds it.

Views (all over the pinned resources):

``release_entity``
    One row per ``entity_key``: type, canonical namespace and identifier, label,
    reference key and structure group. When resources differ on a value, the
    resource that sorts first and has one wins, so a rebuild reads the same.
``release_identifier``
    Every (entity_key, ns, id) the release holds, with the canonical identity
    included.
``release_relation``
    Every relation with the resource it came from.
``release_relation_evidence``
    Every relation with one row per evidence it carries (dataset and row id).
``release_entity_annotation``
    The annotations of each entity, with the resource that stated them.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import duckdb

from omnipath_core.versioning import input_module_revision

# A full release needs more than the default; the environment sets it for a build.
MEMORY_ENVIRONMENT = "OMNIPATH_SUBSETS_MEMORY"
DEFAULT_MEMORY_LIMIT = "4GB"

# The entity types a chemical subset reads.
CHEMICAL_ENTITY_TYPES = ("chemical_entity", "small_molecule")

# A relation of these kinds is an edge of an ontology.
ONTOLOGY_STATEMENT_KINDS = ("ontology", "ontology_axiom")


def resource_code_revision(data_root: str | Path, source: str, version: str) -> str | None:
    """The code revision a resource was built with, from its build manifest."""
    path = Path(data_root) / "resources" / source / version / "build_manifest.json"
    if not path.is_file():
        return None
    return input_module_revision(json.loads(path.read_text()))[0]


def pinned_resources(data_root: str | Path, release: str) -> dict[str, str]:
    path = Path(data_root) / "releases" / f"{release}.json"
    return dict(json.loads(path.read_text())["resources"])


def _files(data_root: Path, resources: dict[str, str], table: str) -> list[str]:
    return [
        str(data_root / "resources" / name / version / f"{table}.parquet")
        for name, version in sorted(resources.items())
    ]


def _scan(files: list[str], columns: str) -> str:
    quoted = ", ".join("'" + path.replace("'", "''") + "'" for path in files)
    return (
        f"SELECT regexp_extract(filename, '/resources/([^/]+)/[^/]+/[^/]+$', 1) AS resource, "
        f"{columns} FROM read_parquet([{quoted}], filename = true, union_by_name = true)"
    )


def open_release(
    data_root: str | Path,
    release: str,
    *,
    memory_limit: str | None = None,
    threads: int | None = None,
) -> duckdb.DuckDBPyConnection:
    """A DuckDB connection with the release's views; close it when done.

    The memory limit is ``memory_limit``, else ``OMNIPATH_SUBSETS_MEMORY``, else 4GB.
    """
    root = Path(data_root).resolve()
    resources = pinned_resources(root, release)
    con = duckdb.connect()
    con.execute(
        "SET memory_limit = ?",
        [memory_limit or os.environ.get(MEMORY_ENVIRONMENT) or DEFAULT_MEMORY_LIMIT],
    )
    if threads is not None:
        con.execute("SET threads = ?", [int(threads)])
    con.execute(
        "CREATE VIEW resource_entity AS "
        + _scan(
            _files(root, resources, "entity"),
            "entity_id, entity_key, entity_type, namespace, identifier, taxon, label, "
            "reference_entity_key, group_connectivity",
        )
    )
    con.execute(
        "CREATE VIEW resource_identifier AS "
        + _scan(_files(root, resources, "entity_identifier"), "entity_id, ns, id, source")
    )
    con.execute(
        "CREATE VIEW release_relation AS "
        + _scan(
            _files(root, resources, "relation"),
            "relation_id, relation_key, statement_kind, subject_entity_key, predicate, "
            "object_entity_key",
        )
    )
    con.execute(
        "CREATE VIEW resource_relation_evidence AS "
        + _scan(
            _files(root, resources, "relation_evidence"), "relation_id, ordinal, dataset, row_id"
        )
    )
    con.execute(
        "CREATE VIEW resource_entity_annotation AS "
        + _scan(_files(root, resources, "entity_annotation"), "entity_id, term, value")
    )
    # A value is taken from the first resource that has it, in name order.
    con.execute(
        """
        CREATE VIEW release_entity AS
        SELECT entity_key,
               arg_min(entity_type, resource) FILTER (WHERE entity_type IS NOT NULL) AS entity_type,
               arg_min(namespace, resource) FILTER (WHERE namespace IS NOT NULL) AS namespace,
               arg_min(identifier, resource) FILTER (WHERE identifier IS NOT NULL) AS identifier,
               arg_min(label, resource) FILTER (WHERE label IS NOT NULL) AS label,
               arg_min(reference_entity_key, resource)
                 FILTER (WHERE reference_entity_key IS NOT NULL) AS reference_entity_key,
               arg_min(group_connectivity, resource)
                 FILTER (WHERE group_connectivity IS NOT NULL) AS group_connectivity
        FROM resource_entity
        GROUP BY entity_key
        """
    )
    con.execute(
        """
        CREATE VIEW release_relation_evidence AS
        SELECT r.resource, r.relation_key, r.statement_kind, r.subject_entity_key, r.predicate,
               r.object_entity_key, e.ordinal, e.dataset, e.row_id
        FROM release_relation r
        LEFT JOIN resource_relation_evidence e
          ON e.resource = r.resource AND e.relation_id = r.relation_id
        """
    )
    con.execute(
        """
        CREATE VIEW release_entity_annotation AS
        SELECT a.resource, e.entity_key, a.term, a.value
        FROM resource_entity_annotation a
        JOIN resource_entity e ON e.resource = a.resource AND e.entity_id = a.entity_id
        """
    )
    con.execute(
        """
        CREATE VIEW release_identifier AS
        SELECT DISTINCT e.entity_key, i.ns, i."id" AS id
        FROM resource_identifier i
        JOIN resource_entity e ON e.resource = i.resource AND e.entity_id = i.entity_id
        UNION
        SELECT DISTINCT entity_key, namespace, identifier
        FROM resource_entity
        WHERE namespace IS NOT NULL AND identifier IS NOT NULL
        """
    )
    return con
