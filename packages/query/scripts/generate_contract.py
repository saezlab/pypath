"""Write the query layer's contract and copy the SQL into both clients.

    uv run python packages/query/scripts/generate_contract.py

``packages/query/sql/`` is the one place the query layer is edited: the macros
(``relations.sql``, ``search.sql``, ``entity.sql``) and ``contract.json``, which this
script generates from ``omnipath_core`` (table names, the serving tables' queries,
Biolink names) and the API's export limit. It then copies ``sql/`` into the Python
package (``python/omnipath_client/sql/``) and the R package (``r/inst/sql/``), which
read only these files at run time. Run it after changing the published or serving
tables, Biolink or a macro; ``--check`` fails when a file is out of date (the Python
client's tests run it).
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
import shutil
import sys

from omnipath_core import biolink
from omnipath_core.schema import (
    PUBLISHED_TABLES,
    RELATION_FACET_COLUMNS,
    RELATION_QUALIFIERS,
    ROW_GROUP_SIZE,
    SERVING_TABLES,
)
from omnipath_core.serving_tables import serving_queries
from omnipath_api.queries.services.export import RELATION_EXPORT_LIMIT

QUERY = Path(__file__).resolve().parents[1]
SQL = QUERY / "sql"
COPIES = [QUERY / "python" / "omnipath_client" / "sql", QUERY / "r" / "inst" / "sql"]

# The resources' tables in the build's scratch database carry a global id (resource index
# * 2^32 + the resource's id), so that their rows find their entity or relation copy by
# one column.
GLOBAL_IDS = {
    "relation": ["relation_id", "_rid"],
    "relation_annotation": ["relation_id", "_rid"],
    "relation_evidence": ["relation_id", "_rid"],
    "entity": ["entity_id", "_eid"],
    "entity_identifier": ["entity_id", "_eid"],
    "entity_annotation": ["entity_id", "_eid"],
    "entity_evidence": ["entity_id", "_eid"],
    "entity_term": ["entity_id", "_eid"],
    "entity_group": ["entity_id", "_eid"],
}

def hierarchy_predicates() -> dict[str, bool]:
    """{predicate: whether its edges point parent -> child} of the hierarchy's predicates."""
    names = {name for kind, _, name in biolink_names() if kind == "predicate"} | {"is_a"}
    return {
        name: reverse
        for name in sorted(names)
        if (reverse := biolink.hierarchy_direction(name)) is not None
    }


def derived_tables() -> dict:
    """The database's tables, derived at build from the resources' tables (written to the
    build's scratch database, ``work``), in build order: every entity and relation merged
    over the resources into one row, numbered in key order (``entity_n``, ``relation_n``),
    with their identifiers, annotations and evidence, and the lookup tables sorted for the
    reads they serve (sql/build/<name>.sql). ``scratch`` tables are stages of the build,
    in its scratch database; ``{{partition column}}`` statements run in slices of keys by
    their first character, so that none needs the memory of all of them."""
    qualifiers = ", ".join(RELATION_QUALIFIERS)
    directions = hierarchy_predicates()
    reversed_ = ", ".join(f"'{name}'" for name, reverse in directions.items() if reverse) or "NULL"
    values = {
        # A relation's filter and facet columns.
        "relation_columns": (
            "category, predicate, subject_type, object_type, taxon, sign, is_directed, sources, "
            f"resources, evidence_count, {qualifiers}"
        ),
        "facet_columns": ", ".join(RELATION_FACET_COLUMNS),
        "flip": f"r.predicate IN ({reversed_})",
        "hierarchy_predicates": ", ".join(f"'{name}'" for name in directions) or "NULL",
    }
    tables = {}
    for name, order, scratch in BUILD:
        select = build_sql(name)
        for key, value in values.items():
            select = select.replace("{{" + key + "}}", value)
        tables[name] = {
            "select": select,
            **({"order": order} if order else {}),
            **({"scratch": True} if scratch else {}),
        }
    return tables


# (table, order, scratch) in build order. The evidence tables keep the resources' order
# (their rows are read through the index by number, and sorted by each query).
BUILD = [
    ("endpoint_stubs", "entity_key, _resource", True),
    ("entity_numbers", "entity_key", True),
    ("relation_numbers", "relation_key", True),
    ("entity_copies", "_eid", True),
    ("relation_copies", "_rid", True),
    ("relation_merged", "relation_n", True),
    ("entity_relation_counts", "entity_n", True),
    ("entity_identifier", "entity_n, ordinal", False),
    ("entity_annotation", "entity_n, ordinal", False),
    ("entity_list_counts", "entity_n", True),
    ("entity", "entity_n", False),
    ("entity_evidence", None, False),
    ("relation_rank", "relation_n", True),
    ("relation_annotation", "relation_n, ordinal", False),
    ("relation_evidence", None, False),
    ("relation", "relation_n", False),
    ("relation_endpoint", "entity_n, page_rank, side", False),
    ("relation_refs", "relation_n", True),
    ("relation_reference", "reference_key, page_rank, side", False),
    ("entity_term", "term, entity_n, kind", False),
    ("entity_group", "group_key, kind, entity_n", False),
    ("ontology_edge", "child_key, parent_key", False),
    ("term_entity", "term_n, entity_n, relation_n", False),
    ("relation_facet", "category, predicate", False),
]


def build_sql(name: str) -> str:
    """A build statement with the fragments it includes (the queries' own merge and
    summary SQL), as the builders run it."""
    from omnipath_client import steps

    return steps.expand(
        (SQL / "build" / f"{name}.sql").read_text(),
        lambda fragment: (SQL / "fragments" / f"{fragment}.sql").read_text(),
    )


# The resources' tables the build reads in one statement only: views over their files
# in the scratch database, not copies (written and read again).
WORK_VIEWS = [
    "entity_identifier", "entity_annotation", "entity_evidence", "relation_annotation",
    "relation_evidence", "entity_term", "entity_group",
]


def work_tables() -> list[str]:
    """The resources' tables the build reads (from its scratch database)."""
    scratch = {name for name, _, is_scratch in BUILD if is_scratch}
    read = set()
    for name, _, _ in BUILD:
        read |= set(re.findall(r"\bwork\.(\w+)", (SQL / "build" / f"{name}.sql").read_text()))
    return [t for t in [*PUBLISHED_TABLES, *SERVING_TABLES] if t in read - scratch]


# ART indexes: rows by a list of numbers or keys written into the SQL as a literal list
# (DuckDB uses an index only then): search pages, relation endpoints and groups read
# dozens to hundreds of entities at a time, relation pages their relations and evidence.
INDEXES = [
    ["entity_n", "entity", "entity_n"],
    ["entity_key", "entity", "entity_key"],
    ["entity_evidence_n", "entity_evidence", "entity_n"],
    ["entity_annotation_n", "entity_annotation", "entity_n"],
    ["relation_n", "relation", "relation_n"],
    ["relation_evidence_n", "relation_evidence", "relation_n"],
]
# Changes when the tables, derived tables or indexes a database holds change. 10: op_taxonomy;
# 11: merged entities; 12: merged entities and relations only, numbered.
DATABASE_FORMAT = 12


def biolink_names() -> list[list[str]]:
    """[kind, lowercase alias, storage name] of Biolink classes and predicates."""
    rows = set()
    for alias, name in biolink._index("class").items():
        rows.add(("class", alias, biolink._storage_name(name)))
    for alias, name in biolink._index("slot").items():
        try:
            rows.add(("predicate", alias, biolink._predicate_name(name)))
        except ValueError:
            pass
    return [list(row) for row in sorted(rows)]


def query_tables() -> dict:
    """The rules of display names and entity groups (omnipath_core) as tables."""
    from omnipath_core.display_names import NAME_NAMESPACES, NAME_RULES
    from omnipath_api.queries.services.groups import chemical_types

    return {
        "op_name_rules": {
            "columns": ["pattern VARCHAR", "penalty INTEGER"],
            "rows": [[pattern, penalty] for pattern, penalty in NAME_RULES],
        },
        "op_name_namespaces": {
            "columns": ["namespace VARCHAR"],
            "rows": [[name] for name in NAME_NAMESPACES],
        },
        "op_chemical_types": {
            "columns": ["entity_type VARCHAR"],
            "rows": [[name] for name in chemical_types()],
        },
    }


def contract() -> dict:
    queries = serving_queries("{relation}", "{entity}", "{entity_identifier}")
    return {
        "about": "Generated by packages/query/scripts/generate_contract.py from packages/core; do not edit.",
        "version": 1,
        "published_tables": list(PUBLISHED_TABLES),
        "serving_tables": list(SERVING_TABLES),
        "evidence_table": "evidence_payloads",
        "serving_queries": {
            name: {"select": select, "order": order} for name, (select, order) in queries.items()
        },
        "serving_sources": ["relation", "entity", "entity_identifier"],
        "parquet_row_group_size": ROW_GROUP_SIZE,
        # Databases are written in the storage format of the oldest supported DuckDB, so
        # clients on different DuckDB versions can share a data directory.
        "storage_version": "v1.5.0",
        # Every connection that runs the queries sets these (API and clients alike):
        # literal id lists of up to 100,000 ids read through the indexes, and no
        # materializing of a whole table read twice in one statement (13 s instead of
        # 0.3 s for a gene's entity keys in a full release).
        "settings": {
            "dynamic_or_filter_threshold": 20000,
            "index_scan_max_count": 100000,
            "index_scan_percentage": 0.5,
            "disabled_optimizers": "common_subplan",
        },
        "global_ids": GLOBAL_IDS,
        "work_tables": work_tables(),
        "work_views": WORK_VIEWS,
        "database_format": DATABASE_FORMAT,
        "derived_tables": derived_tables(),
        "indexes": INDEXES,
        "macro_files": ["shape.sql"],
        "relation_export_limit": RELATION_EXPORT_LIMIT,
        "biolink_names": biolink_names(),
        # Tables of rules the queries read, created on each connection from here (so
        # they match the code, whatever the database): {table: {columns, rows}}.
        "query_tables": query_tables(),
    }


def main(check: bool = False) -> int:
    text = json.dumps(contract(), indent=1) + "\n"
    stale = []
    target = SQL / "contract.json"
    if not target.is_file() or target.read_text() != text:
        stale.append(target)
        if not check:
            target.write_text(text)
    sources = sorted(p.relative_to(SQL) for p in SQL.rglob("*") if p.is_file())
    for copies in COPIES:
        for relative in sources:
            source, copy = SQL / relative, copies / relative
            if not copy.is_file() or copy.read_bytes() != source.read_bytes():
                stale.append(copy)
                if not check:
                    copy.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(source, copy)
        if copies.is_dir():
            for extra in sorted(p for p in copies.rglob("*") if p.is_file()):
                if extra.relative_to(copies) not in sources:
                    stale.append(extra)
                    if not check:
                        extra.unlink()
    for path in stale:
        print(("out of date: " if check else "wrote: ") + str(path.relative_to(QUERY.parent)))
    return 1 if check and stale else 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail when a file is out of date")
    sys.exit(main(parser.parse_args().check))
