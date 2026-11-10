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
from pathlib import Path
import shutil
import sys

from omnipath_core import biolink
from omnipath_core.schema import PUBLISHED_TABLES, RELATION_QUALIFIERS, ROW_GROUP_SIZE, SERVING_TABLES
from omnipath_core.serving_tables import serving_queries
from omnipath_api.queries.services.export import RELATION_EXPORT_LIMIT

QUERY = Path(__file__).resolve().parents[1]
SQL = QUERY / "sql"
COPIES = [QUERY / "python" / "omnipath_client" / "sql", QUERY / "r" / "inst" / "sql"]

# Rows of one table read by id or key (lookups): (table, column, kind). "at" tables are
# read by a resource's own id, "with" tables by a key or term meaning the same in every
# resource.
LOOKUPS = [
    ["relation", "relation_id", "at"],
    ["relation_annotation", "relation_id", "at"],
    ["relation_evidence", "relation_id", "at"],
    ["entity", "entity_id", "at"],
    ["entity_identifier", "entity_id", "at"],
    ["entity_annotation", "entity_id", "at"],
    ["entity_evidence", "entity_id", "at"],
    ["entity", "entity_key", "with"],
    ["relation_endpoint", "key", "with"],
    ["entity_term", "term", "with"],
    ["entity_group", "group_key", "with"],
]
# Native tables carry a global id (resource index * 2^32 + the resource's id), so a
# lookup over every resource is a filter on one column.
GLOBAL_IDS = {
    "relation": ["relation_id", "_rid"],
    "relation_annotation": ["relation_id", "_rid"],
    "relation_evidence": ["relation_id", "_rid"],
    "relation_endpoint": ["relation_id", "_rid"],
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
    """The database's own tables, derived from the tables above once they are written:
    the same statements merged or regrouped across every resource the database holds,
    each sorted for the reads it serves. They keep each row's resource, so a query over
    some resources only (or without the ones a query policy disables) filters them."""
    resource = "_resource || '/' || _resource_version AS resource"
    qualifiers = ", ".join(RELATION_QUALIFIERS)
    # A relation's filter and facet columns.
    columns = (
        "category, predicate, subject_type, object_type, taxon, "
        f"sign, is_directed, sources, evidence_count, {qualifiers}"
    )
    directions = hierarchy_predicates()
    listed = ", ".join(f"'{name}'" for name in directions)
    reversed_ = ", ".join(f"'{name}'" for name, reverse in directions.items() if reverse) or "NULL"
    flip = f"predicate IN ({reversed_})"
    return {
        # Relation keys are sorted only within each resource: one table over all of
        # them, by the key's first 60 bits (the keys are sha256 hex: 8 bytes instead of
        # 64 characters). A reader confirms the full key on the relation rows it finds.
        "relation_key_lookup": {
            "select": "SELECT try_cast('0x' || substr(relation_key, 1, 15) AS BIGINT) AS key_prefix, _rid FROM relation",
            "order": "key_prefix",
        },
        # Every resource's entity rows sorted by key: an entity's copies are adjacent, so
        # reading them by key reads a few row groups. key_rank numbers the keys in that
        # order, so that other tables can name an entity by 4 bytes instead of its
        # 64-character key. (The entity table itself is in global id order, for reads by
        # id.)
        "entity_by_key": {
            "select": "SELECT *, dense_rank() OVER (ORDER BY entity_key)::INTEGER AS key_rank FROM entity",
            "order": "entity_key, _resource, entity_id",
        },
        # The lowercase labels and identifiers, rewritten sorted by term over every
        # resource and with the rank of the entity's key ("replaces": the table itself,
        # not a copy): a prefix is one range, and its entities are read by key rank.
        "entity_term": {
            "select": "SELECT t.*, e.key_rank FROM entity_term t JOIN entity_by_key e USING (_eid)",
            "order": "term, _resource, entity_id",
            "replaces": True,
        },
        # The endpoint table rewritten over every resource, sorted by key and then page
        # order, each row with its relation's filter and facet columns and its place in
        # the page order of every relation (page_rank: category, predicate, labels,
        # resource, id), instead of the labels themselves: the relations of an entity (or
        # of a reference key) are one range, in page order. Its columns are the written
        # table's and more, so readers of it are unchanged. One scan of relation gives
        # each relation's four keys (two ends, entity and reference).
        "relation_endpoint": {
            "select": f"""SELECT * FROM (
                SELECT unnest([subject_entity_key, subject_reference_entity_key,
                        object_entity_key, object_reference_entity_key]) AS key,
                    unnest(['entity', 'reference', 'entity', 'reference']) AS key_kind,
                    unnest(['subject', 'subject', 'object', 'object']) AS side,
                    relation_id, _resource, _resource_version, _rid, {resource},
                    row_number() OVER (ORDER BY category, predicate, subject_label, object_label,
                        _resource || '/' || _resource_version, relation_id)::INTEGER AS page_rank,
                    {columns}
                FROM relation) WHERE key IS NOT NULL""",
            "order": "key, key_kind, page_rank",
            "replaces": True,
        },
        # The ontologies' hierarchies as child -> parent edges.
        "ontology_edge": {
            "select": f"""SELECT
                CASE WHEN {flip} THEN object_entity_key ELSE subject_entity_key END AS child_key,
                CASE WHEN {flip} THEN object_label ELSE subject_label END AS child_label,
                predicate,
                CASE WHEN {flip} THEN subject_entity_key ELSE object_entity_key END AS parent_key,
                CASE WHEN {flip} THEN subject_label ELSE object_label END AS parent_label,
                {resource}
                FROM relation WHERE statement_kind = 'ontology' AND predicate IN ({listed or 'NULL'})""",
            "order": "child_key, parent_key",
        },
        # Each relation touching an ontology term, from the term's side: what the term
        # annotates. The annotated entity and the relation are also numbered (dense
        # ranks of their keys), so that counting distinct ones counts integers.
        "term_entity": {
            "select": f"""WITH terms AS (SELECT DISTINCT entity_key FROM entity WHERE entity_type = 'ontology_class'),
                edges AS (
                    SELECT object_entity_key AS term_key, subject_entity_key AS entity_key, relation_key, {resource}
                    FROM relation WHERE object_entity_key IN (SELECT entity_key FROM terms)
                    UNION ALL SELECT subject_entity_key, object_entity_key, relation_key, {resource}
                    FROM relation WHERE subject_entity_key IN (SELECT entity_key FROM terms))
                SELECT term_key, entity_key, resource,
                    dense_rank() OVER (ORDER BY entity_key)::INTEGER AS entity_n,
                    dense_rank() OVER (ORDER BY relation_key)::INTEGER AS relation_n
                FROM edges""",
            "order": "term_key",
        },
        # Entities merged once, at build, instead of on every request: one row per entity
        # with its summary precomputed and its merged identifier and annotation lists
        # (sql/build/).
        # A request read an entity's copies by long key lists, which scan the tables, and
        # merged them (TP53's 64 copies: about 600 ms of CPU).
        **{name: {"select": build_sql(name), "order": order} for name, order in MERGED_TABLES.items()},
    }


# The merged entity tables (sql/build/<name>.sql) and their order.
MERGED_TABLES = {"entity_merged": "entity_key"}


def build_sql(name: str) -> str:
    """A build statement with the fragments it includes (the queries' own merge and
    summary SQL), as the builders run it."""
    from omnipath_client import steps

    return steps.expand(
        (SQL / "build" / f"{name}.sql").read_text(),
        lambda fragment: (SQL / "fragments" / f"{fragment}.sql").read_text(),
    )


# ART indexes: rows by a list of global ids (DuckDB uses an index only for ids written into
# the SQL as a literal list). Measured on a full release: 7,582 relations in 0.3 s instead
# of 5.6 s; 50 entities from across the release in 0.06-0.09 s instead of 0.4 s (a list
# filter reads every column of the row groups it cannot skip), their evidence and
# annotations 2-4 times faster. Together the entity indexes take 0.66 GB. Not on
# entity_identifier (about 2 GB for 106 million rows) or entity_key.
INDEXES = [
    # The merged entities by key: search pages, relation endpoints and groups read dozens
    # to hundreds of keys at a time.
    ["entity_merged_key", "entity_merged", "entity_key"],
    ["relation_rid", "relation", "_rid"],
    ["relation_evidence_rid", "relation_evidence", "_rid"],
    ["entity_eid", "entity", "_eid"],
    ["entity_evidence_eid", "entity_evidence", "_eid"],
    ["entity_annotation_eid", "entity_annotation", "_eid"],
]
# Changes when the tables, derived tables or indexes a database holds change. 10: op_taxonomy;
# 11: merged entities.
DATABASE_FORMAT = 11


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
        "lookups": LOOKUPS,
        "global_ids": GLOBAL_IDS,
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
