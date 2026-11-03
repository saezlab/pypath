"""Read a resource version's normalized tables for assertions."""

import json
from pathlib import Path

import pyarrow.parquet as pq
from omnipath_core.schema import PUBLISHED_TABLES, SERVING_TABLES

TABLES = (*PUBLISHED_TABLES, *SERVING_TABLES)


class Tables(dict):
    """{table name: rows in file order}."""

    def children(self, table, parent):
        """Rows of an entity_* or relation_* child table for one parent row, by ordinal."""
        key = "relation_id" if table.startswith("relation") else "entity_id"
        rows = [row for row in self[table] if row[key] == parent[key]]
        return sorted(rows, key=lambda row: row.get("ordinal", 0))

    def comparable(self):
        """Every table as a sorted list of JSON rows, for invariance checks."""
        return {
            name: sorted(json.dumps(row, sort_keys=True, default=str) for row in rows)
            for name, rows in self.items()
        }


def read_tables(files) -> Tables:
    """``files`` is the writer's {table: path} result or a resource version directory."""
    if not isinstance(files, dict):
        files = {name: Path(files) / f"{name}.parquet" for name in TABLES}
    return Tables({name: pq.read_table(files[name]).to_pylist() for name in TABLES})
