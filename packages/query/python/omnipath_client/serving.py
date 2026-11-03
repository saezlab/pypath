"""A resource version's serving tables, derived from its published tables.

The queries are the build's (``omnipath_core.serving_tables``), shipped in
``sql/contract.json`` so that every client derives the same tables.
"""

from __future__ import annotations

import os
from pathlib import Path

from .contracts import CONTRACT, SERVING


def _literal(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def write_serving_tables(db, directory, work, tables=SERVING) -> None:
    """Write the serving ``tables`` into ``directory`` from its ``relation``, ``entity``
    and ``entity_identifier`` tables, each sorted as the build sorts it.

    ``db`` is a DuckDB connection (its temporary directory takes what a sort does not
    fit in memory); ``work`` an existing directory for the files being written.
    """
    directory, work = Path(directory), Path(work)
    sources = {
        name: f"read_parquet({_literal(str(directory / f'{name}.parquet'))})"
        for name in CONTRACT["serving_sources"]
    }
    options = (
        f"(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE {CONTRACT['parquet_row_group_size']})"
    )
    for name in tables:
        query = CONTRACT["serving_queries"][name]
        select = query["select"]
        for source, sql in sources.items():
            select = select.replace("{" + source + "}", sql)
        staging = work / f"{name}.parquet"
        db.execute(
            f"COPY (SELECT * FROM ({select}) ORDER BY {query['order']}) "
            f"TO {_literal(str(staging))} {options}"
        )
        os.replace(staging, directory / f"{name}.parquet")
