"""Write a resource version's serving tables from its published tables.

The build writes them after the published tables; a client derives them the same way
after downloading the published tables, which is why they are not downloads.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq

from omnipath_core.schema import RELATION_FACET_COLUMNS, ROW_GROUP_SIZE, SERVING_TABLES

_COPY_OPTIONS = f"(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE {ROW_GROUP_SIZE})"


def _sql_path(path) -> str:
    return str(path).replace("'", "''")


def serving_queries(relations: str, entities: str, identifiers: str) -> dict[str, tuple[str, str]]:
    """{serving table: (SELECT, sort order)} over the given SQL sources of a resource
    version's ``relation``, ``entity`` and ``entity_identifier`` tables.

    The clients ship these queries (``omnipath_client/sql/contract.json``) and run them
    themselves; ``write_serving_tables`` runs them for the build.
    """
    endpoints = " UNION ALL ".join(
        f"SELECT {key} AS key, '{kind}' AS key_kind, '{side}' AS side, relation_id, predicate, category FROM {relations}"
        for side in ("subject", "object")
        for kind, key in (
            ("entity", f"{side}_entity_key"),
            ("reference", f"{side}_reference_entity_key"),
        )
    )
    facet_columns = ", ".join(RELATION_FACET_COLUMNS)
    return {
        "relation_endpoint": (
            f"SELECT * FROM ({endpoints}) WHERE key IS NOT NULL",
            "key, relation_id",
        ),
        "entity_group": (
            f"""SELECT * FROM (
            SELECT reference_entity_key AS group_key, 'reference' AS kind, entity_id FROM {entities}
            UNION ALL SELECT unnest(gene_reference_keys), 'gene', entity_id FROM {entities}
            UNION ALL SELECT group_connectivity, 'connectivity', entity_id FROM {entities}
        ) WHERE group_key IS NOT NULL""",
            "group_key, entity_id",
        ),
        "entity_term": (
            f"""SELECT DISTINCT term, kind, entity_id FROM (
            SELECT lower(label) AS term, 'label' AS kind, entity_id FROM {entities}
            UNION ALL SELECT lower(identifier), 'identifier', entity_id FROM {entities}
            UNION ALL SELECT lower(id), 'identifier', entity_id FROM {identifiers}
        ) WHERE term IS NOT NULL AND term <> ''""",
            "term, entity_id",
        ),
        "relation_facet": (
            f"""SELECT {facet_columns}, max(evidence_count) AS evidence_count,
            count(*) AS relation_count FROM {relations} GROUP BY ALL""",
            "category, predicate",
        ),
    }


def write_serving_tables(db, directory, work, tables=tuple(SERVING_TABLES)) -> None:
    """Write the serving ``tables`` (by default all of ``SERVING_TABLES``) into
    ``directory`` from its ``relation``, ``entity`` and ``entity_identifier`` tables.

    ``db`` is a DuckDB connection; ``work`` an existing directory for temporary files.
    """
    directory, work = Path(directory), Path(work)
    queries = serving_queries(
        *(
            f"read_parquet('{_sql_path(directory / name)}.parquet')"
            for name in ("relation", "entity", "entity_identifier")
        )
    )
    for name in tables:
        query, order = queries[name]
        _sorted_copy(db, query, order, directory / f"{name}.parquet", work)


def _sorted_copy(db, query, order, target, work, max_rows=2_000_000):
    """Write ``query`` sorted by ``order``, in key ranges of at most ``max_rows`` rows
    (boundaries from a sample) appended in order: one global sort of tens of millions
    of rows with long string keys does not fit a build's memory limit."""
    unsorted = work / f"{Path(target).stem}.unsorted.parquet"
    db.execute(f"COPY ({query}) TO '{_sql_path(unsorted)}' (FORMAT PARQUET)")
    key = order.split(",")[0].strip()
    source = f"read_parquet('{_sql_path(unsorted)}')"
    total = pq.read_metadata(unsorted).num_rows
    parts = -(-total // max_rows)
    if parts <= 1:
        db.execute(
            f"COPY (SELECT * FROM {source} ORDER BY {order}) TO '{_sql_path(target)}' {_COPY_OPTIONS}"
        )
    else:
        sample = sorted(
            row[0]
            for row in db.execute(
                f"SELECT {key} FROM (SELECT {key} FROM {source}) USING SAMPLE 200000 ROWS"
            ).fetchall()
        )
        bounds = sorted({sample[len(sample) * i // parts] for i in range(1, parts)})
        writer = None
        try:
            for low, high in zip([None, *bounds], [*bounds, None]):
                where = " AND ".join(
                    [f"{key} >= ?"] * (low is not None) + [f"{key} < ?"] * (high is not None)
                )
                batches = db.execute(
                    f"SELECT * FROM {source} WHERE {where} ORDER BY {order}",
                    [v for v in (low, high) if v is not None],
                ).fetch_record_batch(ROW_GROUP_SIZE)
                for batch in batches:
                    if writer is None:
                        writer = pq.ParquetWriter(target, batch.schema, compression="zstd")
                    writer.write_batch(batch, row_group_size=ROW_GROUP_SIZE)
        finally:
            if writer is not None:
                writer.close()
    unsorted.unlink()
