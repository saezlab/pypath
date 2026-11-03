"""Persistent layout for extracted observations and resolution results.

Observations, one directory per resource::

    <root>/<resource>/queries.parquet      one row per distinct observation
    <root>/<resource>/votes.parquet        every vote row observation_bundle emits
    <root>/<resource>/occurrences.parquet  (input_id, dataset, n) provenance
    <root>/<resource>/extract.json         parameters, counts, per-dataset status

Both ``queries`` and ``votes`` are sorted by ``input_id`` so a resolve run streams
them side by side in bounded memory. ``input_id`` is the observation fingerprint:
a hash of the resolution input (target plus all vote rows), so identical
observations from different datasets or records collapse to one row.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

FINGERPRINT_VERSION = 1
LIBRARY_BY_TARGET = {1: "chemical", 2: "gene_protein", 3: "reaction"}

QUERY_SCHEMA = pa.schema(
    [
        ("input_id", pa.string()),
        ("target", pa.int32()),
        ("library", pa.string()),
        ("entity_type", pa.string()),
        ("namespace", pa.string()),
        ("identifier", pa.string()),
        ("taxon", pa.string()),
    ]
)
VOTE_SCHEMA = pa.schema(
    [
        ("input_id", pa.string()),
        ("ns", pa.string()),
        ("identifier", pa.string()),
        ("scope", pa.string()),
        ("anchor", pa.string()),
        ("target", pa.int32()),
        ("route", pa.int32()),
        ("ordinal", pa.int32()),
        ("primary", pa.bool_()),
        ("gene_only", pa.bool_()),
        ("lookup_key", pa.binary()),
    ]
)
OCCURRENCE_SCHEMA = pa.schema(
    [("input_id", pa.string()), ("dataset", pa.string()), ("n", pa.int64())]
)
RESULT_SCHEMA = pa.schema(
    [
        ("input_id", pa.string()),
        ("outcome", pa.string()),
        ("entities", pa.list_(pa.string())),
        ("gene_mapping_status", pa.string()),
        ("gene_candidates", pa.list_(pa.string())),
        ("protein_entity_id", pa.string()),
        ("candidate_count", pa.int64()),
        ("entity_labels", pa.list_(pa.string())),
    ]
)
_SPECS = {
    "queries": QUERY_SCHEMA,
    "votes": VOTE_SCHEMA,
    "occurrences": OCCURRENCE_SCHEMA,
}


def observation_fingerprint(target: int, rows: list[dict[str, Any]]) -> str:
    """Stable id of one resolution input: target plus every vote row.

    ``input_id``, ``lookup_key`` (derived from the other fields) and everything
    descriptive about the source record are excluded; ``primary`` and
    ``gene_only`` are included so a change in either can never be merged away.
    """
    signature = sorted(
        (
            r["ns"],
            r["identifier"],
            r.get("scope", "") or "",
            int(r["route"]),
            r.get("anchor", "") or "",
            bool(r.get("primary", False)),
            bool(r.get("gene_only", False)),
        )
        for r in rows
    )
    payload = json.dumps([FINGERPRINT_VERSION, int(target), signature], separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()[:32]


def vote_row(input_id: str, row: dict[str, Any]) -> dict[str, Any]:
    """Normalize one bundle row to the stored vote schema."""
    return dict(
        input_id=input_id,
        ns=row["ns"],
        identifier=row["identifier"],
        scope=row.get("scope", "") or "",
        anchor=row.get("anchor", "") or "",
        target=int(row["target"]),
        route=int(row["route"]),
        ordinal=int(row["ordinal"]),
        primary=bool(row.get("primary", False)),
        gene_only=bool(row.get("gene_only", False)),
        lookup_key=bytes(row["lookup_key"]),
    )


class ObservationWriter:
    """Accumulate observations into shard files, then merge them sorted."""

    def __init__(self, directory: str | Path, *, flush_rows: int = 200_000):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.parts = self.directory / ".parts"
        if self.parts.exists():
            shutil.rmtree(self.parts)
        self.parts.mkdir(parents=True)
        self.flush_rows = flush_rows
        self._buffers: dict[str, list[dict[str, Any]]] = {k: [] for k in _SPECS}
        self._shards = 0
        for kind, schema in _SPECS.items():
            pq.write_table(
                pa.Table.from_pylist([], schema=schema), self.parts / f"{kind}-empty.parquet"
            )

    def add(self, input_id: str, query: dict[str, Any], rows: list[dict[str, Any]]) -> None:
        target = int(query["target"])
        self._buffers["queries"].append(
            dict(
                input_id=input_id,
                target=target,
                library=LIBRARY_BY_TARGET.get(target, str(target)),
                entity_type=str(query.get("entity_type") or ""),
                namespace=str(query.get("namespace") or ""),
                identifier=str(query.get("identifier") or ""),
                taxon=str(query.get("taxon") or ""),
            )
        )
        self._buffers["votes"].extend(vote_row(input_id, r) for r in rows)
        self._maybe_flush()

    def add_occurrences(self, dataset: str, counts: dict[str, int]) -> None:
        self._buffers["occurrences"].extend(
            dict(input_id=i, dataset=dataset, n=n) for i, n in counts.items()
        )
        self._maybe_flush()

    def _maybe_flush(self, force: bool = False) -> None:
        if not force and max(len(b) for b in self._buffers.values()) < self.flush_rows:
            return
        for kind, rows in self._buffers.items():
            if rows:
                pq.write_table(
                    pa.Table.from_pylist(rows, schema=_SPECS[kind]),
                    self.parts / f"{kind}-{self._shards}.parquet",
                    compression="zstd",
                )
                rows.clear()
        self._shards += 1

    def finalize(self, *, memory_limit: str = "2GB", threads: int = 2) -> dict[str, int]:
        """Merge shards into the sorted final files; returns row counts."""
        self._maybe_flush(force=True)
        con = duckdb.connect()
        con.execute(f"SET memory_limit='{memory_limit}'")
        con.execute(f"SET threads={int(threads)}")
        con.execute("SET preserve_insertion_order=false")
        spill = self.directory / ".spill"
        spill.mkdir(exist_ok=True)
        con.execute(f"SET temp_directory='{spill}'")
        p = self.parts
        copy = "(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)"
        out = self.directory

        def to_final(name: str, sql: str) -> None:
            tmp = out / f".{name}.parquet.tmp"
            con.execute(f"COPY ({sql}) TO '{tmp}' {copy}")
            tmp.replace(out / f"{name}.parquet")

        to_final(
            "queries",
            f"SELECT * FROM read_parquet('{p}/queries-*.parquet') ORDER BY input_id",
        )
        to_final(
            "votes",
            f"SELECT * FROM read_parquet('{p}/votes-*.parquet') ORDER BY input_id, ordinal",
        )
        to_final(
            "occurrences",
            "SELECT input_id, dataset, sum(n)::BIGINT AS n "
            f"FROM read_parquet('{p}/occurrences-*.parquet') "
            "GROUP BY input_id, dataset ORDER BY input_id, dataset",
        )
        con.close()
        shutil.rmtree(self.parts, ignore_errors=True)
        shutil.rmtree(spill, ignore_errors=True)
        return {kind: pq.ParquetFile(out / f"{kind}.parquet").metadata.num_rows for kind in _SPECS}


def resource_dir(root: str | Path, resource: str) -> Path:
    return Path(root) / resource


def list_resources(root: str | Path) -> list[str]:
    """Resources that have a finished extraction (queries and votes present)."""
    root = Path(root)
    if not root.is_dir():
        return []
    return sorted(
        p.name
        for p in root.iterdir()
        if (p / "queries.parquet").is_file() and (p / "votes.parquet").is_file()
    )


def read_extract_info(root: str | Path, resource: str) -> dict[str, Any]:
    path = resource_dir(root, resource) / "extract.json"
    return json.loads(path.read_text()) if path.is_file() else {}


DEFAULT_SAMPLE = 5000


def sample_select_sql(queries_path: str | Path, per_library: int) -> str:
    """SQL selecting the sampled input_ids: ``per_library`` per library.

    Deterministic: observations are ranked by md5 of their fingerprint, so the same
    observations are picked on every run and for every runtime; chemical and
    gene_protein are sampled independently, hence both appear.
    """
    path = str(queries_path).replace("'", "''")
    return (
        "SELECT input_id FROM (SELECT input_id, row_number() OVER ("
        "PARTITION BY library ORDER BY md5(input_id), input_id) AS rn "
        f"FROM read_parquet('{path}')) WHERE rn <= {int(per_library)}"
    )


def _iter_sample_batches(directory: Path, batch_size: int, sample: int):
    con = duckdb.connect()
    try:
        con.execute(
            "CREATE TEMP TABLE keep AS " + sample_select_sql(directory / "queries.parquet", sample)
        )
        rows = (
            con.execute(
                f"SELECT q.* FROM read_parquet('{directory}/queries.parquet') q "
                "JOIN keep USING (input_id) ORDER BY q.input_id"
            )
            .fetch_arrow_table()
            .to_pylist()
        )
        table = con.execute(
            f"SELECT v.* FROM read_parquet('{directory}/votes.parquet') v "
            "JOIN keep USING (input_id) ORDER BY v.input_id, v.ordinal"
        ).fetch_arrow_table()
        votes: dict[str, list[dict[str, Any]]] = {}
        for vote in table.to_pylist():
            votes.setdefault(vote["input_id"], []).append(vote)
    finally:
        con.close()
    for offset in range(0, len(rows), batch_size):
        chunk = rows[offset : offset + batch_size]
        yield chunk, [v for q in chunk for v in votes.get(q["input_id"], [])]


def iter_resolve_batches(
    directory: str | Path,
    batch_size: int = 5000,
    *,
    vote_chunk: int = 50_000,
    sample: int | None = None,
) -> Iterator[tuple[list[dict[str, Any]], list[dict[str, Any]]]]:
    """Yield (queries, votes) batches; votes of exactly the batch's queries.

    Both files are sorted by ``input_id``, so one forward pass over each suffices.
    With ``sample`` only that many observations per library are yielded (see
    ``sample_select_sql``); a sample is small enough to be read at once.
    """
    directory = Path(directory)
    if sample:
        yield from _iter_sample_batches(directory, batch_size, sample)
        return
    query_file = pq.ParquetFile(directory / "queries.parquet")
    vote_batches = pq.ParquetFile(directory / "votes.parquet").iter_batches(batch_size=vote_chunk)
    pending: list[dict[str, Any]] = []
    exhausted = False

    def refill() -> bool:
        nonlocal exhausted
        if exhausted:
            return False
        try:
            pending.extend(next(vote_batches).to_pylist())
            return True
        except StopIteration:
            exhausted = True
            return False

    for record_batch in query_file.iter_batches(batch_size=batch_size):
        queries = record_batch.to_pylist()
        last = queries[-1]["input_id"]
        # Votes are sorted: keep reading until a vote beyond the batch shows up.
        while (not pending or pending[-1]["input_id"] <= last) and refill():
            pass
        cut = len(pending)
        for i, vote in enumerate(pending):
            if vote["input_id"] > last:
                cut = i
                break
        votes, pending[:] = pending[:cut], pending[cut:]
        yield queries, votes
