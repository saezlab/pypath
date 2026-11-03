"""LMDB point-lookup store of one hub index (spec 3a): ``<index dir>/kv/``.

The Parquet index answers a batch of point lookups by scanning most of its 256 partitions; the
runtime reads this store instead. The Parquet files stay the build artifact (untouched); the
layout is described in ``omnipath_resolver.identity_kv``:

    id.<s>/   db ``id``:  part + ns + 0x00 + identifier [+ 0x00 + taxon] -> "local_id\x1ftag\x1e..."
    rec.<s>/  db ``rec``: part + local_id -> "taxon\x1fanchor\x1fcount\x1freviewed[\x1etype\x1fvalue...]"
    manifest.json

The taxon is part of the key for the symbol namespaces only (format v3): a gene symbol is the same
word in hundreds of species and is always looked up for one taxon, so its rows are split by the
taxon of their record (empty when the record has none) and a lookup reads only its own.

``part`` is the two-hex-digit partition the Parquet index already uses (md5 of the identifier, or
of ``<hub>:<local_id>``), and each partition is one file already sorted by key. Prefixing ``part``
makes the keys of partition 00, 01, ... ff globally ascending, so the build is: read a partition,
group equal keys and join their rows with Arrow (a handful of vectorized calls per partition), and
append. There is no sort, no join and no per-row Python. The only other step is splitting
records.parquet into the same partitions once (DuckDB, for md5). Large hubs are split into 16
shards by the first hex digit of ``part``, built in parallel processes. This module does not touch
``hubindex.py``, whose code hash names the indexes.
"""

from __future__ import annotations

import json
import multiprocessing
import shutil
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import duckdb
import lmdb
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

from omnipath_resolver.identity_kv import (
    FIELD,
    HUB_KV_FORMAT,
    KV_DIR,
    MANIFEST,
    MAX_KEY,
    ROW,
    SYMBOL_NAMESPACES,
    sha256_file,
)

from .common import log

# A hub whose by_id (by_record) holds more rows than this gets 16 id (rec) shards.
SHARD_ROWS = 30_000_000
PARTS = [f"{n:02x}" for n in range(256)]
MAP_SIZE = 64 * 1024**3  # virtual; the file grows with the data, and doubles on MapFullError
PUT_BATCH = 100_000


def quote(value) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def find_index(hub_index_root, hub) -> Path:
    """The newest ``<root>/<hub>/<sha12>/`` directory holding a finished index."""
    candidates = sorted(
        (Path(hub_index_root) / hub).glob("*/manifest.json"), key=lambda p: p.stat().st_mtime
    )
    if not candidates:
        raise FileNotFoundError(f"No finished hub index for {hub} under {hub_index_root}")
    return candidates[-1].parent


def _files(root: Path, part: str) -> list[str]:
    return sorted(str(p) for p in (root / f"part={part}").glob("*.parquet"))


def _parts(shard: int, shards: int) -> list[str]:
    return PARTS if shards == 1 else [p for p in PARTS if int(p[0], 16) == shard]


def _row_count(directory: Path, name: str, default: int) -> int:
    try:
        counts = json.loads((directory / "manifest.json").read_text()).get("counts", {})
    except (OSError, ValueError):
        return default
    if name == "by_id":
        return int((counts.get("by_id") or {}).get("rows", default))
    return int((counts.get("final") or {}).get("by_record", default))


def _connect(spill: Path, memory: str, threads: int):
    spill.mkdir(parents=True, exist_ok=True)
    c = duckdb.connect()
    c.execute(f"SET memory_limit={quote(memory)}")
    c.execute(f"SET threads={int(threads)}")
    c.execute(f"SET temp_directory={quote(spill)}")
    return c


# ------------------------------------------------------------------------ grouping


def _sorted(table: pa.Table, keys: list[str]) -> pa.Table:
    """The partition files are written sorted; sort only if one is not (cheap, in memory)."""
    key = (
        table.column(keys[0])
        if len(keys) == 1
        else pc.binary_join_element_wise(*[table.column(k) for k in keys[:2]], "\x00")
    )
    key = key.combine_chunks() if hasattr(key, "combine_chunks") else key
    if len(key) < 2 or pc.all(pc.less_equal(key[:-1], key[1:])).as_py():
        return table
    return table.sort_by([(k, "ascending") for k in keys])


def join(*columns, sep: str) -> pa.Array:
    """Element-wise ``sep.join`` of string columns; nulls become ''."""
    return pc.binary_join_element_wise(*columns, sep, null_handling="replace", null_replacement="")


def group(keys: pa.Array, frags: pa.Array):
    """Rows sorted by key -> (key per run of equal keys, its fragments joined with ROW), in bulk."""
    n = len(keys)
    if n == 0:
        return keys, frags
    change = pc.not_equal(keys[1:], keys[:-1]).to_numpy(zero_copy_only=False)
    starts = np.concatenate([[0], np.flatnonzero(change) + 1])
    offsets = pa.array(np.concatenate([starts, [n]]).astype(np.int32))
    joined = pc.binary_join(pa.ListArray.from_arrays(offsets, frags), ROW)
    return pc.take(keys, pa.array(starts)), joined


# ------------------------------------------------------------------------- writing


class Writer:
    """Append-only writer of one LMDB database; batches must arrive in increasing key order."""

    def __init__(self, path: Path, name: str, min_free: int):
        path.mkdir(parents=True, exist_ok=True)
        self.path, self.min_free = path, min_free
        self.env = lmdb.open(
            str(path), map_size=MAP_SIZE, max_dbs=2, writemap=False, sync=False, metasync=False
        )
        self.db = self.env.open_db(name.encode())
        self.max_key = min(MAX_KEY, self.env.max_key_size())
        self.last = b""
        self.keys = self.value_bytes = self.skipped = 0

    def add(self, keys: pa.Array, values: pa.Array) -> None:
        """Append a sorted batch (Arrow string arrays of equal length)."""
        if len(keys) == 0:
            return
        keys, values = pc.cast(keys, pa.binary()), pc.cast(values, pa.binary())
        fits = pc.less_equal(pc.binary_length(keys), self.max_key)
        self.skipped += len(keys) - pc.sum(fits).as_py()
        keys, values = pc.filter(keys, fits), pc.filter(values, fits)
        self.value_bytes += pc.sum(pc.binary_length(values)).as_py() or 0
        for start in range(0, len(keys), PUT_BATCH):
            pairs = list(
                zip(
                    keys.slice(start, PUT_BATCH).to_pylist(),
                    values.slice(start, PUT_BATCH).to_pylist(),
                )
            )
            if pairs[0][0] <= self.last:
                raise AssertionError(
                    f"Keys not ascending at {pairs[0][0][:60]!r} after {self.last[:60]!r}"
                )
            self._put(pairs)
            self.last = pairs[-1][0]
            self.keys += len(pairs)
        # Flush once per partition: dirty pages count against the memory cgroup, and with swap
        # off they otherwise pile up until the kernel kills a worker.
        self.env.sync(True)

    def _put(self, pairs) -> None:
        if shutil.disk_usage(self.path).free < self.min_free:
            raise RuntimeError("Hub kv build stopped at the free-disk reserve")
        while True:
            try:
                with self.env.begin(db=self.db, write=True) as txn:
                    _, added = txn.cursor(self.db).putmulti(pairs, dupdata=False, append=True)
                break
            except lmdb.MapFullError:
                self.env.set_mapsize(self.env.info()["map_size"] * 2)
        if added != len(pairs):  # append mode refuses a key that is not strictly increasing
            raise AssertionError(f"{self.path}: stored {added} of {len(pairs)} keys (order?)")

    def close(self) -> dict:
        self.env.sync(True)
        with self.env.begin() as txn:
            entries = txn.stat(self.db)["entries"]
        self.env.close()
        if entries != self.keys:
            raise AssertionError(f"{self.path} holds {entries} keys, {self.keys} were written")
        return dict(
            keys=self.keys,
            value_bytes=self.value_bytes,
            skipped_long_keys=self.skipped,
            file_bytes=(self.path / "data.mdb").stat().st_size,
        )


# -------------------------------------------------------------------------- shards


def build_rec_shard(directory, building, target, shard, shards, hub, min_free):
    """Per partition: every record (head fields) with its by_record rows, one value per record."""
    writer = Writer(target, "rec", min_free)
    rows = orphans = 0
    for part in _parts(shard, shards):
        records = _files(building / "records_parts", part)
        if not records:
            continue
        rec = _sorted(pq.read_table(records), ["local_id"])
        ids = rec.column("local_id").combine_chunks()
        head = join(
            rec.column("taxon"),
            rec.column("anchor"),
            pc.cast(rec.column("anchor_count"), pa.string()),
            pc.if_else(rec.column("reviewed"), "1", "0"),
            sep=FIELD,
        ).combine_chunks()
        files = _files(directory / "by_record", part)
        if files:
            body = _sorted(
                pq.read_table(files, columns=["local_id", "source_type", "value"]), ["local_id"]
            )
            rows += body.num_rows
            keys, joined = group(
                body.column("local_id").combine_chunks(),
                join(body.column("source_type"), body.column("value"), sep=FIELD).combine_chunks(),
            )
            position = pc.index_in(ids, value_set=keys)  # each record's row group, or null
            orphans += len(keys) - pc.count(position).as_py()
            value = pc.binary_join_element_wise(
                head, pc.take(joined, position), ROW, null_handling="skip"
            )
        else:
            value = head
        writer.add(join(pa.scalar(part), ids, sep=""), value)
    stats = writer.close()
    if orphans:
        log("hub_kv_orphan_rows", env=target.name, records_without_head=orphans)
    return dict(stats, rows=rows, orphan_row_groups=orphans)


def taxon_of(values: pa.Array) -> pa.Array:
    """Taxon as the runtime reads it (``identity_runtime._taxon``): '' for none, '0' or 'None'."""
    text = pc.utf8_trim_whitespace(pc.cast(values, pa.string()))
    return pc.if_else(pc.is_in(text, value_set=pa.array(["0", "None"])), "", pc.fill_null(text, ""))


def build_id_shard(directory, building, target, shard, shards, hub, min_free):
    """Per partition: by_id rows grouped by key, one value per key.

    The key is (ns, identifier), and for the symbol namespaces (ns, identifier, taxon).
    """
    writer = Writer(target, "id", min_free)
    rows = 0
    symbols = pa.array(sorted(SYMBOL_NAMESPACES))
    for part in _parts(shard, shards):
        files = _files(directory / "by_id", part)
        if not files:
            continue
        table = pq.read_table(files, columns=["ns", "identifier", "local_id", "tag", "taxon"])
        rows += table.num_rows
        # ns + 0x00 + identifier sorts exactly like (ns, identifier): 0x00 is the smallest byte,
        # and so does the key with a taxon after another 0x00.
        ns, identifier = table.column("ns"), table.column("identifier")
        key = pc.if_else(
            pc.is_in(ns, value_set=symbols),
            join(ns, identifier, taxon_of(table.column("taxon")), sep="\x00"),
            join(ns, identifier, sep="\x00"),
        )
        table = _sorted(
            pa.table(dict(key=key, local_id=table.column("local_id"), tag=table.column("tag"))),
            ["key", "local_id", "tag"],
        )
        keys, joined = group(
            table.column("key").combine_chunks(),
            join(table.column("local_id"), table.column("tag"), sep=FIELD).combine_chunks(),
        )
        writer.add(join(pa.scalar(part), keys, sep=""), joined)
    return dict(writer.close(), rows=rows)


def _run_shard(args):
    kind, directory, building, shard, shards, hub, min_free = args
    directory, building = Path(directory), Path(building)
    name = f"{kind}.{shard}"
    final = building / name
    stats_file = final / "stats.json"
    if stats_file.is_file():  # finished by an earlier run
        return name, json.loads(stats_file.read_text())
    temp = building / (name + ".tmp")
    shutil.rmtree(temp, ignore_errors=True)
    started = time.monotonic()
    build = build_id_shard if kind == "id" else build_rec_shard
    stats = build(directory, building, temp, shard, shards, hub, min_free)
    stats["seconds"] = round(time.monotonic() - started, 2)
    (temp / "stats.json").write_text(json.dumps(stats))
    shutil.rmtree(building / "spill" / temp.name, ignore_errors=True)
    temp.rename(final)
    log("hub_kv_shard", env=name, **stats)
    return name, stats


def partition_records(
    directory: Path, building: Path, hub: str, memory: str, threads: int
) -> float:
    """Rewrite records.parquet once into the by_record partitions (md5 of '<hub>:<local_id>')."""
    target = building / "records_parts"
    done = building / "records_parts.done"
    if done.is_file():
        return float(done.read_text())
    shutil.rmtree(target, ignore_errors=True)
    started = time.monotonic()
    c = _connect(building / "spill" / "records", memory, threads)
    c.execute("SET preserve_insertion_order=false")
    c.execute("SET partitioned_write_max_open_files=256")
    c.execute(
        f"""COPY (SELECT local_id,taxon,anchor,anchor_count,reviewed,
              substr(md5({quote(hub + ":")} || local_id),1,2) part
            FROM read_parquet({quote(directory / "records.parquet")}))
        TO {quote(target)} (FORMAT PARQUET,COMPRESSION ZSTD,PARTITION_BY(part))"""
    )
    c.close()
    seconds = round(time.monotonic() - started, 2)
    done.write_text(str(seconds))
    return seconds


def build_hub_kv_dir(
    directory,
    *,
    id_shards: int | None = None,
    rec_shards: int | None = None,
    memory: str = "4GB",
    threads: int = 2,
    workers: int = 1,
    min_free_gib: float = 20,
) -> dict:
    """Write ``<directory>/kv/`` for one hub index directory (idempotent, resumable per shard)."""
    directory = Path(directory)
    final = directory / KV_DIR
    if (final / MANIFEST).is_file():
        existing = json.loads((final / MANIFEST).read_text())
        if existing.get("format") == HUB_KV_FORMAT:
            log("hub_kv_exists", path=str(final))
            return existing
        # An older store keeps serving until the new one is complete (see below).
        log("hub_kv_outdated", path=str(final), format=existing.get("format"))
    index_manifest = json.loads((directory / "manifest.json").read_text())
    hub = index_manifest["hub"]
    id_shards = id_shards or (16 if _row_count(directory, "by_id", 0) > SHARD_ROWS else 1)
    rec_shards = rec_shards or (16 if _row_count(directory, "by_record", 0) > SHARD_ROWS else 1)
    if id_shards not in (1, 16) or rec_shards not in (1, 16):
        raise ValueError("Shard counts must be 1 or 16")
    min_free = int(min_free_gib * 1024**3)
    building = directory / (KV_DIR + ".building")
    marker = building / "format"
    if building.is_dir() and not (marker.is_file() and marker.read_text() == HUB_KV_FORMAT):
        shutil.rmtree(building)  # a partial build of an older layout
    building.mkdir(parents=True, exist_ok=True)
    marker.write_text(HUB_KV_FORMAT)
    if shutil.disk_usage(building).free < min_free:
        raise RuntimeError("Hub kv build stopped at the free-disk reserve")
    started = time.monotonic()
    workers = max(1, workers)
    log("hub_kv_start", hub=hub, id_shards=id_shards, rec_shards=rec_shards, workers=workers)
    records_seconds = partition_records(directory, building, hub, memory, threads)
    tasks = [
        (kind, str(directory), str(building), shard, count, hub, min_free)
        for kind, count in (("rec", rec_shards), ("id", id_shards))
        for shard in range(count)
    ]
    if workers == 1:
        results = [_run_shard(t) for t in tasks]
    else:
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
            results = list(pool.map(_run_shard, tasks))
    envs = dict(results)
    shutil.rmtree(building / "spill", ignore_errors=True)
    shutil.rmtree(building / "records_parts", ignore_errors=True)
    totals = {
        kind: {
            field: sum(s[field] for n, s in envs.items() if n.startswith(kind + "."))
            for field in ("keys", "rows", "file_bytes", "skipped_long_keys")
        }
        for kind in ("id", "rec")
    }
    manifest = dict(
        format=HUB_KV_FORMAT,
        hub=hub,
        source_index=str(directory),
        source_sha256=sha256_file(directory / "manifest.json"),
        input=index_manifest.get("input"),
        id_shards=id_shards,
        rec_shards=rec_shards,
        taxon_keyed_namespaces=sorted(SYMBOL_NAMESPACES),
        totals=totals,
        envs=dict(sorted(envs.items())),
        records_partition_seconds=records_seconds,
        seconds=round(time.monotonic() - started, 2),
    )
    (building / MANIFEST).write_text(json.dumps(manifest, indent=2))
    marker.unlink()
    shutil.rmtree(final, ignore_errors=True)  # the older store, if any
    building.rename(final)
    log("hub_kv_done", hub=hub, path=str(final), seconds=manifest["seconds"], totals=totals)
    return manifest


def build_hub_kv(hub, hub_index_root, **options) -> dict:
    """Write the kv store of the newest finished index of ``hub`` under ``hub_index_root``."""
    return build_hub_kv_dir(find_index(hub_index_root, hub), **options)
