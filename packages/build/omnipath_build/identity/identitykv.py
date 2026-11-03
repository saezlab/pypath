"""LMDB point-lookup store of the identity decisions (spec 3b): ``<identity dir>/kv/``.

One environment with seven small dbs written from the decision Parquet files (see
``omnipath_resolver.identity_kv`` for the layout): ``exc``, ``exc_members``, ``extra``,
``gp_protein``, ``gp_gene``, ``lipid`` and ``cand`` (from ``record_candidates.parquet``; empty
for snapshots built before that file existed). It is a separate step so existing decisions need not be
rebuilt; ``decisions.py`` (whose code hash names the decisions) is not touched.
"""

from __future__ import annotations

import json
import shutil
import time
from pathlib import Path

import lmdb

from omnipath_resolver.identity_kv import (
    DECISION_DBS,
    DECISIONS_KV_FORMAT,
    KV_DIR,
    MANIFEST,
    MAX_KEY,
    pack,
)

from .common import log

MAP_SIZE = 64 * 1024**3


# db name -> (file, key column, value columns, grouped). Every db is streamed from DuckDB in key
# order (bytewise VARCHAR order = LMDB's memcmp order) and appended in batches, so memory stays
# small whatever the size of the decision tables.
SOURCES = {
    "exc": ("exceptions.parquet", "record_id", ("entity_id", "decision", "quarantined"), False),
    "exc_members": ("exception_members.parquet", "entity_id", ("record_id",), True),
    "extra": (
        "entities_extra.parquet",
        "entity_id",
        ("kind", "taxon", "quarantined", "preferred_record"),
        False,
    ),
    "gp_protein": (
        "gene_products_by_protein.parquet",
        "protein_entity_id",
        ("entrez_id", "taxon"),
        True,
    ),
    "gp_gene": ("gene_products_by_gene.parquet", "entrez_id", ("protein_entity_id", "taxon"), True),
    "lipid": ("lipid_structures.parquet", "goslin", ("inchikey",), False),
    "cand": ("record_candidates.parquet", "record_id", ("entity_id",), True),
}
OPTIONAL = {"cand"}  # absent from snapshots built before the rule
BATCH_ITEMS = 50_000


def _value(name, values):
    if name == "exc":
        entity_id, decision, quarantined = values
        return [entity_id, decision, bool(quarantined)]
    if name == "extra":
        kind, taxon, quarantined, preferred = values
        return [kind, taxon, bool(quarantined), preferred]
    return values[0] if len(values) == 1 else list(values)


def stream(identity_dir: Path, name: str, memory: str = "1GB"):
    """Yield (key, value) for one db in ascending key order, grouped where the db is a list."""
    import duckdb

    file, key_col, value_cols, grouped = SOURCES[name]
    path = Path(identity_dir) / file
    if not path.is_file():
        if name in OPTIONAL:
            return
        raise FileNotFoundError(f"Identity decisions are incomplete: {path} is missing")
    c = duckdb.connect()
    c.execute(f"SET memory_limit='{memory}'")
    c.execute("SET threads=2")
    c.execute("SET preserve_insertion_order=false")
    columns = ",".join([key_col, *value_cols])
    order = ",".join([key_col, *value_cols])
    reader = c.execute(
        f"SELECT DISTINCT {columns} FROM read_parquet(?) WHERE {key_col} IS NOT NULL ORDER BY {order}",
        [str(path)],
    ).fetch_record_batch(100_000)
    current, group = None, []
    for batch in reader:
        cols = [col.to_pylist() for col in batch.columns]
        for row in zip(*cols):
            key, values = row[0], row[1:]
            if not grouped:
                if key == current:
                    raise ValueError(f"Duplicate key {key!r} in {file}")
                current = key
                yield key, _value(name, values)
                continue
            if key != current and current is not None:
                yield current, group
                group = []
            current = key
            group.append(_value(name, values))
    if grouped and current is not None:
        yield current, group
    c.close()


def build_identity_kv(identity_dir, *, force: bool = False, min_free_gib: float = 1) -> dict:
    """Write ``<identity_dir>/kv/`` from the decision Parquet files (idempotent)."""
    identity_dir = Path(identity_dir)
    final = identity_dir / KV_DIR
    if (final / MANIFEST).is_file() and not force:
        log("identity_kv_exists", path=str(final))
        return json.loads((final / MANIFEST).read_text())
    identity = json.loads((identity_dir / "manifest.json").read_text())
    if shutil.disk_usage(identity_dir).free < min_free_gib * 1024**3:
        raise RuntimeError("Identity kv build stopped at the free-disk reserve")
    started = time.monotonic()
    building = identity_dir / (KV_DIR + ".building")
    shutil.rmtree(building, ignore_errors=True)
    building.mkdir(parents=True)
    env = lmdb.open(
        str(building), map_size=MAP_SIZE, max_dbs=len(DECISION_DBS) + 1, sync=False, metasync=False
    )
    counts, skipped = {}, {}
    for name in DECISION_DBS:
        db = env.open_db(name.encode())
        added_total = long_keys = 0
        previous = None
        batch = []

        def flush(batch):
            while True:
                try:
                    with env.begin(db=db, write=True) as txn:
                        _, added = txn.cursor(db).putmulti(batch, dupdata=False, append=True)
                    break
                except lmdb.MapFullError:
                    env.set_mapsize(env.info()["map_size"] * 2)
            if added != len(batch):
                raise AssertionError(f"{name}: stored {added} of {len(batch)} keys")
            return added

        for key, value in stream(identity_dir, name):
            raw = key.encode()
            if previous is not None and raw <= previous:
                raise AssertionError(f"{name}: keys not in ascending byte order at {key!r}")
            previous = raw
            if len(raw) > MAX_KEY:
                long_keys += 1
                continue
            batch.append((raw, pack(value)))
            if len(batch) >= BATCH_ITEMS:
                added_total += flush(batch)
                batch = []
        if batch:
            added_total += flush(batch)
        counts[name], skipped[name] = added_total, long_keys
    env.sync(True)
    env.close()
    manifest = dict(
        format=DECISIONS_KV_FORMAT,
        fingerprint=identity.get("fingerprint"),
        counts=counts,
        skipped_long_keys=skipped,
        file_bytes=(building / "data.mdb").stat().st_size,
        seconds=round(time.monotonic() - started, 2),
    )
    (building / MANIFEST).write_text(json.dumps(manifest, indent=2))
    if final.exists():
        shutil.rmtree(final)
    building.rename(final)
    log("identity_kv_done", path=str(final), **{k: v for k, v in manifest.items() if k != "format"})
    return manifest
