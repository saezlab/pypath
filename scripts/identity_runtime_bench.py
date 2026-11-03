"""Micro-benchmark for IdentityRuntime: lookup_many at 1k/5k/20k keys, record_many cold/warm.

Usage:
    uv run python scripts/identity_runtime_bench.py SNAPSHOT [--keys FILE] [--batches 1000 5000 20000]

SNAPSHOT is an ``<root>/identity/<fingerprint>/`` directory (manifest format omnipath-identity-v2).
Keys come from ``--keys FILE`` or, by default, a random sample drawn from the hub indexes' by_id tables.
FILE has one lookup key per line, either hex of the raw key bytes (``010101...``) or five
tab-separated fields ``target route ns scope identifier`` (scope empty for unscoped).

Each batch size is timed cold (first call of that batch) and warm (immediate repeat, page cache
and DuckDB metadata hot). record_many is timed cold against an empty throwaway cache, then warm
(served from SQLite); use ``--cache-dir`` to keep and reuse a real cache instead.
"""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path
import random
import shutil
import statistics
import tempfile
import time

from omnipath_resolver.identity_runtime import IdentityRuntime
from omnipath_resolver.observations import CODES, key


def read_key_file(path: Path) -> list[bytes]:
    keys = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) == 5:
            target, route, ns, scope, identifier = fields
            keys.append(key(int(target), int(route), ns, scope, identifier))
        else:
            keys.append(bytes.fromhex(line))
    return keys


def sample_keys(
    snapshot: Path, count: int, partitions: int, scoped: bool, seed: int
) -> list[bytes]:
    """Draw ``count`` keys from random ``by_id`` partition files of the hub indexes."""
    rng = random.Random(seed)
    import duckdb

    manifest = json.loads((snapshot / "manifest.json").read_text())
    files = []
    for hub, directory in manifest["hub_indexes"].items():
        directory = Path(directory)
        directory = directory if directory.is_absolute() else snapshot / directory
        files += [(hub, f) for f in sorted(glob.glob(str(directory / "by_id/part=*/*.parquet")))]
    if not files:
        raise SystemExit(f"no by_id files under the hub indexes of {snapshot}")
    chosen = rng.sample(files, min(partitions, len(files)))
    per_file = -(-count * 2 // len(chosen))  # oversample, then trim after shuffling
    con = duckdb.connect()
    keys = set()
    for hub, path in chosen:
        target = 2 if hub in ("uniprot", "entrez", "ramp_gene") else 1
        rows = con.execute(
            f"SELECT ns, identifier, CAST(taxon AS VARCHAR) FROM read_parquet(?) "
            f"USING SAMPLE {per_file} ROWS (reservoir, {rng.randrange(1 << 30)})",
            [path],
        ).fetchall()
        for ns, identifier, taxon in rows:
            if ns not in CODES:
                continue
            route = 2 if target == 2 and ns in ("entrez", "ramp_gene", "kegg_gene") else 1
            known = bool(taxon and taxon.isdigit() and taxon != "0")
            use_scope = known and (scoped or ns in ("genesymbol", "genesymbol-syn"))
            keys.add(key(target, route, ns, taxon if use_scope else "", identifier))
    keys = sorted(keys)
    rng.shuffle(keys)
    if len(keys) < count:
        print(f"warning: only {len(keys)} distinct sampled keys (< {count}); raise --partitions")
    return keys[:count]


def timed(function, *args):
    start = time.perf_counter()
    result = function(*args)
    return result, time.perf_counter() - start


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--keys", type=Path, help="file of lookup keys instead of sampling")
    parser.add_argument("--batches", type=int, nargs="+", default=[1000, 5000, 20000])
    parser.add_argument("--records", type=int, default=5000, help="entities for record_many")
    parser.add_argument("--partitions", type=int, default=64, help="partitions to sample from")
    parser.add_argument("--scoped", action="store_true", help="sampled gene keys taxon-scoped")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--cache-dir", type=Path, help="default: a throwaway directory")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--memory-limit")
    parser.add_argument("--json", type=Path, help="also write the results here")
    args = parser.parse_args()

    need = max(args.batches)
    started = time.perf_counter()
    keys = (
        read_key_file(args.keys)
        if args.keys
        else sample_keys(args.snapshot, need, args.partitions, args.scoped, args.seed)
    )
    print(f"{len(keys)} keys ready in {time.perf_counter() - started:.1f}s")
    if len(keys) < need:
        args.batches = [b for b in args.batches if b <= len(keys)] or [len(keys)]

    cache = args.cache_dir or Path(tempfile.mkdtemp(prefix="identity-bench-"))
    runtime = IdentityRuntime(
        args.snapshot, cache_dir=cache, threads=args.threads, memory_limit=args.memory_limit
    )
    report = {
        "snapshot": str(args.snapshot),
        "fingerprint": runtime.fingerprint,
        "lookup": [],
        "record": [],
    }
    last = {}
    try:
        for size in args.batches:
            batch = keys[:size]
            last, cold = timed(runtime.lookup_many, batch)
            _, warm = timed(runtime.lookup_many, batch)
            candidates = sum(len(v["candidates"]) for v in last.values())
            hits = sum(1 for v in last.values() if v["candidates"])
            widest = max((len(v["candidates"]) for v in last.values()), default=0)
            row = dict(
                keys=size,
                cold_s=cold,
                warm_s=warm,
                cold_keys_per_s=size / cold,
                warm_keys_per_s=size / warm,
                hits=hits,
                candidates=candidates,
                widest=widest,
            )
            report["lookup"].append(row)
            print(
                f"lookup_many {size:>6}: cold {cold:7.2f}s ({size / cold:9.0f}/s)  "
                f"warm {warm:7.2f}s ({size / warm:9.0f}/s)  hits {hits}  "
                f"candidates {candidates}  widest posting {widest}"
            )

        entities = sorted({c[1] for v in last.values() for c in v["candidates"]})
        random.Random(args.seed).shuffle(entities)
        entities = entities[: args.records]
        if entities:
            for label, size in (("cold", len(entities)), ("warm", len(entities))):
                _, seconds = timed(runtime.record_many, entities)
                report["record"].append(dict(phase=label, entities=size, seconds=seconds))
                print(f"record_many {size:>6} {label}: {seconds:7.2f}s ({size / seconds:9.0f}/s)")
            singles = entities[:200]
            per_call = statistics.mean(timed(runtime.record, e)[1] for e in singles)
            print(f"record (single, warm) mean {per_call * 1000:.3f} ms over {len(singles)}")
        else:
            print("no candidates found: nothing to build records for")
    finally:
        runtime.close()
        if args.cache_dir is None:
            shutil.rmtree(cache, ignore_errors=True)
    if args.json:
        args.json.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
