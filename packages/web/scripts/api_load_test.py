#!/usr/bin/env python3
"""Closed-loop load test of the API for programmatic clients: many small calls.

N clients each loop without think time over small calls (entity details, the first
relations page of an entity, a mix with small resolves and searches) on a persistent
connection, as a Python or R client or an AI tool does. Entities are drawn without
repeats from thousands of sampled keys, so the API's answer cache mostly misses;
``--warm K`` cycles K parameter sets instead, so it mostly hits. Each level reports
requests per second, latency percentiles, errors (503 is the slot wait timeout),
response sizes and, with ``--cgroup``, the API's CPU milliseconds per request read
from its cgroup's ``cpu.stat``.

Sample parameters from a release database (needs duckdb), then run levels:

    uv run python packages/web/scripts/api_load_test.py sample release.duckdb params.json
    uv run python packages/web/scripts/api_load_test.py run --params params.json \
        --workload entity --clients 1,8,32,64 --duration 30 \
        --cgroup /sys/fs/cgroup/.../jschaul-api-8421.scope --out entity.jsonl

Run it where it does not take the API's CPUs (its own cgroup, a quota of a few CPUs).
"""

import argparse
import http.client
import json
import multiprocessing as mp
import random
import statistics
import sys
import threading
import time
import urllib.parse

WORKLOADS = ("entity", "relations", "mix")
# The mix: share of calls by kind.
MIX = (("entity", 30), ("relations", 30), ("resolve", 15), ("search", 15), ("relationships", 10))


# ---------------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------------


def sample(database, out, n=20000):
    import duckdb

    con = duckdb.connect(database, read_only=True)
    con.execute("SET threads = 8")
    rng = random.Random(7)
    rows = con.execute(
        "SELECT entity_key, entity_n, entity_type, relation_count, label FROM entity "
        "WHERE relation_count > 0 ORDER BY entity_n"
    ).fetchall()
    by_degree = sorted(rows, key=lambda r: -r[3])

    def pack(rs):
        return [{"key": r[0], "type": r[2], "relations": r[3], "label": r[4]} for r in rs]

    genes = con.execute(
        "SELECT label FROM entity WHERE entity_type = 'protein' AND taxon = '9606' "
        "AND relation_count >= 5 AND label IS NOT NULL ORDER BY entity_n"
    ).fetchall()
    acs = con.execute(
        "SELECT i.id FROM entity_identifier i JOIN entity e USING (entity_n) "
        "WHERE e.relation_count > 0 AND i.ns = 'uniprot' ORDER BY i.id"
    ).fetchall()
    params = {
        # database-wide uniform draw (long tail: the median entity has 2 relations)
        "uniform": pack(rng.sample(rows, n)),
        # the top 5 % by relation count (what clients and tools mostly ask about)
        "popular": pack(rng.sample(by_degree[: len(rows) // 20], n)),
        "genes": [g[0] for g in genes],
        "uniprot": [a[0] for a in acs],
        "terms": [r[4] for r in rng.sample(by_degree[: len(rows) // 5], 30000) if r[4]],
    }
    json.dump(params, open(out, "w"))
    for name in ("uniform", "popular"):
        counts = sorted(e["relations"] for e in params[name])
        print(name, "relation_count median", counts[len(counts) // 2], "p95", counts[int(len(counts) * 0.95)])
    print({k: len(v) for k, v in params.items()})


# ---------------------------------------------------------------------------------
# Calls
# ---------------------------------------------------------------------------------


def entity_call(e):
    return "entity", "GET", "/api/entities/" + urllib.parse.quote(e["key"], safe=""), None


def relations_call(e):
    body = {"filters": {"entity_ids": [e["key"]]}, "limit": 20, "offset": 0}
    return "relations", "POST", "/api/relations/search", body


def relationships_call(e):
    return "relationships", "GET", "/api/entities/" + urllib.parse.quote(e["key"], safe="") + "/relationships?limit=20", None


def resolve_call(e, rng, P):
    k = rng.randint(5, 20)
    ids = []
    for _ in range(k):
        r = rng.random()
        ids.append(rng.choice(P["uniprot"]) if r < 0.4 else rng.choice(P["genes"]) if r < 0.8 else e["key"])
    return "resolve", "POST", "/api/entities/resolve", {"identifiers": ids}


def search_call(e, rng, P):
    term = rng.choice(P["terms"])
    term = term[: rng.randint(3, max(3, len(term)))]  # a typed prefix
    return "search", "POST", "/api/entities/search", {"query": term, "limit": 20}


def next_call(workload, e, rng, P):
    if workload == "entity":
        return entity_call(e)
    if workload == "relations":
        return relations_call(e)
    kind = rng.choices([m[0] for m in MIX], [m[1] for m in MIX])[0]
    if kind == "entity":
        return entity_call(e)
    if kind == "relations":
        return relations_call(e)
    if kind == "relationships":
        return relationships_call(e)
    if kind == "resolve":
        return resolve_call(e, rng, P)
    return search_call(e, rng, P)


# ---------------------------------------------------------------------------------
# Load generation
# ---------------------------------------------------------------------------------


def client_loop(host, port, workload, entities, seed, P, stop_at, record_after, out):
    rng = random.Random(seed)
    conn = http.client.HTTPConnection(host, port, timeout=120)
    i = 0
    while time.perf_counter() < stop_at:
        e = entities[i % len(entities)]
        i += 1
        kind, method, path, body = next_call(workload, e, rng, P)
        data = json.dumps(body).encode() if body is not None else None
        headers = {"content-type": "application/json"} if body is not None else {}
        t0 = time.perf_counter()
        try:
            conn.request(method, path, data, headers)
            r = conn.getresponse()
            payload = r.read()
            status, size = r.status, len(payload)
        except (OSError, http.client.HTTPException) as exc:
            status, size = type(exc).__name__, 0
            conn.close()
            conn = http.client.HTTPConnection(host, port, timeout=120)
        t1 = time.perf_counter()
        if t0 >= record_after:
            out.append((kind, t1 - t0, status, size))


def process_main(args, share, queue):
    host, port, workload, P, entities_for, seed0, stop_at, record_after = args
    out = []
    threads = []
    for ci in share:
        t = threading.Thread(
            target=client_loop,
            args=(host, port, workload, entities_for(ci), seed0 + ci, P, stop_at, record_after, out),
        )
        threads.append(t)
        t.start()
    for t in threads:
        t.join()
    queue.put(out)


def cgroup_usec(path):
    with open(path.rstrip("/") + "/cpu.stat") as f:
        return int(next(line for line in f if line.startswith("usage_usec")).split()[1])


def pct(values, p):
    if not values:
        return float("nan")
    values = sorted(values)
    k = (len(values) - 1) * p
    lo, hi = int(k), min(int(k) + 1, len(values) - 1)
    return values[lo] + (values[hi] - values[lo]) * (k - lo)


def pool_entities(P, pool):
    if pool == "mixed":  # alternate the long tail and the popular entities
        return [e for pair in zip(P["uniform"], P["popular"]) for e in pair]
    return P[pool]


def run_level(a, P, clients, level_index):
    pool = pool_entities(P, a.pool)
    # Each level starts further into the pool, and every client takes its own slice,
    # so no entity repeats within a level and levels do not reuse each other's entities.
    # With --warm K every client cycles the same K entities instead.
    if a.warm:
        base = pool[: a.warm]

        def entities_for(ci):
            order = base[ci % len(base):] + base[: ci % len(base)]
            return order
    else:
        start = (level_index * 4001) % max(1, len(pool) // 2)
        per = max(1, (len(pool) - start) // clients)

        def entities_for(ci):
            return pool[start + ci * per: start + (ci + 1) * per]

    procs = min(a.processes, clients)
    shares = [list(range(p, clients, procs)) for p in range(procs)]
    ctx = mp.get_context("fork")
    queue = ctx.Queue()
    now = time.perf_counter()
    stop_at = now + a.ramp + a.duration
    record_after = now + a.ramp
    # fork: perf_counter is system-wide on Linux, so the deadlines carry over.
    ps = []
    for share in shares:
        args = (a.host, a.port, a.workload, P, entities_for, 1000 * level_index, stop_at, record_after)
        p = ctx.Process(target=process_main, args=(args, share, queue))
        p.start()
        ps.append(p)
    # CPU is read over the measured window only.
    time.sleep(max(0, record_after - time.perf_counter()))
    cpu_start = cgroup_usec(a.cgroup) if a.cgroup else None
    t_start = time.perf_counter()
    time.sleep(max(0, stop_at - time.perf_counter()))
    cpu_end = cgroup_usec(a.cgroup) if a.cgroup else None
    window = time.perf_counter() - t_start
    results = []
    for _ in ps:
        results.extend(queue.get())
    for p in ps:
        p.join()
    return summarize(a, clients, results, window, cpu_start, cpu_end)


def summarize(a, clients, results, window, cpu_start, cpu_end):
    ok = [r for r in results if r[2] == 200]
    codes = {}
    for r in results:
        codes[str(r[2])] = codes.get(str(r[2]), 0) + 1
    lat = [r[1] * 1000 for r in ok]
    row = {
        "workload": a.workload,
        "pool": a.pool,
        "warm": a.warm or 0,
        "clients": clients,
        "window_s": round(window, 1),
        "requests": len(results),
        "req_per_s": round(len(results) / window, 1),
        "ok_per_s": round(len(ok) / window, 1),
        "codes": codes,
        "p50_ms": round(pct(lat, 0.5), 1),
        "p95_ms": round(pct(lat, 0.95), 1),
        "p99_ms": round(pct(lat, 0.99), 1),
        "mean_ms": round(statistics.mean(lat), 1) if lat else None,
        "bytes_mean": round(statistics.mean(r[3] for r in ok)) if ok else None,
        "bytes_p95": round(pct([r[3] for r in ok], 0.95)) if ok else None,
    }
    if cpu_start is not None and results:
        # The CPU of the window covers the requests that finished in it, near enough.
        cpu_ms = (cpu_end - cpu_start) / 1000
        row["api_cpu_cores"] = round(cpu_ms / 1000 / window, 2)
        row["cpu_ms_per_req"] = round(cpu_ms / len(results), 2)
    kinds = {}
    for r in ok:
        kinds.setdefault(r[0], []).append(r)
    row["by_kind"] = {
        k: {
            "n": len(v),
            "p50_ms": round(pct([x[1] * 1000 for x in v], 0.5), 1),
            "p95_ms": round(pct([x[1] * 1000 for x in v], 0.95), 1),
            "bytes_mean": round(statistics.mean(x[3] for x in v)),
        }
        for k, v in sorted(kinds.items())
    }
    return row


def run(a):
    P = json.load(open(a.params))
    if a.warm:
        # Fill the answer cache with the warm set before measuring.
        a2 = argparse.Namespace(**{**vars(a), "duration": 5, "ramp": 0})
        run_level(a2, P, min(8, a.warm), 99)
    out = open(a.out, "a") if a.out else None
    for index, clients in enumerate(int(c) for c in a.clients.split(",")):
        row = run_level(a, P, clients, index + (a.level_offset or 0))
        print(json.dumps(row), flush=True)
        if out:
            out.write(json.dumps(row) + "\n")
            out.flush()
        time.sleep(a.pause)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("sample")
    s.add_argument("database")
    s.add_argument("out")
    r = sub.add_parser("run")
    r.add_argument("--host", default="127.0.0.1")
    r.add_argument("--port", type=int, default=8421)
    r.add_argument("--params", required=True)
    r.add_argument("--workload", choices=WORKLOADS, default="entity")
    r.add_argument("--pool", choices=("mixed", "uniform", "popular"), default="mixed")
    r.add_argument("--clients", default="1,8,32,64")
    r.add_argument("--duration", type=float, default=30)
    r.add_argument("--ramp", type=float, default=3, help="seconds before measuring")
    r.add_argument("--pause", type=float, default=3)
    r.add_argument("--warm", type=int, default=0, help="cycle this many parameter sets (cache hits)")
    r.add_argument("--processes", type=int, default=4, help="client processes (threads share a GIL)")
    r.add_argument("--cgroup", help="the API's cgroup directory, for CPU per request")
    r.add_argument("--level-offset", type=int, default=0, help="shifts the sampled entities")
    r.add_argument("--out")
    a = ap.parse_args()
    if a.cmd == "sample":
        sample(a.database, a.out)
    else:
        run(a)


if __name__ == "__main__":
    main()
