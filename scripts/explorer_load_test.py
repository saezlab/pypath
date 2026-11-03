#!/usr/bin/env python3
"""Closed-loop load test of the explorer, replaying a recorded browser session.

Each virtual user repeats a journey (open the explorer, search, open a result,
relations tab, filter by a source, next page) with think time between steps.
Requests in one step are sent together, as the browser sends them; a step's
latency is the time until its slowest request returns, which is what the user
waits for. Each journey searches a new term, so the API's result caches mostly
miss. Stages run one after another and stop when the error rate or the p95 step
latency passes its limit. Results are in docs/reports/nicesrv-explorer-capacity-20261009.md.

    uv run python scripts/explorer_load_test.py --users 8,12,16,24,32 \
        --abort-step-p95 15 --out load.json
"""

import argparse
import asyncio
import json
import random
import statistics
import time

import httpx

BASE = "https://omnipath-metabo-dev.schaul.click"
API = BASE + "/app-api"
TIMEOUT = 60.0

TERMS = """glucose insulin TP53 EGFR cholesterol ATP pyruvate lactate citrate glutamate
dopamine serotonin histamine adrenaline cortisol estradiol testosterone progesterone
retinol folate cobalamin thiamine riboflavin niacin biotin ascorbate tocopherol
caffeine nicotine ethanol acetate butyrate propionate succinate fumarate malate
oxaloacetate glycine alanine serine threonine valine leucine isoleucine proline
tryptophan tyrosine phenylalanine methionine cysteine lysine arginine histidine
asparagine glutamine aspartate ornithine citrulline creatine carnitine taurine
glutathione NADH FAD coenzyme heme bilirubin urea uric acid ammonia
palmitate oleate linoleate arachidonate prostaglandin leukotriene sphingosine
ceramide phosphatidylcholine phosphatidylserine cardiolipin triglyceride
AKT1 MTOR PIK3CA KRAS BRAF MYC JUN FOS STAT3 NFKB1 TNF IL6 IL1B IFNG TGFB1
VEGFA HIF1A PPARG SREBF1 INSR IGF1R GCK HK2 PFKM PKM LDHA SLC2A1 SLC2A4
ACACA FASN HMGCR LDLR APOE APOB PCSK9 CYP3A4 CYP2D6 ABCB1 ABCG2 SLC22A1
GAPDH ENO1 ALDOA TPI1 PGK1 CS IDH1 IDH2 OGDH SDHA FH MDH2 GLS GLUL ASS1
ARG1 OTC NOS2 NOS3 COX2 PTGS2 ALOX5 MAPK1 MAPK3 MAP2K1 RAF1 SRC ABL1 JAK2
CDK2 CDK4 CCND1 RB1 MDM2 BCL2 BAX CASP3 CASP8 FAS PTEN TSC2 AMPK PRKAA1
LKB1 STK11 SIRT1 FOXO1 PGC1A PPARGC1A NRF2 NFE2L2 KEAP1 HMOX1 SOD1 SOD2
CAT GPX1 TXN PRDX1 glycolysis gluconeogenesis lipogenesis apoptosis autophagy
inflammation diabetes obesity cancer Alzheimer Parkinson asthma hypertension
atherosclerosis steatosis insulin resistance oxidative stress TCA cycle
urea cycle beta oxidation pentose phosphate fatty acid synthesis
sphingolipid metabolism bile acid cholate taurocholate glycocholate
melatonin acetylcholine GABA noradrenaline kynurenine quinolinate NAD
adenosine guanosine inosine hypoxanthine xanthine cAMP cGMP IP3
diacylglycerol calcium magnesium zinc iron copper selenium
metformin aspirin ibuprofen atorvastatin simvastatin rapamycin imatinib
gefitinib erlotinib tamoxifen doxorubicin cisplatin paclitaxel methotrexate
warfarin heparin insulin glargine dexamethasone prednisolone
quercetin resveratrol curcumin genistein catechin anthocyanin lycopene
beta-carotene lutein capsaicin piperine sulforaphane
""".split()


def pct(values, p):
    if not values:
        return float("nan")
    values = sorted(values)
    k = (len(values) - 1) * p
    lo, hi = int(k), min(int(k) + 1, len(values) - 1)
    return values[lo] + (values[hi] - values[lo]) * (k - lo)


class Stats:
    def __init__(self):
        self.requests = []  # (endpoint, seconds, status)
        self.steps = []  # (step, seconds, ok)
        self.journeys = 0
        self.status = []

    def summary(self, elapsed):
        lat = [s for _, s, st in self.requests if st == 200]
        step_lat = [s for _, s, ok in self.steps if ok]
        codes = {}
        for _, _, st in self.requests:
            codes[st] = codes.get(st, 0) + 1
        by_endpoint = {}
        for ep, s, st in self.requests:
            if st == 200:
                by_endpoint.setdefault(ep, []).append(s)
        running = [x["queries"]["running"] for x in self.status if "queries" in x]
        waiting = [x["queries"]["waiting"] for x in self.status if "queries" in x]
        cpu = [x.get("cpuLoad") for x in self.status if x.get("cpuLoad") is not None]
        mem = [x.get("memory") for x in self.status if x.get("memory") is not None]
        failed = {}
        for ep, _, st in self.requests:
            if st != 200:
                failed[f"{ep} {st}"] = failed.get(f"{ep} {st}", 0) + 1
        return {
            "failed": failed,
            "elapsed_s": round(elapsed, 1),
            "requests": len(self.requests),
            "req_per_s": round(len(self.requests) / elapsed, 2),
            "journeys": self.journeys,
            "steps": len(self.steps),
            "steps_per_min": round(len(self.steps) / elapsed * 60, 1),
            "codes": codes,
            "error_rate": round(1 - codes.get(200, 0) / max(1, len(self.requests)), 4),
            "req_p50": round(pct(lat, 0.5), 2),
            "req_p95": round(pct(lat, 0.95), 2),
            "step_p50": round(pct(step_lat, 0.5), 2),
            "step_p95": round(pct(step_lat, 0.95), 2),
            "step_p99": round(pct(step_lat, 0.99), 2),
            "step_max": round(max(step_lat), 2) if step_lat else None,
            "server_running_mean": round(statistics.mean(running), 2) if running else None,
            "server_waiting_mean": round(statistics.mean(waiting), 2) if waiting else None,
            "server_waiting_max": max(waiting) if waiting else None,
            "server_cpu_max": max(cpu) if cpu else None,
            "server_mem_max": max(mem) if mem else None,
            "endpoint_p50": {k: round(pct(v, 0.5), 2) for k, v in sorted(by_endpoint.items())},
            "endpoint_p95": {k: round(pct(v, 0.95), 2) for k, v in sorted(by_endpoint.items())},
        }


class TermPool:
    """Hands out search terms; fresh terms first so result caches mostly miss."""

    def __init__(self, terms, seed, reuse):
        self.rng = random.Random(seed)
        self.terms = list(dict.fromkeys(terms))
        self.rng.shuffle(self.terms)
        self.i = 0
        self.reuse = reuse  # probability of picking a popular (already used) term

    def next(self):
        if self.i and self.rng.random() < self.reuse:
            return self.terms[self.rng.randrange(min(self.i, 15))]
        term = self.terms[self.i % len(self.terms)]
        self.i += 1
        return term


async def call(client, stats, method, path, body=None):
    t0 = time.perf_counter()
    status = -1
    data = None
    try:
        if method == "GET":
            r = await client.get(path)
        else:
            r = await client.post(path, json=body)
        status = r.status_code
        if status == 200 and "json" in r.headers.get("content-type", ""):
            data = r.json()
        else:
            await r.aread()
    except httpx.HTTPError as exc:
        status = type(exc).__name__
    dt = time.perf_counter() - t0
    ep = path.split("?")[0].replace(BASE, "")
    stats.requests.append((ep, dt, status))
    return status, data


async def step(client, stats, name, calls):
    t0 = time.perf_counter()
    results = await asyncio.gather(*(call(client, stats, *c) for c in calls))
    ok = all(st == 200 for st, _ in results)
    stats.steps.append((name, time.perf_counter() - t0, ok))
    return results


async def journey(client, stats, rng, term, think):
    async def pause():
        await asyncio.sleep(rng.uniform(*think))

    # 1. Open the explorer.
    await step(
        client,
        stats,
        "open",
        [
            ("GET", BASE + "/explore?release=latest"),
            ("GET", API + "/status"),
            ("GET", API + "/entities/examples"),
            ("POST", API + "/entities/scoped-facets", {"sources": [], "facetLimit": 100}),
        ],
    )
    await pause()

    # 2. Search.
    res = await step(
        client,
        stats,
        "search",
        [
            (
                "POST",
                API + "/entities/scoped-facets",
                {"sources": [], "query": term, "facetLimit": 100},
            ),
            (
                "POST",
                API + "/entities/groups",
                {"strategy": "auto", "query": term, "filters": {"sources": []}, "member_limit": 1},
            ),
        ],
    )
    groups = (res[1][1] or {}).get("groups") or []
    if not groups:
        return
    await pause()

    # 3. Open one of the top results.
    group = rng.choice(groups[:5])
    key = group["group_key"]
    strategy = (group.get("entity") or {}).get("groupStrategy") or "auto"
    res = await step(
        client,
        stats,
        "open_entity",
        [
            (
                "POST",
                API + "/entities/groups",
                {
                    "strategy": strategy,
                    "group_key": key,
                    "query": term,
                    "filters": {},
                    "resources": [],
                    "include_details": True,
                    "detail_limit": 50,
                    "detail_offset": 0,
                },
            ),
        ],
    )
    detail = (res[0][1] or {}).get("groups") or []
    pks = []
    for member in (detail[0].get("members") if detail else None) or []:
        pks.extend(member.get("sourceEntityPks") or [])
    pks = list(dict.fromkeys(pks))[:8]
    if pks:
        await step(
            client,
            stats,
            "entity_counts",
            [
                ("POST", API + "/entities/by-public-ids", {"public_ids": pks}),
                *[
                    (
                        "POST",
                        API + "/relations/search",
                        {"filters": {"scope_entity_ids": [pk]}, "limit": 1, "offset": 0},
                    )
                    for pk in pks[:4]
                ],
                (
                    "POST",
                    API + "/relations/search",
                    {"filters": {"scope_entity_ids": pks}, "limit": 1, "offset": 0},
                ),
            ],
        )
    await pause()

    # 4. Relations tab for the group.
    res = await step(
        client,
        stats,
        "relations",
        [
            (
                "POST",
                API + "/entities/groups",
                {
                    "strategy": "auto",
                    "query": term,
                    "filters": {"sources": []},
                    "limit": 1,
                    "member_limit": 1,
                },
            ),
            (
                "POST",
                API + "/relations/scoped-facets",
                {
                    "entityIds": [key],
                    "endpointMode": "any",
                    "mode": "union",
                    "sources": [],
                    "taxonomyLimit": 16,
                    "taxonomyQuery": "",
                },
            ),
            (
                "POST",
                API + "/relations/search",
                {"filters": {"sources": [], "entity_ids": [key]}, "limit": 20, "offset": 0},
            ),
        ],
    )
    facets = res[1][1] or []
    sources = [
        f["facetValue"] for f in facets if isinstance(f, dict) and f.get("facetName") == "source"
    ]
    if not sources:
        return
    await pause()

    # 5. Filter by a source.
    src = rng.choice(sources[:6])
    await step(
        client,
        stats,
        "filter_source",
        [
            (
                "POST",
                API + "/entities/groups",
                {
                    "strategy": "auto",
                    "query": term,
                    "filters": {"sources": [src]},
                    "limit": 1,
                    "member_limit": 1,
                },
            ),
            (
                "POST",
                API + "/relations/scoped-facets",
                {
                    "entityIds": [key],
                    "endpointMode": "any",
                    "mode": "union",
                    "sources": [src],
                    "taxonomyLimit": 16,
                    "taxonomyQuery": "",
                },
            ),
            (
                "POST",
                API + "/relations/search",
                {"filters": {"sources": [src], "entity_ids": [key]}, "limit": 20, "offset": 0},
            ),
        ],
    )
    await pause()

    # 6. Next page of relations.
    await step(
        client,
        stats,
        "next_page",
        [
            (
                "POST",
                API + "/relations/search",
                {"filters": {"sources": [src], "entity_ids": [key]}, "limit": 20, "offset": 20},
            ),
        ],
    )
    stats.journeys += 1
    await pause()


async def user(uid, client, stats, pool, think, deadline, seed):
    rng = random.Random(seed * 1000 + uid)
    await asyncio.sleep(rng.uniform(0, max(think[1], 1)))  # stagger starts
    while time.monotonic() < deadline:
        await journey(client, stats, rng, pool.next(), think)


async def poll_status(client, stats, deadline):
    while time.monotonic() < deadline:
        try:
            r = await client.get(API + "/status")
            stats.status.append(r.json())
        except Exception:
            pass
        await asyncio.sleep(2)


async def run_stage(users, seconds, think, pool, seed):
    stats = Stats()
    limits = httpx.Limits(max_connections=users * 8 + 4, max_keepalive_connections=users * 8 + 4)
    async with httpx.AsyncClient(
        timeout=TIMEOUT, limits=limits, http2=False, headers={"user-agent": "omnipath-loadtest/1"}
    ) as client:
        t0 = time.monotonic()
        deadline = t0 + seconds
        tasks = [
            asyncio.create_task(user(i, client, stats, pool, think, deadline, seed))
            for i in range(users)
        ]
        tasks.append(asyncio.create_task(poll_status(client, stats, deadline)))
        # Hard stop: let in-flight steps finish for up to 60 s past the deadline.
        done, pending = await asyncio.wait(tasks, timeout=seconds + 60)
        for t in pending:
            t.cancel()
        elapsed = time.monotonic() - t0
    out = stats.summary(elapsed)
    out.update(users=users, think=list(think))
    return out, stats


async def main():
    global BASE, API, TIMEOUT
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--base", default=BASE, help="explorer origin")
    ap.add_argument("--timeout", type=float, default=TIMEOUT, help="request timeout in seconds")
    ap.add_argument("--users", default="1,2,4,8")
    ap.add_argument("--seconds", type=int, default=90)
    ap.add_argument("--think", default="3,10")
    ap.add_argument("--reuse", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--abort-error-rate", type=float, default=0.10)
    ap.add_argument("--abort-step-p95", type=float, default=30.0)
    args = ap.parse_args()
    BASE = args.base.rstrip("/")
    API = BASE + "/app-api"
    TIMEOUT = args.timeout
    think = tuple(float(x) for x in args.think.split(","))
    pool = TermPool(TERMS, args.seed, args.reuse)
    results = []
    for n in [int(x) for x in args.users.split(",")]:
        summary, stats = await run_stage(n, args.seconds, think, pool, args.seed)
        results.append(summary)
        print(
            json.dumps(
                {
                    k: summary[k]
                    for k in (
                        "users",
                        "req_per_s",
                        "steps_per_min",
                        "journeys",
                        "error_rate",
                        "codes",
                        "step_p50",
                        "step_p95",
                        "step_p99",
                        "failed",
                        "server_running_mean",
                        "server_waiting_mean",
                        "server_waiting_max",
                        "server_cpu_max",
                        "server_mem_max",
                    )
                }
            ),
            flush=True,
        )
        with open(args.out, "w") as fh:
            json.dump({"args": vars(args), "stages": results}, fh, indent=1)
        if (
            summary["error_rate"] > args.abort_error_rate
            or summary["step_p95"] > args.abort_step_p95
        ):
            print(
                f"abort: error_rate={summary['error_rate']} step_p95={summary['step_p95']}",
                flush=True,
            )
            break
        await asyncio.sleep(15)  # let queues drain between stages


if __name__ == "__main__":
    asyncio.run(main())
