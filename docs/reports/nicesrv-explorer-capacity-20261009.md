# Explorer capacity on nicesrv — 9 October 2026

How many people can use https://omnipath-metabo-dev.schaul.click/explore at the
same time, and which API limits to run with. Release `2026.10.9`, image
`omnipath-api:serving-e12b0fb`. Raw stage results are in
`nicesrv-explorer-capacity-20261009.json`; the load generator is
`scripts/explorer_load_test.py`.

## Decision

The API container gets 7 of nicesrv's 8 CPUs and 13 GB of its 16 GB, leaving
1 CPU and about 2.6 GB for the other apps on the host (they use about 1.3 GB).
Queries run in three slots of two DuckDB threads and 3,500 MB each:

```
API_CPUS=7
API_MEMORY_LIMIT=13g
API_QUERY_CONCURRENCY=3
API_QUERY_THREADS=2
API_QUERY_MEMORY=3500MB
API_QUERY_WAIT_SECONDS=30
```

The previous settings (4 CPUs, 6 GB, two slots of two threads and 2,500 MB) are
kept on nicesrv as `serving.env.bak-before-capacity-20261009`.

With these limits the explorer serves about 105 interactions a minute with the
slowest 5% under about 5 s, and about 125 a minute at most. With 15–30 s between
a person's clicks, that is roughly **30–55 people using the explorer at the same
time**, and 35–65 at the most. The old limits supported roughly 20–35.

## Method

The load generator replays a session recorded in the explorer: open the page,
search, open one of the top results (group details, members and their relation
counts), switch to the Relations tab (filter counts and the first 20
relations), filter by one of the result's sources and load the next page. That
is about 18 API requests in six steps. Requests of one step are sent together, as
the browser sends them; a step's latency is the time until its slowest request
returns.

Each simulated person repeats the session with 3–10 s between steps, searching a
new term each time, so the API's result caches mostly miss. The test therefore
measures close to the worst case for caching. Stages of 2 minutes add people
until more than 10% of requests fail or the slowest 5% of steps take more than
15 s (30 s for the first run). Every run used the same search terms in the same
order, and each configuration started from a fresh API container.

The client ran from a laptop, so times include the network round trip, which is
small next to the query times. Exports and ontology browsing were not part of the
session.

## Results

Interactions are steps; the p95 column is the slowest 5% of steps.

**Old limits** — 4 CPUs, 2 slots × 2 threads:

| People | Steps/min | Median step | p95 step | Requests waiting (mean / max) |
| ---: | ---: | ---: | ---: | --- |
| 4 | 33 | 0.8 s | 3.5 s | 0.4 / 4 |
| 8 | 66 | 0.7 s | 2.2 s | 1 / 9 |
| 12 | 81 | 1.5 s | 6.5 s | 6 / 26 |
| 16 | 89 | 2.6 s | 5.8 s | 10 / 29 |
| 24 | 72 | 7.8 s | 23.1 s | 25 / 38 |

The API container reached about 360% of its 400% CPU cap; the host's load
average stayed at or below 3.3.

**Chosen limits (A)** — 7 CPUs, 3 slots × 2 threads:

| People | Steps/min | Median step | p95 step | Requests waiting (mean / max) |
| ---: | ---: | ---: | ---: | --- |
| 8 | 54 | 1.0 s | 8.0 s | 2.6 / 14 |
| 12 | 92 | 1.1 s | 3.6 s | 1.8 / 12 |
| 16 | 111 | 1.4 s | 4.5 s | 5.7 / 23 |
| 24 | 75 | 7.1 s | 14.8 s | 18 / 36 |
| 32 | 126 | 8.0 s | 18.7 s | 23 / 37 |

At 24 people, 53 of 580 requests did not return within the client's 60 s
timeout. Another session's single-threaded job ran on the host from 10:45 to
10:47 CEST, overlapping the start of that stage. The 8-person stage started right
after the container restart, with empty caches. API CPU peaked at about 540% of
700%; host load stayed at or below 3.1.

**Alternative (B)** — 7 CPUs, 6 slots × 1 thread:

| People | Steps/min | Median step | p95 step | Requests waiting (mean / max) |
| ---: | ---: | ---: | ---: | --- |
| 8 | 48 | 1.6 s | 8.8 s | 1.1 / 9 |
| 12 | 92 | 1.1 s | 3.4 s | 0.7 / 12 |
| 16 | 95 | 2.2 s | 6.2 s | 4.8 / 23 |
| 24 | 118 | 4.2 s | 12.3 s | 15 / 33 |
| 32 | 117 | 5.8 s | 25.2 s | 20 / 34 |

At 32 people, 8 requests returned HTTP 503 after waiting 30 s for a query slot.
API CPU peaked at about 580% of 700%; host load stayed at or below 4.3. API memory
stayed under 3.7 GB in every run.

A was chosen over B: it answered faster at moderate load (p95 4.5 s against
6.2 s at 16 people), reached a slightly higher ceiling (126 against 118 steps a
minute) and returned no 503s. One person clicking quickly keeps about 0.7 query
slots busy, with a median step of 0.7 s and a p95 of 2 s.

## From steps to people

Little's law gives the number of people a throughput supports: people = steps
per second × (time per step + time between steps). At 105 steps a minute
(1.75 a second) and about 2 s per step, 15 s between clicks gives about 30 people
and 30 s gives about 55. These are people in the explorer at the same moment;
daily visitors can be many times more.

## Findings

- More CPU helped less than proportionally: 1.75 times the CPUs gave about
  1.6 times the throughput at acceptable latency and 1.4 times the ceiling.
  Neither 7-CPU layout used the full CPU cap.
- The slowest requests under load are the relation filter counts
  (`/relations/scoped-facets`), relation search and entity group details, at
  10–25 s p95 near the ceiling. Making these cheaper would raise capacity as much
  as more hardware.
- Overload shows up as waiting, not errors: requests queue for up to 30 s before
  a 503, so pages get slow well before anything fails.
- Queued queries keep running after the client has gone. When the first test was
  stopped, 33 queries were still waiting and took about 20 s to drain. Cancelling
  a query when its client disconnects would let an overloaded server recover
  faster.
- Memory is not a constraint: the API used under 4 GB of its 13 GB.
