#!/usr/bin/env python3
"""Time the explorer's API requests on a data root: first and repeated call.

Caches are cleared before each call, so both times include the queries; the second
reuses DuckDB's Parquet metadata. Run it with the serving memory limit, for example:

    systemd-run --user --scope -p MemoryMax=3G uv run --frozen python scripts/time_api.py /data
"""

from __future__ import annotations

import argparse
import time

from omnipath_api.queries.parquet import ParquetQueries


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_root")
    parser.add_argument("--only", default="", help="time only requests whose name contains this")
    args = parser.parse_args()
    engine = ParquetQueries(args.data_root)

    def timed(name, call):
        if args.only not in name:
            return
        seconds = []
        for _ in range(2):
            for cache in (engine._detail_cache, engine._relationship_cache, engine._facet_cache):
                cache.clear()
            started = time.monotonic()
            with engine.query_scope():
                call()
            seconds.append(time.monotonic() - started)
        print(f"{seconds[0]:7.3f}s {seconds[1]:7.3f}s  {name}", flush=True)

    def first_key(query, kind=None):
        filters = {"entity_types": [kind]} if kind else None
        with engine.query_scope():
            found = engine.search_entities_api(query, filters=filters, limit=1)["entities"]
        return found[0]["entityPk"] if found else None

    tp53, glucose = first_key("TP53", "protein"), first_key("glucose")
    timed("search TP53", lambda: engine.search_entities_api("TP53"))
    timed("search 'a'", lambda: engine.search_entities_api("a"))
    timed("search entrez:7157", lambda: engine.search_entities_api("entrez:7157"))
    timed("search substring", lambda: engine.search_entities_api("phospholamban"))
    timed("details TP53", lambda: engine.get_entity_details(tp53))
    timed("details glucose", lambda: engine.get_entity_details(glucose))
    timed("evidence page glucose", lambda: engine.get_entity_evidence(glucose))
    timed("relations TP53", lambda: engine.search_relations_api({"entity_pks": [tp53]}))
    timed("relation facets TP53", lambda: engine.get_scoped_relation_facets({"entityIds": [tp53]}))
    timed("entity facets glucose", lambda: engine.get_scoped_entity_facets({"query": "glucose"}))
    timed(
        "gene group PLN",
        lambda: engine.search_entity_groups(
            strategy="auto", group_key="gene:entrez:5350", include_details=True
        ),
    )
    timed("molecular context PLN", lambda: engine.get_molecular_context("gene:entrez:5350"))
    timed("molecular context TP53", lambda: engine.get_molecular_context("gene:entrez:7157"))
    timed(
        "groups glucose",
        lambda: engine.search_entity_groups(strategy="auto", query="glucose", member_limit=1),
    )
    timed(
        "groups default page", lambda: engine.search_entity_groups(strategy="auto", member_limit=1)
    )
    timed("relations browse", lambda: engine.search_relations_api({}))


if __name__ == "__main__":
    main()
