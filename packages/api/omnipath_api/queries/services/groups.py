"""Grouped entity views: chemicals by structure, gene products by their gene.

A query-time view (the step file entity_groups); canonical identities and relation
endpoints stay intact.
"""

from __future__ import annotations

from functools import lru_cache
import json
import time

from omnipath_api.queries.filters import normalize_filters


@lru_cache(maxsize=1)
def chemical_types():
    from omnipath_core.biolink import entity_type, schema

    return tuple(
        entity_type(name) for name in schema().class_descendants("chemical entity", reflexive=True)
    )


def search_groups(queries, *, strategy="chemical_connectivity", **kwargs):
    """The grouped page (the step query caches answers), with the cards' details paged."""
    from omnipath_api.routers.paging import page_details

    start = time.perf_counter()
    detail_limit = kwargs.pop("detail_limit", 20)
    detail_offset = kwargs.pop("detail_offset", 0)
    detail_field = kwargs.pop("detail_field", None)
    result = _search_groups(queries, strategy, **kwargs)
    groups = [dict(group) for group in result["groups"]]
    if kwargs.get("include_details"):
        for group in groups:
            group["entity"] = page_details(
                group["entity"], detail_limit, detail_offset, detail_field
            )
    return dict(result, groups=groups, elapsed_ms=round((time.perf_counter() - start) * 1000, 2))


def _search_groups(
    queries,
    strategy,
    *,
    query="",
    filters=None,
    resources=None,
    limit=20,
    cursor=None,
    group_key=None,
    member_cursor=None,
    member_limit=5,
    include_details=False,
    include_member_keys=False,
    **kwargs,
):
    """The grouped page: the query entity_groups, which the clients run too."""
    filters = normalize_filters(filters)
    scope = resources or filters["sources"] or None
    tokens = list(filters["entity_ids"]) + list(filters["entity_pks"])
    params = {
        "query": query or "",
        "strategy": strategy,
        "sources": list(filters["sources"]) or None,
        "entity_types": list(filters["entity_types"]) or None,
        "taxonomy_ids": [str(t) for t in filters["taxonomy_ids"]] or None,
        "reference_entity_keys": list(filters["reference_entity_keys"]) or None,
        "entity_ids": tokens or None,
        "limit": int(limit),
        "cursor": cursor,
        "group_key": group_key,
        "member_cursor": member_cursor,
        "member_limit": int(member_limit),
        "include_details": bool(include_details),
        "include_member_keys": bool(include_member_keys),
        "group_filters": filters,
        "group_resources": list(scope) if scope else None,
    }
    return json.loads(queries.run_query("entity_groups", params, resources))
