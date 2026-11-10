"""TEMPORARY: the old Python groups against entity_groups (removed after)."""
import json

from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.queries.services import groups
from omnipath_api.release_database import serve_newest
from test_client import service  # noqa: F401


def norm(x):
    return json.loads(json.dumps(x, default=str))


def test_old_and_new_groups(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    cases = []
    for strategy in ("auto", "chemical_connectivity", "gene_reference"):
        for kwargs in ({}, {"query": "TP53"}, {"query": "t"}, {"query": "caf"}, {"limit": 1},
                       {"filters": {"sources": ["alpha"]}}, {"member_limit": 1}, {"include_member_keys": True},
                       {"include_details": True}, {"group_key": "gene:entrez:7157"}, {"member_cursor": "0"}):
            cases.append((strategy, kwargs))
    n = 0
    for strategy, kwargs in cases:
        with api.query_scope():
            old = groups._old_search_groups(api, strategy, **kwargs)
            new = groups._search_groups(api, strategy, **kwargs)
        assert norm(new) == norm(old), (strategy, kwargs)
        if old.get("nextCursor"):
            with api.query_scope():
                old2 = groups._old_search_groups(api, strategy, **{**kwargs, "cursor": old["nextCursor"]})
                new2 = groups._search_groups(api, strategy, **{**kwargs, "cursor": new["nextCursor"]})
            assert norm(new2) == norm(old2), (strategy, kwargs, "page 2")
        n += 1
    print(n, "cases")
