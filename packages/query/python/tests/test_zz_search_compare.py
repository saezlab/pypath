"""TEMPORARY: the old Python search against the entity_search step file (removed after)."""
from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.release_database import serve_newest
from test_client import service  # noqa: F401


def strip(result):
    result = dict(result)
    result.pop("elapsed_ms", None)
    if result.get("nextCursor"):
        result["nextCursor"] = {k: v for k, v in result["nextCursor"].items() if k != "queryKey"}
    return result


def test_old_and_new_search(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    with api.query_scope():
        alias = api._served().alias
        labels = [r[0] for r in api._db.execute(f"SELECT DISTINCT label FROM {alias}.main.entity WHERE label IS NOT NULL").fetchall()]
        keys = [r[0] for r in api._db.execute(f"SELECT DISTINCT entity_key FROM {alias}.main.entity").fetchall()]
        ids = [f"{r[0]}:{r[1]}" for r in api._db.execute(f"SELECT DISTINCT ns, id FROM {alias}.main.entity_identifier").fetchall()]
    queries = ["", "t", "tp", "TP53", "caf", "a", "o'b", "zzz", "53", "entrez:7157", "uniprot|P04637"] + labels + keys + ids
    filter_sets = [None, {"sources": ["alpha"]}, {"entity_types": ["protein"]}, {"taxonomy_ids": ["9606"]},
                   {"reference_entity_keys": ["entrez:7157"]}, {"entity_ids": ["TP53"]}, {"entity_types": ["small_molecule"]}]
    n = 0
    for q in queries:
        for filters in filter_sets:
            for limit in (1, 2, 20):
                with api.query_scope():
                    old = api._old_search_entities_api(q, filters, None, limit)
                    new = api.search_entities_api(q, filters, None, limit)
                assert strip(new) == strip(old), (q, filters, limit)
                n += 1
                if old.get("nextCursor"):
                    with api.query_scope():
                        old2 = api._old_search_entities_api(q, filters, None, limit, old["nextCursor"])
                        new2 = api.search_entities_api(q, filters, None, limit, new["nextCursor"])
                    assert strip(new2) == strip(old2), (q, filters, limit, "page 2")
    print(n, "cases")
