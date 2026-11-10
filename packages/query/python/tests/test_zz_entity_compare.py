"""TEMPORARY: the old Python entity pages against the step files (removed after)."""
import json

from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.release_database import serve_newest
from test_client import service  # noqa: F401


def norm(x):
    return json.loads(json.dumps(x, default=str))


def test_old_and_new_entity_pages(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    with api.query_scope():
        keys = [r[0] for r in api._db.execute(f"SELECT DISTINCT entity_key FROM {api._served().alias}.main.entity").fetchall()]
    n = 0
    for key in keys + ["missing"]:
        with api.query_scope():
            assert norm(api.get_entity_core(key)) == norm(api._old_get_entity_core(key)), key
            from omnipath_api.queries.entity_details import details
            assert norm(api.get_entity_details(key)) == norm(details(api, key, None)), key
            for limit, offset in ((20, 0), (1, 0), (1, 1), (2, 1), (5, 100)):
                assert norm(api.get_entity_evidence(key, None, limit, offset)) == norm(api._old_get_entity_evidence(key, None, limit, offset)), (key, limit, offset)
                assert norm(api.get_entity_relationships(key, None, limit, offset)) == norm(api._old_get_entity_relationships(key, None, limit, offset)), (key, limit, offset)
            n += 1
    with api.query_scope():
        assert norm(api.get_entity_relationships(keys, None, 10, 0)) == norm(api._old_get_entity_relationships(keys, None, 10, 0))
    print(n, "entities")


def test_old_and_new_selection_and_lookups(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    with api.query_scope():
        keys = [r[0] for r in api._db.execute(f"SELECT DISTINCT entity_key FROM {api._served().alias}.main.entity").fetchall()]
    for payload in [{}, {"entityPks": ["TP53"]}, {"entityPks": ["TP53", "caffeine"]}, {"entityPks": keys},
                    {"entityPks": ["TP53"], "includeAssociatedEntities": False}, {"entityPks": ["nope"]},
                    {"entityPks": ["TP53"], "annotationTermIds": ["GO:1"]}]:
        with api.query_scope():
            assert norm(api.resolve_selection_scope(payload)) == norm(api._old_resolve_selection_scope(payload)), payload
    for tokens in [["TP53"], ["caffeine", "TP53"], keys, ["nope"], ["gene:entrez:7157"]]:
        for slim in (False, True):
            with api.query_scope():
                assert norm(api.get_entities_by_pks(tokens, slim=slim)) == norm(api._old_get_entities_by_pks(tokens, slim=slim)), (tokens, slim)
