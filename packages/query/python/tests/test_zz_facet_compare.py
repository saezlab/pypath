"""TEMPORARY: the old Python facets against relation_facets (removed after)."""
import json

from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.release_database import serve_newest
from test_client import service  # noqa: F401


def norm(x):
    return json.loads(json.dumps(x, default=str))


def test_old_and_new_facets(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    cases = [
        {}, {"entityIds": ["TP53"]}, {"entity_ids": ["TP53"], "gene_mode": "source"},
        {"entityIds": ["TP53"], "predicates": ["biolink:interacts_with"]},
        {"entityIds": ["TP53"], "sources": ["alpha"]}, {"sources": ["beta"]},
        {"reference_entity_keys": ["entrez:7157"]}, {"taxonomyIds": ["9606"]},
        {"interaction_types": ["protein"], "partner_types": ["small_molecule"]},
        {"entity_types": ["protein"]}, {"signs": [1]}, {"sign": 1}, {"is_directed": True},
        {"protein_entity_keys": ["x"]}, {"entityIds": ["TP53"], "scope_entity_ids": ["caffeine"]},
        {"taxonomyLimit": 1}, {"taxonomyLimit": 1, "taxonomyQuery": "96"},
        {"taxonomyLimit": 2, "taxonomyIds": ["1234"]}, {"filters": {"entity_ids": ["TP53"]}},
        {"relation_categories": ["interaction"], "min_evidence_count": 2},
    ]
    for payload in cases:
        with api.query_scope():
            old = api._old_get_scoped_relation_facets(dict(payload), None)
            new = api.get_scoped_relation_facets(dict(payload), None)
        assert norm(new) == norm(old), payload
    with api.query_scope():
        old = api._old_get_relation_filter_options(None)
        new = api.get_relation_filter_options(None)
    assert norm(new) == norm(old)
    print(len(cases), "cases")


def test_old_and_new_entity_facets(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    with api.query_scope():
        keys = [r[0] for r in api._db.execute(f"SELECT DISTINCT entity_key FROM {api._served().alias}.main.entity").fetchall()]
    cases = [{}, {"query": "TP53"}, {"query": "t"}, {"query": "caf"}, {"query": "zzz"}, {"query": "uniprot:P04637"},
             {"query": keys[0]}, {"query": "TP53", "sources": ["alpha"]}, {"query": "", "entity_types": ["protein"]},
             {"query": "t", "taxonomy_ids": ["9606"]}, {"query": "", "entity_ids": ["TP53"]}, {"facetLimit": 1},
             {"query": "t", "resources": ["beta"]}, {"query": "53"}]
    for payload in cases:
        with api.query_scope():
            old = api._old_get_scoped_entity_facets(dict(payload), None)
            new = api.get_scoped_entity_facets(dict(payload), None)
        assert norm(new) == norm(old), payload


def test_old_and_new_counts(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    for filters in [{}, {"entity_ids": ["TP53"]}, {"sources": ["alpha"]}, {"predicates": ["biolink:interacts_with"]},
                    {"entity_ids": ["TP53"], "scope_entity_ids": ["caffeine"]}, {"signs": [1]}]:
        with api.query_scope():
            old = api._old_get_facets(filters); new = api.get_facets(filters)
        old.pop("elapsed_ms"); new.pop("elapsed_ms")
        assert norm(new) == norm(old), filters
