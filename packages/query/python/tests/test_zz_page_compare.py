"""TEMPORARY: the old Python relation pages against relation_page (removed after)."""
import json

from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.release_database import serve_newest
from test_client import service  # noqa: F401
from test_zz_export_compare import test_old_and_new_export  # noqa: F401  (its cases)


def strip(r):
    r = json.loads(json.dumps(r, default=str)); r.pop("elapsed_ms", None)
    return r


def test_old_and_new_pages(service):
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    cases = [
        {"entity_ids": ["TP53"]}, {"entity_ids": ["TP53"], "sources": ["alpha"]},
        {"entity_ids": ["TP53"], "gene_mode": "source"}, {"entity_ids": ["TP53"], "gene_mode": "both"},
        {"entity_ids": ["caffeine", "O'Brien"], "gene_mode": "both"}, {"entity_ids": ["gene:entrez:7157"]},
        {"entity_ids": ["uniprot:P04637"]}, {"reference_entity_keys": ["entrez:7157"]},
        {"entity_ids": ["TP53"], "reference_entity_keys": ["entrez:7157"]},
        {"entity_ids": ["TP53"], "scope_entity_ids": ["caffeine"]},
        {"predicates": ["biolink:interacts_with"]},
        {"interaction_types": ["protein"], "partner_types": ["small_molecule"]},
        {"taxonomy_ids": ["9606"], "relation_categories": ["interaction"]}, {"sources": ["beta"]},
        {"signs": [1]}, {"is_directed": True}, {}, {"protein_entity_keys": ["x"]},
        {"entity_ids": ["TP53"], "protein_entity_keys": ["x"]},
        {"entity_ids": ["TP53"], "relation_categories": ["interaction"]},
    ]
    n = 0
    for filters in cases:
        for limit, offset in ((20, 0), (1, 0), (1, 1), (2, 3)):
            with api.query_scope():
                old = api._old_search_relations_api(filters, None, limit, offset)
                new = api.search_relations_api(filters, None, limit, offset)
            assert strip(new) == strip(old), (filters, limit, offset)
            n += 1
    print(n, "pages")
