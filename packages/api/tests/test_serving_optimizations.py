import json
from concurrent.futures import ThreadPoolExecutor
from threading import Event

import pytest
from fastapi.testclient import TestClient

from served import served_queries
from omnipath_api.queries.store.query_cache import QueryCache
from omnipath_api.server import create_app
from omnipath_core.fixtures import write_resource


@pytest.fixture
def engine(tmp_path):
    for source in ["a", "b"]:
        entities = [
            dict(
                entity_key=f"{i:064x}",
                namespace="uniprot",
                identifier=f"P{i}",
                entity_type="protein",
                taxon=str(i),
                label=f"Protein {i:02}",
                identifiers=[dict(ns="alias", id=f"Alias{i}")],
                annotations=[],
            )
            for i in range(30)
        ]
        if source == "b":
            entities[2]["label"] = "Protein 02 alternate"
            entities[2]["has_hierarchy"] = True
            entities[2]["parent_count"] = 1
        rows = [
            dict(
                relation_key=f"r{i}",
                subject_entity_key=f"{i:064x}",
                object_entity_key=f"{0:064x}",
                subject_type="protein",
                object_type="protein",
                subject_label=f"Protein {i:02}",
                object_label="Protein 00",
                predicate="affects",
                category="interaction",
                taxon=str(i),
                sources=[source],
                sign=0,
                evidence_count=1,
                annotations=[
                    dict(term="object_aspect_qualifier", value="activity", scope="relation"),
                    dict(term="object_aspect_qualifier", value="activity", scope="relation"),
                    dict(
                        term="object_direction_qualifier",
                        value="increased" if i % 2 else "decreased",
                        scope="relation",
                    ),
                    dict(term="object_aspect_qualifier", value="expression", scope="subject"),
                    dict(term="description", value="large evidence annotation", scope="relation"),
                ],
            )
            for i in range(30)
        ]
        write_resource(tmp_path / "resources" / source / "1", entities, rows)
    return served_queries(tmp_path)


def test_live_query_policy_preserves_catalog_and_invalidates_cached_queries(engine):
    # Warm caches before changing policy.
    before = engine.get_scoped_entity_facets({})
    assert len(engine._selected_resource_infos()) == 2
    policy = engine.data_root / "query_policy.json"
    policy.write_text(json.dumps({"disabled_resources": ["b"]}))
    assert [info["resource"] for info in engine._selected_resource_infos()] == ["a"]
    for selection in [["b"], ["b/1"], ["a", "b"]]:
        with pytest.raises(ValueError, match="disabled for queries"):
            engine.search_entities_api(resources=selection)
    assert engine.search_entities_api(query="alternate")["total"] == 0
    assert engine.get_scoped_entity_facets({}) != before
    assert [item["slug"] for item in engine.get_stats_sources()] == ["a"]
    catalog = {item["resource_id"]: item for item in engine.list_resource_catalog()}
    assert catalog["b"]["queryable"] is False
    assert engine.get_resource_file_path("b", "entity.parquet").is_file()
    client = TestClient(create_app(queries=engine))
    assert client.post("/entities/search", json={"resources": ["b"]}).status_code == 400
    assert any(
        item["resource_id"] == "b" and not item["queryable"]
        for item in client.get("/resources/catalog").json()["resources"]
    )
    policy.unlink()
    assert engine.search_entities_api(query="alternate")["total"] == 1


def test_paging_search_selected_taxonomy_and_http_gzip(engine):
    all_rows = engine.get_scoped_relation_facets({})
    page = engine.get_scoped_relation_facets({"taxonomyLimit": 10})
    assert len([r for r in page if r["facetName"] == "taxonomy_id"]) == 10
    assert all(r in all_rows for r in page)
    searched = engine.get_scoped_relation_facets({"taxonomyLimit": 10, "taxonomyQuery": "29"})
    assert [r["facetValue"] for r in searched if r["facetName"] == "taxonomy_id"] == ["29"]
    # Selected values are included even when the taxonomy search excludes them.
    selected = engine.get_scoped_relation_facets(
        {"taxonomyLimit": 1, "taxonomyQuery": "absent", "taxonomyIds": ["29"]}
    )
    assert any(r["facetValue"] == "29" for r in selected)
    absent = engine.get_scoped_relation_facets({"taxonomyLimit": 1, "taxonomyIds": ["unknown"]})
    assert any(r["facetValue"] == "unknown" for r in absent)
    client = TestClient(create_app(queries=engine))
    response = client.post(
        "/relations/scoped-facets", json={"taxonomyLimit": 10}, headers={"Accept-Encoding": "gzip"}
    )
    assert response.status_code == 200 and response.headers["content-encoding"] == "gzip"
    assert len([r for r in response.json() if r["facetName"] == "taxonomy_id"]) == 10
    assert client.post("/relations/scoped-facets", json={"taxonomyLimit": 0}).status_code == 422


def _matches(engine, payload, **override):
    """Every entity a search with the payload's query and filters finds."""
    filters = {k: v for k, v in {**payload, **override}.items() if k != "query"}
    found = engine.search_entities_api(query=payload.get("query", ""), filters=filters, limit=1000)
    assert found["nextCursor"] is None
    return found["entities"]


def test_facet_batching_matches_independent_counts(engine):
    # The facets count the resource copies of what the search finds: by type and taxon
    # with every filter, by source without the source filter.
    for payload in [
        {},
        {"query": "02"},
        {"sources": ["a"]},
        {"query": "Protein", "ncbi_tax_id": ["2", "3"]},
    ]:
        rows = engine.get_scoped_entity_facets(payload)
        copies = [
            entity
            for source in payload.get("sources") or ("a", "b")
            for entity in _matches(engine, payload, sources=[source])
        ]
        for facet, field in [("entity_type", "entityType"), ("taxonomy_id", "taxonomyId")]:
            expected = {}
            for entity in copies:
                expected[str(entity[field])] = expected.get(str(entity[field]), 0) + 1
            counted = {r["facetValue"]: r["scopedCount"] for r in rows if r["facetName"] == facet}
            assert counted == expected
        for source in ("a", "b"):
            expected = len(_matches(engine, payload, sources=[source]))
            counted = [r["scopedCount"] for r in rows if r["facetName"] == "source" and r["facetValue"] == source]
            assert counted == ([expected] if expected else [])


def test_shared_cache_bounds_failures_and_mutation_isolation():
    cache = QueryCache(max_entries=1, max_bytes=100)
    entered, finish = Event(), Event()
    calls = []

    def compute():
        calls.append(1)
        entered.set()
        assert finish.wait(5)
        return [{"count": 1}]

    with ThreadPoolExecutor(4) as pool:
        futures = [pool.submit(cache.get, "a", compute) for _ in range(4)]
        assert entered.wait(5)
        finish.set()
        results = [f.result() for f in futures]
    assert len(calls) == 1
    results[0][0]["count"] = 99
    assert cache.get("a", lambda: None) == [{"count": 1}]
    cache.get("b", lambda: ["other"])
    assert len(cache.entries) == 1 and "a" not in cache.entries

    def fail():
        raise ValueError("failure")

    with pytest.raises(ValueError):
        cache.get("bad", fail)
    assert cache.get("bad", lambda: "recovered") == "recovered"
    cache.get("large", lambda: "x" * 101)
    assert "large" not in cache.entries


def test_duplicate_label_choice_is_independent_of_scan_order(tmp_path):
    labels = ["aspirin-triggered resolvin D4", "AT-RvD4"]
    chosen = []
    for root, order in ((tmp_path / "one", labels), (tmp_path / "two", labels[::-1])):
        for source, label in zip("ab", order):
            entity = dict(entity_key="a", entity_type="small_molecule", label=label,
                          identifier="X", namespace="inchikey")
            write_resource(root / "resources" / source / "1", [entity])
        queries = served_queries(root)
        with queries.release_scope():
            entity = queries.get_entity_core("a")["entity"]
        chosen.append((entity["label"], entity["displayName"]))
    assert chosen[0] == chosen[1]


def test_scoped_pages_totals_details_and_shared_scope_cache(engine):
    from collections import Counter

    keys = [f"{i:064x}" for i in range(1, 15)]
    filters = {"scope_entity_ids": keys}
    page = engine.search_relations(filters, limit=6, offset=2)
    assert page["total"] == 28
    assert sum(Counter(row["relation_key"] for row in page["rows"]).values()) == 6
    assert engine.search_relations(filters, limit=6, offset=100)["total"] == 28
    assert (
        engine.search_relations({"scope_entity_ids": keys, "scope_endpoint_mode": "both"})["total"]
        == 0
    )
    detailed = engine.search_relations(filters, limit=2, include_details=True)
    assert any(a["term"] == "description" for a in detailed["rows"][0]["annotations"])
    scope = engine.resolve_selection_scope({"entityPks": [f"{0:064x}"]})
    assert set(scope["entityPks"]) == {f"{i:064x}" for i in range(30)}
    # Cached responses are copied, so caller mutation cannot corrupt the scope.
    scope["entityPks"].clear()
    assert len(engine.resolve_selection_scope({"entityPks": [f"{0:064x}"]})["entityPks"]) == 30
    page = engine.search_relations_api(filters, limit=6)
    page["rows"].clear()
    assert len(engine.search_relations_api(filters, limit=6)["rows"]) == 6
