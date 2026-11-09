"""Paging through an endpoint selection sorts it once and keeps its order."""

import pytest

from served import republish, served_queries
from omnipath_api.queries.services import relations
from omnipath_core.fixtures import rewrite_resource, write_resource

HUB = "0" * 64


def entity(key, identifier, label):
    return dict(
        entity_key=key,
        namespace="entrez",
        identifier=identifier,
        entity_type="protein",
        label=label,
        taxon="9606",
        reference_entity_key=f"entrez:{identifier}",
        gene_reference_keys=[f"entrez:{identifier}"],
    )


def partners(source):
    return [
        entity(f"{source}{i:02d}".ljust(64, "x"), f"{source}{i}", f"P{(i * 7) % 13:02d}")
        for i in range(12)
    ]


def relation(partner, predicate, source):
    return dict(
        relation_key=f"{partner['entity_key'][:8]}-{predicate}",
        subject_entity_key=partner["entity_key"],
        object_entity_key=HUB,
        subject_reference_entity_key=partner["reference_entity_key"],
        object_reference_entity_key="entrez:1",
        subject_type="protein",
        object_type="protein",
        subject_label=partner["label"],
        object_label="HUB",
        predicate=predicate,
        category="interaction",
        taxon="9606",
        sources=[source],
        evidence=[dict(source=source, row_id="1")],
    )


def write(root, source, predicates=("biolink:affects", "biolink:interacts_with")):
    others = partners(source)
    write_resource(
        root / "resources" / source / "1",
        [entity(HUB, "1", "HUB"), *others],
        [relation(p, predicate, source) for p in others for predicate in predicates],
    )


@pytest.fixture
def engine(tmp_path):
    write(tmp_path, "a")
    write(tmp_path, "b", predicates=("biolink:affects",))
    queries = served_queries(tmp_path)
    yield queries
    queries.close()


def page_through(engine, filters, limit=7):
    rows, offset, totals = [], 0, set()
    while True:
        page = engine.search_relations(filters, limit=limit, offset=offset)
        rows += [(r["resource"], r["relation_id"], r["relation_key"]) for r in page["rows"]]
        totals.add(page["total"])
        offset += limit
        if offset >= page["total"]:
            return rows, totals


def test_pages_match_the_uncached_order(engine, monkeypatch):
    filters = {"entity_ids": ["entrez:1"], "relation_categories": ["interaction"]}
    cached, totals = page_through(engine, filters)
    assert totals == {36} and len(cached) == len(set(cached)) == 36
    monkeypatch.setattr(relations, "ORDER_CACHE_LIMIT", 0)
    engine._order_cache.clear()
    uncached, totals = page_through(engine, filters)
    assert totals == {36}
    assert cached == uncached
    assert engine.search_relations(filters, limit=7, offset=100)["rows"] == []


def test_later_pages_skip_the_candidate_search(engine, monkeypatch):
    calls = []
    # A gene's relations are found in the entity_relation table of the release database.
    original = engine._entity_relation_sql
    monkeypatch.setattr(
        engine, "_entity_relation_sql", lambda *a, **k: calls.append(1) or original(*a, **k)
    )
    filters = {"entity_ids": ["entrez:1"]}
    page_through(engine, filters)
    assert len(calls) == 1
    # Other filters are another selection.
    engine.search_relations({"entity_ids": ["entrez:1"], "predicates": ["affects"]}, limit=5)
    assert len(calls) == 2


def test_selections_without_endpoint_filters_are_not_kept(engine):
    first = engine.search_relations({"relation_categories": ["interaction"]}, limit=5)
    assert first["total"] == 36
    assert [value for value, _ in engine._order_cache.entries.values()] == [None]


def test_reloaded_resources_reset_the_order(engine, tmp_path):
    filters = {"entity_ids": ["entrez:1"]}
    assert engine.search_relations(filters, limit=5)["total"] == 36
    path = tmp_path / "resources" / "b" / "1"
    others = partners("b")
    rewrite_resource(path, relations=[relation(others[0], "biolink:affects", "b")])
    republish(engine.data_root, engine)
    assert engine.search_relations(filters, limit=5)["total"] == 25
