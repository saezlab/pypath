"""Paging through an endpoint selection sorts it once and keeps its order."""

import pytest

from served import republish, served_queries
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
        rows += [r["relation_key"] for r in page["rows"]]
        totals.add(page["total"])
        offset += limit
        if offset >= page["total"]:
            return rows, totals


def test_pages_match_the_uncached_order(engine):
    filters = {"entity_ids": ["entrez:1"], "relation_categories": ["interaction"]}
    cached, totals = page_through(engine, filters)
    assert totals == {36} and len(cached) == len(set(cached)) == 36
    engine._query_cache.clear()
    engine._step_cache.clear()
    uncached, totals = page_through(engine, filters)
    assert totals == {36}
    assert cached == uncached
    assert engine.search_relations(filters, limit=7, offset=100)["rows"] == []


def test_reloaded_resources_reset_the_order(engine, tmp_path):
    filters = {"entity_ids": ["entrez:1"]}
    assert engine.search_relations(filters, limit=5)["total"] == 36
    path = tmp_path / "resources" / "b" / "1"
    others = partners("b")
    rewrite_resource(path, relations=[relation(others[0], "biolink:affects", "b")])
    republish(engine.data_root, engine)
    assert engine.search_relations(filters, limit=5)["total"] == 25
