"""A resource selection reads the release as if only those resources were in it.

Three resources state overlapping relations (one relation key is one statement: ends,
predicate and qualifiers, so its sign and direction are the same whoever states it):
the relation fields that depend on who states it (sources, evidence count, annotations)
come from the selected resources only, on a page, in an export, as a single relation
and in the entity's relationships, and the evidence count filter reads them too.
"""

import pytest

from omnipath_core.fixtures import fixture_key as key, write_resource
from served import served_queries

P, Q, C = key("P"), key("Q"), key("C")


def entity(entity_key, identifier, label, entity_type="protein"):
    return dict(
        entity_key=entity_key,
        namespace="test",
        identifier=identifier,
        entity_type=entity_type,
        label=label,
        taxon="9606",
    )


def relation(rkey, resource, evidence, predicate="biolink:affects", direction=None, pmid=None,
             subject=P, object=Q, sign=0, directed=True, category="interaction", subject_type="protein"):
    annotations = [dict(term="publications", value=f"PMID:{pmid}", scope="relation", source=resource)]
    if direction:
        annotations.append(
            dict(term="object_direction_qualifier", value=direction, scope="relation", source=resource)
        )
    return dict(
        relation_key=rkey,
        subject_entity_key=subject,
        object_entity_key=object,
        subject_type=subject_type,
        object_type="protein",
        predicate=predicate,
        category=category,
        taxon="9606",
        sign=sign,
        is_directed=directed,
        sources=[resource],
        annotations=annotations if pmid else annotations[1:],
        evidence=[dict(source=resource, row_id=f"{resource}-{rkey}-{i}") for i in range(evidence)],
    )


@pytest.fixture
def engine(tmp_path):
    entities = [entity(P, "P", "P"), entity(Q, "Q", "Q"), entity(C, "C", "C", "macromolecular_complex")]
    # activation: all three state it (2, 1 and 3 records); inhibition: b and c (1 and 2);
    # binding (undirected): a only; the complex has P as a member: a and b.
    def member(resource):
        return relation("member", resource, 1, predicate="has_part", category="membership",
                        subject=C, subject_type="macromolecular_complex")

    write_resource(
        tmp_path / "resources" / "a" / "1",
        entities,
        [
            relation("activates", "a", 2, direction="increased", sign=1, pmid=1),
            relation("binds", "a", 1, predicate="biolink:interacts_with", directed=False, pmid=5),
            member("a"),
        ],
    )
    write_resource(
        tmp_path / "resources" / "b" / "1",
        entities,
        [
            relation("activates", "b", 1, direction="increased", sign=1, pmid=2),
            relation("inhibits", "b", 1, direction="decreased", sign=-1, pmid=3),
            member("b"),
        ],
    )
    write_resource(
        tmp_path / "resources" / "c" / "1",
        entities,
        [
            relation("activates", "c", 3, direction="increased", sign=1, pmid=4),
            relation("inhibits", "c", 2, direction="decreased", sign=-1, pmid=6),
        ],
    )
    engine = served_queries(tmp_path)
    yield engine
    engine.close()


def page(engine, resources=None, **filters):
    result = engine.search_relations_api({"entity_pks": [P], **filters}, resources=resources, limit=50)
    return {row["relation"]["relationPk"]: row["relation"] for row in result["rows"]}, result["total"]


def facts(relation_json):
    return (
        relation_json["sources"],
        relation_json["evidenceCount"],
        relation_json["sign"],
        relation_json["isDirected"],
    )


def test_page_without_a_selection_is_the_merged_release(engine):
    rows, total = page(engine)
    assert total == 3
    assert facts(rows["activates"]) == (["a", "b", "c"], 6, 1, True)
    assert facts(rows["inhibits"]) == (["b", "c"], 3, -1, True)
    assert facts(rows["binds"]) == (["a"], 1, 0, False)


def test_page_reads_the_selected_resources_only(engine):
    rows, total = page(engine, ["a", "b"])
    assert total == 3
    assert facts(rows["activates"]) == (["a", "b"], 3, 1, True)
    assert facts(rows["inhibits"]) == (["b"], 1, -1, True)
    assert facts(rows["binds"]) == (["a"], 1, 0, False)

    rows, total = page(engine, ["c"])
    assert total == 2
    assert facts(rows["activates"]) == (["c"], 3, 1, True)
    assert facts(rows["inhibits"]) == (["c"], 2, -1, True)

    rows, total = page(engine, ["b", "c"])
    assert total == 2
    assert facts(rows["activates"]) == (["b", "c"], 4, 1, True)
    assert facts(rows["inhibits"]) == (["b", "c"], 3, -1, True)


def test_single_relation_reads_the_selected_resources_only(engine):
    assert facts(engine.get_relation_by_pk("activates")["relation"]) == (["a", "b", "c"], 6, 1, True)
    assert facts(engine.get_relation_by_pk("activates", ["b"])["relation"]) == (["b"], 1, 1, True)
    assert facts(engine.get_relation_by_pk("activates", ["a", "c"])["relation"]) == (["a", "c"], 5, 1, True)
    assert engine.get_relation_by_pk("binds", ["b", "c"]) is None


def annotation_values(row, term):
    return sorted(a["value"] for a in row["annotations"] if a["term"] == term)


def test_export_and_details_read_the_selected_resources_only(engine):
    detailed = engine.search_relations({"entity_pks": [P]}, resources=["a", "b"], limit=50, include_details=True)
    rows = {row["relation_key"]: row for row in detailed["rows"]}
    assert rows["activates"]["sources"] == ["a", "b"]
    assert rows["activates"]["evidence_count"] == 3
    assert len(rows["activates"]["evidence"]) == 3
    assert annotation_values(rows["activates"], "publications") == ["PMID:1", "PMID:2"]
    assert annotation_values(engine.search_relations({"entity_pks": [P]}, resources=["c"], limit=50,
                                                     include_details=True)["rows"][0], "publications") != []
    everything = engine.search_relations({"entity_pks": [P]}, limit=50, include_details=True)
    rows = {row["relation_key"]: row for row in everything["rows"]}
    assert rows["activates"]["sources"] == ["a", "b", "c"]
    assert rows["activates"]["evidence_count"] == 6
    assert annotation_values(rows["activates"], "publications") == ["PMID:1", "PMID:2", "PMID:4"]

    table = engine.export_slice({"entity_pks": [P]}, resources=["b", "c"])
    import io

    import pyarrow.parquet as pq

    exported = {r["relation_key"]: r for r in pq.read_table(io.BytesIO(table[0])).to_pylist()}
    assert exported["activates"]["sources"] == ["b", "c"]
    assert exported["activates"]["evidence_count"] == 4
    assert exported["inhibits"]["evidence_count"] == 3


def test_entity_relationships_read_the_selected_resources_only(engine):
    def members(resources):
        found = engine.get_entity_relationships(Q, resources)
        return {r["relation"]["relationPk"]: facts(r["relation"]) for r in found["relationships"]}

    assert members(None) == {"member": (["a", "b"], 2, 0, True)}
    assert members(["b"]) == {"member": (["b"], 1, 0, True)}
    assert members(["c"]) == {}


def test_the_evidence_count_filter_reads_the_selected_resources_only(engine):
    # activates has 6 records in all, 3 in a and b, 3 in c alone; inhibits 3 and 2.
    rows, total = page(engine, ["a", "b"], min_evidence_count=3)
    assert set(rows) == {"activates"} and total == 1
    rows, total = page(engine, ["b"], min_evidence_count=2)
    assert rows == {} and total == 0
    rows, total = page(engine, ["a", "b", "c"], min_evidence_count=4)
    assert set(rows) == {"activates"} and total == 1
    rows, total = page(engine, ["b", "c"], min_evidence_count=4)
    assert set(rows) == {"activates"} and total == 1
    rows, total = page(engine, ["a", "c"], min_evidence_count=3)
    assert set(rows) == {"activates"}
    rows, total = page(engine, ["a", "c"], min_evidence_count=2)
    assert set(rows) == {"activates", "inhibits"} and total == 2


def test_source_options_list_the_selected_resources_only(engine):
    def sources(resources):
        options = engine.get_scoped_relation_facets({"entity_pks": [P]}, resources=resources)
        return sorted(o["facetValue"] for o in options if o["facetName"] == "source")

    assert sources(None) == ["a", "b", "c"]
    assert sources(["a", "b"]) == ["a", "b"]
    assert sources(["c"]) == ["c"]


def test_the_evidence_count_filter_on_every_selection_path(engine):
    # by entity (a range of the endpoint table), by scope entities (candidates) and
    # without endpoint filters (every relation): the same relations
    for filters in ({"entity_pks": [P]}, {"scope_entity_ids": [P]}, {}):
        found = engine.search_relations_api(
            dict(filters, min_evidence_count=3), resources=["a", "b"], limit=50
        )
        keys = {row["relation"]["relationPk"] for row in found["rows"]}
        assert keys == {"activates"}, filters
        assert found["total"] == 1, filters
        counts = engine.get_facets(dict(filters, min_evidence_count=3), resources=["a", "b"])
        assert counts["signs"] == {"1": 1}, filters
        counts = engine.get_facets(dict(filters, min_evidence_count=3), resources=["b", "c"])
        assert counts["signs"] == {"-1": 1, "1": 1}, filters
