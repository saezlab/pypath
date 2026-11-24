"""The relation export: filters apply before the limit, and evidence and annotations
come in their stored order."""

from __future__ import annotations

import io

import pyarrow.parquet as pq
import pytest

from served import served_queries
from omnipath_core.fixtures import write_resource


def _entity(key, label, identifier):
    return dict(
        entity_key=key, entity_type="protein", namespace="uniprot", identifier=identifier,
        label=label, annotations=[],
        identifiers=[dict(ns="uniprot", id=identifier, source="test", is_canonical=True)],
    )


@pytest.fixture
def engine(tmp_path):
    folder = tmp_path / "resources" / "test" / "1"
    folder.mkdir(parents=True)
    hub = _entity("a" * 64, "HUB", "P00001")
    partners = [_entity(f"{i:x}" * 64, f"PARTNER{i}", f"P1000{i}") for i in range(1, 7)]
    relations = []
    for i, partner in enumerate(partners):
        relations.append(
            dict(
                relation_key=f"{i + 1:x}" * 63 + "f",
                subject_entity_key=hub["entity_key"],
                object_entity_key=partner["entity_key"],
                subject_type="protein", object_type="protein",
                # every second relation has another predicate
                predicate="affects" if i % 2 == 0 else "interacts_with",
                category="interaction", interaction_class="protein_protein",
                taxon="9606", is_directed=True, sign=1, sources=["test"],
                annotations=[
                    dict(term=f"term{j}", value=f"value{i}{j}", source="test", dataset="d")
                    for j in (3, 1, 2)
                ],
                evidence=[
                    dict(source=f"source{j}", dataset="d", row_id=f"{i}-{j}", annotations=[])
                    for j in (3, 1, 2)
                ],
            )
        )
    write_resource(folder, [hub, *partners], relations)
    result = served_queries(data_root=tmp_path)
    yield result
    result.close()


def _export(engine, filters, limit=100_000):
    data, _ = engine.export_slice(filters, limit=limit)
    return pq.read_table(io.BytesIO(data))


def test_export_applies_filters_before_the_limit(engine):
    everything = _export(engine, {"entity_ids": ["HUB"]})
    assert everything.num_rows == 6
    affects = _export(engine, {"entity_ids": ["HUB"], "predicates": ["affects"]})
    assert affects.num_rows == 3
    assert set(affects["predicate"].to_pylist()) == {"affects"}
    # the limit takes the first relations that match, not the first candidates
    limited = _export(engine, {"entity_ids": ["HUB"], "predicates": ["interacts_with"]}, limit=2)
    assert limited.num_rows == 2
    assert set(limited["predicate"].to_pylist()) == {"interacts_with"}
    assert limited["relation_key"].to_pylist() == [
        k for k, p in zip(everything["relation_key"].to_pylist(), everything["predicate"].to_pylist())
        if p == "interacts_with"
    ][:2]
    assert _export(engine, {"entity_ids": ["HUB"]}, limit=4)["relation_key"].to_pylist() == (
        everything["relation_key"].to_pylist()[:4])


def test_export_nests_annotations_and_evidence_in_stored_order(engine):
    table = _export(engine, {"entity_ids": ["HUB"]}, limit=3)
    assert table.num_rows == 3
    for row in table.to_pylist():
        assert [a["term"] for a in row["annotations"]] == ["term3", "term1", "term2"]
        assert [e["source"] for e in row["evidence"]] == ["source3", "source1", "source2"]
        assert row["evidence_count"] == 3


def test_export_without_endpoint_filters_applies_filters_before_the_limit(engine):
    table = _export(engine, {"predicates": ["affects"]}, limit=2)
    assert table.num_rows == 2
    assert set(table["predicate"].to_pylist()) == {"affects"}
