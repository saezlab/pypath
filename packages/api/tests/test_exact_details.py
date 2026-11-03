from fastapi.testclient import TestClient
from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_core.fixtures import rewrite_resource, write_resource


def test_details_never_resolve_accessions_aliases_or_text(tmp_path, monkeypatch):
    folder = tmp_path / "resources" / "test" / "1"
    folder.mkdir(parents=True)
    rows = [
        dict(
            entity_key=key,
            entity_type=kind,
            namespace=ns,
            identifier="1",
            label=label,
            identifiers=[
                dict(ns=ns, id="1", source="test", is_canonical=True),
                dict(ns="name", id="shared-alias", source="test", is_canonical=False),
            ],
            annotations=[],
        )
        for key, kind, ns, label in [
            ("complex-key", "macromolecular_complex", "corum", "BCL6-HDAC4"),
            ("trait-key", "ontology_class", "macdb_trait", "Acinar carcinoma"),
        ]
    ]
    write_resource(folder, rows, [])
    engine = ParquetQueries(data_root=tmp_path)

    def unexpected_search(*args, **kwargs):
        raise AssertionError("Details must never invoke identifier resolution or search")

    monkeypatch.setattr(engine, "get_entities_by_pks", unexpected_search)
    monkeypatch.setattr(engine, "resolve_entity_keys", unexpected_search)
    monkeypatch.setattr(engine, "search_entities_api", unexpected_search)
    client = TestClient(create_app(queries=engine))
    for key, label in [("complex-key", "BCL6-HDAC4"), ("trait-key", "Acinar carcinoma")]:
        response = client.get("/entities/" + key)
        assert response.status_code == 200
        assert response.json()["entity"]["label"] == label
    for token in [
        "1",
        "corum|1",
        "BCL6-HDAC4",
        "Acinar",
        "shared-alias",
        "COMPLEX-KEY",
        "missing-key",
    ]:
        assert client.get("/entities/" + token).status_code == 404


def test_detail_composition_uses_exact_endpoint_keys_and_keeps_quantities(tmp_path):
    folder = tmp_path / "resources" / "test" / "1"
    folder.mkdir(parents=True)
    entities = [
        dict(
            entity_key=key,
            entity_type=kind,
            namespace="test",
            identifier="1",
            label=label,
            identifiers=[],
            annotations=[],
        )
        for key, kind, label in [
            ("food-key", "food", "Food"),
            ("chemical-key", "chemical_entity", "Constituent"),
        ]
    ]
    rewrite_resource(folder / "entity.parquet", entities=entities)
    relation = dict(
        relation_key="composition",
        subject_entity_key="food-key",
        object_entity_key="chemical-key",
        predicate="has_part",
        category="membership",
        sources=["test"],
        evidence_count=1,
        annotations=[
            dict(
                term="PATO:0000033",
                source="test",
                quantity=dict(
                    has_numeric_value=0.0, has_unit="mg/g", source_field="member_content"
                ),
            )
        ],
    )
    rewrite_resource(folder / "relation.parquet", relations=[relation])
    engine = ParquetQueries(data_root=tmp_path)
    detail = TestClient(create_app(queries=engine)).get("/entities/food-key").json()
    assert detail["relationshipsTotal"] == 1
    row = detail["relationships"][0]
    assert row["subjectEntity"]["entityPk"] == "food-key"
    assert row["objectEntity"]["entityPk"] == "chemical-key"
    assert row["objectEntity"]["label"] == "Constituent"
    assert row["annotations"][0]["quantity"]["has_numeric_value"] == 0.0
    assert engine.get_entity_details("1") is None


def test_integrated_entities_retain_distinct_typed_quantities(tmp_path):
    engine = ParquetQueries(data_root=tmp_path)
    annotations = [
        dict(
            term="has_quantitative_value",
            value=None,
            source="model",
            dataset="reactions",
            quantity=dict(has_numeric_value=value, source_field=field),
        )
        for value, field in [(-1000.0, "lower_bound"), (1000.0, "upper_bound")]
    ]
    base = dict(
        entity_key="reaction-key",
        namespace="rhea",
        identifier="10000",
        label="Reaction",
        identifiers=[],
        annotations=annotations,
    )
    merged = engine._merge_entity_row_group([base, dict(base)])
    assert len(merged["annotations"]) == 2
    assert {a["quantity"]["source_field"] for a in merged["annotations"]} == {
        "lower_bound",
        "upper_bound",
    }
