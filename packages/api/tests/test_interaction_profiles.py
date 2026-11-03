from fastapi.testclient import TestClient
from omnipath_api.server import create_app
from omnipath_api.queries.parquet import ParquetQueries
from omnipath_core.interaction_profiles import (
    BINDING_QUALIFIERS,
    TRANSPORT_QUALIFIERS,
    interaction_label,
    relation_qualifiers,
)
from omnipath_core import relation_key
from omnipath_core.fixtures import rewrite_resource


def test_qualified_labels_survive_parquet_list_and_details(tmp_path):
    folder = tmp_path / "resources" / "test" / "1"
    folder.mkdir(parents=True)
    entities = [
        dict(
            entity_key=k,
            entity_type=t,
            namespace="test",
            identifier=k,
            label=k,
            identifiers=[],
            annotations=[],
        )
        for k, t in [("protein", "protein"), ("chemical", "chemical_entity")]
    ]
    rewrite_resource(folder / "entity.parquet", entities=entities)
    rows = []
    for key, predicate, qualifiers in [
        ("binding", "interacts_with", BINDING_QUALIFIERS),
        ("transport", "affects", TRANSPORT_QUALIFIERS),
        ("generic", "interacts_with", ()),
    ]:
        rows.append(
            dict(
                relation_key=key,
                subject_entity_key="protein",
                object_entity_key="chemical",
                subject_type="protein",
                object_type="chemical_entity",
                predicate=predicate,
                category="interaction",
                sources=["test"],
                evidence_count=1,
                annotations=[dict(term=t, value=v, scope="relation") for t, v in qualifiers],
            )
        )
    rewrite_resource(folder / "relation.parquet", relations=rows)
    engine = ParquetQueries(data_root=tmp_path)
    client = TestClient(create_app(queries=engine))
    response = client.post("/relations/search", json={})
    assert response.status_code == 200
    assert {r["relationPk"]: r["displayLabel"] for r in response.json()["relations"]} == {
        "binding": "Binds",
        "transport": "Transports",
        "generic": None,
    }
    assert {r["relationPk"]: r["qualifiers"] for r in response.json()["relations"]} == {
        "binding": dict(BINDING_QUALIFIERS),
        "transport": dict(TRANSPORT_QUALIFIERS),
        "generic": {},
    }
    filtered = client.post(
        "/relations/search", json={"filters": {"causal_mechanism_qualifier": ["binding"]}}
    ).json()
    assert [r["relationPk"] for r in filtered["relations"]] == ["binding"]
    for key, label in [("binding", "Binds"), ("transport", "Transports")]:
        assert engine.get_relation_by_pk(key)["relation"]["displayLabel"] == label
    assert relation_key("p", "interacts_with", "c", rows[0]["annotations"]) != relation_key(
        "p", "interacts_with", "c"
    )
    assert relation_key("p", "affects", "c", rows[1]["annotations"]) != relation_key(
        "c", "affects", "p", rows[1]["annotations"]
    )


def test_relation_qualifiers_keep_effect_terms_only():
    annotations = [
        dict(term="object_direction_qualifier", value="increased"),
        dict(term="object_aspect_qualifier", value="activity"),
        dict(term="object_aspect_qualifier", value="expression"),
        dict(term="causal_mechanism_qualifier", value="binding", scope="evidence"),
        dict(term="stoichiometry", value="1"),
    ]
    assert relation_qualifiers(annotations) == {
        "object_direction_qualifier": "increased",
        "object_aspect_qualifier": "activity",
    }


def test_transport_inhibition_is_not_labelled_transports():
    annotations = [dict(term=t, value=v) for t, v in TRANSPORT_QUALIFIERS]
    annotations.append(dict(term="object_direction_qualifier", value="decreased"))
    assert interaction_label("affects", annotations) is None
