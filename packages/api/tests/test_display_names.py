from fastapi.testclient import TestClient
from served import served_queries
from omnipath_api.server import create_app
from omnipath_core.fixtures import write_resource


def _engine(tmp_path, **entity):
    write_resource(
        tmp_path / "resources/test/1",
        [dict(entity_key="a", entity_type="chemical_entity", namespace="inchikey", **entity)],
    )
    return TestClient(create_app(queries=served_queries(tmp_path)))


def test_resolver_label_is_the_display_name_on_every_page(tmp_path):
    client = _engine(
        tmp_path,
        identifier="QNAYBMKLOCPYGJ-UQEXSWPGSA-N",
        label="alanine-d7",
        identifiers=[dict(ns="name", id=f"Long chemical name {i}") for i in range(25)]
        + [dict(ns="name", id="Alanine")],
    )
    first = client.get("/entities/a?includeRelationships=false").json()["entity"]
    assert first["displayName"] == "alanine-d7"
    assert first["identifiersNextCursor"] == "20"
    later = client.get("/entities/a?includeRelationships=false&detail_offset=20").json()["entity"]
    assert later["displayName"] == "alanine-d7"


def test_numeric_chemical_labels_fall_back_to_a_named_identifier(tmp_path):
    client = _engine(
        tmp_path,
        identifier="1234",
        label="1234",
        identifiers=[dict(ns="chebi", id="CHEBI:15377")],
    )
    entity = client.get("/entities/a?includeRelationships=false").json()["entity"]
    assert entity["displayName"] == "CHEBI:15377"
