from fastapi.testclient import TestClient

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_api.queries.store.inventory import ReleaseStore
from omnipath_core.fixtures import write_resource


def write(root, version, label, key, curie, kind="chemical_entity", taxon=None):
    entity = dict(
        entity_key=key,
        label=label,
        entity_type=kind,
        namespace="inchikey" if kind == "chemical_entity" else "uniprot",
        identifier=key,
        taxon=taxon,
        identifiers=[dict(ns=curie.split(":")[0].lower(), id=curie)],
    )
    write_resource(root / "resources/test" / version, [entity])


def test_examples_are_cached_release_specific_and_do_not_replace_search(tmp_path, monkeypatch):
    write(tmp_path, "1", "Aspirin", "a" * 64, "CHEBI:15365")
    store = ReleaseStore(tmp_path)
    store.publish(dict(schema_version=1, version="1", resources={"test": "1"}))
    engine = ParquetQueries(tmp_path)
    first = engine.get_entity_examples()
    assert len(first["entities"]) == 1 and first["entities"][0]["label"] == "Aspirin"
    assert first["nextCursor"] is None
    # the grouped landing page shows the same examples as group cards
    assert [g["entity"]["displayName"] for g in first["groups"]] == ["Aspirin"]
    assert len(list((tmp_path / ".presentation/examples-v3").glob("*.json"))) == 1
    # A fresh engine reads the prepared page without selecting or hydrating again.
    cached = ParquetQueries(tmp_path)
    monkeypatch.setattr(
        cached,
        "_fetch_entities_by_keys",
        lambda *a, **kw: (_ for _ in ()).throw(AssertionError("cache miss")),
    )
    assert cached.get_entity_examples() == first
    write(tmp_path, "2", "Caffeine", "b" * 64, "CHEBI:27732")
    assert engine.get_entity_examples()["entities"][0]["label"] == "Caffeine"
    with engine.release_scope("1"):
        assert engine.get_entity_examples()["entities"][0]["label"] == "Aspirin"
    assert engine.search_entities_api("Aspirin")["entities"] == []
    client = TestClient(create_app(queries=engine))
    result = client.get("/entities/examples", headers={"X-OmniPath-Release": "1"})
    assert result.status_code == 200 and result.json()["kind"] == "examples"
    assert result.json()["entities"][0]["label"] == "Aspirin"


def test_examples_are_chosen_by_identifier_not_label(tmp_path):
    # 'lactate' also labels a wax ester: a label does not identify an example
    write(tmp_path, "1", "Glucose", "a" * 64, "CHEBI:99999")
    engine = ParquetQueries(tmp_path)
    assert engine.get_entity_examples()["entities"] == []
    write(tmp_path, "2", "D-Glucose", "b" * 64, "CHEBI:17234")
    assert [e["label"] for e in engine.get_entity_examples()["entities"]] == ["D-Glucose"]
