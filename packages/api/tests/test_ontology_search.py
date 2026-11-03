import tempfile
from pathlib import Path

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from fastapi.testclient import TestClient
from omnipath_core.fixtures import write_resource


def test_scoped_ontology_search_deduplicates_and_never_leaks_global_terms():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for source in ("a", "b"):
            write_resource(
                root / "resources" / source / "v1",
                [
                    dict(
                        entity_key="p",
                        identifier="P00533",
                        label="EGFR",
                        namespace="uniprot",
                        entity_type="protein",
                    ),
                    dict(
                        entity_key="t",
                        identifier="GO:1",
                        label="Signaling",
                        namespace="go",
                        entity_type="ontology_class",
                    ),
                    dict(
                        entity_key="u",
                        identifier="GO:2",
                        label="Unrelated",
                        namespace="go",
                        entity_type="ontology_class",
                    ),
                ],
                [dict(relation_key="r", subject_entity_key="p", object_entity_key="t")],
            )
        engine = ParquetQueries(root)
        payload = dict(entityPks=["p"], scoped=True)
        items = engine.search_ontology(payload)
        assert [i["termId"] for i in items] == ["GO:1"]
        assert items[0]["annotatedEntityCount"] == 1
        assert items[0]["annotatedRelationCount"] == 1
        assert engine.search_ontology(dict(scoped=True)) == []
        assert engine.search_ontology(dict(payload, query="unrelated")) == []
        assert engine.search_ontology(dict(payload, ontologyIds=["other"])) == []
        assert engine.search_ontology(dict(payload, offset=1)) == []
        assert engine.search_ontology(payload, counts=True) == [
            dict(ontologyId="go", scopedCount=1)
        ]
        assert [
            i["termId"] for i in engine.search_ontology(dict(scoped=True, termIds=["GO:2"]))
        ] == ["GO:2"]
        assert len(engine.search_ontology()) == 2

        client = TestClient(create_app(queries=engine))
        response = client.post("/ontology/scoped-search", json={"entityPks": ["p"]})
        assert response.status_code == 200
        assert response.json()[0]["termId"] == "GO:1"
        assert client.post("/ontology/scoped-search", json={}).json() == []
        assert client.post("/ontology/ontology-id-counts", json={"entityPks": ["p"]}).json() == [
            dict(ontologyId="go", scopedCount=1)
        ]
