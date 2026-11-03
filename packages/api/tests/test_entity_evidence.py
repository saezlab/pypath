from fastapi.testclient import TestClient

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.server import create_app
from omnipath_core.fixtures import write_resource

KEY = "a" * 64


def resource(root, name, counts):
    entities = [
        dict(
            entity_key=key,
            entity_type="small_molecule",
            namespace="test",
            identifier=key[:4],
            label=key[:4],
            evidence=[dict(source=name, dataset="d", row_id=f"{key[:1]}{i}") for i in range(n)],
        )
        for key, n in counts.items()
    ]
    write_resource(root / "resources" / name / "1", entities)


def test_evidence_is_counted_in_details_and_paged_from_its_endpoint(tmp_path):
    resource(tmp_path, "one", {KEY: 3, "c" * 64: 0})
    resource(tmp_path, "two", {KEY: 2})
    engine = ParquetQueries(tmp_path)
    client = TestClient(create_app(queries=engine))
    details = client.get(f"/entities/{KEY}").json()["entity"]
    pages = [client.get(f"/entities/{KEY}/evidence?limit=2&offset={o}").json() for o in (0, 2, 4)]
    assert details["molecularEvidence"] == [] and details["molecularEvidenceTotal"] == 5
    assert [p["nextCursor"] for p in pages] == ["2", "4", None]
    items = [(e["source"], e["row_id"]) for p in pages for e in p["evidence"]]
    assert items == [("one", "a0"), ("one", "a1"), ("one", "a2"), ("two", "a0"), ("two", "a1")]
    assert client.get(f"/entities/{'d' * 64}/evidence").status_code == 404
    assert client.get(f"/entities/{'c' * 64}/evidence").json()["evidenceTotal"] == 0
