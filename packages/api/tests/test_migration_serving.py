"""The migrated schema-3 resource contract serves without a build or resolver."""

import hashlib
import importlib.util
import io
import json
from pathlib import Path

from fastapi.testclient import TestClient
import pyarrow.parquet as pq
import pytest
from omnipath_core import SERVING_SCHEMA_VERSION

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.releases import ReleaseStore
from omnipath_api.server import create_app
from omnipath_api.settings import Settings
from omnipath_core.schema import PUBLISHED_TABLES


def fixture_helpers():
    """Use the workspace's schema fixtures without modifying the import path."""
    path = Path(__file__).with_name("release_fixture.py")
    spec = importlib.util.spec_from_file_location("migration_release_fixture", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def resource_version(root, source="signor", version="1.0.0", *, marker="original"):
    """Write two qualified statements from two original source fixture rows."""
    fixture = fixture_helpers()
    subject_id, target_id = ("P04637", "P38398") if source == "signor" else ("P05005", "P06493")
    subject = fixture.entity(
        subject_id,
        "protein",
        "uniprot",
        taxon="9606",
        label=f"TP53 {marker}",
        aliases=[("uniprot", "Q15086"), ("genesymbol", "TP53")] if source == "signor" else (),
    )
    target = fixture.entity(target_id, "protein", "uniprot", taxon="9606", label="BRCA1")
    measurement = dict(has_numeric_value=12.5, has_unit="nM", comparator="<=", source_field="IC50")
    relations = []
    payloads = []
    for index, direction in enumerate(("increased", "decreased"), 1):
        annotations = [
            fixture.annotation(
                "object_aspect_qualifier", "activity", source=source, dataset="interactions"
            ),
            fixture.annotation(
                "object_direction_qualifier", direction, source=source, dataset="interactions"
            ),
        ]
        row = fixture.relation(
            subject,
            "affects",
            target,
            source=source,
            dataset="interactions",
            row_id=f"000{index}",
            upstream_id=f"SOURCE-{index}",
            annotations=annotations,
            sign=1 if index == 1 else -1,
        )
        row.update(taxon="9606", category="interaction", interaction_class="protein_protein")
        row["evidence"][0]["annotations"].append(
            fixture.annotation(
                "has_quantitative_value",
                quantity=measurement,
                source=source,
                dataset="interactions",
            )
        )
        relations.append(row)
        payloads.append(
            fixture.payload(
                row,
                {"version_marker": marker, "effect": direction, "IC50": "<=12.5 nM"},
                source=source,
                row_id=f"000{index}",
            )
        )
    directory = fixture.write_resource(
        root, source, [subject, target], relations, payloads, version
    )
    return dict(
        directory=directory, subject=subject, target=target, relations=relations, marker=marker
    )


def publish_snapshot(root, resources):
    store = ReleaseStore(root)
    with pytest.warns(UserWarning, match="taxonomy reference"):
        store.publish(dict(schema_version=1, version="2026.09", resources=resources))
    return root / "releases/2026.09.json"


def catalog(client, release=None):
    response = client.get(
        "/api/resources", params={"shape": "svelte", **({"release": release} if release else {})}
    )
    assert response.status_code == 200, response.text
    return {row["resource_id"]: row["version"] for row in response.json()["resources"]}


def payload_marker(client, relation_key, release=None):
    response = client.get(
        f"/api/relations/{relation_key}/payloads", params={"release": release} if release else {}
    )
    assert response.status_code == 200, response.text
    return response.json()["payloads"][0]["payload"]["version_marker"]


def test_migrated_schema_search_details_evidence_filters_and_exports(tmp_path):
    published = resource_version(tmp_path)
    directory = published["directory"]
    manifest = json.loads((directory / "build_manifest.json").read_text())
    assert manifest["serving_schema_version"] == SERVING_SCHEMA_VERSION
    before = {}
    for table, schema in PUBLISHED_TABLES.items():
        filename = f"{table}.parquet"
        path = directory / filename
        assert pq.read_schema(path).equals(schema)
        before[filename] = hashlib.sha256(path.read_bytes()).hexdigest()
        assert before[filename] == manifest["files"][filename]["sha256"]
    engine = ParquetQueries(tmp_path)
    with TestClient(create_app(queries=engine, settings=Settings(read_only=True))) as client:
        found = client.get("/api/entities/search", params={"q": "Q15086"})
        assert found.status_code == 200, found.text
        assert [row["entityPk"] for row in found.json()["entities"]] == [
            published["subject"]["entity_key"]
        ]
        details = client.get(
            f"/api/entities/{published['subject']['entity_key']}",
            params={"includeRelationships": "false"},
        )
        assert details.status_code == 200, details.text
        identifiers = details.json()["entity"]["identifiers"]
        assert any(row["identifier"] == "Q15086" for row in identifiers)
        assert details.json()["entity"]["canonicalIdentifier"] == "P04637"
        # Detail lookup uses the exact published key; alias search does not change identity.
        assert client.get("/api/entities/P04637").status_code == 404
        expected = published["relations"][0]
        relation_key = expected["relation_key"]
        filters = {
            "sources": ["signor"],
            "taxonomy_ids": ["9606"],
            "object_aspect_qualifier": ["activity"],
            "object_direction_qualifier": ["increased"],
        }
        searched = client.post("/api/relations/search", json={"filters": filters})
        assert searched.status_code == 200, searched.text
        assert searched.json()["total"] == 1
        assert [row["relationPk"] for row in searched.json()["relations"]] == [relation_key]
        detail = client.get(f"/api/relations/{relation_key}")
        assert detail.status_code == 200, detail.text
        assert detail.json()["relation"]["subjectEntityPk"] == published["subject"]["entity_key"]
        assert detail.json()["relation"]["objectEntityPk"] == published["target"]["entity_key"]
        evidence = client.get(f"/api/relations/{relation_key}/evidence")
        assert evidence.status_code == 200, evidence.text
        record = evidence.json()["evidence"][0]
        assert record["source"] == "signor" and record["dataset"] == "interactions"
        assert record["upstreamId"] == "SOURCE-1"
        quantity = next(
            row["quantity"]
            for row in record["recordAttributes"]
            if row["term"] == "has_quantitative_value"
        )
        assert {
            key: quantity[key]
            for key in ("has_numeric_value", "has_unit", "comparator", "source_field")
        } == {
            "has_numeric_value": 12.5,
            "has_unit": "nM",
            "comparator": "<=",
            "source_field": "IC50",
        }
        raw = client.get(f"/api/relations/{relation_key}/payloads").json()["payloads"]
        assert raw == [
            {
                "source": "signor",
                "rowId": "0001",
                "payload": {
                    "version_marker": "original",
                    "effect": "increased",
                    "IC50": "<=12.5 nM",
                },
            }
        ]
        exported = client.post("/api/export", json={"filters": filters})
        assert exported.status_code == 200, exported.text
        rows = pq.read_table(io.BytesIO(exported.content)).to_pylist()
        assert [row["relation_key"] for row in rows] == [relation_key]
        assert rows[0]["evidence"][0]["upstream_id"] == "SOURCE-1"
        assert rows[0]["evidence"][0]["annotations"][-1]["quantity"]["comparator"] == "<="
        empty_filters = {**filters, "entity_ids": ["no-such-identifier"]}
        assert (
            client.post("/api/relations/search", json={"filters": empty_filters}).json()["total"]
            == 0
        )
    assert {
        name: hashlib.sha256((directory / name).read_bytes()).hexdigest() for name in before
    } == before


def test_independent_resource_updates_advance_default_latest_and_keep_named_snapshot(tmp_path):
    old = resource_version(tmp_path, version="1.9.0", marker="old")
    resource_version(tmp_path, source="reactome", version="3.0.0", marker="other old")
    snapshot = publish_snapshot(tmp_path, {"signor": "1.9.0", "reactome": "3.0.0"})
    pinned_text = snapshot.read_bytes()
    engine = ParquetQueries(tmp_path)
    with TestClient(create_app(queries=engine, settings=Settings(read_only=True))) as client:
        assert catalog(client) == {"signor": "1.9.0", "reactome": "3.0.0"}
        relation_key = old["relations"][0]["relation_key"]
        # Populate current caches before publishing another immutable resource version.
        assert payload_marker(client, relation_key) == "old"
        old_detail = client.get(
            f"/entities/{old['subject']['entity_key']}", params={"includeRelationships": "false"}
        ).json()
        resource_version(tmp_path, version="1.10.0", marker="new")
        for release in (None, "latest"):
            assert catalog(client, release) == {"signor": "1.10.0", "reactome": "3.0.0"}
            assert payload_marker(client, relation_key, release) == "new"
        current_detail = client.get(
            f"/entities/{old['subject']['entity_key']}", params={"includeRelationships": "false"}
        ).json()
        assert current_detail["entity"]["label"] == "TP53 new"
        assert old_detail["entity"]["label"] == "TP53 old"
        # A second resource advances independently, without publishing a global release.
        resource_version(tmp_path, source="reactome", version="3.1.0", marker="other new")
        assert catalog(client) == {"signor": "1.10.0", "reactome": "3.1.0"}
        assert catalog(client, "2026.09") == {"signor": "1.9.0", "reactome": "3.0.0"}
        assert payload_marker(client, relation_key, "2026.09") == "old"
        pinned_detail = client.get(
            f"/entities/{old['subject']['entity_key']}",
            headers={"X-OmniPath-Release": "2026.09"},
            params={"includeRelationships": "false"},
        ).json()
        assert pinned_detail["entity"]["label"] == "TP53 old"
        assert client.get("/releases").json()["default"] == "latest"
    assert snapshot.read_bytes() == pinned_text
