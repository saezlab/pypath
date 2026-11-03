"""Malformed metadata fails before anything is installed or queried."""

from copy import deepcopy
import json

import httpx
import pytest

from omnipath_client import ClientError, available, connect, init
from omnipath_client.contracts import validate_catalog, validate_installation


def resource():
    return {
        "resource_id": "alpha",
        "version": "1.0",
        "description": "Additive metadata",
        "files": [
            {"name": f"{table}.parquet", "size_bytes": 100} for table in ("entity", "relation")
        ],
    }


@pytest.mark.parametrize(
    "field,value",
    [
        ("resource_id", 123),
        ("version", 1),
        ("files", None),
        ("files", {}),
        ("files", []),
        ("files", [None]),
    ],
)
def test_invalid_catalog_rejected_by_init(field, value, tmp_path):
    record = resource()
    record[field] = value
    transport = httpx.MockTransport(lambda _: httpx.Response(200, json={"resources": [record]}))
    with httpx.Client(transport=transport) as http:
        with pytest.raises(ClientError):
            init(data_dir=tmp_path / "invalid", http_client=http)
    assert not (tmp_path / "invalid").exists()


@pytest.mark.parametrize(
    "field,value",
    [
        ("size_bytes", True),
        ("size_bytes", -1),
        ("size_bytes", "100"),
        ("rows", False),
        ("rows", -1),
        ("url", 123),
        ("url", "file:///tmp/test"),
        ("name", "../entity.parquet"),
    ],
)
def test_invalid_artifact_metadata_rejected(field, value):
    record = resource()
    record["files"][0][field] = value
    with pytest.raises(ClientError):
        validate_catalog([record])


def test_required_artifacts_and_duplicates():
    record = resource()
    record["files"].pop()
    with pytest.raises(ClientError, match="must contain"):
        validate_catalog([record])
    record = resource()
    record["files"].append(deepcopy(record["files"][0]))
    with pytest.raises(ClientError, match="duplicate artifact"):
        validate_catalog([record])
    with pytest.raises(ClientError, match="Duplicate resource"):
        validate_catalog([resource(), resource()])


def test_installation_manifest_is_validated_before_open(tmp_path):
    manifest = {
        "format": 1,
        "api_url": "https://api.example/api",
        "release": "1.0",
        "resources": [resource()],
    }
    assert validate_installation(manifest)["release"] == "1.0"
    (tmp_path / "installed.json").write_text(json.dumps({**manifest, "release": ""}))
    with pytest.raises(ClientError, match="release"):
        connect(tmp_path)
    (tmp_path / "installed.json").write_text("{")
    with pytest.raises(ClientError, match="Invalid installed.json"):
        connect(tmp_path)
    assert validate_installation(manifest)["materialize"] is False
    with pytest.raises(ClientError, match="materialize"):
        validate_installation({**manifest, "materialize": "yes"})
    manifest["format"] = True
    with pytest.raises(ClientError, match="installation format"):
        validate_installation(manifest)


def test_additive_metadata_is_preserved_without_aliasing():
    record = resource()
    catalog = validate_catalog([record])
    assert catalog["alpha"]["description"] == "Additive metadata"
    record["files"].clear()
    assert len(catalog["alpha"]["files"]) == 2


@pytest.mark.parametrize("payload", [None, [], {"resources": None}, {"resources": "alpha"}])
def test_invalid_catalog_envelopes(payload):
    with httpx.Client(
        transport=httpx.MockTransport(lambda _: httpx.Response(200, json=payload))
    ) as http:
        with pytest.raises(ClientError):
            available(http_client=http)
