import json

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_api.models import ResourcesResponse
from omnipath_core.fixtures import write_resource_version as _write_resource


def test_resource_details_preserve_stats_and_snapshot(tmp_path):
    folder = _write_resource(tmp_path, "cellchat", "1")
    stats = dict(
        input_entities=12, resolved_entities=7, unresolved_entities=3, not_applicable_entities=2
    )
    (folder / "resolution_stats.json").write_text(json.dumps(stats))
    metadata = dict(
        website="https://example.org", license="MIT", downloaded_at="2026-05-01 12:00:00"
    )
    (folder / "build_manifest.json").write_text(json.dumps({"resource_metadata": metadata}))
    items = ParquetQueries(data_root=tmp_path).list_resource_catalog()
    item = ResourcesResponse(resources=items).model_dump()["resources"][0]
    assert item["resolution_stats"] == stats
    assert item["website"] == metadata["website"]
    assert item["downloaded_at"] == metadata["downloaded_at"]
    assert item["metadata_source"] == "build snapshot"


def test_legacy_metadata_does_not_invent_download_date_or_stats(tmp_path):
    folder = _write_resource(tmp_path, "cellchat", "1")
    (folder / "build_manifest.json").write_text(json.dumps({"created_at": "2026-09-05"}))
    item = ParquetQueries(data_root=tmp_path).list_resource_catalog()[0]
    assert item["downloaded_at"] is None
    assert item["resolution_stats"] is None
    assert item["website"] == "https://github.com/jinworks/CellChat"
    assert item["input_module_url"].endswith("/pypath/inputs_v2/cellchat.py")
    assert item["license"] == "GNU General Public License v3.0"


def test_resource_browsing_metadata_survives_api_serialization(tmp_path):
    _write_resource(tmp_path, "reactome", "1")
    items = ParquetQueries(data_root=tmp_path).list_resource_catalog()
    item = ResourcesResponse(resources=items).model_dump()["resources"][0]
    assert item["description"]
    tags = {tag["id"]: tag for tag in item["tags"]}
    assert tags["pathways"] == {
        "id": "pathways",
        "label": "Pathways",
        "dimension": "content",
        "description": "Biological pathways and their constituent entities or activities.",
    }
    assert item["license_use"] == {"academic": "allowed", "commercial": "allowed"}
