import json

from served import served_queries
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
    items = served_queries(data_root=tmp_path).list_resource_catalog()
    item = ResourcesResponse(resources=items).model_dump()["resources"][0]
    assert item["resolution_stats"] == stats
    assert item["website"] == metadata["website"]
    assert item["downloaded_at"] == metadata["downloaded_at"]
    assert item["metadata_source"] == "build snapshot"


def test_legacy_metadata_does_not_invent_download_date_or_stats(tmp_path):
    folder = _write_resource(tmp_path, "cellchat", "1")
    (folder / "build_manifest.json").write_text(json.dumps({"created_at": "2026-09-05"}))
    item = served_queries(data_root=tmp_path).list_resource_catalog()[0]
    assert item["downloaded_at"] is None
    assert item["resolution_stats"] is None
    assert item["website"] == "https://github.com/jinworks/CellChat"
    assert item["input_module_url"].endswith(
        "/packages/resources/omnipath_resources/inputs_v2/cellchat.py"
    )
    assert item["license"] == "GNU General Public License v3.0"


def test_resource_browsing_metadata_survives_api_serialization(tmp_path):
    _write_resource(tmp_path, "reactome", "1")
    items = served_queries(data_root=tmp_path).list_resource_catalog()
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


def test_a_resource_record_carries_its_three_names(tmp_path):
    _write_resource(tmp_path, "reactome", "1")
    item = served_queries(data_root=tmp_path).list_resource_catalog()[0]
    record = ResourcesResponse(resources=[item]).model_dump()["resources"][0]
    assert record["short"] and record["full"] and isinstance(record["synonyms"], list)
    assert record["full"] == item["full_name"]


def test_a_resource_is_selected_by_its_slug_short_name_or_synonym_in_any_case(tmp_path):
    _write_resource(tmp_path, "reactome", "1")
    engine = served_queries(data_root=tmp_path)
    from omnipath_core import resource_metadata

    entry = resource_metadata._catalog()["reactome"]
    for spelling in ("reactome", "REACTOME", entry["short"].upper(), *entry.get("synonyms", [])):
        infos = engine._selected_resource_infos([spelling])
        assert [info["resource"] for info in infos] == ["reactome"], spelling


def test_a_name_that_is_nobodys_is_still_refused(tmp_path):
    import pytest

    _write_resource(tmp_path, "reactome", "1")
    engine = served_queries(data_root=tmp_path)
    with pytest.raises(ValueError, match="unavailable"):
        engine._selected_resource_infos(["nonexistentresource"])


def test_the_sources_of_a_search_are_resolved_to_resource_slugs():
    from omnipath_api.queries.filters import normalize_filters

    assert normalize_filters({"sources": ["REACTOME", "Reactome", "nonexistentresource"]})[
        "sources"
    ] == ["reactome", "nonexistentresource"]
    assert normalize_filters({"resources": ["Reactome"]})["sources"] == ["reactome"]
