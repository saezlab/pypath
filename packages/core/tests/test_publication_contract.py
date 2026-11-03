"""The public contract matches immutable publications, not an older prototype."""

from copy import deepcopy
import pytest

from omnipath_core import (
    BUILD_SCHEMA_VERSION,
    RELEASE_SCHEMA_VERSION,
    RESOURCE_FILES,
    SERVING_FILES,
    SERVING_SCHEMA_VERSION,
    BuildManifest,
    ManifestFile,
    ReleaseManifest,
    validate_build_manifest,
    validate_release_manifest,
)


def publication():
    return {
        "schema_version": 1,
        "serving_schema_version": SERVING_SCHEMA_VERSION,
        "resource": "signor",
        "version": "2.7",
        "created_at": "2026-10-03T00:00:00+00:00",
        "provenance": {"software": {"sha256": "b" * 64}},
        "files": {
            name: {"rows": 0, "sha256": "a" * 64, "size_bytes": 123}
            for name in RESOURCE_FILES + SERVING_FILES
        },
    }


def test_public_models_round_trip_current_filename_keyed_manifest():
    data = publication()
    parsed = validate_build_manifest(data)
    assert isinstance(parsed, BuildManifest)
    assert parsed.schema_version == BUILD_SCHEMA_VERSION == 1
    assert parsed.serving_schema_version == SERVING_SCHEMA_VERSION == 6
    assert parsed.files["entity.parquet"] == ManifestFile(size_bytes=123, rows=0, sha256="a" * 64)
    assert parsed.to_dict() == data
    assert BuildManifest.from_dict(data, strict_fields=True) == parsed
    assert "metadata" not in parsed.to_dict()
    assert "path" not in parsed.files["entity.parquet"].to_dict()


@pytest.mark.parametrize(
    "field,value",
    [
        ("schema_version", "2.0"),
        ("schema_version", True),
        ("serving_schema_version", "3"),
        ("serving_schema_version", False),
        ("resource", "../../signor"),
        ("version", "latest"),
        ("created_at", None),
        ("files", []),
        ("files", {}),
    ],
)
def test_build_rejects_old_or_ambiguous_shapes(field, value):
    data = publication()
    data[field] = value
    with pytest.raises(ValueError):
        validate_build_manifest(data)


@pytest.mark.parametrize(
    "field,value",
    [
        ("rows", True),
        ("rows", -1),
        ("size_bytes", 1.5),
        ("size_bytes", "123"),
        ("sha256", "A" * 64),
        ("sha256", "abc"),
        ("bytes", 123),
        ("num_rows", 0),
        ("path", "entity.parquet"),
    ],
)
def test_file_contract_rejects_stale_fields_and_nonliteral_metadata(field, value):
    data = publication()
    data["files"]["entity.parquet"][field] = value
    with pytest.raises(ValueError):
        validate_build_manifest(data)


def test_reader_metadata_extensions_are_preserved_and_strictness_is_a_policy():
    data = publication()
    data.pop("created_at")  # existing pinned manifests need not contain this field
    data["future_inspection"] = {"counts": [1, 2]}
    assert validate_build_manifest(data).to_dict() == data
    with pytest.raises(ValueError, match="Unknown"):
        validate_build_manifest(data, strict_fields=True)


def test_release_model_keeps_independent_resource_and_release_versions():
    data = {"schema_version": 1, "version": "2026.10", "resources": {"signor": "2.7", "chebi": "9"}}
    parsed = validate_release_manifest(data, strict_fields=True)
    assert isinstance(parsed, ReleaseManifest)
    assert parsed.schema_version == RELEASE_SCHEMA_VERSION == 1
    assert parsed.to_dict() == data
    assert parsed.resources["signor"] == "2.7"


@pytest.mark.parametrize("value", [True, "1", 2, None])
def test_release_schema_version_is_a_literal_integer(value):
    with pytest.raises(ValueError):
        validate_release_manifest(
            {"schema_version": value, "version": "1", "resources": {"go": "2"}}
        )


def test_inventory_name_policy_is_explicit_and_does_not_change_pg_defaults():
    data = {"schema_version": 1, "version": "2026.10", "resources": {"Historical.Name": "working"}}

    def inventory_name(value):
        if not isinstance(value, str) or "/" in value:
            raise ValueError("Unsafe inventory name")
        return value

    assert (
        validate_release_manifest(
            data,
            resource_name_validator=inventory_name,
            resource_version_validator=inventory_name,
        ).to_dict()
        == data
    )
    with pytest.raises(ValueError):
        validate_release_manifest(data)
    assert (
        validate_release_manifest(
            {"schema_version": 1, "version": "1", "resources": {}},
            require_resources=False,
        ).resources
        == {}
    )


def test_taxonomy_structure_validated_without_reading_a_reference():
    data = {
        "schema_version": 1,
        "version": "1",
        "resources": {"go": "2"},
        "references": {
            "taxonomy": {"version": "c" * 64, "taxon_count": 0, "missing_taxon_ids": ["9606"]},
        },
    }
    assert validate_release_manifest(data, strict_fields=True).to_dict() == data
    for field, value in (
        ("version", "latest"),
        ("taxon_count", True),
        ("missing_taxon_ids", [9606]),
    ):
        bad = deepcopy(data)
        bad["references"]["taxonomy"][field] = value
        with pytest.raises(ValueError):
            validate_release_manifest(bad)
