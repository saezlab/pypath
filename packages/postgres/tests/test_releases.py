"""Pinned release validation uses immutable contracts and metadata only."""

from dataclasses import FrozenInstanceError
import hashlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from omnipath_core.fixtures import write_manifest
from omnipath_core.fixtures import write_resource as write_tables
from omnipath_core.versioning import SERVING_SCHEMA_VERSION
from omnipath_postgres.releases import (
    FILE_SCHEMAS,
    ReleaseValidationError,
    load_release,
    verify_release,
)


def write_resource(root, source="signor", version="1.0.0"):
    directory = root / "resources" / source / version
    write_tables(directory)
    manifest = write_manifest(
        directory, source, version, provenance={"dependencies": {"pypath-omnipath": "fixture"}}
    )
    write_json(directory / "build_manifest.json", manifest)
    return directory, manifest


def write_json(path: Path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


@pytest.fixture
def release(tmp_path):
    root = tmp_path / "data"
    directory, build_manifest = write_resource(root)
    manifest = {"schema_version": 1, "version": "2026.09", "resources": {"signor": "1.0.0"}}
    path = tmp_path / "release.json"
    write_json(path, manifest)
    return root, path, directory, manifest, build_manifest


def test_selection_is_pinned_deeply_immutable_and_preserves_build_manifest(release):
    root, path, directory, manifest, _ = release
    result = load_release(root, path)
    assert result.version == "2026.09"
    assert result.schema_version == 1
    assert result.canonical_json == json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    assert result.sha256 == hashlib.sha256(result.canonical_json.encode()).hexdigest()
    (resource,) = result.resources
    assert resource.source == "signor"
    assert resource.version == "1.0.0"
    assert resource.schema_version == 1
    assert resource.serving_schema_version == SERVING_SCHEMA_VERSION
    assert resource.directory == directory
    assert resource.manifest_json == (directory / "build_manifest.json").read_text()
    assert resource.manifest_sha256 == hashlib.sha256(resource.manifest_json.encode()).hexdigest()
    assert resource.table_path("entity") == directory / "entity.parquet"
    assert resource.table_path("relation") == directory / "relation.parquet"
    assert resource.table_path("evidence_payloads") == directory / "evidence_payloads.parquet"
    assert all(artifact.rows == 0 for artifact in resource.files.values())
    with pytest.raises(FrozenInstanceError):
        resource.version = "2"
    with pytest.raises(TypeError):
        result.manifest["resources"]["signor"] = "2"
    with pytest.raises(TypeError):
        resource.manifest["provenance"]["dependencies"]["pypath-omnipath"] = "changed"
    with pytest.raises(TypeError):
        resource.files["entity.parquet"] = None


def test_canonical_digest_ignores_json_order_and_whitespace(release):
    root, path, _, manifest, _ = release
    write_resource(root, "chebi", "2")
    manifest["resources"] = {"signor": "1.0.0", "chebi": "2"}
    write_json(path, manifest)
    first = load_release(root, path)
    path.write_text(
        '{"resources":{"chebi":"2","signor":"1.0.0"},"version":"2026.09","schema_version":1}'
    )
    second = load_release(root, path)
    assert first.sha256 == second.sha256
    assert [resource.source for resource in first.resources] == ["chebi", "signor"]


@pytest.mark.parametrize(
    "resources",
    [
        {},
        [],
        {"signor": "latest"},
        {"signor": "working"},
        {"signor": "current"},
        {"signor": "../1"},
        {"signor": "/1"},
        {"signor": 1},
        {"../signor": "1"},
        {"SIGNOR": "1"},
        {"signor.name": "1"},
    ],
)
def test_rejects_mutable_and_invalid_resource_selections(release, resources):
    root, path, _, manifest, _ = release
    manifest["resources"] = resources
    write_json(path, manifest)
    with pytest.raises(ReleaseValidationError):
        load_release(root, path)


@pytest.mark.parametrize(
    "field,value",
    [
        ("version", "latest"),
        ("version", "working"),
        ("schema_version", 2),
        ("schema_version", True),
        ("extra_selection", {}),
        ("created_at", 0),
    ],
)
def test_rejects_unsupported_release_fields(release, field, value):
    root, path, _, manifest, _ = release
    manifest[field] = value
    write_json(path, manifest)
    with pytest.raises(ReleaseValidationError):
        load_release(root, path)


@pytest.mark.parametrize(
    "field,value",
    [
        ("resource_manifests", {}),
        ("resource_manifests", {"signor": "not-a-digest"}),
        ("partial_resources", "yes"),
    ],
)
def test_rejects_malformed_content_pins(release, field, value):
    root, path, _, manifest, _ = release
    manifest[field] = value
    write_json(path, manifest)
    with pytest.raises(ReleaseValidationError):
        load_release(root, path)


def test_content_pin_must_match_build_manifest_bytes(release):
    root, path, directory, manifest, _ = release
    digest = hashlib.sha256((directory / "build_manifest.json").read_bytes()).hexdigest()
    write_json(path, {**manifest, "resource_manifests": {"signor": digest}})
    assert load_release(root, path).resources[0].manifest_sha256 == digest
    write_json(path, {**manifest, "resource_manifests": {"signor": "0" * 64}})
    with pytest.raises(ReleaseValidationError, match="differs from the release pin"):
        load_release(root, path)


@pytest.mark.parametrize("fields", [{"max_records": 20}, {"datasets": ["interactions"]}])
def test_sample_builds_require_an_explicit_partial_release(release, fields):
    root, path, directory, manifest, build_manifest = release
    write_json(directory / "build_manifest.json", {**build_manifest, **fields})
    with pytest.raises(ReleaseValidationError, match="sample build"):
        load_release(root, path)
    write_json(path, {**manifest, "partial_resources": True})
    assert load_release(root, path).resources[0].source == "signor"


@pytest.mark.parametrize(
    "content",
    [
        '{"schema_version":1,"version":"2026.09","resources":{"signor":"1.0.0","signor":"2"}}',
        '{"schema_version":1,"version":"2026.09","resources":{"signor":NaN}}',
        "[]",
        '{"version":',
    ],
)
def test_rejects_ambiguous_or_invalid_json(release, content):
    root, path, *_ = release
    path.write_text(content)
    with pytest.raises(ReleaseValidationError):
        load_release(root, path)


@pytest.mark.parametrize(
    "field,value",
    [
        ("resource", "chebi"),
        ("version", "2"),
        ("schema_version", 2),
        ("schema_version", True),
        ("serving_schema_version", SERVING_SCHEMA_VERSION - 1),
        ("serving_schema_version", None),
    ],
)
def test_rejects_wrong_resource_or_schema_versions(release, field, value):
    root, path, directory, _, manifest = release
    manifest[field] = value
    write_json(directory / "build_manifest.json", manifest)
    with pytest.raises(ReleaseValidationError):
        load_release(root, path)


@pytest.mark.parametrize(
    "field,value,match",
    [
        ("sha256", "0" * 64, "checksum differs"),
        ("sha256", "invalid", "SHA-256"),
        ("size_bytes", 0, "size differs"),
        ("rows", 1, "row count differs"),
        ("rows", True, "nonnegative integer"),
        ("rows", -1, "nonnegative integer"),
        ("path", "../../outside.parquet", "Unknown artifact metadata"),
    ],
)
def test_rejects_unverified_artifact_metadata(release, field, value, match):
    root, path, directory, _, manifest = release
    manifest["files"]["entity.parquet"][field] = value
    write_json(directory / "build_manifest.json", manifest)
    with pytest.raises(ReleaseValidationError, match=match):
        load_release(root, path)


def test_rejects_schema_mismatch_even_with_updated_checksum(release):
    root, path, directory, _, manifest = release
    artifact = directory / "entity.parquet"
    pq.write_table(pa.table({"wrong_column": ["wrong"]}), artifact)
    manifest["files"][artifact.name] = {
        "rows": 1,
        "size_bytes": artifact.stat().st_size,
        "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
    }
    write_json(directory / "build_manifest.json", manifest)
    with pytest.raises(ReleaseValidationError, match="schema does not match"):
        load_release(root, path)


@pytest.mark.parametrize(
    "name",
    ["entity.parquet", "relation.parquet", "evidence_payloads.parquet", "build_manifest.json"],
)
def test_rejects_missing_artifacts(release, name):
    root, path, directory, *_ = release
    (directory / name).unlink()
    with pytest.raises(ReleaseValidationError, match="Missing"):
        load_release(root, path)


def test_rejects_missing_and_extra_file_entries(release):
    root, path, directory, _, manifest = release
    manifest["files"]["../entity.parquet"] = manifest["files"].pop("entity.parquet")
    write_json(directory / "build_manifest.json", manifest)
    with pytest.raises(ReleaseValidationError, match="exactly the resource and serving tables"):
        load_release(root, path)


@pytest.mark.parametrize("outside", [False, True])
def test_rejects_symlinked_artifacts_inside_or_outside_root(release, outside):
    root, path, directory, *_ = release
    artifact = directory / "entity.parquet"
    moved = (root.parent if outside else root) / "linked.parquet"
    artifact.rename(moved)
    artifact.symlink_to(moved)
    with pytest.raises(ReleaseValidationError, match="escapes|symlinks"):
        load_release(root, path)


def test_rejects_mutable_directory_and_manifest_links(release):
    root, path, directory, *_ = release
    moved = directory.with_name("2.0.0")
    directory.rename(moved)
    directory.symlink_to(moved, target_is_directory=True)
    with pytest.raises(ReleaseValidationError, match="symlinks"):
        load_release(root, path)
    alias = path.with_name("latest.json")
    alias.symlink_to(path)
    with pytest.raises(ReleaseValidationError, match="mutable symlink"):
        load_release(root, alias)


def test_validation_reads_metadata_and_hash_chunks_without_materializing_rows(release, monkeypatch):
    root, path, *_ = release

    def fail(*args, **kwargs):
        raise AssertionError("Release validation must not materialize row data")

    monkeypatch.setattr(pq, "read_table", fail)
    monkeypatch.setattr(pq.ParquetFile, "read", fail)
    monkeypatch.setattr(pq.ParquetFile, "iter_batches", fail)
    assert load_release(root, path).resources[0].files["entity.parquet"].rows == 0


def test_prototype_taxonomy_reference_is_pinned_and_verified(release):
    root, path, _, manifest, _ = release
    temporary = root / "taxonomy.parquet"
    pq.write_table(pa.table({"taxon_id": ["9606"], "scientific_name": ["Homo sapiens"]}), temporary)
    checksum = hashlib.sha256(temporary.read_bytes()).hexdigest()
    target = root / "references" / "taxonomy" / checksum / "taxonomy.parquet"
    target.parent.mkdir(parents=True)
    temporary.rename(target)
    manifest["references"] = {
        "taxonomy": {
            "version": checksum,
            "taxon_count": 1,
            "missing_taxon_ids": [],
            "source_url": "https://example.org/taxonomy",
            "source_sha256": "1" * 64,
        }
    }
    write_json(path, manifest)
    result = load_release(root, path)
    assert result.references["taxonomy"].path == target
    assert result.references["taxonomy"].sha256 == checksum
    assert result.manifest["references"]["taxonomy"]["missing_taxon_ids"] == ()
    target.unlink()
    with pytest.raises(ReleaseValidationError, match="Missing"):
        load_release(root, path)


def test_rejects_unknown_reference_fields(release):
    root, path, _, manifest, _ = release
    manifest["references"] = {"unversioned_reference": "latest"}
    write_json(path, manifest)
    with pytest.raises(ReleaseValidationError, match="Unknown release references"):
        load_release(root, path)


def test_recheck_uses_original_pins_and_rejects_changed_release_manifest(release):
    root, path, _, manifest, _ = release
    pinned = load_release(root, path)
    verify_release(pinned)
    write_resource(root, "signor", "2")
    manifest["resources"]["signor"] = "2"
    write_json(path, manifest)
    with pytest.raises(ReleaseValidationError, match="Pinned manifest changed"):
        verify_release(pinned)
    assert pinned.resources[0].version == "1.0.0"


def test_recheck_rejects_changed_resource_manifest(release):
    root, path, directory, _, manifest = release
    pinned = load_release(root, path)
    manifest["new_metadata"] = "modified during load"
    write_json(directory / "build_manifest.json", manifest)
    with pytest.raises(ReleaseValidationError, match="Pinned manifest changed"):
        verify_release(pinned)


def test_recheck_rejects_rewritten_parquet(release):
    root, path, directory, *_ = release
    pinned = load_release(root, path)
    artifact = directory / "entity.parquet"
    schema = FILE_SCHEMAS[artifact.name]
    pq.write_table(pa.Table.from_pylist([{"identifier": "new"}], schema=schema), artifact)
    with pytest.raises(
        ReleaseValidationError, match="row count differs|size differs|checksum differs"
    ):
        verify_release(pinned)
