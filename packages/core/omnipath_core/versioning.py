"""Structural contracts for immutable resource publications and release selections.

These models describe the existing JSON format. File integrity, storage locations,
publication and inventory policies belong to the callers, not these validators.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any

from .schema import PUBLISHED_TABLES, SERVING_TABLES

BUILD_SCHEMA_VERSION = 1
SERVING_SCHEMA_VERSION = 6
RELEASE_SCHEMA_VERSION = 1
# Published tables are the resource download; serving tables are derived lookup indexes.
RESOURCE_FILES = tuple(f"{table}.parquet" for table in PUBLISHED_TABLES)
SERVING_FILES = tuple(f"{table}.parquet" for table in SERVING_TABLES)


def validate_version(version: str | None) -> str:
    """Require a literal numeric version, never a mutable selection."""
    if not isinstance(version, str) or not re.fullmatch(r"[0-9]+(?:\.[0-9]+)*", version):
        raise ValueError("An explicit numeric resource version is required (for example 1.0.0)")
    return version


def validate_source(source: str) -> str:
    if not isinstance(source, str) or not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", source):
        raise ValueError(f"Invalid resource name: {source!r}")
    return source


def _object(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be an object")
    return value


def _fields(
    value: dict, required: set[str], optional: set[str], context: str, *, strict: bool = True
) -> None:
    missing = required - value.keys()
    unknown = value.keys() - required - optional if strict else set()
    if missing:
        raise ValueError(f"Missing {context} fields: {', '.join(sorted(missing))}")
    if unknown:
        raise ValueError(f"Unknown {context} fields: {', '.join(sorted(unknown))}")


def _integer(value: Any, context: str) -> int:
    if type(value) is not int or value < 0:
        raise ValueError(f"{context} must be a nonnegative integer")
    return value


def _schema_version(value: Any, expected: int, context: str) -> int:
    if type(value) is not int or value != expected:
        raise ValueError(f"{context} must be integer {expected}")
    return value


def _digest(value: Any, context: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"{context} must be a lowercase SHA-256 digest")
    return value


@dataclass(frozen=True)
class ManifestFile:
    """Metadata for one filename-keyed serving Parquet artifact."""

    size_bytes: int
    rows: int
    sha256: str

    def to_dict(self) -> dict[str, Any]:
        return {"rows": self.rows, "sha256": self.sha256, "size_bytes": self.size_bytes}

    @classmethod
    def from_dict(cls, value: Any, *, context: str = "Artifact metadata") -> ManifestFile:
        value = _object(value, context)
        _fields(value, {"rows", "sha256", "size_bytes"}, set(), context)
        return cls(
            size_bytes=_integer(value["size_bytes"], f"{context} size_bytes"),
            rows=_integer(value["rows"], f"{context} rows"),
            sha256=_digest(value["sha256"], f"{context} sha256"),
        )


@dataclass(frozen=True)
class BuildManifest:
    """Current resource publication: schema 1, serving schema 6, one entry per table file.

    Additional publisher metadata is retained verbatim in ``metadata`` and emitted
    at the JSON root. It does not introduce a nested ``metadata`` output field.
    ``created_at`` is optional for legacy reader compatibility.
    """

    resource: str
    version: str
    files: Mapping[str, ManifestFile]
    created_at: str | None = None
    schema_version: int = BUILD_SCHEMA_VERSION
    serving_schema_version: int = SERVING_SCHEMA_VERSION
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        result = dict(self.metadata)
        result.update(
            schema_version=self.schema_version,
            resource=self.resource,
            version=self.version,
            serving_schema_version=self.serving_schema_version,
            files={name: value.to_dict() for name, value in self.files.items()},
        )
        if self.created_at is not None:
            result["created_at"] = self.created_at
        return result

    @classmethod
    def from_dict(cls, value: Any, *, strict_fields: bool = False) -> BuildManifest:
        return validate_build_manifest(value, strict_fields=strict_fields)


_BUILD_FIELDS = {"schema_version", "serving_schema_version", "resource", "version", "files"}
_BUILD_METADATA = {
    "provenance",
    "reference",
    "input_fingerprints",
    "relation_taxon_policy",
    "max_records",
    "datasets",
    "processing_policy",
    "entity_resolution_policy",
    "complex_resolution_policy",
    "batch_limits",
    "build_execution",
    "batch_workers",
    "duckdb_memory_limit",
    "duckdb_threads",
    "payload_origins",
    "phase_metrics",
    "resource_metadata",
}


def validate_build_manifest(value: Any, *, strict_fields: bool = False) -> BuildManifest:
    """Validate JSON structure without opening files or enforcing artifact policy."""
    value = _object(value, "Build manifest")
    _fields(
        value,
        _BUILD_FIELDS,
        _BUILD_METADATA | {"created_at"},
        "Build manifest",
        strict=strict_fields,
    )
    _schema_version(value["schema_version"], BUILD_SCHEMA_VERSION, "Build schema version")
    _schema_version(
        value["serving_schema_version"], SERVING_SCHEMA_VERSION, "Serving schema version"
    )
    resource = validate_source(value["resource"])
    version = validate_version(value["version"])
    if "created_at" in value and not isinstance(value["created_at"], str):
        raise ValueError("Build created_at must be a string")
    entries = _object(value["files"], "Build files")
    if entries.keys() != set(RESOURCE_FILES + SERVING_FILES):
        raise ValueError("Build manifest must describe exactly the resource and serving tables")
    return BuildManifest(
        resource=resource,
        version=version,
        created_at=value.get("created_at"),
        files={
            name: ManifestFile.from_dict(entries[name], context=f"artifact metadata for {name}")
            for name in RESOURCE_FILES + SERVING_FILES
        },
        metadata={
            key: item for key, item in value.items() if key not in _BUILD_FIELDS | {"created_at"}
        },
    )


@dataclass(frozen=True)
class ReleaseManifest:
    """An independent release schedule pinning exact resource versions."""

    version: str
    resources: Mapping[str, str]
    schema_version: int = RELEASE_SCHEMA_VERSION
    created_at: str | None = None
    references: Mapping[str, Any] | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        result = dict(self.metadata)
        result.update(
            schema_version=self.schema_version, version=self.version, resources=dict(self.resources)
        )
        if self.created_at is not None:
            result["created_at"] = self.created_at
        if self.references is not None:
            result["references"] = dict(self.references)
        return result

    @classmethod
    def from_dict(cls, value: Any, **policies: Any) -> ReleaseManifest:
        return validate_release_manifest(value, **policies)


def validate_release_manifest(
    value: Any,
    *,
    strict_fields: bool = False,
    require_resources: bool = True,
    resource_name_validator: Callable[[Any], str] = validate_source,
    resource_version_validator: Callable[[Any], str] = validate_version,
) -> ReleaseManifest:
    """Validate selections, allowing caller-specific inventory naming policies.

    The release version is always numeric. By default, resource pins use the
    publisher's names and numeric versions. An inventory supporting historical
    names can explicitly provide its existing validators. No files are read.
    """
    value = _object(value, "Release manifest")
    required = {"schema_version", "version", "resources"}
    _fields(
        value,
        required,
        {"created_at", "references", "resource_manifests", "partial_resources"},
        "Release",
        strict=strict_fields,
    )
    _schema_version(value["schema_version"], RELEASE_SCHEMA_VERSION, "Release schema version")
    version = validate_version(value["version"])
    resources = _object(value["resources"], "Release resources")
    if require_resources and not resources:
        raise ValueError("Release must pin at least one resource version")
    for source, resource_version in resources.items():
        resource_name_validator(source)
        resource_version_validator(resource_version)
    if "created_at" in value and not isinstance(value["created_at"], str):
        raise ValueError("Release created_at must be a string")
    if "resource_manifests" in value:
        # Content pins: SHA-256 of each pinned resource's build_manifest.json bytes.
        digests = _object(value["resource_manifests"], "Release resource_manifests")
        if digests.keys() != resources.keys():
            raise ValueError("Release resource_manifests must pin exactly the selected resources")
        for source, digest in digests.items():
            _digest(digest, f"Build manifest checksum for {source}")
    if "partial_resources" in value and type(value["partial_resources"]) is not bool:
        raise ValueError("Release partial_resources must be a boolean")
    references = value.get("references")
    if "references" in value:
        references = _object(references, "Release references")
        _fields(references, set(), {"taxonomy"}, "release references", strict=strict_fields)
        if "taxonomy" in references:
            taxonomy = _object(references["taxonomy"], "Taxonomy reference")
            _fields(
                taxonomy,
                {"version"},
                {"source_url", "source_sha256", "taxon_count", "missing_taxon_ids"},
                "Taxonomy reference",
                strict=strict_fields,
            )
            _digest(taxonomy["version"], "Taxonomy version")
            if "source_sha256" in taxonomy:
                _digest(taxonomy["source_sha256"], "Taxonomy source checksum")
            if "taxon_count" in taxonomy:
                _integer(taxonomy["taxon_count"], "Taxonomy row count")
            if "source_url" in taxonomy and not isinstance(taxonomy["source_url"], str):
                raise ValueError("Taxonomy source_url must be a string")
            if "missing_taxon_ids" in taxonomy and (
                not isinstance(taxonomy["missing_taxon_ids"], list)
                or any(
                    not isinstance(taxon, str) or not taxon.isdigit()
                    for taxon in taxonomy["missing_taxon_ids"]
                )
            ):
                raise ValueError("Taxonomy missing_taxon_ids must contain numeric strings")
    return ReleaseManifest(
        version=version,
        resources=dict(resources),
        created_at=value.get("created_at"),
        references=references,
        metadata={
            key: item
            for key, item in value.items()
            if key not in required | {"created_at", "references"}
        },
    )


def partial_build_reason(manifest: Mapping[str, Any]) -> str | None:
    """Explain why a build manifest describes a sample rather than a complete resource."""
    if manifest.get("max_records") is not None:
        return f"capped at max_records={manifest['max_records']}"
    if manifest.get("datasets") is not None:
        return f"restricted to datasets {list(manifest['datasets'])}"
    return None


def resource_dir(data_root: str | Path, source: str, version: str) -> Path:
    return (
        Path(data_root).resolve()
        / "resources"
        / validate_source(source)
        / validate_version(version)
    )


def manifest_path(data_root: str | Path, source: str, version: str) -> Path:
    return resource_dir(data_root, source, version) / "build_manifest.json"


def table_path(data_root: str | Path, source: str, version: str, table: str) -> Path:
    if table not in PUBLISHED_TABLES and table not in SERVING_TABLES:
        raise ValueError(f"Unknown resource table: {table}")
    return resource_dir(data_root, source, version) / f"{table}.parquet"
