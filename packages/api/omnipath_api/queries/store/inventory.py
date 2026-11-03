"""Inventory scanning and immutable OmniPath release manifests."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from types import MappingProxyType
import json
import os
import re
import tempfile
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq

from omnipath_core.versioning import (
    RESOURCE_FILES,
    partial_build_reason,
    validate_release_manifest,
)

FILES = RESOURCE_FILES


def safe_name(value: Any) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise ValueError(f"Invalid version or resource name: {value!r}")
    return value


def validate_manifest(data: Any) -> dict[str, Any]:
    validate_release_manifest(
        data,
        resource_name_validator=safe_name,
        resource_version_validator=safe_name,
    )
    return data


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, filename = tempfile.mkstemp(prefix=".publish-", dir=path.parent)
    tmp = Path(filename)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(data, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(tmp, 0o644)
        os.link(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


class ReleaseStore:
    def __init__(self, data_root: str | Path):
        self.root = Path(data_root)

    def get(self, version: str) -> dict[str, Any]:
        path = self.root / "releases" / f"{safe_name(version)}.json"
        if not path.is_file():
            raise FileNotFoundError(f"Unknown OmniPath release: {version}")
        manifest = validate_manifest(json.loads(path.read_text()))
        if manifest["version"] != version:
            raise ValueError(f"Release filename and version disagree: {version}")
        return manifest

    def list(self) -> list[dict[str, Any]]:
        versions = [self.get(path.stem) for path in (self.root / "releases").glob("*.json")]
        return sorted(
            versions, key=lambda item: tuple(map(int, item["version"].split("."))), reverse=True
        )

    def default(self) -> str:
        return "latest"

    def validate_files(self, manifest: dict[str, Any]) -> None:
        taxonomy = manifest.get("references", {}).get("taxonomy")
        if taxonomy:
            path = self.root / "references/taxonomy" / taxonomy["version"] / "taxonomy.parquet"
            if hashlib.sha256(path.read_bytes()).hexdigest() != taxonomy["version"]:
                raise ValueError("Taxonomy reference checksum mismatch")
            pq.read_metadata(path)
        for source, version in manifest["resources"].items():
            directory = self.root / "resources" / source / version
            for name in FILES:
                path = directory / name
                if not path.is_file():
                    raise ValueError(f"Missing release artifact: {source}/{version}/{name}")
                pq.read_metadata(path)

    def resource_manifests(self, manifest: dict[str, Any]) -> dict[str, str] | None:
        """Pin each resource's build manifest bytes; refuse unmarked sample builds.

        Returns None for legacy stores where some resource has no build manifest.
        """
        digests, partial, legacy = {}, [], False
        for source, version in sorted(manifest["resources"].items()):
            path = self.root / "resources" / source / version / "build_manifest.json"
            if not path.is_file():
                legacy = True
                continue
            content = path.read_bytes()
            digests[source] = hashlib.sha256(content).hexdigest()
            reason = partial_build_reason(json.loads(content))
            if reason:
                partial.append(f"{source}/{version} ({reason})")
        if partial and manifest.get("partial_resources") is not True:
            raise ValueError(
                "Release pins sample builds: "
                + ", ".join(partial)
                + '; rebuild them uncapped or set "partial_resources": true'
            )
        pinned = manifest.get("resource_manifests")
        if legacy:
            if pinned is not None:
                raise ValueError("Cannot verify resource_manifests: a build manifest is missing")
            return None
        if pinned is not None and pinned != digests:
            raise ValueError("Release resource_manifests do not match the local build manifests")
        return digests

    def publish(self, manifest: dict[str, Any]) -> dict[str, Any]:
        manifest = dict(validate_manifest(manifest))
        self.validate_files(manifest)
        digests = self.resource_manifests(manifest)
        if digests is not None:
            manifest["resource_manifests"] = digests
        # Legacy stores can still publish without a reference. Once provisioned,
        # every new release resolves its complete taxon set from the cached dump.
        if (
            not manifest.get("references", {}).get("taxonomy")
            and (self.root / "references/taxonomy/taxdump.tar.gz").is_file()
        ):
            from omnipath_core.taxonomy import prepare_reference

            manifest["references"] = {
                **manifest.get("references", {}),
                "taxonomy": prepare_reference(self.root, manifest),
            }
        if not manifest.get("references", {}).get("taxonomy"):
            warnings.warn(
                "Release has no taxonomy reference; provision references/taxonomy/taxdump.tar.gz for name resolution",
                stacklevel=2,
            )
        manifest.setdefault("created_at", datetime.now(timezone.utc).isoformat())
        _write_json(self.root / "releases" / f"{manifest['version']}.json", manifest)
        return manifest

    def references(self, source: str, version: str) -> list[str]:
        return [item["version"] for item in self.list() if item["resources"].get(source) == version]


@dataclass(frozen=True)
class InventorySnapshot:
    """A complete inventory published with one pointer replacement."""

    resources: dict
    latest: dict
    fingerprint: tuple

    def __post_init__(self):
        object.__setattr__(
            self,
            "resources",
            MappingProxyType(
                {key: MappingProxyType(dict(value)) for key, value in self.resources.items()}
            ),
        )
        object.__setattr__(
            self,
            "latest",
            MappingProxyType(
                {key: MappingProxyType(dict(value)) for key, value in self.latest.items()}
            ),
        )
