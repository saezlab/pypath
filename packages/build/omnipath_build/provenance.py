"""Content fingerprints for reproducible published artifacts and checkpoints."""

from __future__ import annotations

import hashlib
import json
import importlib.metadata
import platform
from pathlib import Path


def file_fingerprint(path: str | Path) -> dict[str, str | int]:
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"sha256": digest.hexdigest(), "size_bytes": path.stat().st_size}


def runtime_provenance() -> dict:
    packages = (
        "omnipath-build",
        "omnipath-core",
        "omnipath-resolver",
        "pypath-omnipath",
        "biolink-model",
        "duckdb",
        "pyarrow",
        "rdkit",
        "pygoslin",
    )
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    from .runtime_files import software_provenance

    software = software_provenance()
    builder = software["packages"]["omnipath-build"]
    return {
        "python": platform.python_version(),
        "dependencies": versions,
        "software": software,
        # Retain existing keys for readers of earlier publication metadata.
        "code_revision": builder["git"]["commit"] if builder["git"] else None,
        "build_code_sha256": builder["content_sha256"],
    }


def reference_provenance(library_dir: str | Path | None) -> dict | None:
    if library_dir is None:
        return None
    root = Path(library_dir)
    manifest = root / "manifest.json"
    if manifest.is_file():
        value = json.loads(manifest.read_text())
        return {
            "path": str(root),
            "fingerprint": value.get("fingerprint"),
            "format": value.get("format"),
            "manifest": file_fingerprint(manifest),
        }
    return {
        "path": str(root),
        "manifest": None,
        "files": {
            str(p.relative_to(root)): file_fingerprint(p)
            for p in sorted(root.rglob("*"))
            if p.is_file() and p.suffix in (".parquet", ".json")
        },
    }
