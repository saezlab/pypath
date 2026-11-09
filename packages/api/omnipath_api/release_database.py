"""The served release databases: build one for a published release, prune old ones.

The API serves a release from ``<data_root>/serving/<version>.duckdb``, the release's
DuckDB database (``omnipath_client.database``, the format local clients write too),
with ``serving/<version>.json`` describing how it was built. The newest databases are
kept; the release manifests (``releases/<version>.json``) and the resources' files are
never removed here, so older releases stay downloadable.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Callable

from omnipath_client.database import build_release_database
from omnipath_core.versioning import RESOURCE_FILES, SERVING_FILES

from omnipath_api.queries.release import SERVING_DIR
from omnipath_api.queries.store.inventory import ReleaseStore

DEFAULT_KEEP = 3
# Free disk to leave after writing a database, on top of its size.
DISK_RESERVE_BYTES = 20 * 1024**3
_NUMERIC = re.compile(r"[0-9]+(?:\.[0-9]+)*")


def _version_key(version: str) -> tuple[int, ...]:
    return tuple(map(int, version.split(".")))


def _code_version() -> str | None:
    try:
        return subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip() or None
    except (OSError, subprocess.CalledProcessError):
        return None


def served_databases(data_root: str | Path) -> list[str]:
    """The versions with a database in ``serving/``, newest first."""
    directory = Path(data_root) / SERVING_DIR
    if not directory.is_dir():
        return []
    versions = [p.stem for p in directory.glob("*.duckdb") if _NUMERIC.fullmatch(p.stem)]
    return sorted(versions, key=_version_key, reverse=True)


def build_served_database(
    data_root: str | Path,
    version: str,
    *,
    keep: int | None = DEFAULT_KEEP,
    memory_limit: str = "4GB",
    threads: int | None = None,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    """Write ``serving/<version>.duckdb`` for a published release, then keep only the
    newest ``keep`` databases (all with ``keep=None``). Returns the build record that
    is also written to ``serving/<version>.json``."""
    root = Path(data_root)
    ReleaseStore(root).get(version)  # a published release
    directory = root / SERVING_DIR
    directory.mkdir(parents=True, exist_ok=True)
    previous = [directory / f"{v}.duckdb" for v in served_databases(root) if v != version]
    needed = max((p.stat().st_size for p in previous if p.exists()), default=0) + DISK_RESERVE_BYTES
    free = shutil.disk_usage(directory).free
    if free < needed:
        raise OSError(
            f"Not enough free disk in {directory}: {free / 1e9:.1f} GB free, "
            f"{needed / 1e9:.1f} GB needed (the last database plus the reserve)"
        )
    target = directory / f"{version}.duckdb"
    result = build_release_database(
        root, version, target, memory_limit=memory_limit, threads=threads, progress=progress
    )
    digest = hashlib.sha256()
    with target.open("rb") as stream:
        for chunk in iter(lambda: stream.read(16 * 1024**2), b""):
            digest.update(chunk)
    record = {
        "release": version,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "code_version": _code_version(),
        "size_bytes": result["size_bytes"],
        "sha256": digest.hexdigest(),
        "build_seconds": result["seconds"],
        "resources": result["resources"],
    }
    staging = directory / f".{version}.json"
    staging.write_text(json.dumps(record, indent=2) + "\n")
    os.replace(staging, directory / f"{version}.json")
    record["pruned"] = prune_served_databases(root, keep) if keep is not None else []
    return record


def prune_served_databases(data_root: str | Path, keep: int = DEFAULT_KEEP) -> list[str]:
    """Remove all but the newest ``keep`` release databases; returns the versions removed."""
    if keep < 1:
        raise ValueError("keep must be at least 1")
    directory = Path(data_root) / SERVING_DIR
    removed = []
    for version in served_databases(data_root)[keep:]:
        for suffix in (".duckdb", ".json"):
            (directory / f"{version}{suffix}").unlink(missing_ok=True)
        removed.append(version)
    return removed


def publish_newest(data_root: str | Path, version: str | None = None, **extra) -> str:
    """Publish a release of the newest complete version of every resource under
    ``data_root/resources`` (for development and tests); ``version`` defaults to one
    more than the newest published release. Returns the version."""
    root = Path(data_root)
    store = ReleaseStore(root)
    newest: dict[str, str] = {}
    for directory in sorted((root / "resources").glob("*/*")):
        if not directory.is_dir() or not _NUMERIC.fullmatch(directory.name):
            continue
        if not all(
            (directory / name).is_file()
            for name in RESOURCE_FILES + SERVING_FILES
            if name != "evidence_payloads.parquet"
        ):
            continue
        current = newest.get(directory.parent.name)
        if current is None or _version_key(directory.name) > _version_key(current):
            newest[directory.parent.name] = directory.name
    if version is None:
        published = [m["version"] for m in store.list()]
        version = str(int(published[0].split(".")[0]) + 1) if published else "1"
    store.publish(dict(schema_version=1, version=version, resources=newest, **extra))
    return version


def serve_newest(data_root: str | Path, version: str | None = None, **options) -> str:
    """Publish the newest resources as a release and build its database; returns the
    release version."""
    version = publish_newest(data_root, version)
    build_served_database(data_root, version, **options)
    return version
