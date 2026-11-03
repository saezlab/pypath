"""The local data directory: choose resources, download them, derive serving tables."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from copy import deepcopy
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any
from urllib.parse import quote

import duckdb
import httpx
import platformdirs

from .contracts import (
    SERVING,
    Installation,
    ResourceRecord,
    validate_catalog,
    validate_installation,
    validate_resource,
)
from .contracts import validate_url as _url
from .database import DATABASE, in_sync, synchronize
from .errors import ClientError
from .serving import write_serving_tables

# Development deployment; override with ``api_url`` or the OMNIPATH_API_URL variable.
DEFAULT_API_URL = "https://omnipath-metabo-dev.schaul.click/api"
MANIFEST = "installed.json"
EVIDENCE = "evidence_payloads.parquet"
# The catalog's tag dimensions; each is a selection argument of init().
DIMENSIONS = ("topic", "content", "molecule", "knowledge_origin")


def data_directory(path: str | Path | None = None) -> Path:
    """The data directory: ``path``, else OMNIPATH_DATA_DIR, else the per-user default."""
    chosen = path or os.environ.get("OMNIPATH_DATA_DIR") or platformdirs.user_data_dir("omnipath")
    # Normalised: the database records its directory, and "a/../b" is "b".
    return Path(os.path.normpath(Path(chosen).expanduser().absolute()))


def read_installation(root: Path) -> Installation | None:
    """The directory's manifest, or None when nothing was installed there."""
    try:
        text = (root / MANIFEST).read_text()
    except FileNotFoundError:
        return None
    try:
        return validate_installation(json.loads(text))
    except (ValueError, UnicodeError) as exc:
        raise ClientError(f"Invalid {MANIFEST} in {root}") from exc


def _write_installation(root: Path, installation: Installation) -> None:
    fd, staging = tempfile.mkstemp(dir=root, prefix=".installed-")
    with os.fdopen(fd, "w") as stream:
        stream.write(json.dumps(installation, indent=2) + "\n")
    os.replace(staging, root / MANIFEST)


def _get(http: httpx.Client, api_url: str, path: str, **params: Any) -> Any:
    response = http.get(f"{api_url}/{path}", params=params)
    response.raise_for_status()
    try:
        return response.json()
    except ValueError as exc:
        raise ClientError(f"Invalid JSON response from {path}") from exc


def _catalog(http: httpx.Client, api_url: str, release: str) -> dict[str, ResourceRecord]:
    response = _get(http, api_url, "resources", shape="svelte", release=release)
    if not isinstance(response, dict):
        raise ClientError("Expected a catalog object")
    return validate_catalog(response.get("resources"))


def _as_list(value: str | Sequence[str]) -> list[str]:
    return [value] if isinstance(value, str) else list(value)


def _tags(record: ResourceRecord, dimension: str) -> set[str]:
    return {
        tag["id"]
        for tag in record.get("tags") or []
        if isinstance(tag, dict) and tag.get("dimension") == dimension and "id" in tag
    }


def select(
    catalog: dict[str, ResourceRecord],
    resources: str | Sequence[str] | None = None,
    **tags: str | Sequence[str] | None,
) -> list[ResourceRecord]:
    """The catalog's resources that are named (when names are given) and carry one of
    the given tags in every dimension that is given. No criteria selects all."""
    chosen = list(catalog.values())
    if resources is not None:
        names = _as_list(resources)
        if not names:
            raise ValueError("Select at least one resource")
        unknown = set(names) - catalog.keys()
        if unknown:
            raise ValueError(f"Unknown resources in selected release: {sorted(unknown)}")
        chosen = [record for record in chosen if record["resource_id"] in names]
    for dimension, values in tags.items():
        if values is None:
            continue
        wanted = set(_as_list(values))
        known = set().union(*(_tags(record, dimension) for record in catalog.values()))
        if not wanted or wanted - known:
            raise ValueError(
                f"Unknown {dimension}: {sorted(wanted - known)}; available: {sorted(known)}"
            )
        chosen = [record for record in chosen if _tags(record, dimension) & wanted]
    if not chosen:
        raise ValueError("No resource matches the selection")
    return chosen


def available(
    *,
    release: str = "latest",
    api_url: str | None = None,
    timeout: float = 30,
    http_client: httpx.Client | None = None,
) -> list[dict[str, Any]]:
    """The resources of a release that init() can install, with their tags and sizes."""
    api_url = _url(api_url or os.environ.get("OMNIPATH_API_URL") or DEFAULT_API_URL).rstrip("/")
    http = http_client or httpx.Client(timeout=timeout, follow_redirects=True)
    try:
        return deepcopy(list(_catalog(http, api_url, release).values()))
    finally:
        if http_client is None:
            http.close()


def _fetch(http: httpx.Client, artifact: dict[str, Any], path: Path) -> None:
    """Download one file; a wrong size or a non-Parquet file is an error."""
    if not artifact.get("url"):
        raise ClientError("Server does not provide a static download URL")
    size = 0
    with http.stream(
        "GET", _url(artifact["url"]), headers={"Accept-Encoding": "identity"}
    ) as response:
        response.raise_for_status()
        with path.open("wb") as stream:
            for chunk in response.iter_bytes():
                stream.write(chunk)
                size += len(chunk)
    if size != artifact["size_bytes"] or size < 8:
        raise ClientError(f"Truncated download: {artifact['name']}")
    with path.open("rb") as stream:
        header = stream.read(4)
        stream.seek(-4, 2)
        if header != b"PAR1" or stream.read(4) != b"PAR1":
            raise ClientError(f"Invalid Parquet download: {artifact['name']}")


def _present(folder: Path, artifacts: list[dict[str, Any]]) -> bool:
    return all(
        (folder / a["name"]).is_file() and (folder / a["name"]).stat().st_size == a["size_bytes"]
        for a in artifacts
    )


def _evidence_artifact(
    http: httpx.Client, api_url: str, release: str, record: ResourceRecord
) -> dict[str, Any]:
    name = record["resource_id"]
    listed = validate_resource(
        _get(
            http,
            api_url,
            f"resources/{quote(name, safe='')}/files",
            release=release,
            include_evidence="true",
        )
    )
    if listed["resource_id"] != name or listed["version"] != record["version"]:
        raise ClientError(f"{name} changed during init; run init again or use a pinned release")
    evidence = [f for f in listed["files"] if f["name"] == EVIDENCE]
    if not evidence:
        raise ClientError(f"No evidence artifact available for {name}")
    return evidence[0]


def _install(http, db, root: Path, record: ResourceRecord) -> bool:
    """Make one resource version complete on disk: its files, then its serving tables.

    A version directory appears only when complete, and then replaces the resource's
    other versions. Returns whether anything was downloaded.
    """
    downloaded = False
    parent = root / "resources" / record["resource_id"]
    target = parent / record["version"]
    published = [f for f in record["files"] if f["name"] != EVIDENCE]
    evidence = [f for f in record["files"] if f["name"] == EVIDENCE]
    if _present(target, published):
        # A serving table added by a newer serving schema is derived in place; otherwise
        # only raw evidence can be missing from a complete version.
        missing = [name for name in SERVING if not (target / f"{name}.parquet").is_file()]
        if missing:
            work = Path(tempfile.mkdtemp(dir=target, prefix=".work-"))
            try:
                write_serving_tables(db, target, work, tables=missing)
            finally:
                shutil.rmtree(work, ignore_errors=True)
        for artifact in evidence:
            if not _present(target, [artifact]):
                fd, staging = tempfile.mkstemp(dir=target, prefix=".download-")
                os.close(fd)
                try:
                    _fetch(http, artifact, Path(staging))
                    os.replace(staging, target / artifact["name"])
                    downloaded = True
                finally:
                    Path(staging).unlink(missing_ok=True)
    else:
        parent.mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.mkdtemp(dir=parent, prefix=".staging-"))
        try:
            for artifact in record["files"]:
                _fetch(http, artifact, stage / artifact["name"])
            work = stage / ".work"
            work.mkdir()
            write_serving_tables(db, stage, work)
            work.rmdir()
            if target.exists():
                shutil.rmtree(target)
            os.rename(stage, target)
            downloaded = True
        finally:
            shutil.rmtree(stage, ignore_errors=True)
    for other in parent.iterdir():
        if other != target and other.is_dir():
            shutil.rmtree(other)
    return downloaded


def init(
    resources: str | Sequence[str] | None = None,
    *,
    topic: str | Sequence[str] | None = None,
    content: str | Sequence[str] | None = None,
    molecule: str | Sequence[str] | None = None,
    knowledge_origin: str | Sequence[str] | None = None,
    release: str = "latest",
    data_dir: str | Path | None = None,
    include_evidence: bool = False,
    materialize: bool | None = None,
    api_url: str | None = None,
    timeout: float = 30,
    memory_limit: str = "2GB",
    http_client: httpx.Client | None = None,
    progress: Callable[[str], None] | None = None,
) -> Path:
    """Download resources into the data directory and prepare them for querying.

    Selects by resource name and by the catalog's tags; with no selection, every
    resource of the release. Resources already installed stay installed and move to
    this release's versions, so running init again updates the directory. Only missing
    or new versions are downloaded. Each resource becomes available when it is
    complete; an interrupted init keeps the resources it finished. Raw evidence
    payloads are downloaded only with ``include_evidence``, and then kept on updates.

    Ends by writing the directory's DuckDB database (``omnipath.duckdb``): native
    tables copied from the downloaded files (the default for a new directory), or with
    ``materialize=False`` views over the files, which need almost no disk but make the
    prebuilt queries slower the more resources are installed. The choice is kept on
    later runs unless given again. Returns the data directory.
    """
    if not isinstance(release, str) or not release.strip():
        raise ValueError("release must not be empty")
    api_url = _url(api_url or os.environ.get("OMNIPATH_API_URL") or DEFAULT_API_URL).rstrip("/")
    root = data_directory(data_dir)
    say = progress or (lambda _message: None)
    http = http_client or httpx.Client(timeout=timeout, follow_redirects=True)
    db = duckdb.connect()
    try:
        db.execute("SET memory_limit = ?", [memory_limit])
        catalog = _catalog(http, api_url, release)
        chosen = select(
            catalog,
            resources,
            topic=topic,
            content=content,
            molecule=molecule,
            knowledge_origin=knowledge_origin,
        )
        root.mkdir(parents=True, exist_ok=True)
        db.execute("SET temp_directory = ?", [str(root / ".duckdb-tmp")])
        installation = read_installation(root) or Installation(
            format=1, api_url=api_url, release=release, materialize=True, resources=[]
        )
        if materialize is None:
            materialize = installation["materialize"]
        installed = {record["resource_id"]: record for record in installation["resources"]}
        wanted = {record["resource_id"] for record in chosen}
        # Installed resources follow the release too: one directory, one release.
        targets = [r for name, r in catalog.items() if name in wanted or name in installed]
        try:
            for number, entry in enumerate(targets, 1):
                record = deepcopy(entry)
                name = record["resource_id"]
                had_evidence = name in installed and any(
                    f["name"] == EVIDENCE for f in installed[name]["files"]
                )
                record["files"] = [f for f in record["files"] if f["name"] != EVIDENCE]
                if (include_evidence and name in wanted) or had_evidence:
                    record["files"].append(_evidence_artifact(http, api_url, release, record))
                size = sum(f["size_bytes"] for f in record["files"]) / 1e6
                downloaded = _install(http, db, root, record)
                state = f"{size:.0f} MB" if downloaded else "up to date"
                say(f"[{number}/{len(targets)}] {name} {record['version']} ({state})")
                installed[name] = record
                installation = Installation(
                    format=1,
                    api_url=api_url,
                    release=release,
                    materialize=materialize,
                    resources=list(installed.values()),
                )
                _write_installation(root, installation)
        finally:
            # The database describes what is installed, also after a failure midway.
            if installation["resources"]:
                installation["materialize"] = materialize
                _write_installation(root, installation)
                if not in_sync(root, installation):
                    say("Writing the database" + (" (native tables)" if materialize else ""))
                    synchronize(root, installation, memory_limit)
        say(f"{len(installed)} resources installed in {root}")
        say(f"Database: {root / DATABASE}")
    finally:
        db.close()
        shutil.rmtree(root / ".duckdb-tmp", ignore_errors=True)
        if http_client is None:
            http.close()
    return root
