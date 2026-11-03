"""Validate catalog and installation records once while preserving additive server metadata."""

from __future__ import annotations

from copy import deepcopy
from importlib.resources import files
import json
import re
from typing import Any, NotRequired, TypedDict, cast
from urllib.parse import urlsplit

from .errors import ClientError

# What both clients share (``sql/contract.json``, generated from omnipath_core): table
# names, the serving tables' queries, lookups and Biolink names.
CONTRACT: dict[str, Any] = json.loads(
    files("omnipath_client").joinpath("sql", "contract.json").read_text()
)
# The published tables of a resource version (entity_* and relation_* rows refer to
# their entity or relation by ``entity_id`` / ``relation_id``) and the lookup tables
# derived from them at installation.
TABLES = tuple(CONTRACT["published_tables"])
SERVING = tuple(CONTRACT["serving_tables"])


class Artifact(TypedDict):
    name: str
    size_bytes: int
    url: NotRequired[str | None]
    rows: NotRequired[int | None]


class ResourceRecord(TypedDict):
    resource_id: str
    version: str
    files: list[Artifact]


class Installation(TypedDict):
    format: int
    api_url: str
    release: str
    # Whether the database holds native tables (True) or views over the files.
    materialize: bool
    resources: list[ResourceRecord]


def validate_url(value: Any) -> str:
    if not isinstance(value, str):
        raise ClientError("Expected a public HTTP(S) URL without embedded credentials")
    try:
        parsed = urlsplit(value)
        valid = parsed.scheme in {"http", "https"} and parsed.hostname and parsed.username is None
        _ = parsed.port
    except ValueError:
        valid = False
    if not valid:
        raise ClientError("Expected a public HTTP(S) URL without embedded credentials")
    return value


def validate_resource(value: Any) -> ResourceRecord:
    if not isinstance(value, dict):
        raise ClientError("Expected a resource object")
    source, version = value.get("resource_id"), value.get("version")
    if not isinstance(source, str) or not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", source):
        raise ClientError("Invalid resource ID in catalog")
    if not isinstance(version, str) or not re.fullmatch(r"[0-9]+(?:\.[0-9]+)*", version):
        raise ClientError("Invalid resource version in catalog")
    artifacts = value.get("files")
    if not isinstance(artifacts, list):
        raise ClientError("Resource files must be a list")
    names: set[str] = set()
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise ClientError("Expected an artifact object")
        name = artifact.get("name")
        if (
            not isinstance(name, str)
            or name not in {f"{t}.parquet" for t in TABLES}
            or name in names
        ):
            raise ClientError("Unknown or duplicate artifact in catalog")
        size = artifact.get("size_bytes")
        if type(size) is not int or size < 8:
            raise ClientError("Invalid artifact size in catalog")
        if artifact.get("url") is not None:
            validate_url(artifact["url"])
        rows = artifact.get("rows")
        if rows is not None and (type(rows) is not int or rows < 0):
            raise ClientError("Invalid artifact row count")
        names.add(name)
    if not {"entity.parquet", "relation.parquet"} <= names:
        raise ClientError("Resource must contain entity.parquet and relation.parquet")
    return cast(ResourceRecord, deepcopy(value))


def validate_catalog(value: Any) -> dict[str, ResourceRecord]:
    if not isinstance(value, list):
        raise ClientError("Resources must be a list")
    catalog: dict[str, ResourceRecord] = {}
    for raw in value:
        record = validate_resource(raw)
        source = record["resource_id"]
        if source in catalog:
            raise ClientError(f"Duplicate resource in catalog: {source}")
        catalog[source] = record
    return catalog


def validate_installation(value: Any) -> Installation:
    if not isinstance(value, dict):
        raise ClientError("Expected an installation object")
    if type(value.get("format")) is not int or value["format"] != 1:
        raise ClientError("Unsupported installation format")
    api_url = validate_url(value.get("api_url"))
    release = value.get("release")
    if not isinstance(release, str) or not release.strip():
        raise ClientError("Installation release must be a non-empty string")
    materialize = value.get("materialize", False)
    if type(materialize) is not bool:
        raise ClientError("Installation materialize must be true or false")
    catalog = validate_catalog(value.get("resources"))
    return Installation(
        format=1,
        api_url=api_url,
        release=release,
        materialize=materialize,
        resources=list(catalog.values()),
    )
