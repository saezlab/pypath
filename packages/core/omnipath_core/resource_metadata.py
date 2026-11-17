"""Source metadata usable by the API without importing pypath input modules."""

from enum import StrEnum
from functools import lru_cache
import json
from pathlib import Path
import re


@lru_cache(maxsize=1)
def _catalog():
    return json.loads((Path(__file__).parent / "vocab/resource_metadata.yaml").read_text())


@lru_cache(maxsize=1)
def _tag_vocabulary():
    return json.loads((Path(__file__).parent / "vocab/resource_tags.yaml").read_text())


def _slug(name: str) -> str:
    """Lowercase letters and digits only, as a resource's slug is spelled."""
    return re.sub(r"[^a-z0-9]", "", (name or "").lower())


@lru_cache(maxsize=1)
def _name_index() -> dict[str, str]:
    """{slug, slugified short name, slugified synonym} -> the resource's slug."""
    index: dict[str, str] = {}
    for slug, entry in _catalog().items():
        index[slug] = slug
    for slug, entry in _catalog().items():
        index.setdefault(_slug(entry.get("short", "")), slug)
        for synonym in entry.get("synonyms") or ():
            index.setdefault(_slug(synonym), slug)
    return index


def resolve_resource_name(name: str) -> str | None:
    """The slug of the resource a name stands for: its slug, its short name or a synonym,
    in any case and with or without punctuation; ``None`` when no resource has the name."""
    key = _slug(name)
    return _name_index().get(key) if key else None


def _dimension(name, dimension):
    return StrEnum(
        name,
        {
            tag.upper(): tag
            for tag, entry in _tag_vocabulary().items()
            if entry["dimension"] == dimension
        },
    )


# The controlled tags a resource may carry, one enum per dimension the resources page
# groups them by, e.g. ``Content.PATHWAYS``. Members come from ``vocab/resource_tags.yaml``,
# which also holds each tag's label and description.
Topic = _dimension("Topic", "topic")
Content = _dimension("Content", "content")
KnowledgeOrigin = _dimension("KnowledgeOrigin", "knowledge_origin")


def resource_metadata(source: str, manifest: dict) -> dict:
    # Legacy builds have no source snapshot. Never substitute build time for
    # download time: cached inputs can predate a build by months.
    snapshot = manifest.get("resource_metadata")
    current = _catalog().get(source, {})
    metadata = dict(snapshot if isinstance(snapshot, dict) else current)
    vocabulary = _tag_vocabulary()
    # The input module's current name, description and tags apply to existing builds,
    # but never overwrite their acquisition provenance. License use is that of the
    # build's own license: a license no input declares stays unknown.
    license_use = next(
        (
            entry["license_use"]
            for entry in _catalog().values()
            if entry["license"] == metadata.get("license")
        ),
        {"academic": "unknown", "commercial": "unknown"},
    )
    return {
        **metadata,
        "display_name": current.get("display_name") or source,
        # The names are the input module's current ones, whatever the build's snapshot says.
        "short": current.get("short") or source,
        "full_name": current.get("full_name") or current.get("display_name") or source,
        "synonyms": list(current.get("synonyms") or ()),
        "description": current.get("description", ""),
        "tags": [{"id": tag, **vocabulary[tag]} for tag in current.get("tags", [])],
        "license_use": dict(license_use),
        "metadata_source": "build snapshot"
        if isinstance(snapshot, dict)
        else "pypath input configuration",
        "downloaded_at": metadata.get("downloaded_at"),
    }
