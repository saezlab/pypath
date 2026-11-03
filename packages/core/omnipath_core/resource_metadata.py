"""Source metadata usable by the API without importing pypath input modules."""

from enum import StrEnum
from functools import lru_cache
import json
from pathlib import Path


@lru_cache(maxsize=1)
def _catalog():
    return json.loads((Path(__file__).parent / "vocab/resource_metadata.yaml").read_text())


@lru_cache(maxsize=1)
def _tag_vocabulary():
    return json.loads((Path(__file__).parent / "vocab/resource_tags.yaml").read_text())


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
        "description": current.get("description", ""),
        "tags": [{"id": tag, **vocabulary[tag]} for tag in current.get("tags", [])],
        "license_use": dict(license_use),
        "metadata_source": "build snapshot"
        if isinstance(snapshot, dict)
        else "pypath input configuration",
        "downloaded_at": metadata.get("downloaded_at"),
    }
