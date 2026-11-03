"""Raw entity and relation observation dataclasses."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
from omnipath_resolver import RawEntityObservation as _ResolverObservation


@dataclass(slots=True)
class RawEntityObservation(_ResolverObservation):
    """The resolver's observation, plus what the build keys an entity on beyond it."""

    # The compartments an event takes place in (``|``-joined GO ids). Not read by
    # resolution: it splits one resolved reaction into its compartment variants.
    compartments: str | None = None


@dataclass(slots=True)
class RawRelationObservation:
    """Extracted raw relation observation before canonicalization."""

    relation_key: str
    subject_entity_key: str
    predicate: str
    object_entity_key: str
    source: str
    dataset: str
    row_id: str
    upstream_id: str
    statement_kind: str = "relation"
    annotations: list[dict[str, Any]] = field(default_factory=list)
    payload_json: str | None = None
