"""Shape raw DuckDB relation rows into API and Svelte domain models."""

from __future__ import annotations

from typing import Any
from omnipath_core.interaction_profiles import interaction_label, relation_qualifiers


def shape_relation_summary(row: dict[str, Any]) -> dict[str, Any]:
    """Convert DuckDB relation row into standard API relation summary."""
    participant_types: list[str] = []
    for val in (row.get("subject_type"), row.get("object_type")):
        if val and val not in participant_types:
            participant_types.append(str(val))

    raw_sources = row.get("sources") or []
    if not isinstance(raw_sources, list):
        raw_sources = [str(raw_sources)] if raw_sources else []
    sources = [str(s) for s in raw_sources if s]

    return {
        "relationPk": str(row.get("relation_key") or ""),
        "subjectEntityPk": str(row.get("subject_entity_key") or ""),
        "subjectReferenceEntityKey": row.get("subject_reference_entity_key"),
        "objectReferenceEntityKey": row.get("object_reference_entity_key"),
        "predicate": str(row.get("predicate") or ""),
        "displayLabel": interaction_label(str(row.get("predicate") or ""), row.get("annotations")),
        "qualifiers": relation_qualifiers(row.get("annotations")),
        "objectEntityPk": str(row.get("object_entity_key") or ""),
        "relationCategory": row.get("category"),
        "interactionClass": row.get("interaction_class"),
        "isDirected": bool(row.get("is_directed", False)),
        "sign": int(row.get("sign") or 0),
        "participantTypes": participant_types,
        "evidenceCount": int(row.get("evidence_count") or 0),
        "sources": sources,
    }
