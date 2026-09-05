"""Silver schema re-exports from omnipath_core.

Canonical definitions live in omnipath_core.silver_schema.
"""
from __future__ import annotations

from omnipath_core.silver_schema import (
    ANNOTATION_FIELDS,
    ASSOCIATION_FIELDS,
    BASE_ENTITY_FIELDS,
    BASE_MEMBERSHIP_FIELDS,
    ENTITY_FIELDS,
    ENTITY_REF_FIELDS,
    ENTITY_SCHEMA,
    IDENTIFIER_FIELDS,
    MEMBERSHIP_FIELDS,
    NESTED_ENTITY_FIELDS,
    ONTOLOGY_RELATION_FIELDS,
    RELATION_FIELDS,
    RELATION_SCHEMA,
    Annotation,
    Association,
    Entity,
    EntityRef,
    Identifier,
    Membership,
    OntologyRelation,
    Relation,
    format_term,
)

__all__ = [
    "Identifier",
    "Annotation",
    "Association",
    "Membership",
    "EntityRef",
    "OntologyRelation",
    "Entity",
    "Relation",
    "format_term",
    "ENTITY_SCHEMA",
    "RELATION_SCHEMA",
    "IDENTIFIER_FIELDS",
    "ANNOTATION_FIELDS",
    "ENTITY_REF_FIELDS",
    "ONTOLOGY_RELATION_FIELDS",
    "ASSOCIATION_FIELDS",
    "BASE_ENTITY_FIELDS",
    "BASE_MEMBERSHIP_FIELDS",
    "NESTED_ENTITY_FIELDS",
    "MEMBERSHIP_FIELDS",
    "ENTITY_FIELDS",
]
