"""Shared normalization helpers for evidence ingest.

This module is intentionally backend-neutral. It converts pypath silver
entities, identifiers, annotations, and memberships into plain dictionaries and
deterministic content keys that can be consumed by COPY-based ingest or any
future persistence backend. It also centralizes the rules that decide when an
entity tree should produce interaction relation evidence, membership relation
evidence, or ontology-annotation relation evidence.
"""

from __future__ import annotations

import json
import uuid
import hashlib
from dataclasses import dataclass
from collections.abc import Iterable

from omnipath_build.cv_terms import (
    normalize_entity_type,
)
from pypath.internals.cv_terms import (
    IdentifierNamespaceCv,
    cv_term_label_accession,
)
from omnipath_build.relation_rules import (
    ASSOCIATION_CATEGORY,
    ASSOCIATION_PREDICATE,
    CONTROL_PREDICATE,
    INTERACTION_LIKE_TYPES,
    PredicateRule,
    string_or_none,
    predicate_for_membership,
    predicate_for_interaction,
    is_unprojectable_transport,
    order_relation_participants,
)
from pypath.internals.silver_schema import Entity

@dataclass(frozen=True)
class IngestStats:
    """Summary counts from an evidence ingest run."""

    source_rows: int = 0
    entity_evidence: int = 0
    relation_evidence: int = 0
    annotations: int = 0
    identifiers: int = 0


@dataclass
class MutableStats:
    """Mutable accumulator for ingest counts."""

    source_rows: int = 0
    entity_evidence: int = 0
    relation_evidence: int = 0
    annotations: int = 0
    identifiers: int = 0

    def freeze(self) -> IngestStats:
        """Return an immutable public stats snapshot."""

        return IngestStats(
            source_rows=self.source_rows,
            entity_evidence=self.entity_evidence,
            relation_evidence=self.relation_evidence,
            annotations=self.annotations,
            identifiers=self.identifiers,
        )


@dataclass(frozen=True)
class AnnotationRelationSpec:
    """Prepared ontology-annotation relation evidence."""

    relation_occurrence_id: str
    subject_occurrence_id: str
    predicate_rule: PredicateRule
    object_entity_type: str
    object_id_type: str
    object_id: str


@dataclass(frozen=True)
class RelationSpec:
    """Prepared relation evidence before backend-specific persistence."""

    relation_occurrence_id: str
    subject_ref: object
    predicate_rule: PredicateRule
    object_ref: object


def unwrap_record(item: object) -> tuple[object, None]:
    """Return the payload and no external provenance."""

    return item, None


def entity_to_row(entity: Entity) -> dict[str, object]:
    """Convert a silver entity object into the row shape used by ingest."""

    return {
        'type': normalize_entity_type(getattr(entity, 'type', None)),
        'identifiers': [
            identifier_to_row(identifier)
            for identifier in getattr(entity, 'identifiers', None) or []
        ],
        'annotations': annotations_to_rows(
            getattr(entity, 'annotations', None) or []
        ),
        'associations': associations_to_rows(
            getattr(entity, 'associations', None) or []
        ),
        'membership': [
            {
                'member': entity_to_row(membership.member),
                'is_parent': getattr(membership, 'is_parent', None),
                'annotations': annotations_to_rows(
                    getattr(membership, 'annotations', None) or []
                ),
                'associations': associations_to_rows(
                    getattr(membership, 'associations', None) or []
                ),
            }
            for membership in getattr(entity, 'membership', None) or []
            if getattr(membership, 'member', None) is not None
        ],
        'ontology_relations': [
            ontology_relation_to_row(relation)
            for relation in getattr(entity, 'ontology_relations', None) or []
        ],
    }


def ontology_relation_to_row(relation: object) -> dict[str, object]:
    """Convert a silver ontology relation into serializable fields."""

    if isinstance(relation, dict):
        object_ref = relation.get('object') or {}
        return {
            'predicate': text_or_none(relation.get('predicate')),
            'ontology_id': text_or_none(relation.get('ontology_id')),
            'object': entity_ref_to_row(object_ref),
        }
    return {
        'predicate': text_or_none(getattr(relation, 'predicate', None)),
        'ontology_id': text_or_none(getattr(relation, 'ontology_id', None)),
        'object': entity_ref_to_row(getattr(relation, 'object', None)),
    }


def association_to_row(association: object) -> dict[str, object]:
    """Convert a silver association into serializable fields."""

    if isinstance(association, dict):
        object_ref = association.get('object') or {}
        return {
            'predicate': text_or_none(association.get('predicate')),
            'object': entity_ref_to_row(object_ref),
        }
    return {
        'predicate': text_or_none(getattr(association, 'predicate', None)),
        'object': entity_ref_to_row(getattr(association, 'object', None)),
    }


def associations_to_rows(
    associations: Iterable[object],
) -> list[dict[str, object]]:
    """Convert association objects into serializable fields."""

    return [association_to_row(association) for association in associations]


def entity_ref_to_row(ref: object) -> dict[str, str | None]:
    """Convert an entity reference endpoint into serializable fields."""

    if isinstance(ref, dict):
        return {
            'type': text_or_none(ref.get('type')),
            'identifier_type': text_or_none(ref.get('identifier_type')),
            'identifier': text_or_none(ref.get('identifier')),
        }
    return {
        'type': text_or_none(getattr(ref, 'type', None)),
        'identifier_type': text_or_none(getattr(ref, 'identifier_type', None)),
        'identifier': text_or_none(getattr(ref, 'identifier', None)),
    }


def identifier_to_row(identifier: object) -> dict[str, str | None]:
    """Convert an identifier object into serializable fields."""

    return {
        'type': text_or_none(getattr(identifier, 'type', None)),
        'value': text_or_none(getattr(identifier, 'value', None)),
    }


def annotations_to_rows(
    annotations: Iterable[object],
) -> list[dict[str, str | None]]:
    """Convert annotation objects into serializable fields."""

    return [annotation_to_row(annotation) for annotation in annotations]


def interaction_relation_annotations(
    row: dict[str, object],
) -> list[dict[str, str | None]]:
    """Return annotations that should be carried by relation evidence."""

    return annotations_to_rows(row.get('annotations') or [])


def annotation_to_row(annotation: object) -> dict[str, str | None]:
    """Convert one annotation object or dict into serializable fields."""

    if isinstance(annotation, dict):
        return {
            'term': text_or_none(annotation.get('term')),
            'value': text_or_none(annotation.get('value')),
            'unit': text_or_none(
                annotation.get('unit', annotation.get('units'))
            ),
        }
    return {
        'term': text_or_none(getattr(annotation, 'term', None)),
        'value': text_or_none(getattr(annotation, 'value', None)),
        'unit': text_or_none(
            getattr(annotation, 'unit', None)
            or getattr(annotation, 'units', None)
        ),
    }


def annotation_key(
    term: str,
    value: str | None,
    unit: str | None,
) -> str:
    """Return a deterministic UUID key for an annotation value."""

    return content_uuid([term, value, unit])


def identifier_key(vocab_identifier_type: str, value: str) -> str:
    """Return a deterministic UUID key for an identifier evidence value."""

    return content_uuid([vocab_identifier_type, value])


def entity_evidence_key(
    source: str,
    dataset: str,
    row_id: int,
    occurrence_id: str,
) -> str:
    """Return a deterministic UUID key for one entity evidence occurrence."""

    return content_uuid([source, dataset, row_id, occurrence_id])


def relation_evidence_key(
    source: str,
    dataset: str,
    row_id: int,
    relation_occurrence_id: str,
) -> str:
    """Return a deterministic UUID key for one relation evidence occurrence."""

    return content_uuid([source, dataset, row_id, relation_occurrence_id])


def content_uuid(parts: list[object]) -> str:
    """Return a deterministic UUID from a canonical JSON payload."""

    payload = json.dumps(
        parts,
        ensure_ascii=False,
        separators=(',', ':'),
    )
    digest = hashlib.sha256(payload.encode('utf-8')).digest()
    return str(uuid.UUID(bytes=digest[:16]))


def text_or_none(value: object) -> str | None:
    """Normalize enum-like values and blank strings to nullable text."""

    if value is None:
        return None
    text = cv_term_label_accession(value)
    if text:
        return text
    if hasattr(value, 'value'):
        value = value.value
    text = str(value).strip()
    return text or None


def entity_type_accession(vocab_entity_type: str | None) -> str | None:
    """Return the accession part of an entity type."""

    if vocab_entity_type is None:
        return None
    parts = vocab_entity_type.split(':', 2)
    if len(parts) >= 2 and parts[1].isdigit():
        return f'{parts[0]}:{parts[1]}'
    if len(parts) == 3 and parts[2].isdigit():
        return f'{parts[1]}:{parts[2]}'
    return vocab_entity_type


def copy_value(value: object) -> str:
    """Render a Python value for PostgreSQL CSV COPY."""

    if value is None:
        return '\\N'
    if isinstance(value, bool):
        return 'true' if value else 'false'
    if isinstance(value, bytes | bytearray | memoryview):
        return '\\x' + bytes(value).hex()
    return str(value)


TAXONOMY_IDENTIFIER_TERM = cv_term_label_accession(
    IdentifierNamespaceCv.NCBI_TAX_ID
)
STANDARD_INCHI_IDENTIFIER_TERM = cv_term_label_accession(
    IdentifierNamespaceCv.STANDARD_INCHI
)


def include_identifier(ident_type: str | None, ident_value: str | None) -> bool:
    """Return whether an identifier belongs in the generic identifier table."""

    return (
        ident_type is not None
        and ident_value is not None
        and ident_type != STANDARD_INCHI_IDENTIFIER_TERM
    )


def extract_taxonomy_id(row: dict[str, object]) -> str | None:
    """Extract NCBI taxonomy from normalized omnipath_build identifiers/annotations."""

    for ident in row.get('identifiers') or []:
        if not isinstance(ident, dict):
            continue
        if (
            text_or_none(ident.get('type')) == TAXONOMY_IDENTIFIER_TERM
            and text_or_none(ident.get('value'))
        ):
            return _taxonomy_int_text(text_or_none(ident.get('value')))
    for annotation in row.get('annotations') or []:
        if not isinstance(annotation, dict):
            continue
        if (
            text_or_none(annotation.get('term')) == TAXONOMY_IDENTIFIER_TERM
            and text_or_none(annotation.get('value'))
        ):
            return _taxonomy_int_text(text_or_none(annotation.get('value')))

    member_tax_ids: set[str] = set()
    for membership in row.get('membership') or []:
        if not isinstance(membership, dict):
            continue
        member = membership.get('member') or {}
        if not isinstance(member, dict):
            continue
        for annotation in member.get('annotations') or []:
            if not isinstance(annotation, dict):
                continue
            if (
                text_or_none(annotation.get('term')) == TAXONOMY_IDENTIFIER_TERM
                and text_or_none(annotation.get('value'))
            ):
                member_tax_ids.add(
                    _taxonomy_int_text(text_or_none(annotation.get('value')))
                    or ''
                )
    member_tax_ids.discard('')
    if len(member_tax_ids) == 1:
        return next(iter(member_tax_ids))
    return None


def _taxonomy_int_text(value: str | None) -> str | None:
    if value is None:
        return None
    text = value.strip()
    return text if text.isdigit() else None


def association_relation(
    association: dict[str, object],
    *,
    subject_occurrence_id: str,
) -> AnnotationRelationSpec | None:
    """Build relation evidence for an explicit association."""

    object_ref = association.get('object') or {}
    if not isinstance(object_ref, dict):
        return None
    predicate = string_or_none(association.get('predicate')) or ASSOCIATION_PREDICATE
    object_entity_type = string_or_none(object_ref.get('type'))
    object_id_type = string_or_none(object_ref.get('identifier_type'))
    object_id = string_or_none(object_ref.get('identifier'))
    if (
        object_entity_type is None
        or object_id_type is None
        or object_id is None
    ):
        return None

    digest = hashlib.md5(
        f'{predicate}\0{object_entity_type}\0{object_id_type}\0{object_id}'.encode()
    ).hexdigest()
    return AnnotationRelationSpec(
        relation_occurrence_id=f'{subject_occurrence_id}:association:{digest}',
        subject_occurrence_id=subject_occurrence_id,
        predicate_rule=PredicateRule(
            predicate,
            ASSOCIATION_CATEGORY,
        ),
        object_entity_type=object_entity_type,
        object_id_type=object_id_type,
        object_id=object_id,
    )


def is_interaction_like(vocab_entity_type: str | None) -> bool:
    """Return whether an entity type should be handled as an interaction."""

    return entity_type_accession(vocab_entity_type) in INTERACTION_LIKE_TYPES


def interaction_relation_spec(
    row: dict[str, object],
    member_refs: list[tuple[object, object]],
    *,
    occurrence_id: str,
) -> RelationSpec | None:
    """Build an interaction relation spec from two member references."""

    participants = [
        {
            'ref': member_ref,
            'entity_type': text_or_none(
                getattr(getattr(membership, 'member', None), 'type', None)
            ),
            'membership_annotations': annotations_to_rows(
                getattr(membership, 'annotations', None) or []
            ),
        }
        for member_ref, membership in member_refs
    ]
    if is_unprojectable_transport(row, participants):
        return None
    ordered = order_relation_participants(row, participants)
    if len(ordered) != 2:
        return None
    return RelationSpec(
        relation_occurrence_id=f'{occurrence_id}:interaction',
        subject_ref=ordered[0]['ref'],
        predicate_rule=predicate_for_interaction(row, ordered),
        object_ref=ordered[1]['ref'],
    )


def membership_relation_spec(
    *,
    parent_ref: object,
    member_ref: object,
    membership: object,
    parent_type: str | None,
    relation_occurrence_id: str,
) -> RelationSpec:
    """Build a membership relation spec from parent/member references."""

    member_is_parent = bool(getattr(membership, 'is_parent', False))
    semantic_parent_ref = member_ref if member_is_parent else parent_ref
    semantic_child_ref = parent_ref if member_is_parent else member_ref
    semantic_parent_type = parent_type
    if member_is_parent:
        semantic_parent_type = normalize_entity_type(
            getattr(getattr(membership, 'member', None), 'type', None)
        )
    membership_row = {
        'annotations': annotations_to_rows(
            getattr(membership, 'annotations', None) or []
        ),
    }
    predicate_rule = predicate_for_membership(semantic_parent_type, membership_row)
    if predicate_rule.predicate == CONTROL_PREDICATE and not member_is_parent:
        semantic_parent_ref, semantic_child_ref = member_ref, parent_ref
    return RelationSpec(
        relation_occurrence_id=relation_occurrence_id,
        subject_ref=semantic_parent_ref,
        predicate_rule=predicate_rule,
        object_ref=semantic_child_ref,
    )
