"""Rules for projecting silver entity structures into graph predicates.

The ingest phase receives heterogeneous pypath entity records: protein
complexes, pathways, reactions, causal interactions, associations, and
ontology-like annotations all encode relationships in slightly different ways.
This module converts those structures into a small relation vocabulary used by
``relation_evidence`` and later by canonical graph relations.

The rules are conservative. Explicit source/target participant annotations
drive ordering when available. Positive and negative effect annotations drive
regulation predicates. Records without enough semantic detail fall back to
association or generic interaction predicates instead of inventing direction.
"""

from __future__ import annotations

from typing import Any
from dataclasses import dataclass

from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    ExperimentalRoleCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    ParticipantMetadataCv,
    cv_term_label_accession,
)
from pypath.internals.cv_terms.entity_types import EntityTypeCv
from omnipath_build.shared_interaction_schema import (
    SOURCE_ROLE_ACCESSIONS,
    TARGET_ROLE_ACCESSIONS,
    NEGATIVE_SIGN_ACCESSIONS,
    POSITIVE_SIGN_ACCESSIONS,
)

def _term_forms(term: object) -> set[str]:
    label_accession = cv_term_label_accession(term)
    raw = str(term)
    return {raw} if label_accession == raw else {raw, label_accession}


INTERACTION_LIKE_TYPES = {
    str(EntityTypeCv.INTERACTION),
    str(EntityTypeCv.ASSOCIATION),
    str(EntityTypeCv.REACTION),
    str(EntityTypeCv.CATALYSIS),
    str(EntityTypeCv.CONTROL),
    str(EntityTypeCv.DEGRADATION),
    str(EntityTypeCv.TRANSPORT),
}

ASSOCIATION_PREDICATE = 'associated_with'
ASSOCIATION_CATEGORY = 'association'
INTERACTION_CATEGORY = 'interaction'
CONTROL_PREDICATE = 'controls'
TRANSPORT_PREDICATE = 'transports'
TRANSPORTER_ENTITY_TYPES = {
    str(EntityTypeCv.PROTEIN),
    str(EntityTypeCv.GENE),
    str(EntityTypeCv.COMPLEX),
    str(EntityTypeCv.PROTEIN_FAMILY),
}

MEMBERSHIP_RULES: dict[str, str] = {
    str(EntityTypeCv.COMPLEX): 'has_member',
    str(EntityTypeCv.PROTEIN_FAMILY): 'has_member',
    str(EntityTypeCv.PATHWAY): 'has_participant',
    str(EntityTypeCv.REACTION): 'has_participant',
    str(EntityTypeCv.TRANSPORT): 'has_participant',
}

ROLE_TERMS = (
    {str(term) for term in BiologicalRoleCv}
    | {str(term) for term in ExperimentalRoleCv}
    | {str(ParticipantMetadataCv.SOURCE), str(ParticipantMetadataCv.TARGET)}
)

SIGN_POSITIVE_TERMS = {str(term) for term in POSITIVE_SIGN_ACCESSIONS}
SIGN_NEGATIVE_TERMS = {str(term) for term in NEGATIVE_SIGN_ACCESSIONS}
ONTOLOGY_IDENTIFIER_TERM = str(IdentifierNamespaceCv.CV_TERM_ACCESSION)
TAXONOMY_IDENTIFIER_TERM = str(IdentifierNamespaceCv.NCBI_TAX_ID)


@dataclass(frozen=True)
class PredicateRule:
    """Chosen relation predicate and its high-level category."""

    predicate: str
    relation_category: str


def predicate_for_membership(
    parent_type: str | None,
    membership: dict[str, Any],
) -> PredicateRule:
    """Return the relation predicate for parent/member structures."""

    parent_accession = entity_type_accession(parent_type)
    membership_terms = annotation_terms(membership.get('annotations') or [])
    catalytic_terms = (
        _term_forms(BiologicalRoleCv.CATALYST)
        | _term_forms(BiologicalRoleCv.ENZYME)
        | _term_forms(BiologicalRoleCv.CONTROLLER)
    )
    if (
        parent_accession in {str(EntityTypeCv.REACTION), str(EntityTypeCv.TRANSPORT)}
        and membership_terms & catalytic_terms
    ):
        return PredicateRule(CONTROL_PREDICATE, INTERACTION_CATEGORY)

    if parent_accession == str(EntityTypeCv.CV_TERM):
        return PredicateRule(ASSOCIATION_PREDICATE, ASSOCIATION_CATEGORY)

    return PredicateRule(
        MEMBERSHIP_RULES.get(parent_accession or parent_type or '', 'has_member'),
        ASSOCIATION_CATEGORY,
    )


def predicate_for_interaction(
    row: dict[str, Any],
    ordered_participants: list[dict[str, Any]],
) -> PredicateRule:
    """Infer an interaction predicate from entity type, roles, and sign."""

    row_type = string_or_none(row.get('type'))
    row_type_accession = entity_type_accession(row_type)
    annotations = row.get('annotations') or []

    participant_annotations: list[dict[str, Any]] = []
    for participant in ordered_participants:
        participant_annotations.extend(
            participant.get('membership_annotations') or []
        )

    sign = interaction_sign(annotations, participant_annotations)

    if row_type_accession == str(EntityTypeCv.INTERACTION):
        if sign > 0:
            return PredicateRule('positively_regulates', INTERACTION_CATEGORY)
        if sign < 0:
            return PredicateRule('negatively_regulates', INTERACTION_CATEGORY)
        return PredicateRule('interacts_with', INTERACTION_CATEGORY)
    if row_type_accession == str(EntityTypeCv.ASSOCIATION):
        return PredicateRule(ASSOCIATION_PREDICATE, ASSOCIATION_CATEGORY)
    if row_type_accession in {
        str(EntityTypeCv.CONTROL),
        str(EntityTypeCv.CATALYSIS),
        str(EntityTypeCv.DEGRADATION),
    }:
        return PredicateRule(CONTROL_PREDICATE, INTERACTION_CATEGORY)
    if row_type_accession == str(EntityTypeCv.REACTION):
        return PredicateRule('interacts_with', INTERACTION_CATEGORY)
    if row_type_accession == str(EntityTypeCv.TRANSPORT):
        return PredicateRule(TRANSPORT_PREDICATE, INTERACTION_CATEGORY)
    if sign > 0:
        return PredicateRule('positively_regulates', INTERACTION_CATEGORY)
    if sign < 0:
        return PredicateRule('negatively_regulates', INTERACTION_CATEGORY)
    return PredicateRule('related_to', INTERACTION_CATEGORY)


def order_relation_participants(
    row: dict[str, Any],
    participants: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Order participants for relation projection."""

    row_type = entity_type_accession(string_or_none(row.get('type')))
    if row_type == str(EntityTypeCv.TRANSPORT):
        return order_transport_participants(participants)
    return order_interaction_participants(participants)


def count_transporter_participants(participants: list[dict[str, Any]]) -> int:
    """Return the number of transporter-like participants."""

    return sum(
        1
        for participant in participants
        if entity_type_accession(string_or_none(participant.get('entity_type')))
        in TRANSPORTER_ENTITY_TYPES
    )


def is_projectable_transport(
    row: dict[str, Any],
    participants: list[dict[str, Any]],
) -> bool:
    """Return whether a transport record can become a direct relation."""

    row_type = entity_type_accession(string_or_none(row.get('type')))
    return (
        row_type == str(EntityTypeCv.TRANSPORT)
        and len(participants) == 2
        and count_transporter_participants(participants) == 1
    )


def is_unprojectable_transport(
    row: dict[str, Any],
    participants: list[dict[str, Any]],
) -> bool:
    """Return whether a transport direct-relation attempt is invalid."""

    row_type = entity_type_accession(string_or_none(row.get('type')))
    return (
        row_type == str(EntityTypeCv.TRANSPORT)
        and len(participants) == 2
        and not is_projectable_transport(row, participants)
    )


def order_transport_participants(
    participants: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Put the transporter-like participant before the transported substrate."""

    if len(participants) != 2:
        return participants
    transporter_indexes = [
        index
        for index, participant in enumerate(participants)
        if entity_type_accession(string_or_none(participant.get('entity_type')))
        in TRANSPORTER_ENTITY_TYPES
    ]
    if len(transporter_indexes) != 1:
        return order_interaction_participants(participants)
    transporter_index = transporter_indexes[0]
    if transporter_index == 0:
        return participants
    return [participants[1], participants[0]]


def annotation_predicate(annotation: dict[str, Any]) -> str:
    """Return the predicate used for an ontology-valued annotation."""

    return ASSOCIATION_PREDICATE


def interaction_sign(
    record_annotations: list[dict[str, Any]],
    participant_annotations: list[dict[str, Any]] | None = None,
) -> int:
    """Return ``1`` for activating, ``-1`` for inhibitory, or ``0`` unknown."""

    all_annotations = [*record_annotations, *(participant_annotations or [])]
    for annotation in all_annotations:
        term = string_or_none(annotation.get('term'))
        value = (string_or_none(annotation.get('value')) or '').upper()
        if term in SIGN_POSITIVE_TERMS:
            return 1
        if term in SIGN_NEGATIVE_TERMS:
            return -1
        if term == str(InteractionMetadataCv.CONTROL_TYPE):
            if 'ACTIV' in value or 'POSITIVE' in value:
                return 1
            if 'INHIB' in value or 'NEGATIVE' in value or 'REPRESS' in value:
                return -1
    return 0


def order_interaction_participants(
    participants: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Order two participants by source/target role annotations when present."""

    if len(participants) != 2:
        return participants
    first, second = participants
    first_terms = annotation_terms(first.get('membership_annotations') or [])
    second_terms = annotation_terms(second.get('membership_annotations') or [])
    if (
        first_terms & SOURCE_ROLE_ACCESSIONS
        and second_terms & TARGET_ROLE_ACCESSIONS
    ):
        return [first, second]
    if (
        second_terms & SOURCE_ROLE_ACCESSIONS
        and first_terms & TARGET_ROLE_ACCESSIONS
    ):
        return [second, first]
    return participants


def annotation_terms(annotations: list[dict[str, Any]]) -> set[str]:
    """Return non-empty annotation term identifiers."""

    return {
        term
        for annotation in annotations
        if (term := string_or_none(annotation.get('term'))) is not None
    }


def entity_type_accession(entity_type: str | None) -> str | None:
    """Return the accession part of an entity type label."""

    if entity_type is None:
        return None
    parts = entity_type.split(':', 2)
    if len(parts) >= 2 and parts[1].isdigit():
        return f'{parts[0]}:{parts[1]}'
    if len(parts) == 3 and parts[2].isdigit():
        return f'{parts[1]}:{parts[2]}'
    return entity_type


def string_or_none(value: object) -> str | None:
    """Normalize blank values to ``None`` and non-blank values to text."""

    if value is None:
        return None
    text = str(value).strip()
    return text or None
