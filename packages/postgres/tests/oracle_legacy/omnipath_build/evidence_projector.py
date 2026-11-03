"""Shared evidence projection rules for the DuckDB/PostgreSQL load pipeline."""

from __future__ import annotations

from typing import Protocol
from dataclasses import dataclass
from collections.abc import Iterable

from omnipath_build.ingest.common import (
    entity_to_row,
    annotation_key,
    identifier_key,
    annotation_to_row,
    association_to_row,
    include_identifier,
    entity_evidence_key,
    extract_taxonomy_id,
    is_interaction_like,
    association_relation,
    relation_evidence_key,
    membership_relation_spec,
    interaction_relation_spec,
    interaction_relation_annotations,
)
from omnipath_build.relation_rules import (
    ASSOCIATION_CATEGORY,
    CONTROL_PREDICATE,
    string_or_none,
    entity_type_accession,
    is_projectable_transport,
)
from omnipath_build.resolver.chemical_normalization import (
    is_chemical_identifier_type,
    normalize_chemical_identifier,
)
from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    ControlEffectCv,
    EntityTypeCv,
    cv_term_label_accession,
    InteractionMetadataCv,
)


def _term_forms(term: object) -> set[str]:
    label_accession = cv_term_label_accession(term)
    raw = str(term)
    return {raw} if label_accession == raw else {raw, label_accession}


CATALYTIC_CONTROL_ROLE_TERMS = (
    _term_forms(BiologicalRoleCv.CATALYST)
    | _term_forms(BiologicalRoleCv.ENZYME)
    | _term_forms(BiologicalRoleCv.CONTROLLER)
)

@dataclass(frozen=True)
class ProjectionStats:
    """Counts produced while projecting source records into evidence rows."""

    source_rows: int
    entity_evidence: int
    relation_evidence: int
    identifiers: int
    annotations: int


class _RowWriter(Protocol):
    def write(self, row: dict[str, object]) -> None: ...


class _EvidenceWriters(Protocol):
    entity: _RowWriter
    identifier: _RowWriter
    entity_annotation: _RowWriter
    relation_annotation: _RowWriter
    annotation: _RowWriter
    relation: _RowWriter
    annotation_relation: _RowWriter
    ontology_relation: _RowWriter


class EvidenceProjectorBase:
    """Flatten silver entity trees into source-shaped evidence rows."""

    def __init__(
        self,
        *,
        chunk_size: int = 100_000,
    ) -> None:
        self.chunk_size = chunk_size

    def _flatten_entity_tree(
        self,
        entity: object,
        *,
        source: str,
        dataset: str,
        row_id: int,
        occurrence_id: str,
        parent_entity_evidence_id: str | None,
        entity_role: str,
        writers: _EvidenceWriters,
        seen_annotations: set[tuple[str, str, str | None, str | None]],
        stats: _MutableProjectionStats,
    ) -> None:
        row = entity_to_row(entity)
        entity_type = string_or_none(row.get('type'))
        memberships = list(getattr(entity, 'membership', None) or [])
        interaction_member_count = sum(
            1
            for membership in memberships
            if getattr(membership, 'member', None) is not None
        )
        interaction_like = is_interaction_like(entity_type)
        if interaction_like and interaction_member_count < 2:
            return
        participants = [
            {
                'entity_type': string_or_none(
                    getattr(getattr(membership, 'member', None), 'type', None)
                )
            }
            for membership in memberships
            if getattr(membership, 'member', None) is not None
        ]
        entity_type_accession_value = entity_type_accession(entity_type)
        relation_only_interaction = (
            interaction_like
            and interaction_member_count == 2
            and entity_type_accession_value != str(EntityTypeCv.REACTION)
            and (
                entity_type_accession_value != str(EntityTypeCv.TRANSPORT)
                or is_projectable_transport(row, participants)
            )
        )

        entity_evidence_id = None
        if not relation_only_interaction:
            entity_evidence_id = entity_evidence_key(
                source,
                dataset,
                row_id,
                occurrence_id,
            )
            writers.entity.write(
                {
                    'source': source,
                    'dataset': dataset,
                    'row_id': row_id,
                    'entity_evidence_id': entity_evidence_id,
                    'parent_entity_evidence_id': parent_entity_evidence_id,
                    'entity_role': entity_role,
                    'entity_type': entity_type,
                    'taxonomy_id': extract_taxonomy_id(row),
                }
            )
            stats.entity_evidence += 1

            for identifier in row.get('identifiers') or []:
                ident_type = string_or_none(identifier.get('type'))
                ident_value = string_or_none(identifier.get('value'))
                if not include_identifier(ident_type, ident_value):
                    continue
                ident_normalized = normalize_chemical_identifier(
                    ident_type, ident_value,
                )
                # A chemical value that fails its own namespace's syntax --
                # InChIKey=none and similar placeholders -- must never
                # become evidence at all. Then it can never become a
                # canonical identity through any path (spec 011 T043,
                # spec.md US1 acceptance scenario 4).
                if (
                    ident_normalized is None
                    and is_chemical_identifier_type(ident_type)
                ):
                    continue
                writers.identifier.write(
                    {
                        'source': source,
                        'entity_evidence_id': entity_evidence_id,
                        'identifier_id': identifier_key(ident_type, ident_value),
                        'identifier_type': ident_type,
                        'identifier': ident_value,
                        'identifier_normalized': ident_normalized,
                    }
                )
                stats.identifiers += 1

            for annotation in row.get('annotations') or []:
                if _write_annotation(
                    writers.entity_annotation,
                    writers.annotation,
                    seen_annotations,
                    source=source,
                    evidence_id=entity_evidence_id,
                    annotation=annotation,
                ):
                    stats.annotations += 1

            for relation in row.get('ontology_relations') or []:
                object_ref = relation.get('object') or {}
                predicate = string_or_none(relation.get('predicate'))
                object_entity_type = string_or_none(object_ref.get('type'))
                object_identifier_type = string_or_none(
                    object_ref.get('identifier_type')
                )
                object_identifier = string_or_none(object_ref.get('identifier'))
                if (
                    predicate is None
                    or object_entity_type is None
                    or object_identifier_type is None
                    or object_identifier is None
                ):
                    continue
                writers.ontology_relation.write(
                    {
                        'source': source,
                        'dataset': dataset,
                        'subject_entity_evidence_id': entity_evidence_id,
                        'ontology_id': string_or_none(
                            relation.get('ontology_id')
                        ),
                        'subject_entity_type': None,
                        'subject_identifier_type': None,
                        'subject_identifier': None,
                        'predicate': predicate,
                        'object_entity_type': object_entity_type,
                        'object_identifier_type': object_identifier_type,
                        'object_identifier': object_identifier,
                    }
                )
            stats.relation_evidence += self._write_association_relations(
                writers,
                source=source,
                dataset=dataset,
                row_id=row_id,
                occurrence_prefix=occurrence_id,
                associations=row.get('associations') or [],
                subject_occurrence_ids=(occurrence_id,),
            )

        member_refs: list[tuple[str, object]] = []
        for member_index, membership in enumerate(memberships):
            member = getattr(membership, 'member', None)
            if member is None:
                continue
            member_occurrence_id = f'{occurrence_id}:member:{member_index}'
            self._flatten_entity_tree(
                member,
                source=source,
                dataset=dataset,
                row_id=row_id,
                occurrence_id=member_occurrence_id,
                parent_entity_evidence_id=(
                    None if relation_only_interaction else entity_evidence_id
                ),
                entity_role='member',
                writers=writers,
                seen_annotations=seen_annotations,
                stats=stats,
            )
            member_refs.append((member_occurrence_id, membership))

        if relation_only_interaction and len(member_refs) == 2:
            spec = interaction_relation_spec(
                row,
                member_refs,
                occurrence_id=occurrence_id,
            )
            if spec is None:
                return
            relation_evidence_id = relation_evidence_key(
                source,
                dataset,
                row_id,
                spec.relation_occurrence_id,
            )
            self._write_relation(
                writers,
                source=source,
                dataset=dataset,
                row_id=row_id,
                relation_evidence_id=relation_evidence_id,
                subject_occurrence_id=str(spec.subject_ref),
                predicate=spec.predicate_rule.predicate,
                object_occurrence_id=str(spec.object_ref),
                relation_category=(
                    spec.predicate_rule.relation_category or ASSOCIATION_CATEGORY
                ),
            )
            stats.relation_evidence += 1
            relation_annotations = interaction_relation_annotations(row)
            stats.annotations += self._write_relation_annotations(
                writers,
                seen_annotations,
                source=source,
                relation_evidence_id=relation_evidence_id,
                annotations=relation_annotations,
                annotation_scope='relation',
            )
            stats.relation_evidence += self._write_association_relations(
                writers,
                source=source,
                dataset=dataset,
                row_id=row_id,
                occurrence_prefix=f'{spec.relation_occurrence_id}:relation',
                associations=row.get('associations') or [],
                subject_occurrence_ids=(
                    str(spec.subject_ref),
                    str(spec.object_ref),
                ),
            )
            membership_by_occurrence_id = dict(member_refs)
            subject_membership = membership_by_occurrence_id.get(
                str(spec.subject_ref)
            )
            if subject_membership is not None:
                subject_annotations = (
                    getattr(subject_membership, 'annotations', None) or []
                )
                stats.annotations += self._write_relation_annotations(
                    writers,
                    seen_annotations,
                    source=source,
                    relation_evidence_id=relation_evidence_id,
                    annotations=subject_annotations,
                    annotation_scope='subject',
                )
                stats.relation_evidence += self._write_association_relations(
                    writers,
                    source=source,
                    dataset=dataset,
                    row_id=row_id,
                    occurrence_prefix=f'{spec.relation_occurrence_id}:subject',
                    associations=getattr(subject_membership, 'associations', None) or [],
                    subject_occurrence_ids=(str(spec.subject_ref),),
                )
            object_membership = membership_by_occurrence_id.get(
                str(spec.object_ref)
            )
            if object_membership is not None:
                object_annotations = (
                    getattr(object_membership, 'annotations', None) or []
                )
                stats.annotations += self._write_relation_annotations(
                    writers,
                    seen_annotations,
                    source=source,
                    relation_evidence_id=relation_evidence_id,
                    annotations=object_annotations,
                    annotation_scope='object',
                )
                stats.relation_evidence += self._write_association_relations(
                    writers,
                    source=source,
                    dataset=dataset,
                    row_id=row_id,
                    occurrence_prefix=f'{spec.relation_occurrence_id}:object',
                    associations=getattr(object_membership, 'associations', None) or [],
                    subject_occurrence_ids=(str(spec.object_ref),),
                )
        elif member_refs:
            for member_index, (member_occurrence_id, membership) in enumerate(
                member_refs
            ):
                relation_occurrence_id = f'{occurrence_id}:membership:{member_index}'
                spec = membership_relation_spec(
                    parent_ref=occurrence_id,
                    member_ref=member_occurrence_id,
                    membership=membership,
                    parent_type=entity_type,
                    relation_occurrence_id=relation_occurrence_id,
                )
                relation_evidence_id = relation_evidence_key(
                    source,
                    dataset,
                    row_id,
                    spec.relation_occurrence_id,
                )
                self._write_relation(
                    writers,
                    source=source,
                    dataset=dataset,
                    row_id=row_id,
                    relation_evidence_id=relation_evidence_id,
                    subject_occurrence_id=str(spec.subject_ref),
                    predicate=spec.predicate_rule.predicate,
                    object_occurrence_id=str(spec.object_ref),
                    relation_category=(
                        spec.predicate_rule.relation_category
                        or ASSOCIATION_CATEGORY
                    ),
                )
                stats.relation_evidence += 1
                stats.annotations += self._write_relation_annotations(
                    writers,
                    seen_annotations,
                    source=source,
                    relation_evidence_id=relation_evidence_id,
                    annotations=_control_effect_annotations(
                        spec.predicate_rule.predicate,
                        membership,
                    ),
                    annotation_scope='relation',
                )
                member_annotations = getattr(membership, 'annotations', None) or []
                stats.annotations += self._write_relation_annotations(
                    writers,
                    seen_annotations,
                    source=source,
                    relation_evidence_id=relation_evidence_id,
                    annotations=member_annotations,
                    annotation_scope='object',
                )
                stats.relation_evidence += self._write_association_relations(
                    writers,
                    source=source,
                    dataset=dataset,
                    row_id=row_id,
                    occurrence_prefix=f'{relation_occurrence_id}:object',
                    associations=getattr(membership, 'associations', None) or [],
                    subject_occurrence_ids=(member_occurrence_id,),
                )

    @staticmethod
    def _write_relation(
        writers: _EvidenceWriters,
        *,
        source: str,
        dataset: str,
        row_id: int,
        relation_evidence_id: str,
        subject_occurrence_id: str,
        predicate: str,
        object_occurrence_id: str,
        relation_category: str,
    ) -> None:
        writers.relation.write(
            {
                'source': source,
                'dataset': dataset,
                'row_id': row_id,
                'relation_evidence_id': relation_evidence_id,
                'subject_entity_evidence_id': entity_evidence_key(
                    source,
                    dataset,
                    row_id,
                    subject_occurrence_id,
                ),
                'predicate': predicate,
                'object_entity_evidence_id': entity_evidence_key(
                    source,
                    dataset,
                    row_id,
                    object_occurrence_id,
                ),
                'relation_category': relation_category,
            }
        )

    @staticmethod
    def _write_relation_annotations(
        writers: _EvidenceWriters,
        seen_annotations: set[tuple[str, str, str | None, str | None]],
        *,
        source: str,
        relation_evidence_id: str,
        annotations: Iterable[object],
        annotation_scope: str = 'relation',
    ) -> int:
        count = 0
        for annotation in annotations:
            if _write_annotation(
                writers.relation_annotation,
                writers.annotation,
                seen_annotations,
                source=source,
                evidence_id=relation_evidence_id,
                annotation=annotation,
                annotation_scope=annotation_scope,
            ):
                count += 1
        return count

    @staticmethod
    def _write_association_relations(
        writers: _EvidenceWriters,
        *,
        source: str,
        dataset: str,
        row_id: int,
        occurrence_prefix: str,
        associations: Iterable[object],
        subject_occurrence_ids: Iterable[str],
    ) -> int:
        count = 0
        seen_relation_occurrences: set[str] = set()
        for subject_occurrence_id in subject_occurrence_ids:
            for association in associations:
                association_row = association_to_row(association)
                relation_spec = association_relation(
                    association_row,
                    subject_occurrence_id=subject_occurrence_id,
                )
                if relation_spec is None:
                    continue
                if relation_spec.relation_occurrence_id in seen_relation_occurrences:
                    continue
                seen_relation_occurrences.add(relation_spec.relation_occurrence_id)
                writers.annotation_relation.write(
                    {
                        'relation_evidence_id': relation_evidence_key(
                            source,
                            dataset,
                            row_id,
                            (
                                f'{occurrence_prefix}:'
                                f'{relation_spec.relation_occurrence_id}'
                            ),
                        ),
                        'source': source,
                        'dataset': dataset,
                        'row_id': row_id,
                        'subject_entity_evidence_id': entity_evidence_key(
                            source,
                            dataset,
                            row_id,
                            relation_spec.subject_occurrence_id,
                        ),
                        'predicate': relation_spec.predicate_rule.predicate,
                        'object_entity_type': relation_spec.object_entity_type,
                        'object_id_type': relation_spec.object_id_type,
                        'object_id': relation_spec.object_id,
                        'relation_category': (
                            relation_spec.predicate_rule.relation_category
                            or ASSOCIATION_CATEGORY
                        ),
                    }
                )
                count += 1
        return count


def _write_annotation(
    target_writer: _RowWriter,
    value_writer: _RowWriter,
    seen_annotations: set[tuple[str, str, str | None, str | None]],
    *,
    source: str,
    evidence_id: str,
    annotation: object,
    annotation_scope: str | None = None,
) -> bool:
    row = annotation_to_row(annotation)
    term = string_or_none(row.get('term'))
    if term is None:
        return False
    value = string_or_none(row.get('value'))
    unit = string_or_none(row.get('unit', row.get('units')))
    key = annotation_key(term, value, unit)
    target_row: dict[str, object] = {
        'source': source,
        'evidence_id': evidence_id,
        'annotation_key': key,
        'term': term,
        'value': value,
        'unit': unit,
    }
    if annotation_scope is not None:
        target_row['annotation_scope'] = annotation_scope
    target_writer.write(target_row)
    value_tuple = (key, term, value, unit)
    if value_tuple not in seen_annotations:
        seen_annotations.add(value_tuple)
        value_writer.write(
            {
                'annotation_key': key,
                'term': term,
                'value': value,
                'unit': unit,
            }
        )
    return True


def _control_effect_annotations(
    predicate: str,
    membership: object,
) -> list[dict[str, object]]:
    if predicate != CONTROL_PREDICATE:
        return []
    terms = {
        string_or_none(annotation_to_row(annotation).get('term'))
        for annotation in getattr(membership, 'annotations', None) or []
    }
    if terms & CATALYTIC_CONTROL_ROLE_TERMS:
        return [
            {
                'term': InteractionMetadataCv.CONTROL_EFFECT,
                'value': ControlEffectCv.CATALYSIS,
            }
        ]
    return []


@dataclass
class _MutableProjectionStats:
    source_rows: int = 0
    entity_evidence: int = 0
    relation_evidence: int = 0
    identifiers: int = 0
    annotations: int = 0

    def freeze(self) -> ProjectionStats:
        return ProjectionStats(
            source_rows=self.source_rows,
            entity_evidence=self.entity_evidence,
            relation_evidence=self.relation_evidence,
            identifiers=self.identifiers,
            annotations=self.annotations,
        )
