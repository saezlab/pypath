"""Stream and project inputs_v2 datasets into Silver Entity and Relation observations.

Carries source, dataset, row identity, identifiers, annotations, and complete
evidence payloads as processing context directly for downstream canonicalization.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from .extract.observations import RawEntityObservation, RawRelationObservation
from .merge_policy import preferred_label
from omnipath_resolver.canonical.policy import PROTEIN_ENTITY_TYPES, RNA_ENTITY_TYPES
from typing import Any

from omnipath_core import Relation
from omnipath_core.naming import normalize_namespace
from omnipath_core.measurements import quantity_dict
from omnipath_core.molecular_forms import (
    molecular_form_from_identifiers,
    normalize_molecular_form,
)


from omnipath_core.keys import canonical_json, stable_hash, entity_key, relation_key
from omnipath_core.reactions import COMPARTMENT, compartment_context
from omnipath_core.biolink import (
    entity_type as biolink_entity_type,
    predicate as biolink_predicate,
    annotation_value,
    annotation_term,
)


_IDENTIFIER_PRIORITY: dict[str, int] = {
    "uniprot": 10,
    "chebi": 10,
    "pubchem": 10,
    "ensembl": 10,
    "entrez": 10,
    "hgnc": 10,
    "signor": 10,
    "genesymbol": 20,
    "uniprot_keyword": 30,
    "name": 80,
    "synonym": 90,
}


def _text(value: Any) -> str:
    if value is None:
        return ""
    enum_val = getattr(value, "value", None)
    if enum_val is not None and not isinstance(value, (str, bytes)):
        value = enum_val
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _field(value: Any, name: str, default: Any = None) -> Any:
    if isinstance(value, Mapping):
        return value.get(name, default)
    return getattr(value, name, default)


def _identifier_parts(identifier: Any) -> tuple[str, str] | None:
    raw_ns = _field(identifier, "type")
    val = _text(_field(identifier, "value"))
    if raw_ns is None or not val:
        return None
    ns = normalize_namespace(raw_ns)
    if not ns:
        return None
    return ns, val


def _is_taxon_term(term: object) -> bool:
    return annotation_term(term) == "in_taxon"


def _is_entity_ref(value: Any) -> bool:
    """True for silver ``EntityRef`` (type + identifier_type + identifier)."""
    if value is None or isinstance(value, (str, bytes)):
        return False
    if _field(value, "identifiers"):
        return False
    ident = _field(value, "identifier")
    ident_type = _field(value, "identifier_type")
    return bool(_text(ident)) and ident_type is not None


def _coerce_entity(
    value: Any,
    *,
    default_type: str = "unknown",
    default_ns: str = "unknown",
) -> Any:
    """Normalize EntityRef or a bare identifier string into an Entity-like mapping.

    Mirrors ``add_entity_ref`` in omnipath-build: trust the ref's type and
    identifier_type, with no source-specific accession rewriting.
    """
    if isinstance(value, str):
        ident = value.strip()
        return {
            "type": default_type,
            "identifiers": [{"type": default_ns, "value": ident}],
        }
    if _is_entity_ref(value):
        ident = _text(_field(value, "identifier")).strip()
        ns = normalize_namespace(_field(value, "identifier_type"))
        ent_type = biolink_entity_type(_field(value, "type"))
        return {
            "type": ent_type,
            "identifiers": [{"type": ns, "value": ident}],
            "molecular_form": _field(value, "molecular_form"),
        }
    return value


def _extract_identifiers(entity: Any) -> list[tuple[str, str]]:
    raw_ids = _field(entity, "identifiers", ())
    if not raw_ids and _field(entity, "identifier"):
        raw_ids = [entity]
    extracted: list[tuple[str, str]] = []
    for item in raw_ids or ():
        parts = _identifier_parts(item)
        if parts:
            extracted.append(parts)

    def _sort_key(pair: tuple[str, str]) -> int:
        ns = pair[0]
        return _IDENTIFIER_PRIORITY.get(ns, 50)

    extracted.sort(key=_sort_key)
    return extracted


def _extract_taxon(entity: Any, ids: Sequence[tuple[str, str]]) -> str | None:
    taxon = _field(entity, "taxon")
    if taxon is not None:
        t_str = _text(taxon)
        if t_str and t_str.isdigit():
            return t_str

    for ann in _field(entity, "annotations", ()) or ():
        term = _field(ann, "term")
        if _is_taxon_term(term):
            val = _text(_field(ann, "value"))
            prefix, _, identifier = val.partition(":")
            if prefix != "NCBITaxon" or not identifier.isdigit():
                raise ValueError(f"in_taxon requires an NCBITaxon CURIE, got {val!r}")
            return identifier

    return None


# Annotation text that stands for a missing value.
_NO_VALUE = frozenset({"", "nan", "none", "null"})
# Events whose identity includes the compartments they take place in.
_COMPARTMENTALIZED_TYPES = frozenset({"molecular_activity", "biological_process"})


class SilverExtractor:
    """Extracts Silver observations directly from inputs_v2 record streams."""

    def __init__(self, source: str, dataset: str) -> None:
        self.source = source
        self.dataset = dataset
        self.entities: dict[str, RawEntityObservation] = {}
        self.relations: list[RawRelationObservation] = []
        self.payloads: list[dict[str, Any]] = []
        self.record_count = 0
        self.estimated_bytes = 0

    def _annotations(self, owner: Any, scope: str | None) -> list[dict[str, Any]]:
        """Annotation records of an entity, membership or relation, with ``scope``
        (None for entity annotations). Annotations without a value are left out:
        no quantity and no text, or a placeholder such as ``nan``."""
        records = []
        for ann in _field(owner, "annotations", ()) or ():
            term = annotation_term(_field(ann, "term"))
            value = _field(ann, "value")
            quantity = quantity_dict(value)
            # A quantity's fields are its own columns; its text value stays empty.
            text = "" if quantity is not None else annotation_value(term, value) or ""
            if quantity is None and text.strip().lower() in _NO_VALUE:
                continue
            record = {
                "term": term,
                "value": text,
                "quantity": quantity,
                "source": self.source,
                "dataset": self.dataset,
            }
            if scope is not None:
                record["scope"] = scope
            records.append(record)
        return records

    def _extract_entity(
        self,
        entity: Any,
        *,
        row_id: str,
        path: tuple[int, ...] = (),
        persist_annotations: bool = True,
    ) -> str:
        entity_type = biolink_entity_type(_field(entity, "type"))
        ids = _extract_identifiers(entity)
        molecular_form = normalize_molecular_form(
            _field(entity, "molecular_form"), allow_resolved=False
        )
        if entity_type in PROTEIN_ENTITY_TYPES | RNA_ENTITY_TYPES:
            # Only the principal source identifier establishes participant
            # specificity. A catalogue's other cross-references are not evidence
            # that each referenced transcript or isoform participated.
            primary_form = molecular_form_from_identifiers(
                list(_field(entity, "identifiers", ()) or ())[:1]
            )
            if primary_form:
                # A supplied PTM/variant description does not erase the primary
                # sequence identity. Explicit identity fields take precedence.
                combined = {
                    **primary_form,
                    **{k: v for k, v in (molecular_form or {}).items() if v is not None},
                }
                sequences = []
                for specific in (primary_form, molecular_form or {}):
                    for identifier in specific.get("sequence_identifiers") or []:
                        if identifier not in sequences:
                            sequences.append(identifier)
                combined["sequence_identifiers"] = sequences
                molecular_form = normalize_molecular_form(combined, allow_resolved=False)
        if ids:
            namespace, identifier = ids[0]
        else:
            namespace = "row"
            identifier = f"{self.source}:{self.dataset}:{row_id}:{'.'.join(map(str, path))}"

        taxon = _extract_taxon(entity, ids)
        # Labels are identifiers, but not globally unique accession systems.
        # Keep unresolved label identities local to their resource and taxon.
        identity_scope = self.source if namespace in {"name", "synonym"} else None
        # A reaction resolves to compartment-free Rhea or MetaNetX ids; the
        # compartments it states keep its variants apart.
        compartments = (
            compartment_context(
                annotation_value(COMPARTMENT, _field(ann, "value"))
                for ann in _field(entity, "annotations", ()) or ()
                if annotation_term(_field(ann, "term")) == COMPARTMENT
            )
            if entity_type in _COMPARTMENTALIZED_TYPES
            else ""
        )
        # Keep complete contextual requests apart until resolution. Combining
        # alternative identifiers first can change consensus with batch size.
        primary_key = entity_key(entity_type, namespace, identifier)
        context = (taxon, sorted(set(ids)), identity_scope, molecular_form)
        if compartments:
            context = (*context, compartments)
        key = (
            stable_hash(primary_key, context)
            if taxon or identity_scope or len(set(ids)) > 1 or molecular_form or compartments
            else primary_key
        )

        if key not in self.entities:
            self.entities[key] = RawEntityObservation(
                entity_key=key,
                entity_type=entity_type,
                namespace=namespace,
                identifier=identifier,
                taxon=taxon,
                identity_scope=identity_scope,
                molecular_form=molecular_form,
                compartments=compartments or None,
            )

        ent_record = self.entities[key]
        if ent_record.taxon is None and taxon is not None:
            ent_record.taxon = taxon

        for ord_idx, (ns, val) in enumerate(ids):
            ent_record.identifiers.append(
                {
                    "ns": ns,
                    "id": val,
                    "is_canonical": ord_idx == 0,
                    "source": self.source,
                }
            )
        for ns, val in ids:
            if ns == "name":
                ent_record.label = preferred_label(identifier, ent_record.label, val)
        for ann in _field(entity, "annotations", ()) or ():
            if annotation_term(_field(ann, "term")) == "name":
                ent_record.label = preferred_label(
                    identifier,
                    ent_record.label,
                    annotation_value("name", _field(ann, "value")),
                )

        # Only standalone entity records persist global entity annotations.
        # Participant contextual annotations attach to the relation scope instead.
        if persist_annotations:
            occurrence_annotations = self._annotations(entity, None)
            ent_record.annotations.extend(occurrence_annotations)
            ent_record.evidence.append(
                {
                    "source": self.source,
                    "dataset": self.dataset,
                    "row_id": row_id,
                    "upstream_id": row_id if not path else f"{row_id}:entity:{path}",
                    "annotations": occurrence_annotations,
                    "molecular_form": molecular_form,
                }
            )

        # Extract nested membership if present (e.g. complex members)
        for m_idx, membership in enumerate(_field(entity, "membership", ()) or ()):
            member = _field(membership, "member")
            if member is None:
                continue
            m_path = (*path, m_idx)
            m_key = self._extract_entity(
                member,
                row_id=row_id,
                path=m_path,
                persist_annotations=persist_annotations,
            )
            is_parent = bool(_field(membership, "is_parent", False))
            parent_k, child_k = (m_key, key) if is_parent else (key, m_key)
            m_anns = self._annotations(membership, "relation")
            m_predicate = biolink_predicate(_field(membership, "predicate") or "has_member")
            m_rel_k = relation_key(parent_k, m_predicate, child_k, m_anns)
            self.relations.append(
                RawRelationObservation(
                    relation_key=m_rel_k,
                    subject_entity_key=parent_k,
                    predicate=m_predicate,
                    object_entity_key=child_k,
                    source=self.source,
                    dataset=self.dataset,
                    row_id=row_id,
                    upstream_id=f"{row_id}:member:{m_idx}",
                    annotations=m_anns,
                )
            )
            for a_idx, association in enumerate(_field(membership, "associations", ()) or ()):
                self._extract_association(
                    parent_k,
                    association,
                    row_id=row_id,
                    path=(*path, m_idx, a_idx),
                    upstream_id=f"{row_id}:member:{m_idx}:assoc:{a_idx}",
                )

        # Extract ontology relations if present (e.g. is_a, category, part_of)
        for o_idx, ont_rel in enumerate(_field(entity, "ontology_relations", ()) or ()):
            obj_target = _field(ont_rel, "object")
            pred = biolink_predicate(_field(ont_rel, "predicate"))
            if obj_target is None:
                continue
            o_path = (*path, "ont", o_idx)
            obj_target = _coerce_entity(obj_target)
            o_key = self._extract_entity(
                obj_target,
                row_id=row_id,
                path=o_path,
                persist_annotations=False,
            )
            o_anns = self._annotations(ont_rel, "relation")
            o_rel_k = relation_key(key, pred, o_key, o_anns, statement_kind="ontology")
            self.relations.append(
                RawRelationObservation(
                    relation_key=o_rel_k,
                    statement_kind="ontology",
                    subject_entity_key=key,
                    predicate=pred,
                    object_entity_key=o_key,
                    source=self.source,
                    dataset=self.dataset,
                    row_id=row_id,
                    upstream_id=f"{row_id}:ont:{o_idx}",
                    annotations=o_anns,
                )
            )

        for a_idx, association in enumerate(_field(entity, "associations", ()) or ()):
            self._extract_association(
                key,
                association,
                row_id=row_id,
                path=(*path, 20_000 + a_idx),
                upstream_id=f"{row_id}:assoc:{a_idx}",
            )

        return key

    def _extract_association(
        self,
        subject_key: str,
        association: Any,
        *,
        row_id: str,
        path: tuple[int, ...],
        upstream_id: str,
    ) -> None:
        """Project a silver ``Association`` (protein→CV term, etc.) as a relation."""
        reference = _field(association, "object")
        if reference is None:
            return
        obj_key = self._extract_entity(
            _coerce_entity(reference),
            row_id=row_id,
            path=path,
            persist_annotations=False,
        )
        raw_pred = _field(association, "predicate")
        predicate = biolink_predicate(raw_pred) if raw_pred else "associated_with"
        rel_k = relation_key(subject_key, predicate, obj_key)
        self.relations.append(
            RawRelationObservation(
                relation_key=rel_k,
                subject_entity_key=subject_key,
                predicate=predicate,
                object_entity_key=obj_key,
                source=self.source,
                dataset=self.dataset,
                row_id=row_id,
                upstream_id=upstream_id,
            )
        )

    def process_record(
        self,
        record: Any,
        raw_payload: Any,
        row_id: str,
        row_number: int,
        payload_json: str | None = None,
    ) -> None:
        """Process one top-level entity or relation record."""
        payload_str = payload_json if payload_json is not None else canonical_json(raw_payload)
        relation_start = len(self.relations)
        self.record_count += 1
        # Heuristic text/record estimate, not a hard process-RSS limit.
        # Record/relation caps complement it for expansion-heavy mappers.
        self.estimated_bytes += 1024 + 4 * len(payload_str)

        # If record is a Relation (has subject, predicate, object)
        is_rel = (
            isinstance(record, Relation)
            or (
                hasattr(record, "subject")
                and hasattr(record, "predicate")
                and hasattr(record, "object")
            )
            or (
                isinstance(record, Mapping)
                and "subject" in record
                and "predicate" in record
                and "object" in record
            )
        )

        if is_rel:
            subj_raw = _field(record, "subject")
            pred_raw = _field(record, "predicate")
            obj_raw = _field(record, "object")

            subj_key = self._extract_entity(
                _coerce_entity(subj_raw, default_type="unknown", default_ns="unknown"),
                row_id=row_id,
                path=(0,),
                persist_annotations=False,
            )
            obj_key = self._extract_entity(
                _coerce_entity(obj_raw, default_type="unknown", default_ns="unknown"),
                row_id=row_id,
                path=(1,),
                persist_annotations=False,
            )
            predicate = biolink_predicate(pred_raw)

            # Extract relation annotations
            # Participant annotations attach to relation scope. An event's
            # compartments are its identity, already in its key, and not a
            # claim of this relation.
            rel_anns = [
                *self._annotations(record, "relation"),
                *(
                    annotation
                    for side, raw in (("subject", subj_raw), ("object", obj_raw))
                    for annotation in self._annotations(raw, side)
                    if annotation["term"] != COMPARTMENT
                ),
            ]

            rel_k = relation_key(subj_key, predicate, obj_key, rel_anns)
            record_ids = _extract_identifiers(record)
            upstream_id = record_ids[0][1] if record_ids else row_id
            self.relations.append(
                RawRelationObservation(
                    relation_key=rel_k,
                    subject_entity_key=subj_key,
                    predicate=predicate,
                    object_entity_key=obj_key,
                    source=self.source,
                    dataset=self.dataset,
                    row_id=row_id,
                    upstream_id=upstream_id,
                    annotations=rel_anns,
                )
            )
        else:
            # Standalone entity observation
            root_entity_key = self._extract_entity(record, row_id=row_id, path=())
            if len(self.relations) == relation_start:
                self.payloads.append(
                    {
                        "relation_key": None,
                        "entity_key": root_entity_key,
                        "source": self.source,
                        "row_id": row_id,
                        "payload_json": payload_str,
                    }
                )
        # Include derived membership/ontology/association statements as well.
        for rel_key in dict.fromkeys(r.relation_key for r in self.relations[relation_start:]):
            self.payloads.append(
                {
                    "relation_key": rel_key,
                    "source": self.source,
                    "row_id": row_id,
                    "payload_json": payload_str,
                }
            )
        for relation in self.relations[relation_start:]:
            self.estimated_bytes += 512 + sum(
                256 + 4 * len(str(ann.get("value", ""))) for ann in relation.annotations
            )
