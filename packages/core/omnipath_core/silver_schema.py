"""Canonical definitions for silver-layer schemas.

Keeping the PyArrow schema objects and the namedtuple helpers in one place
avoids accidental divergence across the pipeline.
"""

from __future__ import annotations

from typing import Any, List, NamedTuple
from biolink_model.datamodel.model import QuantityValue
from .measurements import Measurement
from .molecular_forms import MolecularForm, MOLECULAR_FORM_STRUCT
import pyarrow as pa

__all__ = [
    "Identifier",
    "Annotation",
    "Association",
    "Membership",
    "EntityRef",
    "OntologyRelation",
    "Entity",
    "Relation",
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
    "ENTITY_SCHEMA",
    "RELATION_FIELDS",
    "RELATION_SCHEMA",
    "format_term",
]


def format_term(val: Any) -> str:
    """Serialize native Biolink primitives without losing enum values or CURIEs."""
    from enum import Enum
    from linkml_runtime.linkml_model.meta import PermissibleValue
    from linkml_runtime.utils.enumerations import EnumDefinitionImpl
    from .biolink import entity_type, slot_name, enum_value

    if val is None:
        return ""
    if isinstance(val, (EnumDefinitionImpl, PermissibleValue)):
        return enum_value(val)
    if isinstance(val, type) and hasattr(val, "class_class_curie"):
        return entity_type(val)
    if hasattr(val, "curie"):
        return slot_name(val)
    if isinstance(val, Enum):
        return str(val.value)
    return str(val)


_format_term = format_term


class Identifier(NamedTuple):
    """An identifier with its namespace/type and value."""

    type: str
    value: str

    def __repr__(self) -> str:
        """Compact representation showing namespace and value."""
        type_name = _format_term(self.type)
        return f"{type_name}:{self.value}"


class Annotation(NamedTuple):
    """Annotation for an interaction or entity."""

    term: str
    value: str | float | QuantityValue | Measurement | None = None
    units: str | None = None

    def __repr__(self) -> str:
        """Compact representation of annotation."""
        term_str = _format_term(self.term)
        if self.value is None:
            return f"{term_str}"
        elif self.units:
            return f"{term_str}={self.value}{self.units}"
        else:
            return f"{term_str}={self.value}"


class EntityRef(NamedTuple):
    """Reference to an entity endpoint by one identifier."""

    type: str
    identifier_type: str
    identifier: str
    molecular_form: MolecularForm | dict | None = None

    def __repr__(self) -> str:
        type_name = _format_term(self.type)
        id_type = _format_term(self.identifier_type)
        return f"{type_name}({id_type}:{self.identifier})"


class Association(NamedTuple):
    """Association from the enclosing entity to another first-class entity."""

    object: EntityRef
    predicate: str | None = None


class Membership(NamedTuple):
    """Membership relationship representing entity hierarchies / components."""

    member: "Entity"
    is_parent: bool = False
    annotations: list[Annotation] | None = None
    associations: list[Association] | None = None

    predicate: str | None = None

    def __repr__(self) -> str:
        member_repr = repr(self.member)
        if self.annotations:
            annot_str = f" [membership: {len(self.annotations)} annot.]"
        else:
            annot_str = ""
        return f"Membership(entity={member_repr}{annot_str})"

    def pretty(self, indent: int = 0) -> str:
        prefix = "  " * indent
        lines = [f"{prefix}Membership:"]
        lines.append(f"{prefix}├─ entity:")
        entity_lines = self.member.pretty(0).split("\n")
        for line in entity_lines:
            lines.append(f"{prefix}│  {line}")
        if self.annotations:
            lines.append(f"{prefix}└─ membership_annotations:")
            for i, ann in enumerate(self.annotations):
                connector = "└─" if i == len(self.annotations) - 1 else "├─"
                lines.append(f"{prefix}   {connector} {ann!r}")
        else:
            lines[1] = f"{prefix}└─ entity:"
            lines = [lines[0]] + [lines[1]] + [f"{prefix}   {line}" for line in entity_lines]
        return "\n".join(lines)


class OntologyRelation(NamedTuple):
    """Ontology edge from enclosing entity to an object entity.

    ``annotations`` describe the axiom, e.g. the ontology's own predicate when
    the Biolink predicate is broader than it.
    """

    predicate: str
    object: EntityRef
    ontology_id: str | None = None
    annotations: list[Annotation] | None = None


class Entity(NamedTuple):
    """Entity record matching silver entities schema."""

    type: str
    identifiers: List[Identifier]
    annotations: List[Annotation] | None = None
    membership: List[Membership] | None = None
    ontology_relations: List[OntologyRelation] | None = None
    associations: List[Association] | None = None
    molecular_form: MolecularForm | dict | None = None

    def __repr__(self) -> str:
        return self.pretty()

    def __str__(self) -> str:
        return self.pretty()

    def pretty(self, indent: int = 0) -> str:
        prefix = "  " * indent
        type_name = _format_term(self.type)
        lines = [f"{prefix}{type_name}"]

        # Identifiers
        lines.append(f"{prefix}├─ identifiers:")
        for i, id_ in enumerate(self.identifiers):
            connector = (
                "└─"
                if i == len(self.identifiers) - 1
                and not self.annotations
                and not self.associations
                and not self.membership
                and not self.ontology_relations
                else "├─"
            )
            lines.append(f"{prefix}│  {connector} {id_!r}")

        # Annotations
        if self.annotations:
            is_last = not self.associations and not self.membership and not self.ontology_relations
            connector = "└─" if is_last else "├─"
            lines.append(f"{prefix}{connector} annotations:")
            for i, ann in enumerate(self.annotations):
                ann_connector = "└─" if i == len(self.annotations) - 1 else "├─"
                ann_prefix = "   " if is_last else "│  "
                lines.append(f"{prefix}{ann_prefix}{ann_connector} {ann!r}")

        if self.associations:
            is_last = not self.membership and not self.ontology_relations
            connector = "└─" if is_last else "├─"
            lines.append(f"{prefix}{connector} associations: ({len(self.associations)})")

        # Membership
        if self.membership:
            is_last = not self.ontology_relations
            connector = "└─" if is_last else "├─"
            lines.append(f"{prefix}{connector} membership: ({len(self.membership)})")
            for i, member in enumerate(self.membership):
                member_connector = "└─" if i == len(self.membership) - 1 else "├─"
                member_prefix = "   " if is_last else "│  "

                lines.append(f"{prefix}{member_prefix}{member_connector} Membership:")
                continuation = "   " if i == len(self.membership) - 1 else "│  "
                member_lines = member.member.pretty(0).split("\n")
                lines.append(f"{prefix}{member_prefix}{continuation}  ├─ entity: {member_lines[0]}")
                for line in member_lines[1:]:
                    lines.append(f"{prefix}{member_prefix}{continuation}  │        {line}")

                if member.annotations:
                    lines.append(
                        f"{prefix}{member_prefix}{continuation}  └─ membership_annotations:"
                    )
                    for j, ann in enumerate(member.annotations):
                        ann_connector = "└─" if j == len(member.annotations) - 1 else "├─"
                        lines.append(
                            f"{prefix}{member_prefix}{continuation}     {ann_connector} {ann!r}"
                        )
                else:
                    if lines[-len(member_lines) - 1].endswith("├─ entity: " + member_lines[0]):
                        lines[-len(member_lines) - 1] = (
                            f"{prefix}{member_prefix}{continuation}  └─ entity: {member_lines[0]}"
                        )

        if self.ontology_relations:
            lines.append(f"{prefix}└─ ontology_relations: ({len(self.ontology_relations)})")

        return "\n".join(lines)


class Relation(NamedTuple):
    """Relation record representing a directed edge/interaction between entities."""

    subject: EntityRef | Entity
    predicate: str
    object: EntityRef | Entity
    category: str | None = None
    interaction_class: str | None = None
    identifiers: List[Identifier] | None = None
    annotations: List[Annotation] | None = None

    def __repr__(self) -> str:
        sub_str = repr(self.subject)
        obj_str = repr(self.object)
        pred_str = _format_term(self.predicate)
        extra = []
        if self.category:
            extra.append(f"category={self.category}")
        if self.interaction_class:
            extra.append(f"interaction_class={self.interaction_class}")
        if self.identifiers:
            extra.append(f"identifiers={self.identifiers!r}")
        if self.annotations:
            extra.append(f"annotations={self.annotations!r}")
        extra_str = f" ({', '.join(extra)})" if extra else ""
        return f"Relation({sub_str} -[{pred_str}]-> {obj_str}{extra_str})"

    def __str__(self) -> str:
        return repr(self)

    def pretty(self, indent: int = 0) -> str:
        prefix = "  " * indent
        pred_str = _format_term(self.predicate)
        lines = [f"{prefix}Relation: {pred_str}"]
        lines.append(f"{prefix}├─ subject:")
        sub_lines = (
            self.subject.pretty(0) if hasattr(self.subject, "pretty") else repr(self.subject)
        ).split("\n")
        for line in sub_lines:
            lines.append(f"{prefix}│  {line}")
        lines.append(f"{prefix}├─ object:")
        obj_lines = (
            self.object.pretty(0) if hasattr(self.object, "pretty") else repr(self.object)
        ).split("\n")
        for line in obj_lines:
            lines.append(f"{prefix}│  {line}")
        if self.category:
            lines.append(f"{prefix}├─ category: {self.category}")
        if self.interaction_class:
            lines.append(f"{prefix}├─ interaction_class: {self.interaction_class}")
        if self.identifiers:
            lines.append(f"{prefix}├─ identifiers: ({len(self.identifiers)})")
            for id_ in self.identifiers:
                lines.append(f"{prefix}│  ├─ {id_!r}")
        if self.annotations:
            lines.append(f"{prefix}└─ annotations: ({len(self.annotations)})")
            for ann in self.annotations:
                lines.append(f"{prefix}   └─ {ann!r}")
        return "\n".join(lines)


#### PyArrow Schemas (needed for silver parquet files) ####

IDENTIFIER_FIELDS = [
    pa.field("type", pa.string()),
    pa.field("value", pa.string()),
]

ANNOTATION_FIELDS = [
    pa.field("term", pa.string()),
    pa.field("value", pa.string()),
    pa.field("units", pa.string()),
]

ENTITY_REF_FIELDS = [
    pa.field("type", pa.string()),
    pa.field("identifier_type", pa.string()),
    pa.field("identifier", pa.string()),
    pa.field("molecular_form", MOLECULAR_FORM_STRUCT),
]

ONTOLOGY_RELATION_FIELDS = [
    pa.field("predicate", pa.string()),
    pa.field("object", pa.struct(ENTITY_REF_FIELDS)),
    pa.field("ontology_id", pa.string()),
    pa.field("annotations", pa.list_(pa.struct(ANNOTATION_FIELDS))),
]

ASSOCIATION_FIELDS = [
    pa.field("object", pa.struct(ENTITY_REF_FIELDS)),
    pa.field("predicate", pa.string()),
]

BASE_ENTITY_FIELDS = [
    pa.field("type", pa.string()),
    pa.field("identifiers", pa.list_(pa.struct(IDENTIFIER_FIELDS))),
    pa.field("annotations", pa.list_(pa.struct(ANNOTATION_FIELDS))),
    pa.field("ontology_relations", pa.list_(pa.struct(ONTOLOGY_RELATION_FIELDS))),
    pa.field("associations", pa.list_(pa.struct(ASSOCIATION_FIELDS))),
    pa.field("molecular_form", MOLECULAR_FORM_STRUCT),
]

BASE_MEMBERSHIP_FIELDS = [
    pa.field("predicate", pa.string()),
    pa.field("member", pa.struct(BASE_ENTITY_FIELDS)),
    pa.field("is_parent", pa.bool_()),
    pa.field("annotations", pa.list_(pa.struct(ANNOTATION_FIELDS))),
    pa.field("associations", pa.list_(pa.struct(ASSOCIATION_FIELDS))),
]

NESTED_ENTITY_FIELDS = [
    *BASE_ENTITY_FIELDS,
    pa.field("membership", pa.list_(pa.struct(BASE_MEMBERSHIP_FIELDS))),
]

MEMBERSHIP_FIELDS = [
    pa.field("predicate", pa.string()),
    pa.field("member", pa.struct(NESTED_ENTITY_FIELDS)),
    pa.field("is_parent", pa.bool_()),
    pa.field("annotations", pa.list_(pa.struct(ANNOTATION_FIELDS))),
    pa.field("associations", pa.list_(pa.struct(ASSOCIATION_FIELDS))),
]

ENTITY_FIELDS = [
    pa.field("type", pa.string(), nullable=False),
    pa.field("identifiers", pa.list_(pa.struct(IDENTIFIER_FIELDS)), nullable=False),
    pa.field("annotations", pa.list_(pa.struct(ANNOTATION_FIELDS))),
    pa.field("membership", pa.list_(pa.struct(MEMBERSHIP_FIELDS))),
    pa.field("ontology_relations", pa.list_(pa.struct(ONTOLOGY_RELATION_FIELDS))),
    pa.field("associations", pa.list_(pa.struct(ASSOCIATION_FIELDS))),
    pa.field("molecular_form", MOLECULAR_FORM_STRUCT),
]

ENTITY_SCHEMA = pa.schema(ENTITY_FIELDS)

RELATION_FIELDS = [
    pa.field("subject", pa.struct(NESTED_ENTITY_FIELDS), nullable=False),
    pa.field("predicate", pa.string(), nullable=False),
    pa.field("object", pa.struct(NESTED_ENTITY_FIELDS), nullable=False),
    pa.field("category", pa.string()),
    pa.field("interaction_class", pa.string()),
    pa.field("identifiers", pa.list_(pa.struct(IDENTIFIER_FIELDS))),
    pa.field("annotations", pa.list_(pa.struct(ANNOTATION_FIELDS))),
]

RELATION_SCHEMA = pa.schema(RELATION_FIELDS)
