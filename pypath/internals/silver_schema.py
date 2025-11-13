"""
Canonical definitions for silver-layer schemas.

Keeping the PyArrow schema objects and the namedtuple helpers in one place
avoids accidental divergence across the pipeline.
"""
from typing import List, NamedTuple, Self


import pyarrow as pa

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
)

__all__ = [
    'Identifier',
    'Annotation',
    'Membership',
    'Entity',
    'ENTITY_SCHEMA',
]

class Identifier(NamedTuple):
    """An identifier with its type."""

    type: IdentifierNamespaceCv
    value: str

    def __repr__(self) -> str:
        """Compact representation showing namespace and value."""
        type_name = self.type.name if hasattr(self.type, 'name') else str(self.type)
        return f"{type_name}:{self.value}"

class Annotation(NamedTuple):
    """Annotation for an interaction or entity."""

    term: str
    value: str | float | None = None
    units: str | None = None

    def __repr__(self) -> str:
        """Compact representation of annotation."""
        term_str = self.term.name if hasattr(self.term, 'name') else str(self.term)
        if self.value is None:
            return f"{term_str}"
        elif self.units:
            return f"{term_str}={self.value}{self.units}"
        else:
            return f"{term_str}={self.value}"

class Membership(NamedTuple):
    member: 'Entity'  # Forward reference since Entity is defined below
    is_parent: bool = False  # For an entity self, in each Membership, is_parent=True means member is a parent of self; is_parent=False means self is parent of member
    annotations: list[Annotation] | None = None

    def __repr__(self) -> str:
        """Compact representation of membership."""
        member_repr = repr(self.member)
        if self.annotations:
            annot_str = f" [membership: {len(self.annotations)} annot.]"
        else:
            annot_str = ""
        return f"Membership(entity={member_repr}{annot_str})"

    def pretty(self, indent: int = 0) -> str:
        """Tree-like representation with indentation."""
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
            # Update entity connector to └─ if no annotations
            lines[1] = f"{prefix}└─ entity:"
            lines = [lines[0]] + [lines[1]] + [f"{prefix}   {line}" for line in entity_lines]
        return "\n".join(lines)

class Entity(NamedTuple):
    """ Entity record matching entities schema."""

    type: EntityTypeCv # e.g. EntityTypeCv.INTERACTION, EntityTypeCv.CV_TERM
    identifiers: List[Identifier]  # e.g. IdentifierNamespaceCv.NAME and .SYNONYM)
    annotations: List[Annotation] | None = None

    members: List[Membership] | None = None 

    def __repr__(self) -> str:
        """Tree-like detailed representation."""
        return self.pretty()

    def __str__(self) -> str:
        """Tree-like detailed representation."""
        return self.pretty()

    def pretty(self, indent: int = 0) -> str:
        """Tree-like representation with full details."""
        prefix = "  " * indent
        type_name = self.type.name if hasattr(self.type, 'name') else str(self.type)

        lines = [f"{prefix}{type_name}"]

        # Identifiers
        lines.append(f"{prefix}├─ identifiers:")
        for i, id_ in enumerate(self.identifiers):
            connector = "└─" if i == len(self.identifiers) - 1 and not self.annotations and not self.members else "├─"
            lines.append(f"{prefix}│  {connector} {id_!r}")

        # Annotations
        if self.annotations:
            is_last = not self.members
            connector = "└─" if is_last else "├─"
            lines.append(f"{prefix}{connector} annotations:")
            for i, ann in enumerate(self.annotations):
                ann_connector = "└─" if i == len(self.annotations) - 1 else "├─"
                ann_prefix = "   " if is_last else "│  "
                lines.append(f"{prefix}{ann_prefix}{ann_connector} {ann!r}")

        # Members
        if self.members:
            is_last = True
            connector = "└─" if is_last else "├─"
            lines.append(f"{prefix}{connector} members: ({len(self.members)})")
            for i, member in enumerate(self.members):
                member_connector = "└─" if i == len(self.members) - 1 else "├─"
                member_prefix = "   " if is_last else "│  "

                # Membership wrapper
                lines.append(f"{prefix}{member_prefix}{member_connector} Membership:")

                # Member entity
                continuation = "   " if i == len(self.members) - 1 else "│  "
                member_lines = member.member.pretty(0).split("\n")
                lines.append(f"{prefix}{member_prefix}{continuation}  ├─ entity: {member_lines[0]}")
                for line in member_lines[1:]:
                    lines.append(f"{prefix}{member_prefix}{continuation}  │        {line}")

                # Membership annotations
                if member.annotations:
                    lines.append(f"{prefix}{member_prefix}{continuation}  └─ membership_annotations:")
                    for j, ann in enumerate(member.annotations):
                        ann_connector = "└─" if j == len(member.annotations) - 1 else "├─"
                        lines.append(f"{prefix}{member_prefix}{continuation}     {ann_connector} {ann!r}")
                else:
                    # Close the entity branch if no annotations
                    # Replace the last "├─ entity:" with "└─ entity:"
                    if lines[-len(member_lines)-1].endswith("├─ entity: " + member_lines[0]):
                        lines[-len(member_lines)-1] = f"{prefix}{member_prefix}{continuation}  └─ entity: {member_lines[0]}"

        return "\n".join(lines)

#### PyArrow Schemas (needed for parquet files) ####

# Reusable field definitions
IDENTIFIER_FIELDS = [
    pa.field('type', pa.string()),
    pa.field('value', pa.string()),
]

ANNOTATION_FIELDS = [
    pa.field('term', pa.string()),
    pa.field('value', pa.string()),
    pa.field('units', pa.string()),
]

# Base entity fields (without members/is_member_of to avoid circular reference)
BASE_ENTITY_FIELDS = [
    pa.field('type', pa.string()),
    pa.field('identifiers', pa.list_(pa.struct(IDENTIFIER_FIELDS))),
]

# Membership structure (entity + is_parent + annotations)
MEMBERSHIP_FIELDS = [
    pa.field('member', pa.struct(BASE_ENTITY_FIELDS)),
    pa.field('is_parent', pa.bool_()),
    pa.field('annotations', pa.list_(pa.struct(ANNOTATION_FIELDS))),
]

# Full entity schema
ENTITY_FIELDS = [
    pa.field('type', pa.string(), nullable=False),
    pa.field('identifiers', pa.list_(pa.struct(IDENTIFIER_FIELDS)), nullable=False),
    pa.field('annotations', pa.list_(pa.struct(ANNOTATION_FIELDS))),
    pa.field('members', pa.list_(pa.struct(MEMBERSHIP_FIELDS))),
]

ENTITY_SCHEMA = pa.schema(ENTITY_FIELDS)
