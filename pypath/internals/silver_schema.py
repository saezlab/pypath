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
    AnnotationTypeCv,  # This is the Union of all annotation CV terms
    LicenseCV,
    UpdateCategoryCV,
)

__all__ = [
    'Identifier',
    'Annotation',
    'Membership',
    'Entity',
    'DownloadConfig',
    'Source',
    'ENTITY_SCHEMA',
    'SOURCE_SCHEMA',
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

    term: AnnotationTypeCv
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
    annotations: list[Annotation] | None = None

    def __repr__(self) -> str:
        """Compact representation of membership."""
        member_repr = repr(self.member)
        if self.annotations:
            annot_str = f" [{len(self.annotations)} annot.]"
        else:
            annot_str = ""
        return f"Member({member_repr}{annot_str})"

    def pretty(self, indent: int = 0) -> str:
        """Tree-like representation with indentation."""
        prefix = "  " * indent
        lines = [f"{prefix}Member:"]
        lines.append(self.member.pretty(indent + 1))
        if self.annotations:
            lines.append(f"{prefix}  annotations:")
            for ann in self.annotations:
                lines.append(f"{prefix}    - {ann!r}")
        return "\n".join(lines)

class Entity(NamedTuple):
    """ Entity record matching entities schema."""

    type: EntityTypeCv # e.g. EntityTypeCv.INTERACTION, EntityTypeCv.CV_TERM
    identifiers: List[Identifier]  # e.g. IdentifierNamespaceCv.NAME and .SYNONYM)
    annotations: List[Annotation] | None = None

    members: List[Membership] | None = None # e.g. for complexes and families
    is_member_of: List[Membership] | None = None # e.g. for proteins that are part of complexes or families

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
            connector = "└─" if i == len(self.identifiers) - 1 and not self.annotations and not self.members and not self.is_member_of else "├─"
            lines.append(f"{prefix}│  {connector} {id_!r}")

        # Annotations
        if self.annotations:
            is_last = not self.members and not self.is_member_of
            connector = "└─" if is_last else "├─"
            lines.append(f"{prefix}{connector} annotations:")
            for i, ann in enumerate(self.annotations):
                ann_connector = "└─" if i == len(self.annotations) - 1 else "├─"
                ann_prefix = "   " if is_last else "│  "
                lines.append(f"{prefix}{ann_prefix}{ann_connector} {ann!r}")

        # Members
        if self.members:
            is_last = not self.is_member_of
            connector = "└─" if is_last else "├─"
            lines.append(f"{prefix}{connector} members: ({len(self.members)})")
            for i, member in enumerate(self.members):
                member_connector = "└─" if i == len(self.members) - 1 else "├─"
                member_prefix = "   " if is_last else "│  "

                # Member entity
                member_lines = member.member.pretty(0).split("\n")
                lines.append(f"{prefix}{member_prefix}{member_connector} {member_lines[0]}")
                for line in member_lines[1:]:
                    continuation = "   " if i == len(self.members) - 1 else "│  "
                    lines.append(f"{prefix}{member_prefix}{continuation}  {line}")

                # Member annotations
                if member.annotations:
                    lines.append(f"{prefix}{member_prefix}{'   ' if i == len(self.members) - 1 else '│  '}  └─ annotations:")
                    for ann in member.annotations:
                        lines.append(f"{prefix}{member_prefix}{'   ' if i == len(self.members) - 1 else '│  '}     └─ {ann!r}")

        # Is member of
        if self.is_member_of:
            lines.append(f"{prefix}└─ is_member_of: ({len(self.is_member_of)})")
            for i, membership in enumerate(self.is_member_of):
                member_connector = "└─" if i == len(self.is_member_of) - 1 else "├─"

                # Member entity
                member_lines = membership.member.pretty(0).split("\n")
                lines.append(f"{prefix}   {member_connector} {member_lines[0]}")
                for line in member_lines[1:]:
                    continuation = "   " if i == len(self.is_member_of) - 1 else "│  "
                    lines.append(f"{prefix}   {continuation}  {line}")

                # Member annotations
                if membership.annotations:
                    lines.append(f"{prefix}   {'   ' if i == len(self.is_member_of) - 1 else '│  '}  └─ annotations:")
                    for ann in membership.annotations:
                        lines.append(f"{prefix}   {'   ' if i == len(self.is_member_of) - 1 else '│  '}     └─ {ann!r}")

        return "\n".join(lines)
  
class Resource(NamedTuple):
    """ Source database record."""

    id: str  # e.g. 'uniprot'
    name: str  # e.g. 'UniProt'
    license: LicenseCV
    update_category: UpdateCategoryCV
    
    publication: str | None = None
    url: str | None = None
    description: str | None = None

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
    pa.field('source', pa.string()),
    pa.field('type', pa.string()),
    pa.field('identifiers', pa.list_(pa.struct(IDENTIFIER_FIELDS))),
]

# Membership structure (entity + annotations)
MEMBERSHIP_FIELDS = [
    pa.field('member', pa.struct(BASE_ENTITY_FIELDS)),
    pa.field('annotations', pa.list_(pa.struct(ANNOTATION_FIELDS))),
]

# Full entity schema
ENTITY_FIELDS = [
    pa.field('source', pa.string(), nullable=False),
    pa.field('type', pa.string(), nullable=False),
    pa.field('identifiers', pa.list_(pa.struct(IDENTIFIER_FIELDS)), nullable=False),
    pa.field('annotations', pa.list_(pa.struct(ANNOTATION_FIELDS))),
    pa.field('members', pa.list_(pa.struct(MEMBERSHIP_FIELDS))),
    pa.field('is_member_of', pa.list_(pa.struct(MEMBERSHIP_FIELDS))),
]

ENTITY_SCHEMA = pa.schema(ENTITY_FIELDS)

# Download config fields
DOWNLOAD_CONFIG_FIELDS = [
    pa.field('url', pa.string()),
    pa.field('method', pa.string()),
    pa.field('additional_params', pa.map_(pa.string(), pa.string())),
]

# Source schema
RESOURCE_FIELDS = [
    pa.field('id', pa.string(), nullable=False),
    pa.field('name', pa.string(), nullable=False),
    pa.field('download_config', pa.struct(DOWNLOAD_CONFIG_FIELDS), nullable=False),
    pa.field('license', pa.string(), nullable=False),
    pa.field('update_category', pa.string(), nullable=False),
    pa.field('publication', pa.string()),
    pa.field('url', pa.string()),
    pa.field('description', pa.string()),
]

RESOURCE_SCHEMA = pa.schema(RESOURCE_FIELDS)
