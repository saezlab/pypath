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

class Annotation(NamedTuple):
    """Annotation for an interaction or entity."""

    term: AnnotationTypeCv
    value: str | float | None = None
    units: str | None = None

class Membership(NamedTuple):
    member: 'Entity'  # Forward reference since Entity is defined below
    annotations: list[Annotation] | None = None

class Entity(NamedTuple):
    """ Entity record matching entities schema."""

    # Required fields
    source: str
    type: EntityTypeCv # e.g. EntityTypeCv.INTERACTION, EntityTypeCv.CV_TERM
    identifiers: List[Identifier]  # e.g. IdentifierNamespaceCv.NAME and .SYNONYM)
    annotations: List[Annotation] | None = None

    members: List[Membership] | None = None # e.g. for complexes and families
    is_member_of: List[Membership] | None = None # e.g. for proteins that are part of complexes or families
  
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
