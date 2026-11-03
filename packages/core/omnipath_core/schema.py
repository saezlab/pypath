"""Canonical PyArrow schemas of the OmniPath Parquet tables.

Each resource version publishes the normalized tables in ``PUBLISHED_TABLES`` and
serves lookups from the ``SERVING_TABLES``; see the section below.
"""

from __future__ import annotations

from omnipath_core.measurements import QUANTITY_STRUCT
from omnipath_core.molecular_forms import MOLECULAR_FORM_STRUCT
import pyarrow as pa

ANNOTATION_STRUCT = pa.struct(
    [
        ("term", pa.string()),
        ("value", pa.string()),
        ("quantity", QUANTITY_STRUCT),
        ("source", pa.string()),
        ("dataset", pa.string()),
        ("scope", pa.string()),
    ]
)

ENTITY_ANNOTATION_STRUCT = pa.struct(
    [
        ("term", pa.string()),
        ("value", pa.string()),
        ("quantity", QUANTITY_STRUCT),
        ("source", pa.string()),
        ("dataset", pa.string()),
    ]
)

ENTITY_EVIDENCE_STRUCT = pa.struct(
    [
        ("source", pa.string()),
        ("dataset", pa.string()),
        ("row_id", pa.string()),
        ("upstream_id", pa.string()),
        ("annotations", pa.list_(ENTITY_ANNOTATION_STRUCT)),
        ("molecular_form", MOLECULAR_FORM_STRUCT),
    ]
)

PAYLOAD_SCHEMA = pa.schema(
    [
        ("relation_key", pa.string()),
        ("entity_key", pa.string()),
        ("source", pa.string()),
        ("row_id", pa.string()),
        ("payload_json", pa.string()),
    ]
)


# ---------------------------------------------------------------------------------------
# Normalized resource tables. One file per table and resource version, each sorted by its
# lookup key in 8,192-row groups. ``entity_id`` and ``relation_id`` are a row's position in
# ``entity.parquet`` / ``relation.parquet`` (which are sorted by their hash keys); child
# tables refer to their parent by that id, the 64-character keys stay the cross-resource
# identities. Serving tables are derived lookup indexes and are not part of downloads.

ROW_GROUP_SIZE = 8192

_QUANTITY_FIELDS = [(f"quantity_{field.name}", field.type) for field in QUANTITY_STRUCT]

ENTITY_TABLE = pa.schema(
    [
        ("entity_id", pa.int32()),
        ("entity_key", pa.string()),
        ("entity_type", pa.string()),
        ("namespace", pa.string()),
        ("identifier", pa.string()),
        ("taxon", pa.string()),
        ("label", pa.string()),
        ("has_hierarchy", pa.bool_()),
        ("parent_count", pa.int64()),
        ("child_count", pa.int64()),
        ("reference_entity_key", pa.string()),
        ("gene_reference_keys", pa.list_(pa.string())),
        # First InChIKey block of an InChIKey reference: the chemical structure group.
        ("group_connectivity", pa.string()),
        ("identifier_count", pa.int64()),
        ("annotation_count", pa.int64()),
        ("evidence_count", pa.int64()),
        ("relation_count", pa.int64()),
    ]
)

ENTITY_IDENTIFIER_TABLE = pa.schema(
    [
        ("entity_id", pa.int32()),
        ("ordinal", pa.int32()),
        ("ns", pa.string()),
        ("id", pa.string()),
        ("is_canonical", pa.bool_()),
        ("source", pa.string()),
    ]
)

ENTITY_ANNOTATION_TABLE = pa.schema(
    [
        ("entity_id", pa.int32()),
        ("ordinal", pa.int32()),
        ("term", pa.string()),
        ("value", pa.string()),
        *_QUANTITY_FIELDS,
        ("source", pa.string()),
        ("dataset", pa.string()),
    ]
)

ENTITY_EVIDENCE_TABLE = pa.schema(
    [
        ("entity_id", pa.int32()),
        ("ordinal", pa.int32()),
        ("source", pa.string()),
        ("dataset", pa.string()),
        ("row_id", pa.string()),
        ("upstream_id", pa.string()),
        ("annotations", pa.list_(ENTITY_ANNOTATION_STRUCT)),
        ("molecular_form", MOLECULAR_FORM_STRUCT),
    ]
)

RELATION_QUALIFIERS = (
    "object_aspect_qualifier",
    "object_direction_qualifier",
    "causal_mechanism_qualifier",
)

RELATION_TABLE = pa.schema(
    [
        ("relation_id", pa.int32()),
        ("relation_key", pa.string()),
        ("statement_kind", pa.string()),
        ("subject_entity_key", pa.string()),
        ("subject_reference_entity_key", pa.string()),
        ("subject_label", pa.string()),
        ("subject_type", pa.string()),
        ("predicate", pa.string()),
        ("object_entity_key", pa.string()),
        ("object_reference_entity_key", pa.string()),
        ("object_label", pa.string()),
        ("object_type", pa.string()),
        ("taxon", pa.string()),
        ("is_directed", pa.bool_()),
        ("sign", pa.int32()),
        ("category", pa.string()),
        ("interaction_class", pa.string()),
        ("sources", pa.list_(pa.string())),
        ("evidence_count", pa.int64()),
        ("annotation_count", pa.int64()),
        *[(name, pa.list_(pa.string())) for name in RELATION_QUALIFIERS],
    ]
)

RELATION_ANNOTATION_TABLE = pa.schema(
    [
        ("relation_id", pa.int32()),
        ("ordinal", pa.int32()),
        ("term", pa.string()),
        ("value", pa.string()),
        *_QUANTITY_FIELDS,
        ("source", pa.string()),
        ("dataset", pa.string()),
        ("scope", pa.string()),
    ]
)

RELATION_EVIDENCE_TABLE = pa.schema(
    [
        ("relation_id", pa.int32()),
        ("ordinal", pa.int32()),
        ("source", pa.string()),
        ("dataset", pa.string()),
        ("row_id", pa.string()),
        ("upstream_id", pa.string()),
        ("annotations", pa.list_(ANNOTATION_STRUCT)),
        ("subject_molecular_form", MOLECULAR_FORM_STRUCT),
        ("object_molecular_form", MOLECULAR_FORM_STRUCT),
    ]
)

# Each relation under its endpoints' entity keys and reference keys (key_kind).
RELATION_ENDPOINT_TABLE = pa.schema(
    [
        ("key", pa.string()),
        ("key_kind", pa.string()),
        ("side", pa.string()),
        ("relation_id", pa.int32()),
        ("predicate", pa.string()),
        ("category", pa.string()),
    ]
)

# Entities by group: kind ``reference`` (reference_entity_key), ``gene`` (a product's
# gene reference) or ``connectivity`` (first InChIKey block).
ENTITY_GROUP_TABLE = pa.schema(
    [("group_key", pa.string()), ("kind", pa.string()), ("entity_id", pa.int32())]
)

# Lowercase labels and identifiers for prefix and exact search.
ENTITY_TERM_TABLE = pa.schema(
    [("term", pa.string()), ("kind", pa.string()), ("entity_id", pa.int32())]
)

# Each distinct combination of the relation filter columns, with the largest evidence count
# and the number of its relations: which filter options a selection leaves, without reading
# the relations. ``sources`` is the resource's own name, as in ``relation``.
RELATION_FACET_COLUMNS = (
    "category",
    "predicate",
    "subject_type",
    "object_type",
    "taxon",
    "sign",
    "is_directed",
    "sources",
    *RELATION_QUALIFIERS,
)
RELATION_FACET_TABLE = pa.schema(
    [
        *[RELATION_TABLE.field(name) for name in RELATION_FACET_COLUMNS],
        ("evidence_count", pa.int64()),
        ("relation_count", pa.int64()),
    ]
)

PUBLISHED_TABLES = {
    "entity": ENTITY_TABLE,
    "entity_identifier": ENTITY_IDENTIFIER_TABLE,
    "entity_annotation": ENTITY_ANNOTATION_TABLE,
    "entity_evidence": ENTITY_EVIDENCE_TABLE,
    "relation": RELATION_TABLE,
    "relation_annotation": RELATION_ANNOTATION_TABLE,
    "relation_evidence": RELATION_EVIDENCE_TABLE,
    "evidence_payloads": PAYLOAD_SCHEMA,
}

SERVING_TABLES = {
    "relation_endpoint": RELATION_ENDPOINT_TABLE,
    "entity_group": ENTITY_GROUP_TABLE,
    "entity_term": ENTITY_TERM_TABLE,
    "relation_facet": RELATION_FACET_TABLE,
}
