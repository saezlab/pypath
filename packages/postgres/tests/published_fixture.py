"""Small published Parquet fixtures shared by the PostgreSQL projection tests."""

from copy import deepcopy

from omnipath_core.fixtures import write_resource

ENTITY_A = "published::protein/A:'α"
ENTITY_B = "published::protein/B"
ENTITY_C = "published::chemical/standalone"
RELATION = "published::relation/key:'β"
QUANTITY = {
    "has_numeric_value": 12.5,
    "has_unit": "UO:0000061",
    "has_unit_prefix": "nano",
    "has_binary_relation": "less_than",
    "source_field": "IC50 (nM)",
    "comparator": "<=",
}
RAW_PAYLOAD = '{ "upstream": "SRC-42", "value": "≤ 12.5 nM", "nested": [null, true] }\n'


def annotation(term="has_attribute", *, value=None, quantity=None, scope="evidence"):
    return {
        "term": term,
        "value": value,
        "quantity": quantity,
        "source": "reported-source",
        "dataset": "reported-dataset",
        "scope": scope,
    }


def entity(key, identifier, *, namespace="uniprot", taxon="9606", entity_type="protein"):
    return {
        "entity_key": key,
        "entity_type": entity_type,
        "namespace": namespace,
        "identifier": identifier,
        "taxon": taxon,
        "label": identifier,
        "has_hierarchy": False,
        "parent_count": 0,
        "child_count": 0,
        "identifiers": [],
        "annotations": [],
        "reference_entity_key": None,
        "gene_reference_keys": None,
        "evidence": None,
    }


def fixture_rows():
    a = entity(ENTITY_A, "P04637")
    identifier = {"ns": "uniprot", "id": "P04637", "is_canonical": True, "source": "source"}
    a["identifiers"] = [deepcopy(identifier), deepcopy(identifier)]
    quantity = annotation(quantity=deepcopy(QUANTITY))
    quantity.pop("scope")  # Entity annotations have no scope.
    scalar = annotation("description", value="A source description with α")
    scalar.pop("scope")
    a["annotations"] = [quantity, deepcopy(scalar), deepcopy(scalar)]
    b = entity(ENTITY_B, "P0DP23")
    b["identifiers"] = None
    b["annotations"] = None
    c = entity(ENTITY_C, "CHEBI:15377", namespace="chebi", taxon=None, entity_type="small_molecule")
    evidence = {
        "source": "reported-source",
        "dataset": "reported-dataset",
        "row_id": "00001",
        "upstream_id": "SRC-42",
        "subject_molecular_form": None,
        "object_molecular_form": None,
        "annotations": [
            annotation(quantity=deepcopy(QUANTITY)),
            annotation("publications", value="PMID:123456"),
        ],
    }
    relation = {
        "relation_key": RELATION,
        "statement_kind": "relation",
        "subject_entity_key": ENTITY_A,
        "subject_reference_entity_key": None,
        "object_reference_entity_key": None,
        "subject_label": "Subject α",
        "subject_type": "protein",
        "predicate": "affects",
        "object_entity_key": ENTITY_B,
        "object_label": "Object β",
        "object_type": "protein",
        "taxon": "9606",
        "is_directed": True,
        "sign": -1,
        "category": "interaction",
        "interaction_class": "chemical_modulation",
        "sources": ["reported-source", "other-source"],
        "evidence_count": 2,
        "evidence": [deepcopy(evidence), deepcopy(evidence)],
        "annotations": [
            annotation(quantity=deepcopy(QUANTITY), scope="relation"),
            annotation("object_direction_qualifier", value="decreased", scope="relation"),
        ],
    }
    payloads = [
        {
            "relation_key": RELATION,
            "entity_key": None,
            "source": "reported-source",
            "row_id": "00001",
            "payload_json": RAW_PAYLOAD,
        },
        {
            "relation_key": None,
            "entity_key": ENTITY_C,
            "source": "entity-source",
            "row_id": "entity-only",
            "payload_json": '{"standalone": true}',
        },
    ]
    return [a, b, c], [relation], payloads


def write_fixture(directory, *, entities=None, relations=None, payloads=None):
    defaults = fixture_rows()
    rows = (
        entities if entities is not None else defaults[0],
        relations if relations is not None else defaults[1],
        payloads if payloads is not None else defaults[2],
    )
    write_resource(directory, *rows)
    return rows
