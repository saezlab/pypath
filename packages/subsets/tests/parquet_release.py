"""Write the content of a test release schema as a published Parquet release.

The subset tests describe a small release once, as rows of the PostgreSQL tables a
product reads. This turns the same rows into per-resource Parquet and a release
manifest, so the Parquet products can be compared with the PostgreSQL ones over
identical content. Entity keys are the entity ids as text.
"""

import json
from pathlib import Path

from omnipath_core.fixtures import write_manifest, write_resource

# The namespace slug a Parquet identifier carries for each identifier type.
SLUGS = {
    "Chebi:MI:0474": "chebi",
    "Kegg Compound:MI:2012": "kegg",
    "Hmdb:OM:0004": "hmdb",
    "Pubchem Compound:OM:0002": "pubchem",
    "Lipidmaps:OM:0003": "lipidmaps",
    "Swisslipids:OM:0009": "swisslipids",
    "Chembl Compound:MI:0967": "chembl",
    "Drugbank:MI:2002": "drugbank",
    "Name:OM:0202": "name",
    "Standard Inchi Key:MI:1101": "inchikey",
}


def _rows(cur, query, *args):
    cur.execute(query, args)
    return cur.fetchall()


def _exists(cur, table):
    cur.execute("SELECT to_regclass(%s) IS NOT NULL", [table])
    return cur.fetchone()[0]


def entity_row(key, entity_type, namespace, identifier, label, identifiers, group=None):
    return {
        "entity_key": key,
        "entity_type": entity_type,
        "namespace": namespace,
        "identifier": identifier,
        "taxon": None,
        "label": label,
        "has_hierarchy": False,
        "parent_count": 0,
        "child_count": 0,
        "identifiers": [
            {"ns": ns, "id": value, "is_canonical": canonical, "source": "fixture"}
            for ns, value, canonical in identifiers
        ],
        "annotations": [],
        # The structure group of the entity: its connectivity block is read from here.
        "reference_entity_key": f"inchikey:{group}" if group else None,
        "gene_reference_keys": None,
        "evidence": None,
    }


def relation_row(
    key,
    subject,
    predicate,
    obj,
    *,
    kind,
    source,
    subject_type="small_molecule",
    row=None,
    object_type="small_molecule",
):
    relation = {
        "relation_key": key,
        "statement_kind": kind,
        "subject_entity_key": subject,
        "subject_reference_entity_key": None,
        "object_reference_entity_key": None,
        "subject_label": subject,
        "subject_type": subject_type,
        "predicate": predicate,
        "object_entity_key": obj,
        "object_label": obj,
        "object_type": object_type,
        "taxon": None,
        "is_directed": True,
        "sign": 0,
        "category": "ontology" if kind.startswith("ontology") else "interaction",
        "interaction_class": None,
        "sources": [source],
        "evidence_count": 0,
        "evidence": [],
        "annotations": [],
    }
    if row is not None:
        relation["evidence_count"] = 1
        relation["evidence"] = [
            {
                "source": source,
                "dataset": "fixture",
                "row_id": str(row),
                "upstream_id": None,
                "subject_molecular_form": None,
                "object_molecular_form": None,
                "annotations": [],
            }
        ]
    return relation


def write_resources(root, resources, version="2026.11", *, commit="0123456789abcdef"):
    """A release of the named resources, each given as ``(entities, relations)`` rows."""
    root = Path(root)
    pinned = {}
    for source, (entities, relations) in sorted(resources.items()):
        directory = root / "resources" / source / "1"
        write_resource(directory, entities, relations, [])
        write_manifest(directory, source, "1", provenance={"code_revision": commit})
        pinned[source] = "1"
    path = root / "releases" / f"{version}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({"schema_version": 1, "version": version, "resources": pinned}, indent=2) + "\n"
    )
    return version


def write_release_from_schema(cur, root, version="2026.11"):
    """The release of the schema ``cur`` is set to, written under ``root``."""
    root = Path(root)
    entity_types = dict(_rows(cur, "SELECT entity_type_id, name FROM vocab_entity_type"))
    identifier_types = dict(
        _rows(cur, "SELECT identifier_type_id, name FROM vocab_identifier_type")
    )
    identifiers = {}
    for entity_id, type_id, value in _rows(
        cur,
        "SELECT l.entity_id::text, v.identifier_type_id, v.value FROM entity_identifier_lookup l "
        "JOIN identifier_evidence v ON v.identifier_id = l.identifier_id ORDER BY 1, 2, 3",
    ):
        identifiers.setdefault(entity_id, []).append((SLUGS[identifier_types[type_id]], value))
    groups = {
        entity_id: group_key
        for entity_id, group_key in _rows(
            cur,
            "SELECT entity_id::text, group_key FROM chemical_resolution_group_member "
            "WHERE level_id = 1",
        )
    }
    entities = []
    for entity_id, type_id, canonical_type, canonical, label in _rows(
        cur,
        "SELECT entity_id::text, entity_type_id, canonical_identifier_type_id, "
        "canonical_identifier, label FROM entity ORDER BY 1",
    ):
        namespace = SLUGS[identifier_types[canonical_type]]
        entities.append(
            entity_row(
                entity_id,
                entity_types[type_id],
                namespace,
                canonical,
                label,
                [
                    (ns, value, ns == namespace and value == canonical)
                    for ns, value in identifiers.get(entity_id, [])
                ],
                groups.get(entity_id),
            )
        )

    resources = {}
    if _exists(cur, "entity_ontology_relation"):
        for number, (subject, predicate, obj, source) in enumerate(
            _rows(
                cur,
                "SELECT r.subject_entity_id::text, p.name, r.object_entity_id::text, s.name "
                "FROM entity_ontology_relation r "
                "JOIN vocab_relation_predicate p ON p.relation_predicate_id = r.predicate_id "
                "JOIN data_source s ON s.source_id = r.source_id ORDER BY 4, 1, 2, 3",
            )
        ):
            resources.setdefault(source, []).append(
                relation_row(
                    f"ontology:{number}", subject, predicate, obj, kind="ontology", source=source
                )
            )
    if _exists(cur, "interaction_party"):
        # A reaction is an entity of its own with a has_participant edge to each member.
        reactions = {}
        for interaction, entity_id in _rows(
            cur, "SELECT DISTINCT interaction_id::text, entity_id::text FROM interaction_party"
        ):
            reactions.setdefault(interaction, []).append(entity_id)
        for interaction, members in sorted(reactions.items()):
            key = f"reaction:{interaction}"
            entities.append(entity_row(key, "molecular_activity", "reactome", interaction, key, []))
            for member in members:
                resources.setdefault("reactome", []).append(
                    relation_row(
                        f"{key}>{member}",
                        key,
                        "has_participant",
                        member,
                        kind="relation",
                        source="reactome",
                        subject_type="molecular_activity",
                    )
                )
    resources.setdefault("chebi", [])

    pinned = {}
    for source, relations in sorted(resources.items()):
        directory = root / "resources" / source / "1"
        write_resource(directory, entities, relations, [])
        write_manifest(directory, source, "1")
        pinned[source] = "1"
    path = root / "releases" / f"{version}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({"schema_version": 1, "version": version, "resources": pinned}, indent=2) + "\n"
    )
    return version
