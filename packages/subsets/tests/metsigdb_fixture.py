"""One small MetSigDB release, described once and written two ways.

``write_parquet_release`` writes it as per-resource Parquet and a release manifest.
``load_postgres_schema`` writes the same content into the tables the PostgreSQL build of
the product reads, so both builds can be compared over identical data. Entities are
named by a readable key; the PostgreSQL ids are UUIDs derived from it.
"""

import json
from pathlib import Path
import uuid

from omnipath_core.fixtures import write_manifest, write_resource

from parquet_release import entity_row, relation_row

COMMIT = "0123456789abcdef0123456789abcdef01234567"
STAMP = "parity_metsigdb_v1"
DATASET = "fixture"

# Each resource's source name in the release.
RESOURCES = ("reactome", "wikipathways", "kegg", "macdb", "hmdb", "chemont")

CHEMICAL_A = dict(
    type="chemical_entity",
    canonical=("chebi", "CHEBI:123"),
    label="Chemical A",
    identifiers=[
        ("inchikey", "AAAAAAAAAAAAAA-BBBBBBBBBB-C"),
        ("smiles", "CCO"),
        ("hmdb", "HMDB00008"),
        ("pubchem", "11"),
        ("chebi", "CHEBI:123"),
        ("kegg", "C00001"),
    ],
    group="AAAAAAAAAAAAAA",
)
CHEMICAL_B = dict(
    type="small_molecule",
    canonical=("chebi", "456"),
    identifiers=[("hmdb", "HMDB0000002"), ("chebi", "456")],
)


def set_entity(kind, identifier, **extra):
    return dict(type=kind, canonical=("omnipath", identifier), **extra)


ENTITIES = {
    "chemical_a": CHEMICAL_A,
    "chemical_b": CHEMICAL_B,
    "protein": dict(type="protein", canonical=("uniprot", "P12345")),
    "reactome_set": set_entity(
        "pathway",
        "R-HSA-12345",
        label="Curated Reactome label",
        taxon="9606",
        label_from="reactome",
        names=["Ignored Reactome name"],
    ),
    "wiki_set": set_entity("pathway", "WP123", taxon="10090", names=["Z Wiki name", "A Wiki name"]),
    "kegg_set": set_entity("pathway", "rn00010"),
    "overview_set": set_entity("pathway", "rn01100"),
    "reaction": dict(type="molecular_activity", canonical=("kegg", "kegg_parity_reaction")),
    "trait": set_entity(
        "ontology_class",
        "123",
        label="MACdb trait label",
        taxon="9606",
        label_from="macdb",
        trait_type="cancer_trait",
    ),
    "leaf": set_entity("ontology_class", "CHEMONTID:L", label="Leaf class", label_from="chemont"),
    "middle": set_entity(
        "ontology_class", "CHEMONTID:M", label="Middle class", label_from="chemont"
    ),
    "root": set_entity("ontology_class", "CHEMONTID:R", label="Root class", label_from="chemont"),
    "outside": set_entity("ontology_class", "CHEMONTID:O", label_from="chemont"),
}

# resource, row, subject, predicate, object (the published orientation)
RECORDS = (
    ("reactome", 10, "reactome_set", "has_part", "chemical_a"),
    ("reactome", 11, "chemical_a", "part_of", "reactome_set"),
    ("reactome", 12, "reactome_set", "has_member", "chemical_b"),
    ("reactome", 13, "reactome_set", "has_part", "protein"),
    ("wikipathways", 20, "chemical_a", "participates_in", "wiki_set"),
    ("wikipathways", 21, "wiki_set", "associated_with", "chemical_b"),
    ("kegg", 30, "kegg_set", "associated_with", "reaction"),
    ("kegg", 31, "reaction", "part_of", "overview_set"),
    ("kegg", 32, "reaction", "has_input", "chemical_a"),
    ("kegg", 33, "reaction", "has_output", "chemical_a"),
    ("kegg", 34, "reaction", "has_output", "chemical_b"),
    ("macdb", 40, "chemical_a", "associated_with", "trait"),
    ("macdb", 41, "trait", "associated_with", "chemical_b"),
    ("hmdb", 50, "chemical_a", "associated_with", "leaf"),
    ("hmdb", 51, "chemical_a", "associated_with", "root"),
    ("hmdb", 52, "leaf", "associated_with", "chemical_b"),
)
# resource, subject, predicate, object
ONTOLOGY = (
    ("chemont", "leaf", "subclass_of", "middle"),
    ("chemont", "middle", "subclass_of", "root"),
    ("chemont", "leaf", "part_of", "outside"),
    ("hmdb", "root", "subclass_of", "outside"),
)


def entity_id(key):
    return str(uuid.uuid5(uuid.NAMESPACE_URL, "metsigdb-fixture:" + key))


def key_of(entity_uuid):
    return {entity_id(key): key for key in ENTITIES}[str(entity_uuid)]


def _participants(resource):
    keys = set()
    for source, _, subject, _, obj in RECORDS:
        if source == resource:
            keys |= {subject, obj}
    for source, subject, _, obj in ONTOLOGY:
        if source == resource:
            keys |= {subject, obj}
    return keys


def _parquet_entity(key, resource):
    spec = ENTITIES[key]
    namespace, identifier = spec["canonical"]
    identifiers = [(namespace, identifier, True)]
    identifiers += [(ns, value, False) for ns, value in spec.get("identifiers", ())]
    identifiers += [("name", name, False) for name in spec.get("names", ())]
    row = entity_row(
        key,
        spec["type"],
        namespace,
        identifier,
        spec.get("label") if spec.get("label_from", resource) == resource else None,
        identifiers,
        spec.get("group"),
    )
    # A resource states the organism of its own sets.
    if spec.get("taxon") and _in_record(resource, key):
        row["taxon"] = spec["taxon"]
    if spec.get("trait_type") and resource == "macdb":
        row["annotations"] = [
            {
                "term": "macdb:trait_type",
                "value": spec["trait_type"],
                "quantity": None,
                "source": "macdb",
                "dataset": DATASET,
            }
        ]
    return row


def _in_record(resource, key):
    return any(
        source == resource and key in (subject, obj) for source, _, subject, _, obj in RECORDS
    )


def write_parquet_release(root, version="2026.11", *, commit=COMMIT):
    """The release under ``root``; returns the version."""
    root = Path(root)
    pinned = {}
    for resource in RESOURCES:
        entities = [_parquet_entity(key, resource) for key in sorted(_participants(resource))]
        relations = []
        for source, row, subject, predicate, obj in RECORDS:
            if source != resource:
                continue
            relation = relation_row(
                f"{resource}:{row}",
                subject,
                predicate,
                obj,
                kind="relation",
                source=resource,
                subject_type=ENTITIES[subject]["type"],
            )
            relation["object_type"] = ENTITIES[obj]["type"]
            relation["evidence"] = [
                {
                    "source": resource,
                    "dataset": DATASET,
                    "row_id": str(row),
                    "upstream_id": None,
                    "subject_molecular_form": None,
                    "object_molecular_form": None,
                    "annotations": [],
                }
            ]
            relation["evidence_count"] = 1
            relations.append(relation)
        for number, (source, subject, predicate, obj) in enumerate(ONTOLOGY):
            if source == resource:
                relations.append(
                    relation_row(
                        f"{resource}:ontology:{number}",
                        subject,
                        predicate,
                        obj,
                        kind="ontology",
                        source=resource,
                        subject_type=ENTITIES[subject]["type"],
                    )
                )
        directory = root / "resources" / resource / "1"
        write_resource(directory, entities, relations, [])
        write_manifest(directory, resource, "1", provenance={"code_revision": commit})
        pinned[resource] = "1"
    path = root / "releases" / f"{version}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps({"schema_version": 1, "version": version, "resources": pinned}, indent=2) + "\n"
    )
    return version


IDENTIFIER_TYPES = {
    "inchikey": "Standard Inchi Key:MI:1101",
    "smiles": "Smiles:MI:0239",
    "hmdb": "Hmdb:OM:0004",
    "pubchem": "Pubchem Compound:OM:0002",
    "chebi": "Chebi:MI:0474",
    "kegg": "Kegg Compound:MI:2012",
    "name": "Name:OM:0202",
    "omnipath": "Cv Term Accession:OM:0204",
    "uniprot": "Uniprot:MI:1097",
}


def load_postgres_schema(cur):
    """Create and fill the tables the PostgreSQL MetSigDB build reads, in the current schema."""
    cur.execute(
        """
        CREATE TABLE data_source(source_id bigint PRIMARY KEY, name text);
        CREATE TABLE resources(resource_id text PRIMARY KEY, input_module text,
          input_module_commit text);
        CREATE TABLE vocab_entity_type(entity_type_id serial PRIMARY KEY, name text UNIQUE);
        CREATE TABLE vocab_identifier_type(identifier_type_id bigint PRIMARY KEY, name text);
        CREATE TABLE vocab_relation_predicate(relation_predicate_id serial PRIMARY KEY,
          name text UNIQUE);
        CREATE TABLE entity(entity_id uuid PRIMARY KEY, entity_type_id bigint,
          canonical_identifier_type_id bigint, canonical_identifier text, label text);
        CREATE TABLE identifier_evidence(identifier_id uuid PRIMARY KEY,
          identifier_type_id bigint, value text);
        CREATE TABLE entity_identifier_lookup(entity_id uuid, identifier_id uuid);
        CREATE TABLE relation_evidence(source_id bigint, relation_evidence_id uuid,
          dataset_id bigint, row_id bigint, subject_entity_id uuid, predicate_id bigint,
          object_entity_id uuid);
        CREATE TABLE entity_ontology_relation(source_id bigint, subject_entity_id uuid,
          predicate_id bigint, object_entity_id uuid, ontology_id text);
        CREATE TABLE entity_ontology_term(term_entity_id uuid, term_id text, label text,
          ontology_id text);
        CREATE TABLE chemical_resolution_group_member(level_id smallint, group_key text,
          entity_id uuid, inchikey text);
        CREATE TABLE entity_evidence(source_id bigint, entity_evidence_id uuid, taxonomy_id bigint);
        CREATE TABLE entity_evidence_resolution(source_id bigint, entity_evidence_id uuid,
          entity_id uuid);
        CREATE TABLE annotation(annotation_key uuid PRIMARY KEY, term text, value text);
        CREATE TABLE entity_evidence_annotation(source_id bigint, entity_evidence_id uuid,
          annotation_key uuid);
        """
    )
    source_ids = {name: 500001 + number for number, name in enumerate(RESOURCES)}
    for name, source_id in source_ids.items():
        cur.execute("INSERT INTO data_source VALUES (%s,%s)", [source_id, name])
        cur.execute(
            "INSERT INTO resources VALUES (%s,%s,%s)",
            [name, "omnipath_resources.inputs_v2." + name, COMMIT],
        )
    type_ids = {}
    for kind in sorted({spec["type"] for spec in ENTITIES.values()}):
        cur.execute(
            "INSERT INTO vocab_entity_type(name) VALUES (%s) RETURNING entity_type_id", [kind]
        )
        type_ids[kind] = cur.fetchone()[0]
    namespace_ids = {}
    for number, (slug, name) in enumerate(IDENTIFIER_TYPES.items(), start=1):
        cur.execute("INSERT INTO vocab_identifier_type VALUES (%s,%s)", [number, name])
        namespace_ids[slug] = number
    predicates = {r[3] for r in RECORDS} | {o[2] for o in ONTOLOGY}
    predicate_ids = {}
    for name in sorted(predicates):
        cur.execute(
            "INSERT INTO vocab_relation_predicate(name) VALUES (%s) RETURNING relation_predicate_id",
            [name],
        )
        predicate_ids[name] = cur.fetchone()[0]

    for key, spec in ENTITIES.items():
        namespace, identifier = spec["canonical"]
        cur.execute(
            "INSERT INTO entity VALUES (%s,%s,%s,%s,%s)",
            [
                entity_id(key),
                type_ids[spec["type"]],
                namespace_ids[namespace],
                identifier,
                spec.get("label") if key in ("chemical_a",) else None,
            ],
        )
        if spec.get("label") and key != "chemical_a":
            cur.execute(
                "INSERT INTO entity_ontology_term VALUES (%s,%s,%s,%s)",
                [entity_id(key), identifier, spec["label"], spec["label_from"]],
            )
        items = list(spec.get("identifiers", ())) + [("name", n) for n in spec.get("names", ())]
        for number, (slug, value) in enumerate(items):
            identifier_id = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([key, slug, value])))
            cur.execute(
                "INSERT INTO identifier_evidence VALUES (%s,%s,%s)",
                [identifier_id, namespace_ids[slug], value],
            )
            cur.execute(
                "INSERT INTO entity_identifier_lookup VALUES (%s,%s)",
                [entity_id(key), identifier_id],
            )
        if spec.get("group"):
            cur.execute(
                "INSERT INTO chemical_resolution_group_member VALUES (1,%s,%s,%s)",
                [spec["group"], entity_id(key), dict(spec["identifiers"])["inchikey"]],
            )

    for source, row, subject, predicate, obj in RECORDS:
        source_id = source_ids[source]
        evidence = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps(["record", source, row])))
        cur.execute(
            "INSERT INTO relation_evidence VALUES (%s,%s,%s,%s,%s,%s,%s)",
            [
                source_id,
                evidence,
                source_id,
                row,
                entity_id(subject),
                predicate_ids[predicate],
                entity_id(obj),
            ],
        )
        for key in (subject, obj):
            taxon = ENTITIES[key].get("taxon")
            occurrence = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([source, row, key])))
            cur.execute(
                "INSERT INTO entity_evidence VALUES (%s,%s,%s)",
                [source_id, occurrence, int(taxon) if taxon else None],
            )
            cur.execute(
                "INSERT INTO entity_evidence_resolution VALUES (%s,%s,%s)",
                [source_id, occurrence, entity_id(key)],
            )
            if source == "macdb" and ENTITIES[key].get("trait_type"):
                annotation = str(uuid.uuid5(uuid.NAMESPACE_URL, "macdb:trait_type=" + key))
                cur.execute(
                    "INSERT INTO annotation VALUES (%s,'macdb:trait_type',%s) ON CONFLICT DO NOTHING",
                    [annotation, ENTITIES[key]["trait_type"]],
                )
                cur.execute(
                    "INSERT INTO entity_evidence_annotation VALUES (%s,%s,%s)",
                    [source_id, occurrence, annotation],
                )
    for source, subject, predicate, obj in ONTOLOGY:
        cur.execute(
            "INSERT INTO entity_ontology_relation VALUES (%s,%s,%s,%s,%s)",
            [
                source_ids[source],
                entity_id(subject),
                predicate_ids[predicate],
                entity_id(obj),
                source,
            ],
        )
