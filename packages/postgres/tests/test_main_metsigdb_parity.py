"""Bounded all-five MetSigDB oracle against unchanged frozen main.

Twenty declared source records represent the same resolved chemistry and sets
in the original vocabulary and Biolink. The reference has genuine matched
occurrences; the port has truthful published occurrences. That approved
eligibility difference and the chemical type vocabulary are explicit, while
complete memberships, source records, context, identifiers and DDL must agree.
No parser, resolver, resource build, network mapping or existing schema is used.
"""

from __future__ import annotations

from collections import Counter
import importlib
import json
from pathlib import Path
import uuid
from oracles import read_legacy

# Reuse the independent private importer and random-schema lifecycle, not any
# assertions or data from the completed scientific fixtures.
import test_main_parity_contract as reference_contract

MAIN = reference_contract.MAIN
ORACLE_PREFIX = reference_contract.ORACLE_PREFIX
_catalogue = reference_contract._catalogue
oracle = reference_contract.oracle
port = reference_contract.port
parity_connection = reference_contract.parity_connection
paired_schemas = reference_contract.paired_schemas


FROZEN_DIGESTS = {
    "__init__.py": "d9dd048ee1edd05156d187114161982d376daa9deea16e793626089147b1814d",
    "build.py": "c6ca09f3ce888bbdb216ea2c3557c55486e03c2d451aff7a46a71bbd59398ce7",
    "mapping.py": "68812f9cf15c78fccbc6eb1402307a0031de1877d00706b0426da9969aa757d9",
    "sql/extract_classyfire.sql": "b1cf3eebb7e7d221a7c2ca8b68860e665be226847365600e8fe1f243ea0f4fd9",
    "sql/extract_kegg.sql": "a068e3fabb645e1c668e7f623723808308709c477fb84f0bc653dca5efc0c231",
    "sql/extract_macdb.sql": "aaca63a6a3d9898c61e4c2e3fcdf6a096650585ef82e348cedd3b8b03348d694",
    "sql/extract_onehop.sql": "864af5d7c85d78e60ab698580eb60b0cfc6f4e83d85d90e381ed9c2451fbca17",
    "sql/membership_table.sql": "8efa89ee10d5eaf70a7e5e7347d3e87cc745d42bee7452d124fb125cbe365780",
    "sql/publish_membership.sql": "885a570a74c5a813c12258bfa41b4a83975cb870a18c168dd0dd3556561aa129",
    "sql/upsert_membership.sql": "aa5722304f9deec512bcd91a4d82ac64daa53115f6fdf1ac49528b020d5afebd",
}
SOURCE_IDS = {
    name: 500001 + i
    for i, name in enumerate(("reactome", "wikipathways", "kegg", "macdb", "hmdb", "chemont"))
}
ENTITY_IDS = {
    name: str(uuid.UUID(int=60000 + i))
    for i, name in enumerate(
        (
            "chemical_a",
            "chemical_b",
            "protein",
            "reactome_set",
            "wiki_set",
            "kegg_set",
            "overview_set",
            "reaction",
            "trait",
            "trait_type",
            "leaf",
            "middle",
            "root",
            "outside",
        )
    )
}
COMMIT = "0123456789abcdef0123456789abcdef01234567"
STAMP = "parity_metsigdb_v1"

# Each tuple is one source row; endpoint occurrences share that row identity.
# Main is set -> chemical associated_with. The port exercises equivalent
# pathway membership orientation and both KEGG participant roles.
RECORDS = (
    ("reactome", 10, "reactome_set", "associated_with", "chemical_a", "has_part", False),
    ("reactome", 11, "reactome_set", "associated_with", "chemical_a", "part_of", True),
    ("reactome", 12, "reactome_set", "associated_with", "chemical_b", "has_member", False),
    ("reactome", 13, "reactome_set", "associated_with", "protein", "has_part", False),
    ("wikipathways", 20, "wiki_set", "associated_with", "chemical_a", "participates_in", True),
    ("wikipathways", 21, "wiki_set", "associated_with", "chemical_b", "associated_with", False),
    ("kegg", 30, "kegg_set", "associated_with", "reaction", "has_participant", False),
    ("kegg", 31, "overview_set", "associated_with", "reaction", "part_of", True),
    ("kegg", 32, "reaction", "has_participant", "chemical_a", "has_input", False),
    ("kegg", 33, "reaction", "has_participant", "chemical_a", "has_output", False),
    ("kegg", 34, "reaction", "has_participant", "chemical_b", "has_output", False),
    ("macdb", 40, "trait", "associated_with", "chemical_a", "associated_with", True),
    ("macdb", 41, "trait", "associated_with", "chemical_b", "associated_with", False),
    ("hmdb", 50, "leaf", "associated_with", "chemical_a", "associated_with", True),
    ("hmdb", 51, "root", "associated_with", "chemical_a", "associated_with", False),
    ("hmdb", 52, "leaf", "associated_with", "chemical_b", "associated_with", False),
)
ONTOLOGY_RECORDS = (
    ("chemont", 60, "leaf", "is_a", "middle"),
    ("chemont", 61, "middle", "is_a", "root"),
    ("chemont", 62, "leaf", "part_of", "outside"),
    ("hmdb", 63, "root", "is_a", "outside"),
)


def test_metsigdb_reference_is_exact_main():
    for relative, digest in FROZEN_DIGESTS.items():
        read_legacy("metsigdb/" + relative, digest)
    assert len(RECORDS) + len(ONTOLOGY_RECORDS) == 20
    assert {r[0] for r in RECORDS} == set(SOURCE_IDS) - {"chemont"}


def test_metsigdb_publication_identity_and_table_ddl_are_exact_main(port):
    from omnipath_subsets import metsigdb

    current = Path(metsigdb.__file__).parent / "sql"
    for name in ("membership_table.sql", "upsert_membership.sql"):
        assert (current / name).read_text() == read_legacy("metsigdb/sql/" + name)


def _evidence_id(source, row, entity):
    return str(
        uuid.uuid5(uuid.NAMESPACE_URL, json.dumps(["metsigdb_occurrence", source, row, entity]))
    )


def _substrate(conn, schema, *, biolink):
    from psycopg2 import sql

    def q(text):
        return sql.SQL(text).format(s=sql.Identifier(schema))

    with conn.cursor() as cur:
        for source, sid in SOURCE_IDS.items():
            cur.execute(
                q("INSERT INTO {s}.data_source(source_id,name) VALUES (%s,%s)"), [sid, source]
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.dataset(dataset_id,source_id,name) VALUES (%s,%s,'bounded_main_parity')"
                ),
                [sid, sid],
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.resources(resource_id,input_module,input_module_commit) VALUES (%s,%s,%s)"
                ),
                [source, "pypath.inputs_v2." + source, COMMIT],
            )
        type_names = {
            "chemical_a": "chemical_entity" if biolink else "Chemical:OM:0037",
            "chemical_b": "small_molecule" if biolink else "Chemical:OM:0037",
            "protein": "protein" if biolink else "Protein:MI:0326",
            "pathway": "pathway" if biolink else "Pathway:OM:0014",
            "reaction": "molecular_activity" if biolink else "Reaction:OM:0016",
            "ontology": "ontology_class" if biolink else "Cv Term:OM:0012",
        }
        for name in set(type_names.values()):
            cur.execute(
                q(
                    "INSERT INTO {s}.vocab_entity_type(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
                ),
                [name],
            )
        cur.execute(q("SELECT name,entity_type_id FROM {s}.vocab_entity_type"))
        types = dict(cur.fetchall())
        cur.execute(
            q(
                "INSERT INTO {s}.vocab_identifier_type(identifier_type_id,name) VALUES (590000,'Smiles:MI:0239') ON CONFLICT(name) DO NOTHING"
            )
        )
        cur.execute(q("SELECT name,identifier_type_id FROM {s}.vocab_identifier_type"))
        namespaces = dict(cur.fetchall())
        canonical = {
            "chemical_a": "123",
            "chemical_b": "456",
            "protein": "P12345",
            "reactome_set": "R-HSA-12345",
            "wiki_set": "WP123",
            "kegg_set": "rn00010",
            "overview_set": "rn01100",
            "reaction": "kegg_parity_reaction",
            "trait": "123",
            "trait_type": "cancer_trait",
            "leaf": "CHEMONT:L",
            "middle": "CHEMONT:M",
            "root": "CHEMONT:R",
            "outside": "CHEMONT:O",
        }
        entity_types = {}
        for name, eid in ENTITY_IDS.items():
            kind = (
                name
                if name in ("chemical_a", "chemical_b", "protein", "reaction")
                else "pathway"
                if name.endswith("set")
                else "ontology"
            )
            entity_types[name] = types[type_names[kind]]
            namespace = (
                "Chebi:MI:0474"
                if name.startswith("chemical_")
                else "Uniprot:MI:1097"
                if name == "protein"
                else "omnipath:reaction_member_hash"
                if name == "reaction"
                else "Cv Term Accession:OM:0204"
            )
            label = "Chemical A" if name == "chemical_a" else None
            cur.execute(
                q(
                    "INSERT INTO {s}.entity(entity_id,entity_type_id,canonical_identifier_type_id,canonical_identifier,resolution_status_id,label) VALUES (%s,%s,%s,%s,2,%s)"
                ),
                [eid, entity_types[name], namespaces[namespace], canonical[name], label],
            )

        def identifier(entity, namespace, value, source="reactome"):
            key = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([namespace, value])))
            cur.execute(
                q(
                    "INSERT INTO {s}.identifier_evidence(identifier_id,identifier_type_id,value) VALUES (%s,%s,%s) ON CONFLICT DO NOTHING"
                ),
                [key, namespaces[namespace], value],
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.entity_identifier(source_id,entity_id,identifier_id) VALUES (%s,%s,%s) ON CONFLICT DO NOTHING"
                ),
                [SOURCE_IDS[source], ENTITY_IDS[entity], key],
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.entity_identifier_lookup(entity_id,identifier_id) VALUES (%s,%s) ON CONFLICT DO NOTHING"
                ),
                [ENTITY_IDS[entity], key],
            )

        for namespace, value in (
            ("Standard Inchi Key:MI:1101", "AAAAAAAAAAAAAA-BBBBBBBBBB-C"),
            ("Smiles:MI:0239", "CCO"),
            ("Hmdb:OM:0004", "HMDB00008"),
            ("Pubchem Compound:OM:0002", "11"),
            ("Chebi:MI:0474", "CHEBI:123" if biolink else "123"),
            ("Kegg Compound:MI:2012", "C00001"),
        ):
            identifier("chemical_a", namespace, value)
        # A duplicate source-owned identifier must not multiply projection rows.
        identifier("chemical_a", "Hmdb:OM:0004", "HMDB00008", "hmdb")
        identifier("chemical_b", "Hmdb:OM:0004", "HMDB0000002")
        identifier("chemical_b", "Chebi:MI:0474", "456")
        identifier("reactome_set", "Name:OM:0202", "Ignored Reactome name")
        identifier("wiki_set", "Name:OM:0202", "Z Wiki name", "wikipathways")
        identifier("wiki_set", "Name:OM:0202", "A Wiki name", "wikipathways")
        cur.execute(
            q(
                "INSERT INTO {s}.chemical_resolution_group_member(level_id,group_key,entity_id,inchikey) VALUES (1,'AAAAAAAAAAAAAA',%s,'AAAAAAAAAAAAAA-BBBBBBBBBB-C')"
            ),
            [ENTITY_IDS["chemical_a"]],
        )
        for entity, label, ontology in (
            ("reactome_set", "Curated Reactome label", "reactome"),
            ("trait", "MACdb trait label", "macdb"),
            ("leaf", "Leaf class", "chemont"),
            ("middle", "Middle class", "chemont"),
            ("root", "Root class", "chemont"),
        ):
            cur.execute(
                q(
                    "INSERT INTO {s}.entity_ontology_term(term_entity_id,term_id,label,ontology_id) VALUES (%s,%s,%s,%s)"
                ),
                [ENTITY_IDS[entity], canonical[entity], label, ontology],
            )
        predicate_names = {r[5] if biolink else r[3] for r in RECORDS}
        predicate_names |= {"part_of", "subclass_of" if biolink else "is_a"}
        for name in predicate_names:
            cur.execute(
                q(
                    "INSERT INTO {s}.vocab_relation_predicate(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
                ),
                [name],
            )
        cur.execute(q("SELECT name,relation_predicate_id FROM {s}.vocab_relation_predicate"))
        predicates = dict(cur.fetchall())
        cur.execute(
            q(
                "INSERT INTO {s}.vocab_relation_category(name) VALUES ('interaction') ON CONFLICT(name) DO NOTHING"
            )
        )
        cur.execute(
            q(
                "SELECT relation_category_id FROM {s}.vocab_relation_category WHERE name='interaction'"
            )
        )
        category = cur.fetchone()[0]
        for (
            source,
            row,
            original_subject,
            original_predicate,
            original_object,
            published_predicate,
            reverse,
        ) in RECORDS:
            sid = SOURCE_IDS[source]
            subject, obj = (
                (original_object, original_subject)
                if biolink and reverse
                else (original_subject, original_object)
            )
            predicate = published_predicate if biolink else original_predicate
            for entity in (subject, obj):
                evidence = _evidence_id(source, row, entity)
                taxonomy = (
                    9606
                    if entity in ("reactome_set", "trait")
                    else 10090
                    if entity == "wiki_set"
                    else None
                )
                cur.execute(
                    q(
                        "INSERT INTO {s}.entity_evidence(source_id,entity_evidence_id,dataset_id,row_id,entity_role_id,entity_type_id,taxonomy_id) VALUES (%s,%s,%s,%s,1,%s,%s)"
                    ),
                    [sid, evidence, sid, row, entity_types[entity], taxonomy],
                )
                status = 5 if biolink else 2 if entity == "reaction" else 1
                cur.execute(
                    q(
                        "INSERT INTO {s}.entity_evidence_resolution(source_id,entity_evidence_id,status_id,entity_id) VALUES (%s,%s,%s,%s)"
                    ),
                    [sid, evidence, status, ENTITY_IDS[entity]],
                )
            evidence = str(
                uuid.uuid5(uuid.NAMESPACE_URL, json.dumps(["metsigdb_record", source, row]))
            )
            # Main evidence endpoints have an exclusive canonical-or-occurrence
            # contract. Published claims already use their canonical endpoints.
            if biolink:
                subject_id, subject_evidence = ENTITY_IDS[subject], None
                object_id, object_evidence = ENTITY_IDS[obj], None
            else:
                subject_occurrence = source == "macdb" or source == "kegg" and row >= 32
                subject_id = None if subject_occurrence else ENTITY_IDS[subject]
                subject_evidence = (
                    _evidence_id(source, row, subject) if subject_occurrence else None
                )
                object_id, object_evidence = None, _evidence_id(source, row, obj)
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence(source_id,relation_evidence_id,dataset_id,row_id,subject_entity_id,subject_entity_evidence_id,predicate_id,object_entity_id,object_entity_evidence_id,relation_category_id) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
                ),
                [
                    sid,
                    evidence,
                    sid,
                    row,
                    subject_id,
                    subject_evidence,
                    predicates[predicate],
                    object_id,
                    object_evidence,
                    category,
                ],
            )
            relation = str(
                uuid.uuid5(
                    uuid.NAMESPACE_URL, json.dumps(["metsigdb_triple", subject, predicate, obj])
                )
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.relation(relation_id,subject_entity_id,predicate_id,object_entity_id,relation_category_id) VALUES (%s,%s,%s,%s,%s) ON CONFLICT DO NOTHING"
                ),
                [relation, ENTITY_IDS[subject], predicates[predicate], ENTITY_IDS[obj], category],
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence_relation(source_id,relation_id,relation_evidence_id) VALUES (%s,%s,%s)"
                ),
                [sid, relation, evidence],
            )
        # MACdb's published trait attribute corresponds to its old structured
        # ontology relation; neither side infers a subtype from the label.
        if biolink:
            key = str(uuid.uuid5(uuid.NAMESPACE_URL, "macdb:trait_type=cancer_trait"))
            cur.execute(
                q(
                    "INSERT INTO {s}.annotation(annotation_key,term,value) VALUES (%s,'macdb:trait_type','cancer_trait')"
                ),
                [key],
            )
            cur.execute(
                q("INSERT INTO {s}.entity_evidence_annotation VALUES (%s,%s,%s)"),
                [SOURCE_IDS["macdb"], _evidence_id("macdb", 40, "trait"), key],
            )
        else:
            cur.execute(
                q("INSERT INTO {s}.entity_ontology_relation VALUES (%s,%s,%s,%s,'macdb')"),
                [
                    SOURCE_IDS["macdb"],
                    ENTITY_IDS["trait"],
                    predicates["part_of"],
                    ENTITY_IDS["trait_type"],
                ],
            )
        for source, row, subject, predicate, obj in ONTOLOGY_RECORDS:
            mapped = "subclass_of" if biolink and predicate == "is_a" else predicate
            cur.execute(
                q("INSERT INTO {s}.entity_ontology_relation VALUES (%s,%s,%s,%s,'chemont')"),
                [SOURCE_IDS[source], ENTITY_IDS[subject], predicates[mapped], ENTITY_IDS[obj]],
            )
    conn.commit()


def _rows(conn, schema):
    from psycopg2 import sql

    with conn.cursor() as cur:
        cur.execute(sql.SQL("SELECT * FROM {}.metsigdb_membership").format(sql.Identifier(schema)))
        names = [d.name for d in cur.description]
        rows = []
        for values in cur.fetchall():
            row = dict(zip(names, values, strict=True))
            # The approved policy admits published chemical entities and their
            # concrete small_molecule subtype without inventing matched flags.
            assert row["metabolite_entity_type"] in {
                "Chemical:OM:0037",
                "chemical_entity",
                "small_molecule",
            }
            row["metabolite_entity_type"] = "published_chemical"
            rows.append(json.dumps(row, sort_keys=True, default=str))
        return Counter(rows)


def _assert_rows(expected, actual):
    missing, extra = expected - actual, actual - expected
    assert not missing and not extra, {
        "missing": list(missing.items())[:3],
        "extra": list(extra.items())[:3],
        "missing_total": sum(missing.values()),
        "extra_total": sum(extra.values()),
    }


def _load(conn, schema, module, rules):
    from psycopg2 import sql

    stats = []
    with conn.cursor() as cur:
        cur.execute(sql.SQL("SET search_path = {},public").format(sql.Identifier(schema)))
    for rule in rules:
        stats.append(module.load_resource(conn, rule, stamp=STAMP, max_records=20))
    return [(s.resource, s.rows, s.sets, s.metabolites, s.removed) for s in stats]


def test_all_five_metsigdb_complete_rows_provenance_and_republication_match_main(
    paired_schemas, oracle, port
):
    conn, (main, adapted) = paired_schemas
    _substrate(conn, main, biolink=False)
    _substrate(conn, adapted, biolink=True)
    reference = importlib.import_module(ORACLE_PREFIX + ".metsigdb.build")
    original_rules = importlib.import_module(ORACLE_PREFIX + ".metsigdb.mapping").RESOURCES
    current = importlib.import_module("omnipath_subsets.metsigdb.build")
    current_rules = importlib.import_module("omnipath_subsets.metsigdb.mapping").RESOURCES
    reference.ensure_membership_table(conn, schema=main)
    current.ensure_membership_table(conn, schema=adapted)
    expected_stats = _load(conn, main, reference, original_rules)
    assert expected_stats == _load(conn, adapted, current, current_rules)
    expected, actual = _rows(conn, main), _rows(conn, adapted)
    _assert_rows(expected, actual)
    rows = [json.loads(row) for row in actual]
    from psycopg2 import sql

    with conn.cursor() as cur:
        cur.execute(
            sql.SQL("SELECT DISTINCT metabolite_entity_type FROM {}.metsigdb_membership").format(
                sql.Identifier(adapted)
            )
        )
        assert {r[0] for r in cur.fetchall()} == {"chemical_entity", "small_molecule"}
        cur.execute(
            sql.SQL("SELECT DISTINCT status_id FROM {}.entity_evidence_resolution").format(
                sql.Identifier(adapted)
            )
        )
        assert {r[0] for r in cur.fetchall()} == {5}
    assert {r["resource"] for r in rows} == {
        "Reactome",
        "WikiPathways",
        "KEGG",
        "MACdb",
        "ClassyFire",
    }
    assert all(r["metabolite_entity_id"] != ENTITY_IDS["protein"] for r in rows)
    assert all(
        r["metabolite_entity_id"] in {ENTITY_IDS["chemical_a"], ENTITY_IDS["chemical_b"]}
        for r in rows
    )
    reactome = [r for r in rows if r["resource"] == "Reactome"]
    assert all(
        r["set_label"] == "Curated Reactome label" and r["organism"] == 9606 and r["set_size"] == 2
        for r in reactome
    )
    assert (
        next(r for r in reactome if r["metabolite_entity_id"] == ENTITY_IDS["chemical_a"])[
            "provenance_record"
        ]["row_id"]
        == 10
    )
    assert all(
        r["set_label"] == "A Wiki name" and r["organism"] == 10090
        for r in rows
        if r["resource"] == "WikiPathways"
    )
    assert all(
        r["set_sub_type"] == "cancer_trait" and r["set_label"] == "MACdb trait label"
        for r in rows
        if r["resource"] == "MACdb"
    )
    kegg = [r for r in rows if r["resource"] == "KEGG"]
    assert {r["set_sub_type"] for r in kegg} == {"metabolic_map", "overview_map"}
    assert all(r["provenance_record"]["via_reaction"] == ENTITY_IDS["reaction"] for r in kegg)
    assert all(
        r["provenance_record"]["row_id"] == 32
        for r in kegg
        if r["metabolite_entity_id"] == ENTITY_IDS["chemical_a"]
    )
    classyfire = [r for r in rows if r["resource"] == "ClassyFire"]
    assert all(r["set_source_id"] != "CHEMONT:O" for r in classyfire)
    root_a = next(
        r
        for r in classyfire
        if r["set_source_id"] == "CHEMONT:R"
        and r["metabolite_entity_id"] == ENTITY_IDS["chemical_a"]
    )
    root_b = next(
        r
        for r in classyfire
        if r["set_source_id"] == "CHEMONT:R"
        and r["metabolite_entity_id"] == ENTITY_IDS["chemical_b"]
    )
    assert root_a["set_context"] == {"assignment": "direct", "depth": 0}
    assert root_b["set_context"] == {"assignment": "ancestor", "depth": 2, "via": "CHEMONT:L"}
    assert (
        root_a["provenance_record"]["row_id"] == 51 and root_b["provenance_record"]["row_id"] == 52
    )
    for row in rows:
        source = "hmdb" if row["resource"] == "ClassyFire" else row["resource"].lower()
        assert row["provenance_source"] == "pypath.inputs_v2." + source + "@" + COMMIT[:12]
        assert row["build_id"] == STAMP
        if row["metabolite_entity_id"] == ENTITY_IDS["chemical_a"]:
            assert (
                row["metabolite_label"],
                row["metabolite_structure_key"],
                row["inchikey"],
                row["smiles"],
                row["hmdb"],
                row["pubchem"],
                row["chebi"],
                row["kegg"],
            ) == (
                "Chemical A",
                "AAAAAAAAAAAAAA",
                "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                "CCO",
                "HMDB0000008",
                "11",
                "123",
                "C00001",
            )
        else:
            assert (
                row["metabolite_label"] == "456"
                and row["metabolite_structure_key"] is None
                and row["hmdb"] == "HMDB0000002"
            )
    # Main's actual table, constraints and filter indexes remain unchanged.
    original_catalogue, current_catalogue = _catalogue(conn, main), _catalogue(conn, adapted)
    for grain in ("columns", "constraints", "indexes"):
        assert Counter(
            {r: n for r, n in original_catalogue[grain].items() if r[0] == "metsigdb_membership"}
        ) == Counter(
            {r: n for r, n in current_catalogue[grain].items() if r[0] == "metsigdb_membership"}
        )
    # Republish the same build: full rows and original record identities stay.
    assert _load(conn, main, reference, original_rules) == _load(
        conn, adapted, current, current_rules
    )
    _assert_rows(expected, _rows(conn, main))
    _assert_rows(expected, _rows(conn, adapted))
    # An upstream pair disappearing under an unchanged stamp removes only that
    # resource's old membership, and refreshes its surviving set_size.
    for schema in (main, adapted):
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL(
                    "DELETE FROM {}.relation_evidence_relation WHERE source_id=%s AND relation_evidence_id=%s"
                ).format(sql.Identifier(schema)),
                [
                    SOURCE_IDS["reactome"],
                    str(
                        uuid.uuid5(
                            uuid.NAMESPACE_URL, json.dumps(["metsigdb_record", "reactome", 12])
                        )
                    ),
                ],
            )
            cur.execute(
                sql.SQL("DELETE FROM {}.relation_evidence WHERE source_id=%s AND row_id=12").format(
                    sql.Identifier(schema)
                ),
                [SOURCE_IDS["reactome"]],
            )
        conn.commit()
    assert (
        _load(conn, main, reference, original_rules[:1])
        == _load(conn, adapted, current, current_rules[:1])
        == [("Reactome", 1, 1, 1, 1)]
    )
    expected_after, actual_after = _rows(conn, main), _rows(conn, adapted)
    _assert_rows(expected_after, actual_after)
    surviving = [json.loads(row) for row in actual_after]
    assert next(r for r in surviving if r["resource"] == "Reactome")["set_size"] == 1
    _assert_rows(
        Counter(row for row in expected if json.loads(row)["resource"] != "Reactome"),
        Counter(row for row in actual_after if json.loads(row)["resource"] != "Reactome"),
    )
