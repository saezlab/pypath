"""Metabolic domain from catalogs, lipid names, ChEBI and linked classes."""

from contextlib import closing
import uuid

import psycopg2
from psycopg2 import sql
import pytest

from omnipath_postgres.relational.classify import classify_metabolic_domain
from omnipath_postgres.relational.classify import metabolic_domain


CHEMICALS = {
    # A class decides.
    "palmitate": [("chemont", "Lipids and lipid-like molecules")],
    "sucrose": [("chemont", "Disaccharides"), ("chemont", "Ketals")],
    # Most of its classes: an alternative parent does not outvote the lineage.
    "ser-gln": [
        ("chemont", "Amino acids, peptides, and analogues"),
        ("chemont", "Dipeptides"),
        ("chemont", "Serine and derivatives"),
        ("chemont", "Fatty acids and conjugates"),
    ],
    # As many lipid as amino acid labels: precedence makes it a lipid.
    "phosphatidylserine": [("chemont", "Glycerophosphoserines")],
    # The class is the subject: an `associated_with` link has no direction.
    "lysine": [("chemont<", "Amino acids, peptides, and analogues")],
    # SwissLipids catalogs lipids only: its entries are lipids with no class
    # link, and its curation outranks a computed class.
    "cardiolipin": [("swisslipids", None)],
    "laccer": [("swisslipids", None), ("chemont", "Purine nucleotides")],
    # A name in lipid shorthand, from any source; a drug code or a trivial
    # name that Goslin would read as a lipid is not one.
    "hmdb species": [("name", "PC(16:0/18:1(9Z))"), ("chemont", "Amino acids")],
    "bortezomib": [("name", "PS-341")],
    "ethanol": [("name", "Ethanol")],
    # ChEBI's hierarchy, through an intermediate class, outranks ChemOnt; two
    # roots tie, and precedence makes a glycolipid a lipid.
    "ctp": [("chebi", "CHEBI:36976"), ("chemont", "Pentose phosphates")],
    "glycolipid": [("chebi", "CHEBI:18059"), ("chebi", "CHEBI:16646")],
    # LIPID MAPS decides only what nothing else places: it lists valine as an
    # amino fatty acid.
    "lmsd lipid": [("lipidmaps", None)],
    "valine": [("lipidmaps", None), ("chemont", "Amino acids, peptides, and analogues")],
    # Any other resource says nothing about the domain, and neither do pathways.
    "imatinib": [("chembl", None)],
    "atp": [("pathway", "Purine metabolism"), ("pathway", "Glycerophospholipid metabolism")],
}


SOURCES = {"swisslipids": 1, "chembl": 2, "lipidmaps": 3}


@pytest.fixture
def domain_database(postgres_dsn):
    schema = "metabolic_domain_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as connection:

        def query(text):
            return sql.SQL(text).format(s=sql.Identifier(schema))

        with connection.cursor() as cursor:
            cursor.execute(query("CREATE SCHEMA {s}"))
            cursor.execute(
                query("CREATE TABLE {s}.vocab_entity_type (entity_type_id bigint, name text)")
            )
            cursor.execute(
                query(
                    "CREATE TABLE {s}.vocab_identifier_type (identifier_type_id bigint, name text)"
                )
            )
            cursor.execute(
                query(
                    "CREATE TABLE {s}.vocab_metabolic_domain "
                    "(metabolic_domain_id smallint PRIMARY KEY, name text, precedence integer)"
                )
            )
            cursor.execute(
                query(
                    "CREATE TABLE {s}.entity (entity_id bigint, entity_type_id bigint, "
                    "canonical_identifier_type_id bigint, canonical_identifier text, "
                    "label text, metabolic_domain_id smallint)"
                )
            )
            cursor.execute(
                query(
                    "CREATE TABLE {s}.relation (subject_entity_id bigint, object_entity_id bigint)"
                )
            )
            cursor.execute(query("CREATE TABLE {s}.data_source (source_id bigint, name text)"))
            cursor.execute(
                query(
                    "CREATE TABLE {s}.entity_evidence_resolution (source_id bigint, entity_id bigint)"
                )
            )
            cursor.execute(query("CREATE TABLE {s}.entity_name (entity_id bigint, name text)"))
            cursor.execute(
                query(
                    "CREATE TABLE {s}.vocab_relation_predicate (relation_predicate_id bigint, name text)"
                )
            )
            cursor.execute(
                query("INSERT INTO {s}.vocab_relation_predicate VALUES (48, 'subclass_of')")
            )
            cursor.execute(
                query(
                    "CREATE TABLE {s}.entity_ontology_relation (source_id bigint, "
                    "subject_entity_id bigint, predicate_id bigint, object_entity_id bigint, "
                    "ontology_id text)"
                )
            )
            cursor.executemany(
                query("INSERT INTO {s}.data_source VALUES (%s, %s)"),
                [(source_id, name) for name, source_id in SOURCES.items()],
            )
            cursor.execute(
                query(
                    "INSERT INTO {s}.vocab_entity_type VALUES "
                    "(1, 'small_molecule'), (2, 'named_thing'), (3, 'pathway')"
                )
            )
            cursor.execute(query("INSERT INTO {s}.vocab_identifier_type VALUES (1, 'chemont')"))
            ids = iter(range(1, 1000))
            roots = {}

            def entity(type_id, namespace_id, identifier, label):
                entity_id = next(ids)
                cursor.execute(
                    query("INSERT INTO {s}.entity VALUES (%s, %s, %s, %s, %s, NULL)"),
                    (entity_id, type_id, namespace_id, identifier, label),
                )
                return entity_id

            def subclass_of(child, parent):
                cursor.execute(
                    query(
                        "INSERT INTO {s}.entity_ontology_relation VALUES (4, %s, 48, %s, 'chebi')"
                    ),
                    (child, parent),
                )

            for label, links in CHEMICALS.items():
                chemical = entity(1, None, None, label)
                for kind, target in links:
                    if kind == "chebi":
                        if target not in roots:
                            roots[target] = entity(1, None, target, f"root {target}")
                        intermediate = entity(1, None, None, f"class under {target}")
                        subclass_of(intermediate, roots[target])
                        subclass_of(chemical, intermediate)
                        continue
                    if kind == "name":
                        cursor.execute(
                            query("INSERT INTO {s}.entity_name VALUES (%s, %s)"), (chemical, target)
                        )
                        continue
                    if kind in SOURCES:
                        cursor.execute(
                            query("INSERT INTO {s}.entity_evidence_resolution VALUES (%s, %s)"),
                            (SOURCES[kind], chemical),
                        )
                        continue
                    linked = entity(
                        3 if kind == "pathway" else 2,
                        None if kind == "pathway" else 1,
                        None,
                        target,
                    )
                    pair = (linked, chemical) if kind == "chemont<" else (chemical, linked)
                    cursor.execute(query("INSERT INTO {s}.relation VALUES (%s, %s)"), pair)
        connection.commit()
        try:
            yield connection, schema
        finally:
            with connection.cursor() as cursor:
                cursor.execute(query("DROP SCHEMA {s} CASCADE"))
            connection.commit()


@pytest.mark.parametrize("chunk", [5000, 1], ids=["inline", "process pool"])
def test_sources_decide_in_priority_order(domain_database, monkeypatch, chunk):
    connection, schema = domain_database
    monkeypatch.setattr(metabolic_domain, "_PARSE_CHUNK", chunk)
    classify_metabolic_domain(connection, schema=schema)
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL(
                "SELECT e.label, d.name FROM {s}.entity e "
                "JOIN {s}.vocab_metabolic_domain d USING (metabolic_domain_id)"
            ).format(s=sql.Identifier(schema))
        )
        domains = {label: name for label, name in cursor.fetchall() if label in CHEMICALS}
    assert domains == {
        "palmitate": "lipid",
        "sucrose": "carbohydrate",
        "ser-gln": "amino_acid",
        "phosphatidylserine": "lipid",
        "lysine": "amino_acid",
        "cardiolipin": "lipid",
        "laccer": "lipid",
        "hmdb species": "lipid",
        "bortezomib": "other",
        "ethanol": "other",
        "ctp": "nucleotide",
        "glycolipid": "lipid",
        "lmsd lipid": "lipid",
        "valine": "amino_acid",
        "imatinib": "other",
        "atp": "other",
    }
