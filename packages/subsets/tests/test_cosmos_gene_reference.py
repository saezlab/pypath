"""A GeneID catalyst is a reference, not a selection of a protein product."""

from contextlib import closing
import uuid

import psycopg2
from psycopg2 import sql
import pytest

from omnipath_subsets.cosmos import build, translate
from omnipath_subsets.metsigdb import build as metsigdb


def catalyst(identifier, id_type="Entrez:MI:0477"):
    return translate.LabelEntity(str(uuid.uuid4()), identifier, id_type, 9606, "gene")


class Mappings:
    failed = False

    def __init__(self, count):
        self.count = count
        self.calls = []

    def resolve_proteins(self, source_type, taxonomy, identifiers):
        self.calls.append((source_type, taxonomy, identifiers))
        return {identifier: ("P00001", self.count) for identifier in identifiers}

    def map_to(self, source_type, target_type, taxonomy, identifiers):
        assert not identifiers
        return {}

    def known_chemicals(self, source_type, identifiers):
        assert not identifiers
        return set()


@pytest.mark.parametrize("product_count", [1, 2])
def test_optional_utils_does_not_select_a_product_for_geneid(product_count):
    entity = catalyst("101")
    mappings = Mappings(product_count)
    labels = translate.translate_identifiers([entity], mappings)
    assert labels[entity.entity_id] == translate.TranslatedLabel("101", "entrez", False)
    assert mappings.calls == []


@pytest.mark.parametrize("target_count", [1, 2])
def test_optional_other_catalyst_mapping_requires_one_target(target_count):
    entity = catalyst("HGNC:1", "Hgnc:MI:1095")
    labels = translate.translate_identifiers([entity], Mappings(target_count))
    expected = (
        translate.TranslatedLabel("P00001", "uniprot", True)
        if target_count == 1
        else translate.TranslatedLabel("HGNC:1", "hgnc", False)
    )
    assert labels[entity.entity_id] == expected


def test_canonical_uniprot_passes_through_optional_utils():
    entity = catalyst("P99999", "Uniprot:MI:1097")
    mappings = Mappings(2)
    assert translate.translate_identifiers([entity], mappings)[entity.entity_id] == (
        translate.TranslatedLabel("P99999", "uniprot", False)
    )
    assert mappings.calls == []


# This bounded fixture exposes only the normalized product inputs. It needs no
# published-provenance copies, gene-product catalogue or reference database.
@pytest.fixture
def product_input(postgres_dsn):
    schema = "subset_gene_reference_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            with connection.cursor() as cur:
                cur.execute(sql.SQL("CREATE SCHEMA {}").format(sql.Identifier(schema)))
                cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(schema)))
                cur.execute(
                    """
                    CREATE TABLE vocab_entity_type(entity_type_id int, name text);
                    CREATE TABLE vocab_identifier_type(identifier_type_id int, name text);
                    CREATE TABLE vocab_relation_role(relation_role_id int, name text);
                    CREATE TABLE vocab_relation_predicate(relation_predicate_id int, name text);
                    CREATE TABLE entity(entity_id uuid, entity_type_id int,
                      canonical_identifier_type_id int, canonical_identifier text,
                      taxonomy_id int);
                    CREATE TABLE identifier_evidence(identifier_id int, identifier_type_id int,
                      value text);
                    CREATE TABLE entity_identifier_lookup(entity_id uuid, identifier_id int);
                    CREATE TABLE interaction(interaction_id uuid, sources text[]);
                    CREATE TABLE interaction_party(interaction_id uuid, entity_id uuid,
                      role_id int, compartment text);
                    CREATE TABLE interaction_fact_resource(interaction_id uuid,
                      subject_entity_id uuid);
                    CREATE TABLE relation_evidence(source_id int, relation_evidence_id int,
                      predicate_id int, subject_entity_id uuid, object_entity_id uuid,
                      dataset_id int, row_id text);
                    CREATE TABLE relation_evidence_annotation(source_id int,
                      relation_evidence_id int, annotation_key text);
                    CREATE TABLE annotation(annotation_key text, term text, value text);
                    CREATE TABLE entity_evidence_annotation(source_id int,
                      entity_evidence_id int, annotation_key text);
                    CREATE TABLE entity_evidence_resolution(source_id int,
                      entity_evidence_id int, entity_id uuid);
                    INSERT INTO vocab_entity_type VALUES
                      (1,'gene'),(2,'protein'),(3,'small_molecule'),
                      (4,'molecular_activity'),(5,'pathway');
                    INSERT INTO vocab_identifier_type VALUES
                      (1,'Entrez:MI:0477'),(2,'Uniprot:MI:1097'),(3,'Chebi:MI:0474');
                    INSERT INTO vocab_relation_role VALUES
                      (1,'reactant'),(2,'product'),(3,'enzyme');
                    INSERT INTO vocab_relation_predicate VALUES (1,'part_of');
                    """
                )
                ids = {
                    name: str(uuid.UUID(int=index))
                    for index, name in enumerate(
                        (
                            "gene",
                            "reported_protein",
                            "product",
                            "input",
                            "output",
                            "event",
                            "interaction",
                            "pathway",
                        ),
                        start=1,
                    )
                }
                cur.executemany(
                    "INSERT INTO entity VALUES (%s,%s,%s,%s,9606)",
                    [
                        (ids["gene"], 1, 1, "101"),
                        (ids["reported_protein"], 2, 1, "102"),
                        (ids["product"], 2, 2, "P99999"),
                        (ids["input"], 3, 3, "15377"),
                        (ids["output"], 3, 3, "15422"),
                        (ids["event"], 4, None, "reaction"),
                        (ids["pathway"], 5, None, "pathway"),
                    ],
                )
                # Product aliases sort before the canonical accession. Both
                # reported entity types must retain their GeneID identity.
                cur.execute("INSERT INTO identifier_evidence VALUES (1,2,'P00001'),(2,2,'P99998')")
                cur.executemany(
                    "INSERT INTO entity_identifier_lookup VALUES (%s,%s)",
                    [
                        (ids[name], alias)
                        for name in ("gene", "reported_protein", "product")
                        for alias in (1, 2)
                    ],
                )
                cur.execute(
                    "INSERT INTO interaction VALUES (%s,ARRAY['rhea'])", [ids["interaction"]]
                )
                cur.execute(
                    "INSERT INTO interaction_fact_resource VALUES (%s,%s)",
                    [ids["interaction"], ids["event"]],
                )
                cur.executemany(
                    "INSERT INTO interaction_party VALUES (%s,%s,%s,%s)",
                    [
                        (ids["interaction"], ids[name], role, compartment)
                        for name, role, compartment in (
                            ("gene", 3, None),
                            ("reported_protein", 3, None),
                            ("product", 3, None),
                            ("input", 1, "cytosol"),
                            ("output", 2, "cytosol"),
                        )
                    ],
                )
                cur.execute(
                    "INSERT INTO annotation VALUES ('direction',%s,'reversible')",
                    [build.DIRECTION_TERM],
                )
                cur.execute("INSERT INTO entity_evidence_annotation VALUES (1,1,'direction')")
                cur.execute(
                    "INSERT INTO entity_evidence_resolution VALUES (1,1,%s)", [ids["event"]]
                )
                cur.executemany(
                    "INSERT INTO relation_evidence VALUES (1,%s,1,%s,%s,1,%s)",
                    [
                        (index, ids[name], ids["pathway"], str(index))
                        for index, name in enumerate(
                            ("input", "gene", "reported_protein", "product")
                        )
                    ],
                )
                yield cur, ids
        finally:
            # Roll back the whole synthetic schema; never touch a deployed schema.
            connection.rollback()


@pytest.mark.integration
def test_published_geneid_labels_and_connectors_preserve_reaction_participants(product_input):
    cur, ids = product_input
    stats = translate.stage_label_translation(cur)
    assert (stats.gene_entities, stats.gene_in_namespace, stats.gene_mapped) == (3, 1, 0)
    cur.execute("SELECT entity_id::text,identifier,id_namespace FROM _cos_label")
    labels = {
        entity_id: (identifier, namespace) for entity_id, identifier, namespace in cur.fetchall()
    }
    assert labels[ids["gene"]] == ("101", "entrez")
    assert labels[ids["reported_protein"]] == ("102", "entrez")
    assert labels[ids["product"]] == ("P99999", "uniprot")
    assert labels[ids["input"]] == ("CHEBI:15377", "chebi")

    cur.execute(build._sql_text("cosmos_edge_table.sql"))
    cur.execute(
        build._sql_text("project_edges.sql"),
        {
            "build_id": "fixture",
            "reaction_entity_types": list(build.REACTION_ENTITY_TYPES),
            "reaction_sources": list(build.REACTION_SOURCES),
            "direction_term": build.DIRECTION_TERM,
            "reaction_type_term": build.REACTION_TYPE_TERM,
            "transport_types": build.TRANSPORT_TYPES,
        },
    )
    cur.execute(build._sql_text("project_connectors.sql"), {"build_id": "fixture"})
    cur.execute(
        "SELECT source_label,target_label,source_id_type,target_id_type,reverse,direction "
        "FROM cosmos_edge WHERE interaction_type='connector'"
    )
    connectors = set(cur.fetchall())
    assert connectors == {
        (identifier, f"Gene1__{identifier}" + suffix, namespace, namespace, reverse, "reversible")
        for identifier, namespace in (("101", "entrez"), ("102", "entrez"), ("P99999", "uniprot"))
        for suffix, reverse in (("", False), ("_rev", True))
    }
    cur.execute(
        "SELECT source_label,target_label,source_entity_id::text,target_entity_id::text,"
        "source_compartment,target_compartment FROM cosmos_edge "
        "WHERE interaction_type='catalysis'"
    )
    edges = cur.fetchall()
    assert len(edges) == 12  # Two chemical participants, three catalysts, two directions.
    chemical_nodes = {"Metab__CHEBI:15377_cytosol", "Metab__CHEBI:15422_cytosol"}
    gene_nodes = {row[1] for row in connectors}
    for (
        source,
        target,
        source_entity,
        target_entity,
        source_compartment,
        target_compartment,
    ) in edges:
        assert {source, target} & chemical_nodes
        assert {source, target} & gene_nodes
        chemical = source if source in chemical_nodes else target
        expected_id = ids["input"] if "15377" in chemical else ids["output"]
        assert (source_entity if source == chemical else target_entity) == expected_id
        assert (source_compartment if source == chemical else target_compartment) == "cytosol"
    assert all("P00001" not in source + target for source, target, *_ in edges)


@pytest.mark.integration
def test_metsigdb_memberships_remain_chemical_only_with_gene_reference_catalysts(product_input):
    cur, ids = product_input
    params = {
        "source_id": 1,
        "set_entity_type_id": 5,
        "chemical_entity_type_ids": metsigdb._chemical_entity_type_ids(cur),
        "is_pathway": True,
        "max_records": None,
    }

    def memberships():
        cur.execute(metsigdb._sql_text("extract_onehop.sql"), params)
        cur.execute(
            "SELECT set_source_id,set_entity_id::text,metabolite_entity_id::text "
            "FROM metsigdb_stage"
        )
        return cur.fetchall()

    before = memberships()
    assert before == [("pathway", ids["pathway"], ids["input"])]
    cur.execute("DELETE FROM relation_evidence WHERE subject_entity_id<>%s", [ids["input"]])
    assert memberships() == before
