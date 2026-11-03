"""Small published fixtures exercise main's relational COPY contracts."""

from copy import deepcopy
from types import SimpleNamespace
import json
import uuid

import duckdb
import pytest

from omnipath_postgres.projection import (
    _prepare_dimensions,
    companion_ddl,
    prepare_aligned_release,
)
from published_fixture import ENTITY_A, ENTITY_B, QUANTITY, annotation, fixture_rows, write_fixture


def selected(path, source="fixture", version="1"):
    return SimpleNamespace(directory=path, source=source, version=version)


def rows(con, plan, table):
    query = next(item for item in plan.queries if item.table == table)
    return [
        dict(zip(query.columns, row, strict=True)) for row in con.execute(query.query).fetchall()
    ]


def test_main_base_and_narrow_companions_preserve_published_occurrences(tmp_path):
    write_fixture(tmp_path)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        assert plan.counts["entity"] == 3
        assert plan.counts["identifier_evidence"] == 4
        assert plan.counts["entity_evidence"] == 3
        assert plan.counts["relation"] == 1
        assert plan.counts["relation_evidence"] == 2
        assert plan.counts["parquet_identifier_occurrence"] == 2
        assert plan.counts["parquet_annotation_occurrence"] == 9
        assert (
            len({r["identifier_id"] for r in rows(con, plan, "parquet_identifier_occurrence")}) == 1
        )
        assert len({r["relation_evidence_id"] for r in rows(con, plan, "relation_evidence")}) == 2
        assert all(
            r["subject_entity_id"] and r["object_entity_id"]
            for r in rows(con, plan, "relation_evidence")
        )
        assert all(
            r["subject_entity_evidence_id"] is None for r in rows(con, plan, "relation_evidence")
        )
        assert all(r["resolution_status_id"] == 5 for r in rows(con, plan, "entity"))
        assert all(r["status_id"] == 5 for r in rows(con, plan, "entity_evidence_resolution"))
        assert [r["original_row_id"] for r in rows(con, plan, "parquet_evidence")] == [
            "00001",
            "00001",
        ]
        assert all(r["row_id"] == 1 for r in rows(con, plan, "relation_evidence"))
        assert all(isinstance(r["entity_id"], uuid.UUID) for r in rows(con, plan, "entity"))
        assert not any(
            "record_json" in q.columns or "payload_json" in q.columns for q in plan.queries
        )


def test_safe_bare_aliases_preserve_canonical_values_and_original_occurrences(tmp_path):
    entities, relations, _ = fixture_rows()
    entities[0].update(namespace="chebi", identifier="CHEBI:15377")
    entities[0]["identifiers"] = [
        dict(ns="chebi", id="CHEBI:15377", is_canonical=True, source="fixture"),
        dict(ns="chebi", id="15377", is_canonical=False, source="fixture"),
        dict(ns="kegg_reaction", id="rn:R00001", is_canonical=False, source="fixture"),
    ]
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        occurrence = next(
            r for r in rows(con, plan, "parquet_entity") if r["entity_key"] == ENTITY_A
        )
        canonical = next(
            r for r in rows(con, plan, "entity") if r["entity_id"] == occurrence["entity_id"]
        )
        assert canonical["canonical_identifier"] == occurrence["identifier"] == "CHEBI:15377"
        dictionary = {r["identifier_id"]: r for r in rows(con, plan, "identifier_evidence")}
        links = [
            r
            for r in rows(con, plan, "entity_identifier")
            if r["entity_id"] == occurrence["entity_id"]
        ]
        assert {dictionary[r["identifier_id"]]["value"] for r in links} == {
            "CHEBI:15377",
            "15377",
            "rn:R00001",
            "R00001",
        }
        evidence_links = [
            r
            for r in rows(con, plan, "entity_evidence_identifier")
            if r["entity_evidence_id"] == occurrence["entity_evidence_id"]
        ]
        assert {dictionary[r["identifier_id"]]["value"] for r in evidence_links} == {
            "CHEBI:15377",
            "15377",
            "rn:R00001",
        }
        original = sorted(
            rows(con, plan, "parquet_identifier_occurrence"), key=lambda r: r["ordinal"]
        )
        assert [r["identifier"] for r in original] == ["CHEBI:15377", "15377", "rn:R00001"]
        assert [r["is_canonical"] for r in original] == [True, False, False]
        assert plan.compatibility["safe_bare_identifier_alias_pairs"] == 2
        assert "Kegg Reaction:MI:2013" in {
            name for _, name in plan.dimensions["vocab_identifier_type"]
        }


def test_aliases_only_strip_known_outer_prefixes_with_valid_suffixes(tmp_path):
    entities, relations, _ = fixture_rows()
    inputs = [
        ("hmdb", "hmdb:HMDB0000001", "HMDB0000001"),
        ("chembl", "chembl:CHEMBL1", "CHEMBL1"),
        ("swisslipids", "swisslipids:SLM:000000001", "SLM:000000001"),
        ("kegg", "cpd:C00001", "C00001"),
        ("pubchem", "CID:123", "123"),
        ("chebi", "CHEBI:not-a-number", None),
        ("unknown", "CHEBI:123", None),
        ("name", "CHEBI:123", None),
        ("swisslipids", "SLM:000000002", None),
        ("kegg_reaction", "cpd:C00002", None),
    ]
    entities[0]["identifiers"] = [
        dict(ns=ns, id=value, is_canonical=False, source="fixture") for ns, value, _ in inputs
    ]
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        actual = con.execute("""SELECT original.value,alias.alias_value
            FROM ap_identifier_alias alias JOIN ap_identifier original
            ON original.identifier_id=alias.original_identifier_id""").fetchall()
        assert set(actual) == {(value, bare) for _, value, bare in inputs if bare is not None} | {
            ("CHEBI:15377", "15377"),
        }
        assert plan.compatibility["safe_bare_identifier_alias_pairs"] == 6


def test_quantity_identity_includes_prefix_comparator_and_preserves_original_value(tmp_path):
    entities, relations, _ = fixture_rows()
    relations[0]["annotations"] = []
    relations[0]["evidence"] = []
    quantity1 = dict(QUANTITY)
    quantity2 = dict(QUANTITY, comparator=">", has_unit_prefix="micro")
    value = json.dumps(quantity1)
    a1 = annotation("has_quantitative_value", value=value, quantity=quantity1)
    a2 = annotation("has_quantitative_value", value=value, quantity=quantity2)
    relations[0]["annotations"] = [a1, deepcopy(a1), a2]
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        values = [r for r in rows(con, plan, "annotation") if r["term"] == "has_quantitative_value"]
        assert len(values) == 2
        assert {r["value"] for r in values} == {"12.5"}
        quantities = {r["annotation_key"]: r for r in rows(con, plan, "annotation_quantity")}
        assert {quantities[r["annotation_key"]]["comparator"] for r in values} == {"<=", ">"}
        assert {quantities[r["annotation_key"]]["published_value"] for r in values} == {value}
        occurrences = [
            r
            for r in rows(con, plan, "parquet_annotation_occurrence")
            if r["owner_kind"] == "relation"
        ]
        assert [r["ordinal"] for r in occurrences] == [0, 1, 2]


def test_aggregated_statement_annotations_do_not_contaminate_other_evidence_rows(tmp_path):
    entities, relations, _ = fixture_rows()
    forward = annotation("biopax:conversionDirection", value="LEFT-TO-RIGHT", scope="relation")
    reverse = annotation("biopax:conversionDirection", value="RIGHT-TO-LEFT", scope="relation")
    compartment1 = annotation("biopax:cellularLocation", value="cytosol", scope="object")
    compartment2 = annotation("biopax:cellularLocation", value="extracellular", scope="object")
    relations[0]["evidence"][0].update(row_id="first-event", annotations=[forward, compartment1])
    relations[0]["evidence"][1].update(row_id="second-event", annotations=[reverse, compartment2])
    # Writer's relation annotation array is the union of the occurrence arrays.
    relations[0]["annotations"] = [
        deepcopy(a) for a in [forward, reverse, compartment1, compartment2]
    ]
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        values = {r["annotation_key"]: r["value"] for r in rows(con, plan, "annotation")}
        ownership = {
            r["relation_evidence_id"]: r["ordinal"] for r in rows(con, plan, "parquet_evidence")
        }
        actual = {0: set(), 1: set()}
        for row in rows(con, plan, "relation_evidence_annotation"):
            actual[ownership[row["relation_evidence_id"]]].add(values[row["annotation_key"]])
        assert actual == {0: {"LEFT-TO-RIGHT", "cytosol"}, 1: {"RIGHT-TO-LEFT", "extracellular"}}


def test_true_statement_annotation_respects_declared_source_and_dataset(tmp_path):
    entities, relations, _ = fixture_rows()
    relations[0]["evidence"][0]["annotations"] = []
    relations[0]["evidence"][1].update(
        source="other-source", dataset="other-dataset", annotations=[]
    )
    relations[0]["annotations"] = [
        annotation("description", value="source owned", scope="relation")
    ]
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        links = rows(con, plan, "relation_evidence_annotation")
        assert len(links) == 1
        evidence = {r["relation_evidence_id"]: r for r in rows(con, plan, "parquet_evidence")}
        assert evidence[links[0]["relation_evidence_id"]]["source"] == "reported-source"


def test_shared_entities_and_qualified_statements_collapse_only_main_graph_triple(tmp_path):
    entities, relations, _ = fixture_rows()
    r2 = deepcopy(relations[0])
    r2.update(relation_key="qualified-two", sign=1)
    r2["annotations"] = [
        annotation("object_direction_qualifier", value="increased", scope="object")
    ]
    write_fixture(tmp_path / "one", entities=entities, relations=relations + [r2])
    write_fixture(tmp_path / "two", entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(
            con,
            (selected(tmp_path / "one", "a"), selected(tmp_path / "two", "b")),
            retain_published_provenance=True,
        )
        assert plan.counts["entity"] == 3
        assert plan.counts["parquet_entity"] == 6
        assert plan.counts["relation"] == 1
        assert plan.counts["parquet_statement"] == 3
        assert plan.counts["relation_evidence"] == 6
        assert len({r["relation_id"] for r in rows(con, plan, "parquet_statement")}) == 1
        assert {r["sign"] for r in rows(con, plan, "parquet_statement")} == {-1, 1}
        assert {r["scope"] for r in rows(con, plan, "parquet_annotation_occurrence")} >= {
            "object",
            "relation",
            "evidence",
        }


def test_shared_statement_key_preserves_each_resources_taxonomy_consensus_and_evidence(tmp_path):
    entities, relations, _ = fixture_rows()
    for source, taxon in (("one", "9606"), ("two", "10090")):
        scoped = deepcopy(relations)
        scoped[0]["taxon"] = taxon
        taxon_annotation = annotation("in_taxon", value="NCBITaxon:" + taxon, scope="relation")
        taxon_annotation["source"] = source
        scoped[0]["annotations"].append(deepcopy(taxon_annotation))
        for evidence in scoped[0]["evidence"]:
            evidence["source"] = source
            evidence["annotations"].append(deepcopy(taxon_annotation))
        write_fixture(tmp_path / source, entities=entities, relations=scoped)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(
            con,
            (selected(tmp_path / "one", "one"), selected(tmp_path / "two", "two")),
            retain_published_provenance=True,
        )
        claims = rows(con, plan, "parquet_statement")
        assert {(r["resource"], r["taxon"]) for r in claims} == {("one", "9606"), ("two", "10090")}
        assert len({r["relation_key"] for r in claims}) == 1
        assert len({r["relation_id"] for r in claims}) == 1
        assert plan.counts["relation_evidence"] == 4
        assert {(r["resource"], r["source"]) for r in rows(con, plan, "parquet_evidence")} == {
            ("one", "one"),
            ("two", "two"),
        }
        values = {r["annotation_key"]: r["value"] for r in rows(con, plan, "annotation")}
        assert {
            (r["resource"], values[r["annotation_key"]])
            for r in rows(con, plan, "parquet_annotation_occurrence")
            if r["term"] == "in_taxon"
        } == {
            ("one", "NCBITaxon:9606"),
            ("two", "NCBITaxon:10090"),
        }


@pytest.mark.parametrize(
    "field,value",
    [
        ("subject_entity_key", ENTITY_B),
        ("object_entity_key", ENTITY_A),
        ("predicate", "interacts_with"),
        ("statement_kind", "ontology_axiom"),
    ],
)
def test_shared_statement_key_with_conflicting_identity_core_still_fails(tmp_path, field, value):
    entities, relations, _ = fixture_rows()
    write_fixture(tmp_path / "one", entities=entities, relations=relations)
    relations[0][field] = value
    write_fixture(tmp_path / "two", entities=entities, relations=relations)
    with (
        duckdb.connect() as con,
        pytest.raises(ValueError, match="statement identity has conflicting"),
    ):
        prepare_aligned_release(
            con,
            (selected(tmp_path / "one", "one"), selected(tmp_path / "two", "two")),
            retain_published_provenance=True,
        )


def test_distinct_published_entities_are_not_silently_merged_for_main_unique_index(tmp_path):
    entities, relations, _ = fixture_rows()
    duplicate = deepcopy(entities[0])
    duplicate["entity_key"] = "distinct-published-key-same-main-natural-tuple"
    write_fixture(tmp_path, entities=entities + [duplicate], relations=relations)
    with duckdb.connect() as con, pytest.raises(ValueError, match="natural-key uniqueness"):
        prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)


def test_shared_published_key_with_conflicting_canonical_fields_is_rejected(tmp_path):
    entities, relations, _ = fixture_rows()
    write_fixture(tmp_path / "one", entities=entities, relations=relations)
    entities[0]["identifier"] = "DIFFERENT"
    write_fixture(tmp_path / "two", entities=entities, relations=relations)
    with duckdb.connect() as con, pytest.raises(ValueError, match="conflicting canonical"):
        prepare_aligned_release(
            con,
            (selected(tmp_path / "one", "a"), selected(tmp_path / "two", "b")),
            retain_published_provenance=True,
        )


def test_shared_identity_taxonomy_conflicts_keep_each_source_occurrence(tmp_path):
    entities, relations, _ = fixture_rows()
    write_fixture(tmp_path / "one", entities=entities, relations=relations)
    entities[0]["taxon"] = "10090"
    write_fixture(tmp_path / "two", entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(
            con,
            (selected(tmp_path / "one", "a"), selected(tmp_path / "two", "b")),
            retain_published_provenance=True,
        )
        published = [r for r in rows(con, plan, "parquet_entity") if r["entity_key"] == ENTITY_A]
        canonical = {r["entity_id"]: r for r in rows(con, plan, "entity")}
        assert {r["taxon"] for r in published} == {"9606", "10090"}
        assert len({r["entity_id"] for r in published}) == 1
        assert canonical[published[0]["entity_id"]]["taxonomy_id"] is None
        assert plan.compatibility["shared_identity_taxonomy_conflicts"] == 1


def test_writer_empty_unknown_taxon_is_not_fabricated_or_rejected(tmp_path):
    entities, relations, _ = fixture_rows()
    entities[0]["taxon"] = ""
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        published = next(
            r for r in rows(con, plan, "parquet_entity") if r["entity_key"] == ENTITY_A
        )
        canonical = next(
            r for r in rows(con, plan, "entity") if r["entity_id"] == published["entity_id"]
        )
        assert published["taxon"] == ""
        assert canonical["taxonomy_id"] is None


def test_source_scoped_names_keep_distinct_uuid_and_main_unique_canonical_tuple(tmp_path):
    entities, relations, _ = fixture_rows()
    entities[0].update(namespace="name", identifier="same name")
    duplicate = deepcopy(entities[0])
    duplicate["entity_key"] = "source-scoped-name-two"
    write_fixture(tmp_path, entities=entities + [duplicate], relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        canonical = {r["entity_id"]: r for r in rows(con, plan, "entity")}
        published = [r for r in rows(con, plan, "parquet_entity") if r["namespace"] == "name"]
        assert len({r["entity_id"] for r in published}) == 2
        assert {canonical[r["entity_id"]]["canonical_identifier"] for r in published} == {
            ENTITY_A,
            "source-scoped-name-two",
        }
        assert {r["identifier"] for r in published} == {"same name"}
        assert plan.compatibility["source_scoped_name_fallback_entities"] == 2
        dictionary = {r["identifier_id"]: r for r in rows(con, plan, "identifier_evidence")}
        linked = rows(con, plan, "entity_identifier")
        for occurrence in published:
            aliases = {
                dictionary[r["identifier_id"]]["value"]
                for r in linked
                if r["entity_id"] == occurrence["entity_id"]
            }
            assert "same name" in aliases
            assert occurrence["entity_key"] in aliases


def test_base_ontology_terms_stays_empty_and_relational_inputs_preserve_synonyms(tmp_path):
    entities, relations, _ = fixture_rows()
    entities[0].update(namespace="go", identifier="GO:0000001", entity_type="ontology_class")
    entities[0]["identifiers"] = [
        dict(ns="synonym", id="first synonym", is_canonical=False, source="go")
    ]
    entities[1]["entity_type"] = "ontology_class"  # Protein identity is not a GO term.
    relations[0].update(statement_kind="ontology_axiom", predicate="subclass_of")
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(
            con, (selected(tmp_path, "go"),), retain_published_provenance=True
        )
        assert plan.counts["ontology_terms"] == 0
        assert rows(con, plan, "ontology_terms") == []
        assert "first synonym" in {
            row["identifier"] for row in rows(con, plan, "parquet_identifier_occurrence")
        }
        assert rows(con, plan, "entity_ontology_relation")[0]["ontology_id"] == "gene_ontology"


def test_source_owned_row_context_surrogates_and_original_strings_are_recoverable(tmp_path):
    entities, relations, _ = fixture_rows()
    relations[0]["evidence"][0]["row_id"] = "RXN:source:'α"
    relations[0]["evidence"][1]["row_id"] = "RXN:source:'α"
    relations[0]["evidence"][1]["source"] = "different-source"
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        evidence = rows(con, plan, "relation_evidence")
        provenance = rows(con, plan, "parquet_evidence")
        assert len({r["source_id"] for r in evidence}) == 2
        assert all(isinstance(r["row_id"], int) for r in evidence)
        assert {r["original_row_id"] for r in provenance} == {"RXN:source:'α"}
        assert len({r["row_id"] for r in evidence}) == 2


def test_null_empty_lists_and_all_null_structs_are_retained_without_fabricated_dictionary_values(
    tmp_path,
):
    entities, relations, _ = fixture_rows()
    entities[0]["identifiers"].append(dict(ns=None, id=None, is_canonical=None, source=None))
    entities[0]["annotations"].append(
        dict(term=None, value=None, quantity=dict.fromkeys(QUANTITY), source=None, dataset=None)
    )
    relations[0]["evidence"] = None
    relations[0]["sources"] = ["quoted,'\\\"", None, "", "quoted,'\\\""]
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        assert plan.counts["relation_evidence"] == 1
        assert rows(con, plan, "parquet_evidence")[0]["synthetic"] is True
        occurrence = rows(con, plan, "parquet_identifier_occurrence")[-1]
        assert occurrence["identifier_id"] is None
        assert any(r["term"] is None for r in rows(con, plan, "parquet_annotation_occurrence"))
        provenance = {r["entity_key"]: r for r in rows(con, plan, "parquet_entity")}
        assert provenance[ENTITY_A]["identifiers_present"] is True
        assert provenance[ENTITY_B]["identifiers_present"] is False
        assert rows(con, plan, "parquet_statement")[0]["evidence_present"] is False
        assert "NULL" in rows(con, plan, "parquet_statement")[0]["sources"]


def test_ontology_axioms_do_not_enter_main_graph_relation(tmp_path):
    entities, relations, _ = fixture_rows()
    relations[0].update(statement_kind="ontology_axiom", predicate="subclass_of")
    write_fixture(tmp_path, entities=entities, relations=relations)
    with duckdb.connect() as con:
        plan = prepare_aligned_release(con, (selected(tmp_path),), retain_published_provenance=True)
        assert plan.counts["relation"] == 0
        assert plan.counts["entity_ontology_relation"] == 1
        assert rows(con, plan, "parquet_statement")[0]["relation_id"] is None


def test_dimensions_reuse_existing_database_ids_and_input_order_is_irrelevant(tmp_path):
    write_fixture(tmp_path / "one")
    write_fixture(tmp_path / "two")
    existing = {
        "data_source": ((91, "fixture"), (123, "reported-source")),
        "vocab_relation_predicate": ((1001, "affects"),),
    }
    with duckdb.connect() as con:
        plan = prepare_aligned_release(
            con,
            (selected(tmp_path / "one"),),
            dimension_rows=existing,
            retain_published_provenance=True,
        )
        assert (
            dict((name, identifier) for identifier, name in plan.sources)["reported-source"] == 123
        )
        assert rows(con, plan, "relation")[0]["predicate_id"] == 1001
        first = sorted(rows(con, plan, "entity"), key=lambda row: row["entity_id"])
        second = prepare_aligned_release(
            con,
            (selected(tmp_path / "two"),),
            dimension_rows=existing,
            retain_published_provenance=True,
        )
        assert sorted(rows(con, second, "entity"), key=lambda row: row["entity_id"]) == first


def test_dimension_collection_fetches_only_distinct_names_and_preserves_existing_ids():
    class ObservedConnection:
        prefix = "SELECT DISTINCT dimension_name FROM ("
        suffix = ") required(dimension_name) WHERE dimension_name IS NOT NULL"

        def __init__(self, connection, *, emulate_previous=False):
            self.connection = connection
            self.emulate_previous = emulate_previous
            self.name_fetches = []
            self.dataset_fetches = []

        def execute(self, statement):
            self.name_query = statement.startswith(self.prefix)
            self.dataset_query = statement.startswith("SELECT resource,")
            if self.name_query and self.emulate_previous:
                assert statement.endswith(self.suffix)
                statement = statement[len(self.prefix) : -len(self.suffix)]
            self.connection.execute(statement)
            return self

        def executemany(self, *args):
            return self.connection.executemany(*args)

        def fetchall(self):
            result = self.connection.fetchall()
            if self.name_query:
                self.name_fetches.append(result)
            if self.dataset_query:
                self.dataset_fetches.append(result)
            return result

    # Ten thousand repeated rows per staged input exercise the former Python
    # allocation issue without a resource build or a full publication fixture.
    with duckdb.connect() as con:
        con.execute("""CREATE TABLE ap_entity_raw AS SELECT
            CASE WHEN range%2=0 THEN 'one' ELSE 'two' END resource,
            CASE WHEN range%2=0 THEN 'protein' ELSE 'small_molecule' END entity_type,
            CASE WHEN range%2=0 THEN 'uniprot' ELSE 'chebi' END namespace
            FROM range(10000)""")
        con.execute("""CREATE TABLE ap_statement_raw AS SELECT
            CASE WHEN range%2=0 THEN 'one' ELSE 'two' END resource,
            CASE WHEN range%2=0 THEN 'affects' ELSE 'interacts_with' END predicate,
            CASE WHEN range%2=0 THEN 'interaction' ELSE NULL END category
            FROM range(10000)""")
        con.execute("""CREATE TABLE ap_identifier_raw AS SELECT
            struct_pack(ns:=CASE WHEN range%2=0 THEN 'uniprot' ELSE NULL END) item
            FROM range(10000)""")
        con.execute("""CREATE TABLE ap_evidence_raw AS SELECT
            CASE WHEN range%2=0 THEN 'one' ELSE 'two' END resource,
            struct_pack(source:=CASE WHEN range%2=0 THEN 'declared' ELSE NULL END,
                        dataset:=CASE WHEN range%2=0 THEN 'observations' ELSE NULL END) item
            FROM range(10000)""")
        con.execute("""CREATE TABLE ap_annotation_raw AS SELECT
            CASE range%4 WHEN 0 THEN 'relation' WHEN 1 THEN 'object'
                         WHEN 2 THEN 'evidence' ELSE NULL END AS "scope" FROM range(10000)""")
        existing = {
            "data_source": ((41, "one"), (70, "unobserved-source")),
            "vocab_relation_predicate": ((1001, "affects"),),
            "vocab_entity_type": ((501, "protein"),),
            "dataset": ((123, 41, "omnipath:published_entities"),),
        }
        previous = ObservedConnection(con, emulate_previous=True)
        expected = _prepare_dimensions(previous, existing)
        bounded = ObservedConnection(con)
        actual = _prepare_dimensions(bounded, existing)
        assert actual == expected
        assert len(bounded.name_fetches) == 6
        assert all(
            len(result) <= 3
            and len(set(result)) == len(result)
            and all(row[0] is not None for row in result)
            for result in bounded.name_fetches
        )
        assert sum(map(len, previous.name_fetches)) > 30000
        assert sum(map(len, bounded.name_fetches)) == 14
        # The dataset UNION already has pair grain; source IDs and dataset IDs
        # remain identical, including seeded names absent from the raw input.
        assert len(bounded.dataset_fetches) == len(previous.dataset_fetches) == 1
        assert sorted(bounded.dataset_fetches[0]) == sorted(previous.dataset_fetches[0])
        assert len(bounded.dataset_fetches[0]) == 5
        assert (41, "one") in actual["data_source"]
        assert (70, "unobserved-source") in actual["data_source"]
        assert (1001, "affects") in actual["vocab_relation_predicate"]
        assert (123, 41, "omnipath:published_entities") in actual["dataset"]


def test_companion_ddl_quotes_schema_and_contains_no_full_record_json():
    ddl = companion_ddl('schema"quoted', retain_published_provenance=True)
    assert len(ddl) == 9
    assert all('"schema""quoted"' in statement for statement in ddl)
    assert not any("record_json" in statement for statement in ddl)
    assert sum("jsonb" in statement for statement in ddl) == 2
