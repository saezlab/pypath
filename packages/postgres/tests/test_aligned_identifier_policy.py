"""Main's InChI exclusion is a namespace policy, never a value-length filter."""

import ast
from collections import Counter
from contextlib import closing
from copy import deepcopy
import hashlib
from types import SimpleNamespace

import duckdb
import psycopg2
from psycopg2 import sql
import pytest

from omnipath_postgres import loader as aligned_loader, projection as aligned_projection
from omnipath_subsets.metsigdb import build as metsigdb_build
from omnipath_subsets.metsigdb.mapping import PROJECTION_IDENTIFIERS
import test_main_parity_contract as reference_contract
from published_fixture import QUANTITY, annotation, entity, fixture_rows, write_fixture


oracle = reference_contract.oracle
port = reference_contract.port
paired_schemas = reference_contract.paired_schemas
MAIN_TYPES = {
    "chebi": "Chebi:MI:0474",
    "inchi": "Standard Inchi:MI:2010",
    "smiles": "Smiles:MI:0239",
    "inchikey": "Standard Inchi Key:MI:1101",
    "name": "Name:OM:0202",
    "synonym": "Synonym:OM:0203",
    "hmdb": "Hmdb:OM:0004",
}
INCHI_SHORT = "InChI=1S/CH4/h1H4"
INCHI_LONG = "InChI=1S/" + "/".join(hashlib.sha256(str(i).encode()).hexdigest() for i in range(64))
SMILES_LONG = "C(O)" * 1024
NAME_LONG = "Published systematic name α " * 256
SYNONYM_LONG = "Published systematic synonym β " * 256
UNKNOWN_LONG = INCHI_SHORT * 256


@pytest.fixture(scope="module")
def parity_connection(postgres_dsn):
    # This fixture never accepts a remote/private deployment DSN.
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            yield connection
        finally:
            connection.rollback()


@pytest.fixture(scope="module")
def main_include_identifier():
    path = reference_contract.MAIN / "omnipath_build/ingest/common.py"
    assert (
        hashlib.sha256(path.read_bytes()).hexdigest()
        == "fce46a9d9a375017acf5eff7c7ef33a5ea6bfa1d8c0c932531f4fbe3a4eb9ee1"
    )
    function = next(
        node
        for node in ast.parse(path.read_text()).body
        if isinstance(node, ast.FunctionDef) and node.name == "include_identifier"
    )
    namespace = {"STANDARD_INCHI_IDENTIFIER_TERM": MAIN_TYPES["inchi"]}
    # Execute only the frozen main pure policy, with its declared constant.
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), "exec"), namespace)
    return namespace["include_identifier"]


def _fixture(directory):
    chemical = entity(
        "policy-chemical", "15377", namespace="chebi", entity_type="small_molecule", taxon=None
    )
    chemical["label"] = "Published chemical"
    identifiers = [
        ("chebi", "15377", True, "input"),
        ("inchi", INCHI_SHORT, False, "input"),
        ("inchi", INCHI_LONG, False, "input"),
        ("inchi", INCHI_LONG, False, "reference"),
        (MAIN_TYPES["inchi"], INCHI_SHORT, False, "old-spelling"),
        ("inchi", None, False, "null-occurrence"),
        ("smiles", SMILES_LONG, False, "input"),
        ("synonym", SYNONYM_LONG, False, "input"),
        ("name", NAME_LONG, False, "input"),
        ("inchikey", "AAAAAAAAAAAAAA-BBBBBBBBBB-C", False, "input"),
        ("hmdb", "HMDB0000001", False, "input"),
        # Long InChI-shaped content under another namespace remains untouched.
        ("unregistered-source-id", UNKNOWN_LONG, False, "other-namespace"),
    ]
    chemical["identifiers"] = [
        dict(ns=ns, id=value, is_canonical=canonical, source=source)
        for ns, value, canonical, source in identifiers
    ]
    pathway = entity(
        "policy-pathway", "source-set-42", namespace="pathway-fixture", entity_type="pathway"
    )
    pathway["identifiers"] = [
        dict(ns="name", id="Published pathway", is_canonical=False, source="input")
    ]
    # Canonical identities remain truthful even for a future canonical InChI.
    raw_canonical = entity(
        "policy-canonical-inchi",
        "InChI=1S/H2O/h1H2",
        namespace="inchi",
        entity_type="small_molecule",
        taxon=None,
    )
    cv_canonical = entity(
        "policy-canonical-cv-inchi",
        "InChI=1S/H2/h1H",
        namespace=MAIN_TYPES["inchi"],
        entity_type="small_molecule",
        taxon=None,
    )
    relation = deepcopy(fixture_rows()[1][0])
    relation.update(
        relation_key="policy-membership",
        subject_entity_key=chemical["entity_key"],
        subject_type="small_molecule",
        object_entity_key=pathway["entity_key"],
        object_type="pathway",
        predicate="part_of",
        sign=0,
        category="membership",
        interaction_class="directed",
        sources=["reactome"],
        evidence_count=1,
    )
    measurement = annotation(
        "has_quantitative_value",
        value="measured affinity",
        quantity=deepcopy(QUANTITY),
        scope="object",
    )
    relation["annotations"] = [deepcopy(measurement)]
    relation["evidence"] = [
        dict(
            source="reactome",
            dataset="membership",
            row_id="42",
            upstream_id="source-membership-42",
            annotations=[measurement],
        )
    ]
    entities = [chemical, pathway, raw_canonical, cv_canonical]
    write_fixture(directory, entities=entities, relations=[relation], payloads=[])
    assert len(entities) + 1 == 5
    return (SimpleNamespace(source="reactome", version="1.0.0", directory=directory),), entities


def _rows(connection, plan, table):
    query = next(query for query in plan.queries if query.table == table)
    return [
        dict(zip(query.columns, row, strict=True))
        for row in connection.execute(query.query).fetchall()
    ]


@pytest.mark.parametrize("retain", [False, True])
def test_generic_identifiers_match_frozen_main_policy_without_length_or_content_filter(
    tmp_path,
    main_include_identifier,
    retain,
):
    selected, entities = _fixture(tmp_path)
    expected = set()
    for item in entities:
        candidates = [(item["namespace"], item["identifier"])]
        candidates.extend(
            (identifier["ns"], identifier["id"]) for identifier in item["identifiers"] or ()
        )
        for namespace, value in candidates:
            typed = MAIN_TYPES.get(namespace, namespace)
            if main_include_identifier(typed, value):
                expected.add((typed, value))
    with duckdb.connect() as connection:
        plan = aligned_projection.prepare_aligned_release(
            connection, selected, retain_published_provenance=retain
        )
        dimensions = dict(plan.dimensions["vocab_identifier_type"])
        dictionary = _rows(connection, plan, "identifier_evidence")
        actual = {(dimensions[row["identifier_type_id"]], row["value"]) for row in dictionary}
        assert actual == expected
        assert all(main_include_identifier(namespace, value) for namespace, value in actual)
        assert (MAIN_TYPES["smiles"], SMILES_LONG) in actual
        assert (MAIN_TYPES["name"], NAME_LONG) in actual
        assert (MAIN_TYPES["synonym"], SYNONYM_LONG) in actual
        assert ("unregistered-source-id", UNKNOWN_LONG) in actual
        identifiers = {row["identifier_id"] for row in dictionary}
        for table in ("entity_identifier", "entity_evidence_identifier"):
            assert all(
                row["identifier_id"] in identifiers for row in _rows(connection, plan, table)
            )
        assert plan.compatibility["generic_identifier_policy"] == "main_standard_inchi_excluded"
        assert plan.compatibility["excluded_standard_inchi_identifier_occurrences"] == 4
        assert plan.compatibility["excluded_standard_inchi_dictionary_values"] == 5
        assert (
            _rows(connection, plan, "annotation_quantity")[0]["has_numeric_value"]
            == QUANTITY["has_numeric_value"]
        )
        assert plan.counts["relation_evidence_annotation"]
        actual_entities = {row["entity_id"]: row for row in _rows(connection, plan, "entity")}
        published = connection.execute(
            "SELECT entity_key,entity_id FROM ap_entity_occurrence"
        ).fetchall()
        for key, identifier in published:
            original = next(item for item in entities if item["entity_key"] == key)
            assert actual_entities[identifier]["canonical_identifier"] == original["identifier"]
            # No replacement identifier or recanonicalized entity is invented.
            if original["namespace"] in ("inchi", MAIN_TYPES["inchi"]):
                assert (
                    dimensions[actual_entities[identifier]["canonical_identifier_type_id"]]
                    == original["namespace"]
                )
                assert all(
                    row["entity_id"] != identifier
                    for row in _rows(connection, plan, "entity_identifier")
                )
        if retain:
            original = entities[0]["identifiers"]
            occurrences = sorted(
                (
                    row
                    for row in _rows(connection, plan, "parquet_identifier_occurrence")
                    if row["entity_key"] == "policy-chemical"
                ),
                key=lambda row: row["ordinal"],
            )
            assert [
                (row["ordinal"], row["ns"], row["identifier"], row["is_canonical"], row["source"])
                for row in occurrences
            ] == [
                (ordinal, row["ns"], row["id"], row["is_canonical"], row["source"])
                for ordinal, row in enumerate(original)
            ]
            excluded = [row for row in occurrences if row["ns"] in ("inchi", MAIN_TYPES["inchi"])]
            assert len(excluded) == 5  # Includes the preserved NULL source occurrence.
            assert all(row["identifier_id"] is None for row in excluded)
        else:
            assert "parquet_identifier_occurrence" not in plan.counts


def test_smiles_mapping_reaches_actual_main_metsigdb_publication(
    tmp_path,
    paired_schemas,
    port,
):
    connection, (_main, adapted) = paired_schemas
    selected, _entities = _fixture(tmp_path)
    assert dict(PROJECTION_IDENTIFIERS)["smiles"] == MAIN_TYPES["smiles"]
    with duckdb.connect() as duck:
        plan = aligned_projection.prepare_aligned_release(
            duck, selected, dimension_rows=aligned_loader._read_dimensions(connection, adapted)
        )
        aligned_loader._write_dimensions(connection, adapted, plan.dimensions)
        for _, source in plan.sources:
            port.schema.ensure_source_partitions(connection, schema=adapted, source=source)
        with connection.cursor() as cursor:
            for statement in aligned_projection.companion_ddl(adapted):
                # paired_schemas already creates annotation_quantity.
                cursor.execute(statement.replace("CREATE TABLE ", "CREATE TABLE IF NOT EXISTS ", 1))
        for query in plan.queries:
            aligned_loader._copy(connection, adapted, query, duck, tmp_path)
        with connection.cursor() as cursor:
            port.derive._populate_entity_identifier_lookup(cursor, adapted)
            cursor.execute(sql.SQL("SET search_path = {},public").format(sql.Identifier(adapted)))
            cursor.execute(
                "CREATE TEMP TABLE metsigdb_stage (metabolite_entity_id uuid,set_entity_id uuid)"
            )
            cursor.execute(
                "INSERT INTO metsigdb_stage SELECT chemical.entity_id,pathway.entity_id FROM "
                "entity chemical CROSS JOIN entity pathway WHERE chemical.canonical_identifier='15377' "
                "AND pathway.canonical_identifier='source-set-42'"
            )
            source_id = next(
                identifier for identifier, source in plan.sources if source == "reactome"
            )
            expected = None
            for old_main, text in (
                (
                    True,
                    (
                        reference_contract.MAIN
                        / "omnipath_build/metsigdb/sql/publish_membership.sql"
                    ).read_text(),
                ),
                (False, metsigdb_build._sql_text("publish_membership.sql")),
            ):
                # Frozen main has the original matched flag. Published inputs
                # use truthful status 5, with the approved chemical eligibility
                # adaptation; give the oracle its original status vocabulary.
                cursor.execute(
                    "UPDATE entity_evidence_resolution SET status_id=%s", [1 if old_main else 5]
                )
                cursor.execute(text, {"source_id": source_id})
                cursor.execute("SELECT * FROM metsigdb_projection")
                actual = Counter(tuple(row) for row in cursor.fetchall())
                assert len(actual) == 1
                cursor.execute("SELECT * FROM metsigdb_set_attrs")
                attrs = Counter(tuple(row) for row in cursor.fetchall())
                if expected is None:
                    expected = (actual, attrs)
                else:
                    assert (actual, attrs) == expected
                cursor.execute("SELECT smiles,inchikey,hmdb,chebi FROM metsigdb_projection")
                assert cursor.fetchall() == [
                    (SMILES_LONG, "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "HMDB0000001", "15377")
                ]
                cursor.execute("SELECT set_label,organism FROM metsigdb_set_attrs")
                assert cursor.fetchall() == [("Published pathway", 9606)]
            cursor.execute("SELECT DISTINCT status_id FROM entity_evidence_resolution")
            assert cursor.fetchall() == [(5,)]
        connection.commit()
