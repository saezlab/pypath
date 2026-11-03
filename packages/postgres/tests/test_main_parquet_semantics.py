"""Independent bounded Parquet quantity and frozen-main ontology checks.

These fixtures contain fifteen source entities/statements in total. They never
call a parser, resolver, resource builder, network or server deployment. Main's
legacy definition/comment distinction is absent from published descriptions;
the ontology oracle therefore tests the supported published description only.
"""

from collections import Counter
from contextlib import closing
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import uuid

import duckdb
import psycopg2
from psycopg2 import sql
import pytest

from omnipath_postgres import loader as aligned_loader
from omnipath_postgres.projection import companion_ddl, prepare_aligned_release
import test_main_parity_contract as reference_contract
from published_fixture import annotation, entity, fixture_rows, write_fixture


oracle = reference_contract.oracle
port = reference_contract.port
paired_schemas = reference_contract.paired_schemas


@pytest.fixture(scope="module")
def parity_connection(postgres_dsn):
    # Deliberately use only conftest's disposable local PostgreSQL fixture.
    # This test file does not inspect or connect to an externally supplied DSN.
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            yield connection
        finally:
            connection.rollback()


def _uuid_json(payload):
    text = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    return uuid.UUID(hex=hashlib.md5(text.encode()).hexdigest())


def _entity_id(key):
    return _uuid_json(["published-entity", key])


def _rows(connection, plan, table):
    query = next(item for item in plan.queries if item.table == table)
    return [
        dict(zip(query.columns, row, strict=True))
        for row in connection.execute(query.query).fetchall()
    ]


def test_quantity_parquet_projection_retains_complete_structured_rows_and_identities(tmp_path):
    entities, relations, _ = fixture_rows()
    entities = entities[:2]
    for item in entities:
        item["annotations"] = []
    base = {
        "has_numeric_value": 12.5,
        "has_unit": "UO:0000061",
        "has_unit_prefix": "nano",
        "has_binary_relation": "less_than",
        "source_field": "IC50",
        "comparator": "<=",
    }
    quantities = [
        base,
        dict(base, has_unit_prefix="micro"),
        dict(base, has_unit="UO:0000062"),
        dict(base, comparator=">="),
        dict(base, source_field="KD"),
        dict(base, has_binary_relation="greater_than"),
        dict(base, has_numeric_value=13.5),
    ]
    # Same published text and term: structured fields alone distinguish the
    # quantities. (Published tables store a quantity as columns: all null is none.)
    published_value = "original published measurement"
    measurements = [
        annotation(
            "BAO:0000192", value=published_value, quantity=deepcopy(quantity), scope="relation"
        )
        for quantity in [*quantities, None]
    ]
    measurements.append(deepcopy(measurements[0]))
    relations[0]["annotations"] = deepcopy(measurements)
    relations[0]["evidence"] = [
        {
            "source": "reported-source",
            "dataset": "reported-dataset",
            "row_id": "source:7",
            "upstream_id": "assay-seven",
            "annotations": deepcopy(measurements),
        }
    ]
    relations[0]["evidence_count"] = 1
    write_fixture(tmp_path, entities=entities, relations=relations, payloads=[])
    expected_annotations, expected_quantities = {}, {}
    for quantity in [*quantities, None]:
        key = _uuid_json(
            {
                "kind": "annotation",
                "term": "BAO:0000192",
                "value": published_value,
                "quantity": quantity,
            }
        )
        numeric = quantity["has_numeric_value"] if quantity is not None else None
        expected_annotations[key] = {
            "annotation_key": key,
            "term": "BAO:0000192",
            "value": str(numeric) if numeric is not None else published_value,
            "unit": quantity["has_unit"] if quantity is not None else None,
        }
        if quantity is not None:
            expected_quantities[key] = {
                "annotation_key": key,
                **quantity,
                "published_value": published_value,
            }
    with duckdb.connect() as connection:
        plan = prepare_aligned_release(
            connection,
            (
                SimpleNamespace(
                    directory=tmp_path, source="quantity_fixture", version="quantity-v1"
                ),
            ),
            retain_published_provenance=True,
        )
        assert {
            row["annotation_key"]: row for row in _rows(connection, plan, "annotation")
        } == expected_annotations
        assert {
            row["annotation_key"]: row for row in _rows(connection, plan, "annotation_quantity")
        } == expected_quantities
        evidence = _rows(connection, plan, "parquet_evidence")
        assert len(evidence) == 1
        assert evidence[0]["original_row_id"] == "source:7"
        assert evidence[0]["upstream_id"] == "assay-seven"
        scopes = {
            name: identifier for identifier, name in plan.dimensions["vocab_annotation_scope"]
        }
        assert Counter(
            tuple(row.values()) for row in _rows(connection, plan, "relation_evidence_annotation")
        ) == Counter(
            (evidence[0]["source_id"], evidence[0]["relation_evidence_id"], key, scopes["relation"])
            for key in expected_annotations
        )
        occurrences = _rows(connection, plan, "parquet_annotation_occurrence")
        for kind, ordinal in (("relation", -1), ("evidence", 0)):
            actual = sorted(
                (row for row in occurrences if row["owner_kind"] == kind),
                key=lambda row: row["ordinal"],
            )
            expected_keys = [
                _uuid_json(
                    {
                        "kind": "annotation",
                        "term": "BAO:0000192",
                        "value": published_value,
                        "quantity": item["quantity"],
                    }
                )
                for item in measurements
            ]
            assert [
                (
                    row["ordinal"],
                    row["evidence_ordinal"],
                    row["annotation_key"],
                    row["source"],
                    row["dataset"],
                    row["scope"],
                )
                for row in actual
            ] == [
                (i, ordinal, key, "reported-source", "reported-dataset", "relation")
                for i, key in enumerate(expected_keys)
            ]


TERMS = (
    ("chemont", "a", "CHEMONT:A", "A label", "a published description", ("A synonym",)),
    ("chemont", "b", "CHEMONT:B", "B label", "b published description", ("B synonym",)),
    ("chemont", "c", "CHEMONT:C", "C label", None, ()),
    ("chemont", "d", "CHEMONT:D", "D label", "d published description", ()),
    ("go", "g", "GO:0000001", "G label", "g published description", ("G synonym",)),
    ("go", "h", "GO:0000002", None, None, ()),
)
AXIOMS = (
    ("chemont", "a", "subclass_of", "b"),
    ("chemont", "b", "subclass_of", "c"),
    ("chemont", "a", "subclass_of", "c"),
    ("chemont", "a", "part_of", "d"),
    ("chemont", "b", "related_to", "d"),
    ("go", "g", "subclass_of", "h"),
)
SOURCE_IDS = {"chemont": 9006, "go": 9007}
ONTOLOGY_IDS = {"chemont": "chemont", "go": "gene_ontology"}


def _ontology_fixture(root):
    by_name = {}
    for source, name, canonical, label, description, synonyms in TERMS:
        item = entity(
            "ontology-published:" + name,
            canonical,
            namespace=source,
            taxon=None,
            entity_type="ontology_class",
        )
        item.update(label=label, has_hierarchy=True)
        item["identifiers"] = [dict(ns=source, id=canonical, is_canonical=True, source=source)]
        item["identifiers"].extend(
            dict(ns="synonym", id=value, is_canonical=False, source=source) for value in synonyms
        )
        if label is not None:
            item["identifiers"].append(dict(ns="name", id=label, is_canonical=False, source=source))
        item["annotations"] = (
            []
            if description is None
            else [
                {
                    "term": "description",
                    "value": description,
                    "quantity": None,
                    "source": source,
                    "dataset": "ontology",
                },
            ]
        )
        by_name[name] = item
    selected = []
    for source in SOURCE_IDS:
        statements = []
        for number, (owner, subject, predicate, obj) in enumerate(AXIOMS, 1):
            if owner != source:
                continue
            statement = deepcopy(fixture_rows()[1][0])
            statement.update(
                relation_key=f"ontology-published-axiom:{number}",
                statement_kind="ontology_axiom",
                subject_entity_key=by_name[subject]["entity_key"],
                object_entity_key=by_name[obj]["entity_key"],
                predicate=predicate,
                subject_type="ontology_class",
                object_type="ontology_class",
                subject_label=by_name[subject]["label"],
                object_label=by_name[obj]["label"],
                taxon=None,
                sign=0,
                category="ontology",
                interaction_class="directed",
                sources=[source],
                evidence_count=1,
                annotations=[],
                evidence=[
                    {
                        "source": source,
                        "dataset": "ontology",
                        "row_id": str(number),
                        "upstream_id": f"axiom-{number}",
                        "annotations": [],
                    }
                ],
            )
            statements.append(statement)
        directory = root / source
        write_fixture(
            directory,
            entities=[
                value
                for name, value in by_name.items()
                if next(term[0] for term in TERMS if term[1] == name) == source
            ],
            relations=statements,
            payloads=[],
        )
        selected.append(SimpleNamespace(source=source, version="ontology-v1", directory=directory))
    assert len(TERMS) + len(AXIOMS) == 12
    return tuple(selected), by_name


def _main_ontology_facts(connection, schema, oracle, by_name):
    """Declare the resolved supported facts directly for unchanged frozen main."""

    def query(text):
        return sql.SQL(text).format(s=sql.Identifier(schema))

    with connection.cursor() as cursor:
        cursor.executemany(
            query("INSERT INTO {s}.data_source(source_id,name) VALUES(%s,%s)").as_string(
                connection
            ),
            [(identifier, name) for name, identifier in SOURCE_IDS.items()],
        )
        cursor.executemany(
            query(
                "INSERT INTO {s}.dataset(dataset_id,source_id,name) VALUES(%s,%s,'ontology')"
            ).as_string(connection),
            [(identifier, identifier) for identifier in SOURCE_IDS.values()],
        )
        cursor.execute(
            query(
                "INSERT INTO {s}.vocab_entity_type(name) VALUES('Cv Term:OM:0012') ON CONFLICT(name) DO NOTHING"
            )
        )
        cursor.execute(
            query("SELECT entity_type_id FROM {s}.vocab_entity_type WHERE name='Cv Term:OM:0012'")
        )
        term_type = cursor.fetchone()[0]
        for identifier, namespace in (
            (500001, "chemont"),
            (500002, "go"),
            (500003, "Synonym:OM:0203"),
        ):
            cursor.execute(
                query(
                    "INSERT INTO {s}.vocab_identifier_type(identifier_type_id,name) VALUES(%s,%s)"
                ),
                (identifier, namespace),
            )
        cursor.execute(query("SELECT name,identifier_type_id FROM {s}.vocab_identifier_type"))
        namespaces = dict(cursor.fetchall())
        for source, name, canonical, _label, description, _synonyms in TERMS:
            entity_id = str(_entity_id(by_name[name]["entity_key"]))
            evidence_id = str(uuid.uuid5(uuid.NAMESPACE_URL, "ontology-main-evidence:" + name))
            cursor.execute(
                query(
                    "INSERT INTO {s}.entity(entity_id,entity_type_id,canonical_identifier_type_id,canonical_identifier,resolution_status_id) VALUES(%s,%s,%s,%s,2)"
                ),
                (entity_id, term_type, namespaces[source], canonical),
            )
            cursor.execute(
                query(
                    "INSERT INTO {s}.entity_evidence(source_id,entity_evidence_id,dataset_id,row_id,entity_role_id,entity_type_id) VALUES(%s,%s,%s,%s,1,%s)"
                ),
                (SOURCE_IDS[source], evidence_id, SOURCE_IDS[source], ord(name), term_type),
            )
            cursor.execute(
                query(
                    "INSERT INTO {s}.entity_evidence_resolution(source_id,entity_evidence_id,status_id,entity_id) VALUES(%s,%s,2,%s)"
                ),
                (SOURCE_IDS[source], evidence_id, entity_id),
            )
            for item in by_name[name]["identifiers"]:
                namespace = {"name": "Name:OM:0202", "synonym": "Synonym:OM:0203"}.get(
                    item["ns"], item["ns"]
                )
                identifier_id = str(
                    uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([namespace, item["id"]]))
                )
                cursor.execute(
                    query(
                        "INSERT INTO {s}.identifier_evidence(identifier_id,identifier_type_id,value) VALUES(%s,%s,%s) ON CONFLICT DO NOTHING"
                    ),
                    (identifier_id, namespaces[namespace], item["id"]),
                )
                cursor.execute(
                    query("INSERT INTO {s}.entity_identifier VALUES(%s,%s,%s)"),
                    (SOURCE_IDS[source], entity_id, identifier_id),
                )
                cursor.execute(
                    query("INSERT INTO {s}.entity_evidence_identifier VALUES(%s,%s,%s)"),
                    (SOURCE_IDS[source], evidence_id, identifier_id),
                )
            if description is not None:
                annotation_id = str(
                    uuid.uuid5(uuid.NAMESPACE_URL, "ontology-main-definition:" + name)
                )
                cursor.execute(
                    query("INSERT INTO {s}.annotation(annotation_key,term,value) VALUES(%s,%s,%s)"),
                    (annotation_id, oracle.derive.ONTOLOGY_DEFINITION_TERM, description),
                )
                cursor.execute(
                    query("INSERT INTO {s}.entity_evidence_annotation VALUES(%s,%s,%s)"),
                    (SOURCE_IDS[source], evidence_id, annotation_id),
                )
        for predicate in ("is_a", "part_of", "related_to"):
            cursor.execute(
                query(
                    "INSERT INTO {s}.vocab_relation_predicate(name) VALUES(%s) ON CONFLICT(name) DO NOTHING"
                ),
                (predicate,),
            )
        cursor.execute(query("SELECT name,relation_predicate_id FROM {s}.vocab_relation_predicate"))
        predicates = dict(cursor.fetchall())
        for source, subject, predicate, obj in AXIOMS:
            cursor.execute(
                query("INSERT INTO {s}.entity_ontology_relation VALUES(%s,%s,%s,%s,%s)"),
                (
                    SOURCE_IDS[source],
                    str(_entity_id(by_name[subject]["entity_key"])),
                    predicates["is_a" if predicate == "subclass_of" else predicate],
                    str(_entity_id(by_name[obj]["entity_key"])),
                    ONTOLOGY_IDS[source],
                ),
            )
    connection.commit()


def test_actual_ontology_parquets_copy_and_frozen_main_derivation_agree(
    tmp_path, paired_schemas, oracle, port
):
    connection, (main, adapted) = paired_schemas
    resources, by_name = _ontology_fixture(tmp_path)
    _main_ontology_facts(connection, main, oracle, by_name)
    with connection.cursor() as cursor:
        cursor.executemany(
            sql.SQL("INSERT INTO {}.data_source(source_id,name) VALUES(%s,%s)")
            .format(sql.Identifier(adapted))
            .as_string(connection),
            [(identifier, name) for name, identifier in SOURCE_IDS.items()],
        )
    connection.commit()
    with duckdb.connect() as duck:
        plan = prepare_aligned_release(
            duck, resources, dimension_rows=aligned_loader._read_dimensions(connection, adapted)
        )
        aligned_loader._write_dimensions(connection, adapted, plan.dimensions)
        for _, source in plan.sources:
            port.schema.ensure_source_partitions(connection, schema=adapted, source=source)
        with connection.cursor() as cursor:
            for statement in companion_ddl(adapted):
                if '."annotation_quantity"' not in statement:
                    cursor.execute(statement)
        spool = tmp_path / "copy"
        spool.mkdir()
        for query in plan.queries:
            aligned_loader._copy(connection, adapted, query, duck, spool)
        connection.commit()
        assert plan.compatibility["unknown_ontology_scope_statements"] == 0
        assert plan.counts["relation"] == 0
    # Fresh main leaves the base ontology_terms table empty. Its complete
    # serving catalogue below derives from relational ontology/identifier data.
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("""SELECT source_id,term_entity_id::text,term_id,
            ontology_prefix,label,definition,ontology_id,synonyms,synonyms_text,sources
            FROM {}.ontology_terms""").format(sql.Identifier(adapted))
        )
        copied_terms = Counter(json.dumps(row) for row in cursor.fetchall())
    assert copied_terms == Counter()
    expected = reference_contract._ontology_rows(
        connection, main, oracle.derive, reference_contract.MAIN / "omnipath_build"
    )
    actual = reference_contract._ontology_rows(
        connection, adapted, port.derive, Path(port.derive.__file__).parent.parent
    )
    assert actual == expected
    ids = {name: str(_entity_id(item["entity_key"])) for name, item in by_name.items()}
    assert actual["closure"] == Counter(
        json.dumps(row)
        for row in [
            [ids["a"], ids["b"], 1],
            [ids["b"], ids["c"], 1],
            [ids["a"], ids["c"], 1],
        ]
    )
    terms = {json.loads(row)[0]: json.loads(row) for row in actual["ontology_terms"]}
    assert terms[ids["a"]][3:6] == ["A label", "a published description", ["A synonym"]]
    assert terms[ids["h"]][1:5] == ["GO:0000002", "go", "GO:0000002", None]
    assert terms[ids["g"]][9:11] == ["gene_ontology", ["go"]]
    assert len(actual["ontology_assertions"]) == len(AXIOMS)
