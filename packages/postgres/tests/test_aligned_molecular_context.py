"""Serving-v4 molecular facts survive COPY without creating graph states."""

from copy import deepcopy
import json

import duckdb
import pytest

from omnipath_core.fixtures import nested_rows
from omnipath_core.molecular_forms import normalize_molecular_form
from omnipath_core.versioning import SERVING_SCHEMA_VERSION
from omnipath_postgres.projection import companion_ddl, prepare_aligned_release
from omnipath_postgres.relational.db.schema import CONTENT_TABLES
from test_aligned_projection import selected
from published_fixture import ENTITY_A, ENTITY_B, ENTITY_C, RELATION, fixture_rows, write_fixture


CONTEXT_TABLES = {
    "entity_reference_context",
    "statement_reference_context",
    "molecular_evidence_context",
}


def fixture(directory):
    entities, relations, payloads = fixture_rows()
    gene_a, gene_b = "gene:virtual:A", "gene:virtual:B"
    entities[0].update(reference_entity_key=gene_a, gene_reference_keys=[gene_a, gene_a])
    entities[1].update(reference_entity_key=gene_b, gene_reference_keys=[])
    first = normalize_molecular_form(
        {
            "protein_entity_key": ENTITY_A,
            "isoform_identifier": {"ns": "uniprot", "id": "P04637-2"},
            "modifications": [
                {
                    "term": "MOD:00696",
                    "residue": "S",
                    "position": 15,
                    "coordinate_reference": {"coordinate_system": "unknown"},
                }
            ],
        }
    )
    second = normalize_molecular_form(
        {
            "protein_entity_key": ENTITY_B,
            "transcript_entity_key": "transcript:product:B",
            "sequence_identifiers": [{"ns": "enst", "id": "ENST000001"}],
            "variants": [
                {
                    "reference": "A",
                    "alternate": "V",
                    "position": 12,
                    "coordinate_reference": {
                        "identifier": {"ns": "enst", "id": "ENST000001"},
                        "coordinate_system": "transcript",
                        "position_base": 1,
                    },
                }
            ],
        }
    )
    relation = relations[0]
    relation.update(subject_reference_entity_key=gene_a, object_reference_entity_key=gene_b)
    relation["evidence"][0].update(subject_molecular_form=first, object_molecular_form=second)
    relation["evidence"][1].update(subject_molecular_form=second, object_molecular_form=first)
    # Identical source row identifiers do not merge occurrences or mix endpoints.
    relation["evidence"].append(deepcopy(relation["evidence"][0]))
    relation["evidence_count"] = 3
    entities[2]["evidence"] = [
        {
            "source": "standalone-source",
            "dataset": "standalone",
            "row_id": "000007",
            "upstream_id": "standalone-event",
            "annotations": [],
            "molecular_form": second,
        }
    ]
    write_fixture(directory, entities=entities, relations=relations, payloads=payloads)
    return selected(directory, version="1.0.0"), gene_a, gene_b


def rows(connection, plan, table):
    query = next(item for item in plan.queries if item.table == table)
    return [
        dict(zip(query.columns, row, strict=True))
        for row in connection.execute(query.query).fetchall()
    ]


@pytest.mark.parametrize("retain", [False, True])
def test_molecular_context_preserves_exact_pairings_standalone_and_virtual_keys(tmp_path, retain):
    selected, gene_a, gene_b = fixture(tmp_path)
    originals = nested_rows(tmp_path / "relation.parquet")[0]["evidence"]
    standalone = next(
        e for e in nested_rows(tmp_path / "entity.parquet") if e["entity_key"] == ENTITY_C
    )["evidence"][0]
    with duckdb.connect() as connection:
        plan = prepare_aligned_release(connection, (selected,), retain_published_provenance=retain)
        assert CONTEXT_TABLES <= set(plan.counts)
        context = rows(connection, plan, "molecular_evidence_context")
        paired = sorted(
            (item for item in context if item["owner_kind"] == "relation"),
            key=lambda item: item["ordinal"],
        )
        assert [item["ordinal"] for item in paired] == [0, 1, 2]
        assert [json.loads(item["occurrence_json"]) for item in paired] == originals
        assert len({item["relation_evidence_id"] for item in paired}) == 3
        assert {item["relation_evidence_id"] for item in paired} == {
            item["relation_evidence_id"] for item in rows(connection, plan, "relation_evidence")
        }
        entity_occurrence = next(item for item in context if item["owner_kind"] == "entity")
        assert entity_occurrence["relation_evidence_id"] is None
        assert entity_occurrence["row_id"] == "000007"
        assert json.loads(entity_occurrence["occurrence_json"]) == standalone
        references = {
            item["entity_key"]: item for item in rows(connection, plan, "entity_reference_context")
        }
        assert references[ENTITY_A]["reference_entity_key"] == gene_a
        assert json.loads(references[ENTITY_A]["gene_reference_keys"]) == [gene_a, gene_a]
        assert json.loads(references[ENTITY_B]["gene_reference_keys"]) == []
        assert (
            next(item for key, item in references.items() if key not in {ENTITY_A, ENTITY_B})[
                "gene_reference_keys"
            ]
            is None
        )
        assert rows(connection, plan, "statement_reference_context") == [
            {
                "resource": "fixture",
                "version": "1.0.0",
                "relation_key": RELATION,
                "subject_reference_entity_key": gene_a,
                "object_reference_entity_key": gene_b,
            }
        ]
        assert plan.counts["entity"] == 3
        assert plan.counts["relation"] == 1
        assert "state" not in plan.counts
        assert "gene_protein_representative" not in plan.counts


def test_context_is_mandatory_and_resettable():
    assert SERVING_SCHEMA_VERSION == 6
    assert CONTEXT_TABLES <= set(CONTENT_TABLES)
    ddl = companion_ddl("fixture")
    for table in CONTEXT_TABLES:
        assert any(f'."{table}"' in statement for statement in ddl)
    assert all("REFERENCES" not in statement for statement in ddl)


@pytest.mark.integration
def test_local_context_copy_preserves_json_pairs_and_reset_clears_rows(tmp_path, postgres_dsn):
    from contextlib import closing
    import uuid
    import psycopg2
    from psycopg2 import sql
    from omnipath_postgres.loader import _copy
    from omnipath_postgres.relational.db.schema import reset_content_tables

    selected, _, _ = fixture(tmp_path)
    schema = "molecular_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as postgres, duckdb.connect() as duck:
        plan = prepare_aligned_release(duck, (selected,))
        try:
            with postgres.cursor() as cursor:
                cursor.execute(sql.SQL("CREATE SCHEMA {}").format(sql.Identifier(schema)))
                for statement in companion_ddl(schema):
                    cursor.execute(statement)
            for query in plan.queries:
                if query.table in CONTEXT_TABLES:
                    assert (
                        _copy(postgres, schema, query, duck, tmp_path).rows
                        == plan.counts[query.table]
                    )
            postgres.commit()
            with postgres.cursor() as cursor:
                cursor.execute(
                    sql.SQL("""SELECT owner_kind,ordinal,relation_evidence_id,occurrence_json
                    FROM {}.molecular_evidence_context ORDER BY owner_kind,ordinal""").format(
                        sql.Identifier(schema)
                    )
                )
                copied = cursor.fetchall()
                assert len(copied) == 4
                assert copied[0][0] == "entity"
                assert copied[0][2] is None
                originals = nested_rows(tmp_path / "relation.parquet")[0]["evidence"]
                assert [row[3] for row in copied[1:]] == originals
                cursor.execute(
                    sql.SQL("""SELECT gene_reference_keys FROM {}.entity_reference_context
                    WHERE entity_key=%s""").format(sql.Identifier(schema)),
                    (ENTITY_A,),
                )
                assert cursor.fetchone()[0] == ["gene:virtual:A", "gene:virtual:A"]
            assert CONTEXT_TABLES <= set(reset_content_tables(postgres, schema=schema))
            with postgres.cursor() as cursor:
                for table in CONTEXT_TABLES:
                    cursor.execute(
                        sql.SQL("SELECT count(*) FROM {}.{}").format(
                            sql.Identifier(schema), sql.Identifier(table)
                        )
                    )
                    assert cursor.fetchone() == (0,)
        finally:
            postgres.rollback()
            with postgres.cursor() as cursor:
                cursor.execute(
                    sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(sql.Identifier(schema))
                )
            postgres.commit()
