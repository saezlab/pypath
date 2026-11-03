"""Storage policy changes must not change normalized scientific publication.

Eighteen synthetic entity/statement records exercise the actual Parquet reader,
DuckDB preparation, local PostgreSQL COPY and main derivation. No resolver,
resource builder, external database or deployment is called here.
"""

from collections import Counter
from contextlib import closing
from copy import deepcopy
from datetime import datetime
import json
from types import SimpleNamespace
import uuid

import duckdb
import psycopg2
from psycopg2 import sql
import pyarrow.parquet as pq
import pytest

from omnipath_postgres import loader as aligned_loader, projection as aligned_projection
from omnipath_core.fixtures import nested_rows, write_manifest
import test_main_parity_contract as reference_contract
from test_aligned_ontology_memory import _serving_rows
from test_main_parquet_semantics import _ontology_fixture
from published_fixture import QUANTITY, annotation, fixture_rows, write_fixture


oracle = reference_contract.oracle
port = reference_contract.port
AUDIT_TABLES = {
    "parquet_entity",
    "parquet_statement",
    "parquet_evidence",
    "parquet_identifier_occurrence",
    "parquet_annotation_occurrence",
}
NORMALIZED_TABLES = {
    "entity",
    "identifier_evidence",
    "annotation",
    "entity_evidence",
    "entity_evidence_resolution",
    "entity_identifier",
    "entity_evidence_identifier",
    "entity_evidence_annotation",
    "relation",
    "relation_evidence",
    "relation_evidence_relation",
    "relation_evidence_annotation",
    "entity_ontology_relation",
    "ontology_terms",
    "annotation_quantity",
    "entity_reference_context",
    "statement_reference_context",
    "molecular_evidence_context",
}
RETAINED_COPY_ORDER = (
    "entity",
    "identifier_evidence",
    "annotation",
    "entity_evidence",
    "entity_evidence_resolution",
    "entity_identifier",
    "entity_evidence_identifier",
    "entity_evidence_annotation",
    "relation",
    "relation_evidence",
    "relation_evidence_relation",
    "relation_evidence_annotation",
    "entity_ontology_relation",
    "ontology_terms",
    "entity_reference_context",
    "statement_reference_context",
    "molecular_evidence_context",
    "parquet_entity",
    "parquet_statement",
    "parquet_evidence",
    "parquet_identifier_occurrence",
    "parquet_annotation_occurrence",
    "annotation_quantity",
)


def _freeze(value):
    if isinstance(value, datetime):
        # Main's created/updated timestamps describe the separate publications,
        # not source facts. Every other complete PostgreSQL column is compared.
        return "<publication timestamp>"
    if isinstance(value, dict):
        return tuple(sorted((key, _freeze(item)) for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(map(_freeze, value))
    return value


def _copy_rows(connection, plan):
    return {
        item.table: Counter(_freeze(row) for row in connection.execute(item.query).fetchall())
        for item in plan.queries
    }


def _fixture(root):
    entities, relations, _ = fixture_rows()
    # Presentation/type/source distinctions are occurrence facts, even when
    # normalized identifier and annotation dictionaries deduplicate values.
    entities[0]["identifiers"].append(
        dict(entities[0]["identifiers"][0], is_canonical=False, source="other-identifier-source")
    )
    relations[0]["evidence"][1].update(
        row_id="source:row:001", upstream_id="second-event", source="other-source", dataset=None
    )
    second = deepcopy(relations[0])
    second.update(relation_key="qualified-second", sign=1)
    quantity = dict(QUANTITY, has_unit_prefix="micro", comparator=">", source_field="second-assay")
    scoped = annotation(
        "has_quantitative_value",
        value="same source presentation",
        quantity=quantity,
        scope="object",
    )
    second["annotations"] = [deepcopy(scoped), deepcopy(scoped)]
    second["evidence"] = [
        dict(
            source="reported-source",
            dataset="reported-dataset",
            row_id="7",
            upstream_id="third-event",
            annotations=[scoped],
        )
    ]
    write_fixture(
        root / "record_owner", entities=entities, relations=relations + [second], payloads=[]
    )
    shared = deepcopy(entities[0])
    shared["identifiers"][0]["source"] = "cross-resource-identifier-source"
    write_fixture(root / "shared_owner", entities=[shared], relations=[], payloads=[])
    ontology, _ = _ontology_fixture(root)
    selected = (
        SimpleNamespace(
            source="record_owner", version="retention-v1", directory=root / "record_owner"
        ),
        SimpleNamespace(
            source="shared_owner", version="retention-v2", directory=root / "shared_owner"
        ),
        *ontology,
    )
    assert (
        sum(
            pq.ParquetFile(r.directory / name).metadata.num_rows
            for r in selected
            for name in ("entity.parquet", "relation.parquet")
        )
        == 18
    )
    return selected


def _published_counts(selected):
    counts = dict(
        entity_rows=0,
        statement_rows=0,
        evidence_occurrences=0,
        identifier_occurrences=0,
        annotation_occurrences=0,
    )
    for resource in selected:
        for entity in nested_rows(resource.directory / "entity.parquet"):
            counts["entity_rows"] += 1
            counts["identifier_occurrences"] += len(entity["identifiers"] or ())
            counts["annotation_occurrences"] += len(entity["annotations"] or ())
        for statement in nested_rows(resource.directory / "relation.parquet"):
            counts["statement_rows"] += 1
            counts["annotation_occurrences"] += len(statement["annotations"] or ())
            for evidence in statement["evidence"] or ():
                counts["evidence_occurrences"] += 1
                counts["annotation_occurrences"] += len(evidence["annotations"] or ())
    return counts


def test_default_does_not_prepare_audit_tables_and_keeps_all_normalized_rows(tmp_path):
    selected = _fixture(tmp_path)
    plans, copied = {}, {}
    for retain in (False, True):
        observed = []
        with duckdb.connect() as connection:
            kwargs = {"retain_published_provenance": True} if retain else {}
            plan = aligned_projection.prepare_aligned_release(
                connection, selected, progress=lambda fields: observed.append(fields), **kwargs
            )
            plans[retain] = plan
            assert tuple(item.table for item in plan.queries) == tuple(
                table for table in RETAINED_COPY_ORDER if retain or table not in AUDIT_TABLES
            )
            copied[retain] = _copy_rows(connection, plan)
            present = {row[0] for row in connection.execute("SHOW TABLES").fetchall()}
            assert {table for table in AUDIT_TABLES if "ap_copy_" + table in present} == (
                AUDIT_TABLES if retain else set()
            )
            prepared = {fields["table"] for fields in observed if fields["phase"] == "prepare_copy"}
            assert prepared == NORMALIZED_TABLES | (AUDIT_TABLES if retain else set())
            assert set(plan.counts) == prepared
            assert plan.compatibility["retain_published_provenance"] is retain
            assert plan.compatibility["published_provenance_location"] == (
                "postgresql_and_pinned_parquet" if retain else "pinned_parquet"
            )
            assert plan.compatibility["published_input_counts"] == _published_counts(selected)
    assert plans[False].dimensions == plans[True].dimensions
    assert plans[False].sources == plans[True].sources
    assert {table: copied[True][table] for table in NORMALIZED_TABLES} == copied[False]
    assert all(
        plans[False].counts[table] == plans[True].counts[table] for table in NORMALIZED_TABLES
    )
    assert copied[False]["annotation_quantity"]
    assert copied[False]["entity_ontology_relation"]
    assert (
        len(copied[False]["relation"]) == 1
    )  # Distinct qualified claims share only the graph triple.


def test_opt_in_occurrences_retain_source_type_ordinals_and_duplicate_presentation(tmp_path):
    selected = _fixture(tmp_path)
    with duckdb.connect() as connection:
        plan = aligned_projection.prepare_aligned_release(
            connection, selected, retain_published_provenance=True
        )

        def rows(table):
            item = next(query for query in plan.queries if query.table == table)
            return [
                dict(zip(item.columns, row, strict=True))
                for row in connection.execute(item.query).fetchall()
            ]

        identifiers = sorted(
            (
                row
                for row in rows("parquet_identifier_occurrence")
                if row["resource"] == "record_owner"
            ),
            key=lambda row: row["ordinal"],
        )
        assert [
            (row["ordinal"], row["ns"], row["identifier"], row["is_canonical"], row["source"])
            for row in identifiers
        ] == [
            (0, "uniprot", "P04637", True, "source"),
            (1, "uniprot", "P04637", True, "source"),
            (2, "uniprot", "P04637", False, "other-identifier-source"),
        ]
        evidence = sorted(
            (
                row
                for row in rows("parquet_evidence")
                if row["resource"] == "record_owner" and row["relation_key"] != "qualified-second"
            ),
            key=lambda row: row["ordinal"],
        )
        assert [
            (
                row["ordinal"],
                row["original_row_id"],
                row["upstream_id"],
                row["source"],
                row["dataset"],
            )
            for row in evidence
        ] == [
            (0, "00001", "SRC-42", "reported-source", "reported-dataset"),
            (1, "source:row:001", "second-event", "other-source", None),
        ]
        repeated = [
            row
            for row in rows("parquet_annotation_occurrence")
            if row["resource"] == "record_owner"
            and row["owner_kind"] == "relation"
            and row["owner_key"] == "qualified-second"
        ]
        assert [row["ordinal"] for row in repeated] == [0, 1]
        assert len({row["annotation_key"] for row in repeated}) == 1
        assert {(row["scope"], row["source"], row["dataset"]) for row in repeated} == {
            ("object", "reported-source", "reported-dataset")
        }
        assert (
            Counter(row["entity_key"] for row in rows("parquet_entity")).most_common(1)[0][1] == 2
        )


@pytest.mark.parametrize("value", [None, 0, 1, "true"])
def test_projection_retention_requires_a_boolean_before_artifact_io(monkeypatch, value):
    def forbidden(*args, **kwargs):
        pytest.fail("Invalid retention must be rejected before reading resource paths")

    monkeypatch.setattr(aligned_projection, "_resources", forbidden)
    with pytest.raises(ValueError, match="retain_published_provenance"):
        aligned_projection.prepare_aligned_release(None, (), retain_published_provenance=value)
    with pytest.raises(ValueError, match="retain_published_provenance"):
        aligned_projection.companion_ddl("fixture", retain_published_provenance=value)


def test_default_companion_ddl_creates_required_molecular_and_quantity_tables():
    default = aligned_projection.companion_ddl('schema"quoted')
    retained = aligned_projection.companion_ddl('schema"quoted', retain_published_provenance=True)
    assert len(default) == 4
    assert default[-1] == retained[-1]
    assert '"schema""quoted"."annotation_quantity"' in default[-1]
    assert len(retained) == 9
    assert all(table not in statement for statement in default for table in AUDIT_TABLES)


@pytest.mark.parametrize("retain", [False, True])
def test_cli_passes_explicit_retention_choice(monkeypatch, capsys, retain):
    from omnipath_postgres import cli

    calls = []

    def load(*args, **kwargs):
        calls.append((args, kwargs))
        return aligned_loader.LoadResult("fixture", "fixture", "digest", {}, {}, {}, (), {})

    monkeypatch.setattr(cli, "load_release", load)
    arguments = [
        "fixture.json",
        "--data-root",
        "synthetic",
        "--database-url",
        "unused",
        "--schema",
        "fixture",
        "--base-only",
    ]
    if retain:
        arguments.append("--retain-published-provenance")
    assert cli.main(arguments) == 0
    assert calls[0][1]["retain_published_provenance"] is retain
    assert json.loads(capsys.readouterr().out)["layout"] == "main"


def test_cli_cannot_change_retention_of_an_existing_checkpoint(monkeypatch, capsys):
    from omnipath_postgres import cli

    def forbidden(*args, **kwargs):
        pytest.fail("Changing retention during finish must be rejected before connecting")

    monkeypatch.setattr(cli, "finish_release", forbidden)
    with pytest.raises(SystemExit) as failure:
        cli.main(["--finish", "--database-url", "unused", "--retain-published-provenance"])
    assert failure.value.code == 2
    assert "loading and source-audit options do not apply" in capsys.readouterr().err


def _pinned_fixture(root):
    resources = _fixture(root / "temporary")
    versions = {}
    for index, resource in enumerate(resources, 1):
        version = f"1.0.{index}"
        target = root / "resources" / resource.source / version
        target.parent.mkdir(parents=True, exist_ok=True)
        resource.directory.rename(target)
        versions[resource.source] = version
        write_manifest(
            target,
            resource.source,
            version,
            provenance={"dependencies": {"pypath-omnipath": "fixture"}},
        )
    path = root / "release.json"
    path.write_text(json.dumps(dict(schema_version=1, version="2026.10.2", resources=versions)))
    return path


def _pg_rows(connection, schema, tables):
    result = {}
    with connection.cursor() as cursor:
        for table in tables:
            cursor.execute(
                sql.SQL("SELECT * FROM {}.{}").format(sql.Identifier(schema), sql.Identifier(table))
            )
            result[table] = Counter(_freeze(row) for row in cursor.fetchall())
    connection.commit()
    return result


def test_local_default_loader_restores_main_contract_and_finishes_legacy_checkpoint(
    tmp_path,
    postgres_dsn,
    oracle,
    port,
    monkeypatch,
):
    path = _pinned_fixture(tmp_path)
    schemas = {retain: "retention_" + uuid.uuid4().hex for retain in (False, True)}
    baseline = "retention_main_" + uuid.uuid4().hex
    results = {}
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        with connection.cursor() as cursor:
            cursor.execute("CREATE EXTENSION IF NOT EXISTS pg_trgm")
        connection.commit()
        try:
            for retain, schema in schemas.items():
                kwargs = {"retain_published_provenance": True} if retain else {}
                observed = []
                results[retain] = aligned_loader.load_release(
                    tmp_path,
                    path,
                    postgres_dsn,
                    schema=schema,
                    base_only=True,
                    temp_directory=tmp_path / "spool",
                    observer=lambda event, **fields: observed.append((event, fields)),
                    **kwargs,
                )
                assert set(results[retain].counts) == NORMALIZED_TABLES | (
                    AUDIT_TABLES if retain else set()
                )
                assert {
                    fields["table"] for event, fields in observed if event == "copy_complete"
                } == set(results[retain].counts)
                with connection.cursor() as cursor:
                    cursor.execute("SELECT tablename FROM pg_tables WHERE schemaname=%s", (schema,))
                    tables = {row[0] for row in cursor.fetchall()}
                    assert tables & AUDIT_TABLES == (AUDIT_TABLES if retain else set())
                    assert {
                        "parquet_release",
                        "parquet_resource",
                        "parquet_phase",
                        "annotation_quantity",
                    } <= tables
                    cursor.execute(
                        sql.SQL("SELECT phase,result FROM {}.parquet_phase").format(
                            sql.Identifier(schema)
                        )
                    )
                    assert cursor.fetchall() == [("base", {"counts": results[retain].counts})]
                    cursor.execute(
                        sql.SQL("SELECT compatibility FROM {}.parquet_release").format(
                            sql.Identifier(schema)
                        )
                    )
                    assert cursor.fetchone()[0] == results[retain].compatibility
            assert _pg_rows(connection, schemas[False], NORMALIZED_TABLES) == _pg_rows(
                connection, schemas[True], NORMALIZED_TABLES
            )
            oracle.schema.ensure_schema(connection, schema=baseline)
            # Identical source IDs make actual partition bounds comparable.
            dimensions = aligned_loader._read_dimensions(connection, schemas[False])
            aligned_loader._write_dimensions(connection, baseline, dimensions)
            for _, source in dimensions["data_source"]:
                oracle.schema.ensure_source_partitions(connection, schema=baseline, source=source)
            oracle.indexes.create_secondary_indexes(connection, schema=baseline)
            reference = reference_contract._catalogue(connection, baseline)
            connection.commit()
            current = {}
            for retain, schema in schemas.items():
                current[retain] = reference_contract._catalogue(connection, schema)
                connection.commit()
            serving = {}
            for retain, schema in schemas.items():
                with connection.cursor() as cursor:
                    port.derive._create_derived_tables(cursor, schema)
                    port.derive._create_derived_indexes(cursor, schema)
                connection.commit()
                serving[retain] = _serving_rows(connection, schema, port.derive)
                assert serving[retain]
                aligned_loader._checkpoint(
                    connection,
                    schema,
                    "derived",
                    0.0,
                    {"fixture": True},
                    (results[retain].release, results[retain].manifest_sha256),
                    None,
                )
            assert serving[False] == serving[True]
            # A pre-policy twenty-table checkpoint has no new retention keys;
            # finish must keep its stored counts and never repeat COPY.
            schema = schemas[True]
            with connection.cursor() as cursor:
                cursor.execute(
                    sql.SQL("""UPDATE {}.parquet_release SET compatibility=compatibility
                    - 'retain_published_provenance' - 'published_provenance_location' - 'published_input_counts'""").format(
                        sql.Identifier(schema)
                    )
                )
            connection.commit()

            def forbidden(*args, **kwargs):
                pytest.fail("Finish of an existing derived checkpoint must not project or COPY")

            monkeypatch.setattr(aligned_loader, "prepare_aligned_release", forbidden)
            monkeypatch.setattr(aligned_loader, "_copy", forbidden)
            for retain, schema in schemas.items():
                finished = aligned_loader.finish_release(postgres_dsn, schema=schema, products=())
                assert finished.counts == results[retain].counts
                assert finished.resources == results[retain].resources
                assert finished.committed_products == ()
                assert finished.phase_seconds["derived"] == 0.0
                if not retain:
                    assert finished.compatibility == results[retain].compatibility
                else:
                    assert "retain_published_provenance" not in finished.compatibility
            for retain, catalogue in current.items():
                for kind in reference:
                    assert not reference[kind] - catalogue[kind], (
                        f"Missing frozen-main {kind}, retention={retain}: {reference[kind] - catalogue[kind]}"
                    )
        finally:
            connection.rollback()
            with connection.cursor() as cursor:
                for schema in (*schemas.values(), baseline):
                    cursor.execute(
                        sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(sql.Identifier(schema))
                    )
                    connection.commit()
