"""Independent deployed-SQL parity for COSMOS's endpoint direction lookup."""

from oracles import read_deployed_cosmos

from collections import Counter
from contextlib import closing
from dataclasses import asdict
import json
from pathlib import Path
import time
import uuid

import psycopg2
from psycopg2 import sql
import pytest

from omnipath_subsets.cosmos import build
from omnipath_postgres import projection
from omnipath_postgres.relational.db import derived_tables
from omnipath_postgres.relational.db import schema as main_schema


REPO = Path(__file__).resolve().parents[3]
EVENTS = ("reversible", "forward", "right_to_left", "silent")
IDS = {
    name: str(uuid.UUID(int=6000 + index))
    for index, name in enumerate((*EVENTS, "x", "y", "enzyme"))
}
SOURCES = ((9001, "direction_a"), (9002, "direction_b"))
PARAMS = {
    "reaction_entity_types": list(build.REACTION_ENTITY_TYPES),
    "reaction_sources": [name for _, name in SOURCES],
    "direction_term": build.DIRECTION_TERM,
}


@pytest.fixture(scope="module")
def deployed_sql():
    return read_deployed_cosmos()


@pytest.fixture(scope="module")
def direction_fixture(postgres_dsn):
    """20 source records; all derived headers/parties are declared explicitly."""
    name = "cosmos_direction_" + uuid.uuid4().hex
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            with connection.cursor() as cur:
                cur.execute("CREATE EXTENSION IF NOT EXISTS pg_trgm")
            main_schema.ensure_schema(connection, schema=name)
            with connection.cursor() as cur:
                derived_tables._create_derived_tables(cur, name)
                cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(name)))
                cur.executemany("INSERT INTO data_source(source_id,name) VALUES (%s,%s)", SOURCES)
                cur.executemany(
                    "INSERT INTO dataset(dataset_id,source_id,name) VALUES (%s,%s,'fixture')",
                    [(source, source) for source, _ in SOURCES],
                )
                for typ in ("molecular_activity", "small_molecule", "protein"):
                    cur.execute(
                        "INSERT INTO vocab_entity_type(name) VALUES (%s) ON CONFLICT DO NOTHING",
                        [typ],
                    )
                cur.execute("SELECT name,entity_type_id FROM vocab_entity_type")
                types = dict(cur.fetchall())
                cur.execute("SELECT name,identifier_type_id FROM vocab_identifier_type")
                namespaces = dict(cur.fetchall())
                for entity, identifier in IDS.items():
                    typ = (
                        "molecular_activity"
                        if entity in EVENTS
                        else "protein"
                        if entity == "enzyme"
                        else "small_molecule"
                    )
                    namespace = (
                        namespaces["Uniprot:MI:1097"]
                        if entity == "enzyme"
                        else namespaces["Chebi:MI:0474"]
                        if entity in ("x", "y")
                        else None
                    )
                    value = {"x": "123", "y": "456", "enzyme": "P12345"}.get(entity, entity)
                    cur.execute(
                        "INSERT INTO entity(entity_id,entity_type_id,canonical_identifier_type_id,"
                        "canonical_identifier,resolution_status_id) VALUES (%s,%s,%s,%s,2)",
                        [identifier, types[typ], namespace, value],
                    )
                predicates = ("has_input", "has_output", "enabled_by", "catalyzes", "part_of")
                cur.executemany(
                    "INSERT INTO vocab_relation_predicate(name) VALUES (%s) ON CONFLICT DO NOTHING",
                    [(predicate,) for predicate in predicates],
                )
                cur.execute("SELECT name,relation_predicate_id FROM vocab_relation_predicate")
                predicate_ids = dict(cur.fetchall())
                cur.execute(
                    "INSERT INTO vocab_relation_category(name) VALUES ('interaction') ON CONFLICT DO NOTHING"
                )
                cur.execute(
                    "SELECT relation_category_id FROM vocab_relation_category WHERE name='interaction'"
                )
                category = cur.fetchone()[0]
                cur.execute("SELECT relation_role_id,name FROM vocab_relation_role")
                roles = {role: identifier for identifier, role in cur.fetchall()}
                cur.execute(
                    "INSERT INTO vocab_interaction_class(interaction_class_id,name) VALUES (100,'other') "
                    "ON CONFLICT DO NOTHING"
                )
                cur.execute(
                    "SELECT interaction_class_id FROM vocab_interaction_class WHERE name='other'"
                )
                interaction_class = cur.fetchone()[0]
                for index, event in enumerate(EVENTS):
                    interaction = str(uuid.UUID(int=7000 + index))
                    cur.execute(
                        "INSERT INTO interaction(interaction_id,interaction_class_id,arity,sources) "
                        "VALUES (%s,%s,%s,%s)",
                        [
                            interaction,
                            interaction_class,
                            2 if event == "silent" else 3,
                            [source for _, source in SOURCES],
                        ],
                    )
                    for source, _ in SOURCES:
                        cur.execute(
                            "INSERT INTO interaction_fact_resource(interaction_fact_resource_id,"
                            "subject_entity_id,object_entity_id,interaction_class_id,source_id,interaction_id) "
                            "VALUES (%s,%s,%s,%s,%s,%s)",
                            [
                                str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([event, source]))),
                                IDS[event],
                                IDS["x"],
                                interaction_class,
                                source,
                                interaction,
                            ],
                        )
                    parties = [("x", "reactant", "c"), ("y", "product", "m")]
                    if event != "silent":
                        parties.append(("enzyme", "enzyme", "c"))
                    for ordinal, (member, role, compartment) in enumerate(parties):
                        cur.execute(
                            "INSERT INTO interaction_party(interaction_id,entity_id,role_id,ordinal,compartment) "
                            "VALUES (%s,%s,%s,%s,%s)",
                            [interaction, IDS[member], roles[role], ordinal, compartment],
                        )
                evidence_ids = {}
                for index, member in enumerate(("enzyme", "x", "reversible")):
                    evidence_ids[member] = str(uuid.UUID(int=8000 + index))
                    cur.execute(
                        "INSERT INTO entity_evidence(source_id,entity_evidence_id,dataset_id,row_id,"
                        "entity_role_id,entity_type_id) SELECT 9001,%s,9001,%s,1,entity_type_id "
                        "FROM entity WHERE entity_id=%s",
                        [evidence_ids[member], 8000 + index, IDS[member]],
                    )
                    cur.execute(
                        "INSERT INTO entity_evidence_resolution(source_id,entity_evidence_id,status_id,entity_id) "
                        "VALUES (9001,%s,2,%s)",
                        [evidence_ids[member], IDS[member]],
                    )

                def annotation(value, term=build.DIRECTION_TERM):
                    key = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([term, value])))
                    cur.execute(
                        "INSERT INTO annotation(annotation_key,term,value) VALUES (%s,%s,%s) "
                        "ON CONFLICT DO NOTHING",
                        [key, term, value],
                    )
                    return key

                # Event annotation loses to an explicit reversible membership.
                cur.execute(
                    "INSERT INTO entity_evidence_annotation VALUES (9001,%s,%s)",
                    [evidence_ids["reversible"], annotation("LEFT-TO-RIGHT")],
                )
                claims = [
                    (9001, 8100, "reversible", "x", "has_input", "REVERSIBLE"),
                    # Same evidence UUID, other owner/event: annotations cannot leak.
                    (9002, 8100, "silent", "x", "has_input", "unrecognized"),
                    (9001, 8101, "forward", "x", "has_input", "LEFT-TO-RIGHT"),
                    (9001, 8102, "enzyme", "right_to_left", "catalyzes", "RIGHT-TO-LEFT"),
                    (9001, 8103, "reversible", "reversible", "has_output", "reversible"),
                    (9001, 8104, None, "right_to_left", "enabled_by", "right_to_left"),
                    (9001, 8105, "forward", None, "has_output", "left_to_right"),
                    # One assertion legitimately matches each of two event owners.
                    (9001, 8106, "forward", "right_to_left", "enabled_by", "RIGHT_TO_LEFT"),
                    (9001, 8107, "silent", "x", "part_of", "REVERSIBLE"),
                    (9001, 8108, None, None, "has_input", "REVERSIBLE"),
                ]
                assert len(IDS) + len(evidence_ids) + len(claims) == 20
                for source, number, subject, obj, predicate, value in claims:
                    evidence = str(uuid.UUID(int=number))
                    cur.execute(
                        "INSERT INTO relation_evidence(source_id,relation_evidence_id,dataset_id,row_id,"
                        "subject_entity_evidence_id,subject_entity_id,predicate_id,"
                        "object_entity_evidence_id,object_entity_id,relation_category_id) "
                        "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
                        [
                            source,
                            evidence,
                            source,
                            number,
                            evidence_ids["enzyme"] if subject is None else None,
                            IDS[subject] if subject is not None else None,
                            predicate_ids[predicate],
                            evidence_ids["x"] if obj is None else None,
                            IDS[obj] if obj is not None else None,
                            category,
                        ],
                    )
                    cur.execute(
                        "INSERT INTO relation_evidence_annotation VALUES (%s,%s,%s,1)",
                        [source, evidence, annotation(value)],
                    )
                # Duplicate claims through different annotation scopes survive
                # before bool_or; self-loop matching must still happen once.
                cur.execute(
                    "INSERT INTO relation_evidence_annotation VALUES (9001,%s,%s,3)",
                    [str(uuid.UUID(int=8100)), annotation("REVERSIBLE")],
                )
                cur.execute(
                    "INSERT INTO relation_evidence_annotation VALUES (9002,%s,%s,1)",
                    [str(uuid.UUID(int=8100)), annotation(None)],
                )
                cur.execute(
                    "INSERT INTO relation_evidence_annotation VALUES (9002,%s,%s,1)",
                    [str(uuid.UUID(int=8100)), annotation("REVERSIBLE", "unrelated:term")],
                )
                # The loader creates it beside the main tables; COSMOS reads the
                # source accession of a reaction from it.
                cur.execute(
                    "CREATE TABLE entity_reference_context("
                    + projection.COMPANION_DDL["entity_reference_context"]
                    + ")"
                )
                cur.execute("CREATE TABLE build_manifest(build_id text NOT NULL)")
                cur.execute("INSERT INTO build_manifest VALUES ('direction_fixture_v1')")
                cur.execute("RESET search_path")
            connection.commit()
            yield connection, name
        finally:
            connection.rollback()
            with connection.cursor() as cur:
                cur.execute(
                    sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(sql.Identifier(name))
                )
            connection.commit()


def _claims_query(text):
    begin = text.index("WITH reaction_event AS (")
    cte = text.index("reaction_direction AS (", begin)
    claim_begin = text.index("FROM (", cte) + len("FROM (")
    claim_end = text.index(") claims GROUP BY reaction_entity_id", claim_begin)
    return (
        text[begin:cte] + "exact_claims AS (" + text[claim_begin:claim_end] + ") "
        "SELECT reaction_entity_id::text,value FROM exact_claims"
    )


def test_complete_claim_multiset_preserves_null_self_loop_and_source_ownership(
    direction_fixture, deployed_sql
):
    connection, name = direction_fixture
    expected = Counter(
        [
            (IDS["reversible"], "LEFT-TO-RIGHT"),
            (IDS["reversible"], "REVERSIBLE"),
            (IDS["reversible"], "REVERSIBLE"),
            (IDS["reversible"], "reversible"),
            (IDS["forward"], "LEFT-TO-RIGHT"),
            (IDS["forward"], "left_to_right"),
            (IDS["forward"], "RIGHT_TO_LEFT"),
            (IDS["right_to_left"], "RIGHT-TO-LEFT"),
            (IDS["right_to_left"], "right_to_left"),
            (IDS["right_to_left"], "RIGHT_TO_LEFT"),
            (IDS["silent"], "unrecognized"),
            (IDS["silent"], None),
        ]
    )
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(name)))
        for text in (deployed_sql, build._sql_text("project_edges.sql")):
            cur.execute(_claims_query(text), PARAMS)
            assert Counter(cur.fetchall()) == expected
        cur.execute("RESET search_path")


def _complete_rows(connection, name):
    with connection.cursor() as cur:
        # Sequence and every published field, including insertion surrogate,
        # are compared; UUIDs/labels/reaction numbering are never normalized.
        cur.execute(
            sql.SQL("SELECT * FROM {}.cosmos_edge ORDER BY cosmos_edge_id").format(
                sql.Identifier(name)
            )
        )
        columns = [description.name for description in cur.description]
        rows = cur.fetchall()
        return columns, rows


def test_direction_claims_have_one_narrow_producer_and_two_endpoint_scans(direction_fixture):
    connection, name = direction_fixture
    with connection.cursor() as cur:
        cur.execute(sql.SQL("SET search_path={},public").format(sql.Identifier(name)))
        cur.execute(
            "EXPLAIN (FORMAT JSON,VERBOSE) " + _claims_query(build._sql_text("project_edges.sql")),
            PARAMS,
        )
        plan = cur.fetchone()[0][0]["Plan"]
        cur.execute("RESET search_path")

    def nodes(node):
        yield node
        for child in node.get("Plans", []):
            yield from nodes(child)

    entries = list(nodes(plan))
    producers = [
        entry for entry in entries if entry.get("Subplan Name") == "CTE membership_direction"
    ]
    assert len(producers) == 1
    assert len(producers[0]["Output"]) == 3
    assert len([entry for entry in entries if entry.get("CTE Name") == "membership_direction"]) == 2


def test_complete_publication_rows_order_and_statistics_match_deployed(
    direction_fixture, deployed_sql, monkeypatch
):
    connection, name = direction_fixture
    read_sql = build._sql_text
    observations = []
    for deployed in (True, False):
        monkeypatch.setattr(
            build,
            "_sql_text",
            lambda filename: deployed_sql
            if deployed and filename == "project_edges.sql"
            else read_sql(filename),
        )
        started = time.perf_counter()
        stats = build.build_cosmos_projection(
            connection, schema=name, reaction_sources=PARAMS["reaction_sources"], utils_db_url=None
        )
        columns, rows = _complete_rows(connection, name)
        observations.append(
            (
                columns,
                rows,
                {key: value for key, value in asdict(stats).items() if key != "seconds"},
            )
        )
        print(
            f"direction lookup {'deployed' if deployed else 'candidate'} complete fixture: "
            f"{time.perf_counter() - started:.6f}s"
        )
    assert observations[0] == observations[1]
    columns, rows, stats = observations[1]
    decoded = [dict(zip(columns, row, strict=True)) for row in rows]
    reactions = [row for row in decoded if row["interaction_type"] != "connector"]
    assert stats["reactions"] == 4
    assert stats["orphan_reactions"] == 1
    assert stats["reversible_reactions"] == 1
    assert stats["edges"] == 10 and stats["reverse_edges"] == 2
    assert stats["connectors"] > 0
    directions = {str(row["reaction_entity_id"]): row["direction"] for row in reactions}
    assert directions == {
        IDS["reversible"]: "reversible",
        IDS["forward"]: "left_to_right",
        IDS["right_to_left"]: "left_to_right",
        IDS["silent"]: None,
    }
    for row in reactions:
        assert row["reverse"] == ("_rev" in row["source_label"] or "_rev" in row["target_label"])
        assert not row["reverse"] or row["direction"] == "reversible"
        from_metabolite = row["source_type"] == "metabolite"
        metabolite = row["source_entity_id"] if from_metabolite else row["target_entity_id"]
        assert from_metabolite == ((str(metabolite) == IDS["x"]) != row["reverse"])
