"""Bounded full-result parity for the scalar canonical taxonomy aggregate."""

from collections import Counter

import duckdb

from omnipath_postgres.projection import (
    ENTITY_DATASET,
    PUBLISHED_FALLBACK_NAMESPACE,
    _create_entity_projection,
)


# The former canonical input query is an independent oracle. All later entity
# projection steps are executed unchanged in both connections, so this compares
# complete output rows rather than just a taxon scalar or row count.
OLD_CANONICAL_INPUT = """CREATE OR REPLACE TABLE ap_entity_canonical_input AS SELECT
    entity_key,entity_type,namespace,identifier,
    CASE WHEN count(DISTINCT try_cast(taxon AS BIGINT))=1
         THEN min(try_cast(taxon AS BIGINT)) ELSE NULL END taxonomy_id
    FROM ap_entity_raw GROUP BY entity_key,entity_type,namespace,identifier"""


class OldAggregateConnection:
    def __init__(self, connection):
        self.connection = connection
        self.replacements = 0

    def execute(self, query, *args):
        if query.startswith("CREATE OR REPLACE TABLE ap_entity_canonical_input AS"):
            self.replacements += 1
            query = OLD_CANONICAL_INPUT
        return self.connection.execute(query, *args)

    def __getattr__(self, name):
        return getattr(self.connection, name)


def _inputs(connection):
    connection.execute("CREATE MACRO ap_uuid(value) AS md5(value)::UUID")
    cases = {
        "null_only": (None, None),
        "empty_only": ("", ""),
        "invalid_only": ("invalid", "NCBITaxon:9606"),
        "one_known": ("9606", "", None, "invalid"),
        "same_known": ("09606", "9606", "9606"),
        "conflicting_known": ("9606", "10090", None, ""),
        "bigint_extrema": ("-9223372036854775808", "9223372036854775807"),
        "outside_bigint": ("9223372036854775808",),
    }
    # Twenty published occurrences; repeated keys overlap distinct resources.
    # Invalid values are deliberate expression-level cases: release validation
    # still rejects nonempty taxon strings that cannot be represented as BIGINT.
    rows = [
        (f"source_{index}", "1", "published:" + key, "protein", "uniprot", key, taxon)
        for key, taxa in cases.items()
        for index, taxon in enumerate(taxa)
    ]
    assert len(rows) == 20
    connection.execute("""CREATE TABLE ap_entity_raw (
        resource VARCHAR, version VARCHAR, entity_key VARCHAR, entity_type VARCHAR,
        namespace VARCHAR, identifier VARCHAR, taxon VARCHAR)""")
    connection.executemany("INSERT INTO ap_entity_raw VALUES (?,?,?,?,?,?,?)", rows)
    connection.execute("CREATE TABLE ap_namespace_alias(raw VARCHAR,name VARCHAR)")
    connection.execute("CREATE TABLE ap_vocab_entity_type(id BIGINT,name VARCHAR)")
    connection.execute("INSERT INTO ap_vocab_entity_type VALUES (1,'protein')")
    connection.execute("CREATE TABLE ap_vocab_identifier_type(id BIGINT,name VARCHAR)")
    connection.executemany(
        "INSERT INTO ap_vocab_identifier_type VALUES (?,?)",
        [(1, "uniprot"), (2, PUBLISHED_FALLBACK_NAMESPACE)],
    )
    connection.execute("""CREATE TABLE ap_data_source AS
        SELECT row_number() OVER (ORDER BY resource)::BIGINT AS id,resource AS name
        FROM (SELECT DISTINCT resource FROM ap_entity_raw)""")
    connection.execute("CREATE TABLE ap_dataset(id BIGINT,source_id BIGINT,name VARCHAR)")
    connection.execute(
        "INSERT INTO ap_dataset SELECT id,id,? FROM ap_data_source", [ENTITY_DATASET]
    )
    connection.execute("CREATE TABLE ap_identifier_raw(item STRUCT(ns VARCHAR,id VARCHAR))")
    connection.execute(
        "CREATE TABLE ap_annotation_raw(owner_kind VARCHAR,owner_key VARCHAR,term VARCHAR,value VARCHAR)"
    )


def _complete_rows(connection):
    tables = (
        "ap_entity_canonical_input",
        "ap_entity_natural_conflict",
        "ap_entity",
        "ap_entity_occurrence",
        "ap_identifier_input",
        "ap_identifier_alias",
        "ap_identifier",
        "ap_identifier_occurrence",
    )
    return {
        table: Counter(connection.execute(f"SELECT * FROM {table}").fetchall()) for table in tables
    }


def test_min_max_taxonomy_matches_former_distinct_query_for_complete_entity_outputs():
    with duckdb.connect() as old, duckdb.connect() as updated:
        _inputs(old)
        _inputs(updated)
        oracle = OldAggregateConnection(old)
        _create_entity_projection(oracle)
        _create_entity_projection(updated)
        assert oracle.replacements == 1
        assert _complete_rows(updated) == _complete_rows(old)
        assert dict(
            updated.execute(
                "SELECT identifier,taxonomy_id FROM ap_entity_canonical_input"
            ).fetchall()
        ) == {
            "null_only": None,
            "empty_only": None,
            "invalid_only": None,
            "one_known": 9606,
            "same_known": 9606,
            "conflicting_known": None,
            "bigint_extrema": None,
            "outside_bigint": None,
        }
        assert updated.execute("SELECT count(*) FROM ap_entity_occurrence").fetchone() == (20,)
        assert updated.execute("SELECT count(*) FROM ap_entity").fetchone() == (8,)
