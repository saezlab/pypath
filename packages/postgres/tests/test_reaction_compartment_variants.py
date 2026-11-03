"""A reaction stated in several compartments is one entity per compartment set."""

import duckdb
import pytest

from omnipath_postgres.projection import (
    ENTITY_DATASET,
    PUBLISHED_FALLBACK_NAMESPACE,
    _create_entity_projection,
)


def _projection(entities, compartments):
    """``entities``: (entity_key, entity_type, namespace, identifier) published by one
    resource; ``compartments``: (entity_key, GO term) entity annotations."""
    con = duckdb.connect()
    con.execute("CREATE MACRO ap_uuid(value) AS md5(value)::UUID")
    con.execute("""CREATE TABLE ap_entity_raw (resource VARCHAR, version VARCHAR,
        entity_key VARCHAR, entity_type VARCHAR, namespace VARCHAR, identifier VARCHAR,
        taxon VARCHAR)""")
    con.executemany(
        "INSERT INTO ap_entity_raw VALUES ('recon3d','1',?,?,?,?,NULL)", entities
    )
    con.execute(
        "CREATE TABLE ap_annotation_raw(owner_kind VARCHAR,owner_key VARCHAR,term VARCHAR,value VARCHAR)"
    )
    con.executemany(
        "INSERT INTO ap_annotation_raw VALUES ('entity',?,'anatomical_context_qualifier',?)",
        compartments,
    )
    con.execute("CREATE TABLE ap_namespace_alias(raw VARCHAR,name VARCHAR)")
    con.execute("CREATE TABLE ap_vocab_entity_type(id BIGINT,name VARCHAR)")
    con.execute("INSERT INTO ap_vocab_entity_type VALUES (1,'molecular_activity'),(2,'protein')")
    con.execute("CREATE TABLE ap_vocab_identifier_type(id BIGINT,name VARCHAR)")
    con.execute(
        "INSERT INTO ap_vocab_identifier_type VALUES (1,'metanetx_reaction'),(2,'uniprot'),(3,?)",
        [PUBLISHED_FALLBACK_NAMESPACE],
    )
    con.execute("CREATE TABLE ap_data_source AS SELECT 1::BIGINT AS id, 'recon3d' AS name")
    con.execute("CREATE TABLE ap_dataset(id BIGINT,source_id BIGINT,name VARCHAR)")
    con.execute("INSERT INTO ap_dataset VALUES (1,1,?)", [ENTITY_DATASET])
    con.execute("CREATE TABLE ap_identifier_raw(item STRUCT(ns VARCHAR,id VARCHAR))")
    _create_entity_projection(con)
    return con


def test_compartment_variants_take_their_published_keys():
    con = _projection(
        [
            ("cytosol", "molecular_activity", "metanetx_reaction", "MNXR100050"),
            ("lysosome", "molecular_activity", "metanetx_reaction", "MNXR100050"),
            ("unstated", "molecular_activity", "metanetx_reaction", "MNXR100050"),
            ("alone", "molecular_activity", "metanetx_reaction", "MNXR100058"),
        ],
        [("cytosol", "GO:0005829"), ("lysosome", "GO:0005764")],
    )
    rows = dict(
        con.execute(
            """SELECT e.canonical_identifier, t.name FROM ap_entity e
            JOIN ap_vocab_identifier_type t ON t.id=e.canonical_identifier_type_id"""
        ).fetchall()
    )
    assert rows == {
        "cytosol": PUBLISHED_FALLBACK_NAMESPACE,
        "lysosome": PUBLISHED_FALLBACK_NAMESPACE,
        "unstated": PUBLISHED_FALLBACK_NAMESPACE,
        "MNXR100058": "metanetx_reaction",
    }


@pytest.mark.parametrize(
    "entities,compartments",
    [
        # Two reactions in the same compartment are one identity split in two.
        (
            [
                ("a", "molecular_activity", "metanetx_reaction", "MNXR1"),
                ("b", "molecular_activity", "metanetx_reaction", "MNXR1"),
            ],
            [("a", "GO:0005829"), ("b", "GO:0005829")],
        ),
        # A protein accession is never shared by two entities.
        (
            [("a", "protein", "uniprot", "P00533"), ("b", "protein", "uniprot", "P00533")],
            [("a", "GO:0005829"), ("b", "GO:0005764")],
        ),
    ],
)
def test_other_shared_accessions_still_fail(entities, compartments):
    with pytest.raises(ValueError, match="natural-key uniqueness"):
        _projection(entities, compartments)
