"""When a record cannot be resolved because the reference data it needed was
absent, it is marked with the 'missing_translation_table' reason — distinct from
a record the resolver simply had no match for."""

from __future__ import annotations

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

from omnipath_build.duckdb_load import (  # noqa: E402
    PROTEIN_ENTITY_TYPE,
    ENSEMBL_TYPE,
    UNIPROT_TYPE,
    _canonicalize_loaded_duckdb,
    _create_duckdb_content_uuid_macro,
    _create_duckdb_evidence_tables,
)
from omnipath_build.resolver.identifier_types import (  # noqa: E402
    identifier_type_id,
)

UNIPROT_TID = identifier_type_id(UNIPROT_TYPE)

# Statically numbered resolution reasons (see the resolution-reason vocabulary).
MISSING_TRANSLATION_TABLE = 6
NO_ACCEPTED_RESOLVER_CANDIDATE = 3


def _con():
    con = duckdb.connect(':memory:')
    _create_duckdb_content_uuid_macro(con)
    _create_duckdb_evidence_tables(con)
    con.execute(
        'CREATE TABLE identifier_type (identifier_type_id BIGINT, name VARCHAR)'
    )
    for table in ('resolver_lookup', 'protein_uniprot_fallback_lookup'):
        con.execute(
            f"""
            CREATE TABLE {table} (
              entity_type VARCHAR,
              key_identifier_type_id BIGINT,
              key_value VARCHAR,
              taxonomy_id VARCHAR,
              canonical_identifier_type_id BIGINT,
              canonical_identifier VARCHAR
            )
            """
        )
    return con


def _insert_unresolvable_protein(con, ev_id='p'):
    # No resolver or fallback data is loaded, so this protein reaches the
    # ultimate-unresolved branch.
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('s', 'd', 1, ev_id, None, 'r', PROTEIN_ENTITY_TYPE, '9606'),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        ('s', ev_id, f'i_{ev_id}', ENSEMBL_TYPE, 'ENSP00000269305', None),
    )


def _reason(con, ev_id='p'):
    return con.execute(
        'SELECT reason_id FROM entity_evidence_resolution '
        'WHERE entity_evidence_id = ?',
        [ev_id],
    ).fetchone()[0]


def test_record_blocked_by_missing_table_is_marked():
    con = _con()
    # Simulate the pre-flight finding that no resolver relation covering proteins
    # was available.
    con.execute(
        'CREATE TABLE missing_translation_entity_type (entity_type VARCHAR)'
    )
    con.execute(
        'INSERT INTO missing_translation_entity_type VALUES (?)',
        [PROTEIN_ENTITY_TYPE],
    )
    _insert_unresolvable_protein(con)
    _canonicalize_loaded_duckdb(con)

    assert _reason(con) == MISSING_TRANSLATION_TABLE


def test_genuinely_unresolvable_record_is_not_marked():
    con = _con()
    # The resolver database was fully available (no blocked entity types); this
    # record simply had no match.
    _insert_unresolvable_protein(con)
    _canonicalize_loaded_duckdb(con)

    assert _reason(con) == NO_ACCEPTED_RESOLVER_CANDIDATE


def test_resolved_record_has_no_reason():
    con = _con()
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('s', 'd', 1, 'known', None, 'r', PROTEIN_ENTITY_TYPE, '9606'),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        ('s', 'known', 'i_known', UNIPROT_TYPE, 'P00533', None),
    )
    con.execute(
        'INSERT INTO protein_uniprot_fallback_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (PROTEIN_ENTITY_TYPE, UNIPROT_TID, 'P00533', '9606', UNIPROT_TID, 'P00533'),
    )
    _canonicalize_loaded_duckdb(con)

    assert _reason(con, 'known') is None
