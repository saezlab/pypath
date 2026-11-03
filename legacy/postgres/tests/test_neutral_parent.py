"""The neutral form wins over charged forms of one skeleton, read from the
structure key's last character (spec 011 T049, WP2/User Story 2, acceptance
scenario 1).

Uses two of the four ATP variants from research.md R13 -- same connectivity
skeleton (``ZKHQWZAMYRWXGA``), one neutral (``...-N``), one charged
(``...-K``) -- so the tie-break under test (neutral beats charged) is the
only thing that can decide the winner.
"""

from __future__ import annotations

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

from omnipath_build.duckdb_load import (  # noqa: E402
    CHEBI_TYPE,
    CHEMICAL_ENTITY_TYPE,
    KEGG_COMPOUND_TYPE,
    STANDARD_INCHI_KEY_TYPE,
    _canonicalize_loaded_duckdb,
    _create_duckdb_content_uuid_macro,
    _create_duckdb_evidence_tables,
)
from omnipath_build.resolver.identifier_types import (  # noqa: E402
    identifier_type_id,
)


INCHIKEY_TID = identifier_type_id(STANDARD_INCHI_KEY_TYPE)

ATP_NEUTRAL = 'ZKHQWZAMYRWXGA-KQYNXXCUSA-N'
ATP_CHARGED = 'ZKHQWZAMYRWXGA-KQYNXXCUSA-K'


def _con():
    con = duckdb.connect(':memory:')
    _create_duckdb_content_uuid_macro(con)
    _create_duckdb_evidence_tables(con)
    con.execute(
        'CREATE TABLE identifier_type (identifier_type_id BIGINT, name VARCHAR)'
    )
    con.execute(
        """
        CREATE TABLE resolver_lookup (
          entity_type VARCHAR,
          key_identifier_type_id BIGINT,
          key_value VARCHAR,
          taxonomy_id VARCHAR,
          canonical_identifier_type_id BIGINT,
          canonical_identifier VARCHAR
        )
        """
    )
    con.execute(
        """
        CREATE TABLE protein_uniprot_fallback_lookup (
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


def _resolution(con, ev_id: str):
    return con.execute(
        """
        SELECT canonical_identifier_type_id, canonical_identifier, status
        FROM entity_resolution
        WHERE entity_evidence_id = ?
        """,
        [ev_id],
    ).fetchone()


def test_neutral_form_wins_over_charged_form_of_the_same_skeleton():
    con = _con()
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('reactome', 'reactions', 1, 'atp1', None, 'r', CHEMICAL_ENTITY_TYPE, None),
    )
    con.executemany(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        [
            ('reactome', 'atp1', 'i1', CHEBI_TYPE, '15422', 'CHEBI:15422'),
            ('reactome', 'atp1', 'i2', KEGG_COMPOUND_TYPE, 'C00002', 'C00002'),
        ],
    )
    con.executemany(
        'INSERT INTO resolver_lookup VALUES (?, ?, ?, ?, ?, ?)',
        [
            (
                CHEMICAL_ENTITY_TYPE,
                identifier_type_id(CHEBI_TYPE),
                'CHEBI:15422',
                None,
                INCHIKEY_TID,
                ATP_NEUTRAL,
            ),
            (
                CHEMICAL_ENTITY_TYPE,
                identifier_type_id(KEGG_COMPOUND_TYPE),
                'C00002',
                None,
                INCHIKEY_TID,
                ATP_CHARGED,
            ),
        ],
    )

    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id, status = _resolution(con, 'atp1')
    assert status == 'resolved'
    assert canonical_type == INCHIKEY_TID
    assert canonical_id == ATP_NEUTRAL
