"""Protein fallback canonicalization to primary UniProt."""

from __future__ import annotations

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

from omnipath_build.duckdb_load import (  # noqa: E402
    GENE_ENTITY_TYPE,
    PROTEIN_ENTITY_TYPE,
    UNIPROT_TYPE,
    _canonicalize_loaded_duckdb,
    _create_duckdb_content_uuid_macro,
    _create_duckdb_evidence_tables,
)
from omnipath_build.resolver.identifier_types import (  # noqa: E402
    UNRESOLVED_ID_TYPE,
    identifier_type_id,
)


UNIPROT_TID = identifier_type_id(UNIPROT_TYPE)
UNRESOLVED_TID = identifier_type_id(UNRESOLVED_ID_TYPE)


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


def _insert_protein(con, ev_id: str, accession: str) -> None:
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('s', 'd', 1, ev_id, None, 'r', PROTEIN_ENTITY_TYPE, '9606'),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        ('s', ev_id, f'i_{ev_id}', UNIPROT_TYPE, accession, None),
    )


def _resolution(con, ev_id: str):
    return con.execute(
        """
        SELECT entity_type, taxonomy_id, canonical_identifier_type_id,
               canonical_identifier, status, resolution_mechanism
        FROM entity_resolution
        WHERE entity_evidence_id = ?
        """,
        [ev_id],
    ).fetchone()


def test_protein_without_gene_candidate_falls_back_to_primary_uniprot():
    con = _con()
    _insert_protein(con, 'secondary_only', 'Q15086')
    con.execute(
        'INSERT INTO protein_uniprot_fallback_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (
            PROTEIN_ENTITY_TYPE,
            UNIPROT_TID,
            'Q15086',
            '9606',
            UNIPROT_TID,
            'P04637',
        ),
    )

    _canonicalize_loaded_duckdb(con)

    # An unresolved gene-product known only by a primary UniProt resolves to a
    # synthetic "unknown gene" (Gene, keyed by that UniProt), NOT a Protein —
    # so the base graph stays uniformly gene-typed. The AC is still carried as
    # a `protein` state (molecular_entity_type stays Protein).
    assert _resolution(con, 'secondary_only') == (
        GENE_ENTITY_TYPE,
        '9606',
        UNIPROT_TID,
        'P04637',
        'resolved',
        'unknown_gene',
    )


def test_gene_candidate_wins_over_protein_uniprot_fallback():
    con = _con()
    _insert_protein(con, 'gene_and_fallback', 'P00533')
    con.execute(
        'INSERT INTO resolver_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (GENE_ENTITY_TYPE, UNIPROT_TID, 'P00533', '9606', 3, '1956'),
    )
    con.execute(
        'INSERT INTO protein_uniprot_fallback_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (
            PROTEIN_ENTITY_TYPE,
            UNIPROT_TID,
            'P00533',
            '9606',
            UNIPROT_TID,
            'P00533',
        ),
    )

    _canonicalize_loaded_duckdb(con)

    assert _resolution(con, 'gene_and_fallback') == (
        GENE_ENTITY_TYPE,
        '9606',
        3,
        '1956',
        'resolved',
        'resolver',
    )


def test_ambiguous_protein_uniprot_fallback_stays_unresolved():
    con = _con()
    _insert_protein(con, 'ambiguous', 'PAMBIG')
    con.executemany(
        'INSERT INTO protein_uniprot_fallback_lookup VALUES (?, ?, ?, ?, ?, ?)',
        [
            (
                PROTEIN_ENTITY_TYPE,
                UNIPROT_TID,
                'PAMBIG',
                '9606',
                UNIPROT_TID,
                'P00001',
            ),
            (
                PROTEIN_ENTITY_TYPE,
                UNIPROT_TID,
                'PAMBIG',
                '9606',
                UNIPROT_TID,
                'P00002',
            ),
        ],
    )

    _canonicalize_loaded_duckdb(con)

    entity_type, taxonomy_id, canonical_type, _, status, mechanism = _resolution(
        con, 'ambiguous'
    )
    assert (entity_type, taxonomy_id, canonical_type, status, mechanism) == (
        PROTEIN_ENTITY_TYPE,
        '9606',
        UNRESOLVED_TID,
        'unresolved',
        'unresolved',
    )
