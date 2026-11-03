"""The worked KEGG C03189 example resolves, traceable to the arbitrating
authority (spec 011 T052, WP2/User Story 2, acceptance scenario 4).

From research/chemical-resolution-state.md 2.3 -- KEGG C03189 (DL-glycerol
1-phosphate) carries three raw candidates today and is dropped as ambiguous:

| Source of the candidate | InChIKey |
|---|---|
| PubChem CID 6068 | ``WYOMHOATUARGQV-UHFFFAOYSA-N`` |
| KEGG C03189 through the ChEBI hub | ``AWUCVROLDVIAJX-UHFFFAOYSA-N`` |
| KEGG C03189 through the ChEBI hub | ``AWUCVROLDVIAJX-UHFFFAOYSA-L`` |

Two of the three share one skeleton (a protonation-state variant, not a real
disagreement); the third is a genuine PubChem/ChEBI cross-reference
disagreement. After skeleton collapse there are two real candidates, and the
record's own ChEBI identifier (16890) is the structure authority -- so this
must resolve, not drop, with a mechanism traceable to that authority.
"""

from __future__ import annotations

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

from omnipath_build.duckdb_load import (  # noqa: E402
    CHEBI_TYPE,
    CHEMICAL_ENTITY_TYPE,
    KEGG_COMPOUND_TYPE,
    PUBCHEM_COMPOUND_TYPE,
    STANDARD_INCHI_KEY_TYPE,
    _canonicalize_loaded_duckdb,
    _create_duckdb_content_uuid_macro,
    _create_duckdb_evidence_tables,
)
from omnipath_build.resolver.identifier_types import (  # noqa: E402
    identifier_type_id,
)


INCHIKEY_TID = identifier_type_id(STANDARD_INCHI_KEY_TYPE)
CHEBI_TID = identifier_type_id(CHEBI_TYPE)
KEGG_TID = identifier_type_id(KEGG_COMPOUND_TYPE)
PUBCHEM_TID = identifier_type_id(PUBCHEM_COMPOUND_TYPE)

CHEBI_N = 'AWUCVROLDVIAJX-UHFFFAOYSA-N'
CHEBI_L = 'AWUCVROLDVIAJX-UHFFFAOYSA-L'
PUBCHEM_STRUCTURE = 'WYOMHOATUARGQV-UHFFFAOYSA-N'

CHEBI_SOURCE_ID = 1


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
    con.execute(
        """
        CREATE TABLE identifier_authority (
          identifier_type_id BIGINT,
          source_id BIGINT,
          is_structure_authority BOOLEAN
        )
        """
    )
    return con


def _resolution(con, ev_id: str):
    return con.execute(
        """
        SELECT canonical_identifier_type_id, canonical_identifier, status,
               resolution_mechanism
        FROM entity_resolution
        WHERE entity_evidence_id = ?
        """,
        [ev_id],
    ).fetchone()


def test_kegg_c03189_resolves_via_the_chebi_authority():
    con = _con()
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('kegg', 'reactions', 1, 'c03189', None, 'r', CHEMICAL_ENTITY_TYPE, None),
    )
    con.executemany(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        [
            ('kegg', 'c03189', 'i1', KEGG_COMPOUND_TYPE, 'C03189', 'C03189'),
            ('kegg', 'c03189', 'i2', CHEBI_TYPE, '16890', 'CHEBI:16890'),
            ('kegg', 'c03189', 'i3', PUBCHEM_COMPOUND_TYPE, '6068', '6068'),
        ],
    )
    con.executemany(
        'INSERT INTO resolver_lookup VALUES (?, ?, ?, ?, ?, ?)',
        [
            # ChEBI 16890 -- two protonation states of the same skeleton.
            (
                CHEMICAL_ENTITY_TYPE,
                CHEBI_TID,
                'CHEBI:16890',
                None,
                INCHIKEY_TID,
                CHEBI_N,
            ),
            (
                CHEMICAL_ENTITY_TYPE,
                CHEBI_TID,
                'CHEBI:16890',
                None,
                INCHIKEY_TID,
                CHEBI_L,
            ),
            # PubChem CID 6068 -- a different skeleton entirely.
            (
                CHEMICAL_ENTITY_TYPE,
                PUBCHEM_TID,
                '6068',
                None,
                INCHIKEY_TID,
                PUBCHEM_STRUCTURE,
            ),
        ],
    )
    con.execute(
        'INSERT INTO identifier_authority VALUES (?, ?, ?)',
        (CHEBI_TID, CHEBI_SOURCE_ID, True),
    )

    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id, status, mechanism = _resolution(con, 'c03189')
    assert status == 'resolved'
    assert canonical_type == INCHIKEY_TID
    assert canonical_id in (CHEBI_N, CHEBI_L)
    assert mechanism == 'authority_arbitrated'
