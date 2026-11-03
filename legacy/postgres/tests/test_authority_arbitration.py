"""The minting authority's structure wins when candidate skeletons genuinely
differ, with no resource named in the code (spec 011 T050, WP2/User Story 2,
acceptance scenario 2).

A mention citing both a ChEBI id and a KEGG id, whose candidate structures
have *different* connectivity skeletons (a real cross-reference
disagreement, not a protonation/stereo variant of one molecule) -- ChEBI is
declared the structure authority for its own identifier type
(``identifier_authority.is_structure_authority``), so its candidate must
win. Nothing in the resolution code may name "chebi" to make this happen --
the win has to come from the authority table alone.
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
CHEBI_TID = identifier_type_id(CHEBI_TYPE)
KEGG_TID = identifier_type_id(KEGG_COMPOUND_TYPE)

# Two genuinely different skeletons -- not variants of one molecule.
CHEBI_STRUCTURE = 'AWUCVROLDVIAJX-UHFFFAOYSA-N'
KEGG_STRUCTURE = 'WYOMHOATUARGQV-UHFFFAOYSA-N'

# A source_id for the fixture's data_source rows -- not read by the
# resolution logic under test (which must never name a resource), only
# needed to satisfy identifier_authority's foreign key.
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


def test_structure_authority_wins_a_genuine_cross_skeleton_disagreement():
    con = _con()
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('kegg', 'reactions', 1, 'm1', None, 'r', CHEMICAL_ENTITY_TYPE, None),
    )
    con.executemany(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        [
            ('kegg', 'm1', 'i1', CHEBI_TYPE, '16890', 'CHEBI:16890'),
            ('kegg', 'm1', 'i2', KEGG_COMPOUND_TYPE, 'C03189', 'C03189'),
        ],
    )
    con.executemany(
        'INSERT INTO resolver_lookup VALUES (?, ?, ?, ?, ?, ?)',
        [
            (
                CHEMICAL_ENTITY_TYPE,
                CHEBI_TID,
                'CHEBI:16890',
                None,
                INCHIKEY_TID,
                CHEBI_STRUCTURE,
            ),
            (
                CHEMICAL_ENTITY_TYPE,
                KEGG_TID,
                'C03189',
                None,
                INCHIKEY_TID,
                KEGG_STRUCTURE,
            ),
        ],
    )
    # ChEBI is the declared structure authority for its own identifier type;
    # KEGG is not declared for either type here.
    con.execute(
        'INSERT INTO identifier_authority VALUES (?, ?, ?)',
        (CHEBI_TID, CHEBI_SOURCE_ID, True),
    )

    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id, status, mechanism = _resolution(con, 'm1')
    assert status == 'resolved'
    assert canonical_type == INCHIKEY_TID
    assert canonical_id == CHEBI_STRUCTURE
    assert mechanism == 'authority_arbitrated'
