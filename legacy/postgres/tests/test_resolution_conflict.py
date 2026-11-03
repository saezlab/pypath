"""An unarbitrated disagreement leaves the mention unresolved and writes a
conflict record naming sources, identifiers and competing structures (spec
011 T051, WP2/User Story 2, acceptance scenario 3).

Same shape as ``test_authority_arbitration.py`` -- two genuinely different
skeletons on one mention -- but with **no** declared structure authority for
either candidate's identifier type. Nothing can arbitrate, so the mention
must stay unresolved, and the disagreement must not be silently dropped: a
queryable record is expected to name the sources, the identifiers and the
competing structures involved.
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
    UNRESOLVED_ID_TYPE,
    identifier_type_id,
)


INCHIKEY_TID = identifier_type_id(STANDARD_INCHI_KEY_TYPE)
CHEBI_TID = identifier_type_id(CHEBI_TYPE)
KEGG_TID = identifier_type_id(KEGG_COMPOUND_TYPE)
UNRESOLVED_TID = identifier_type_id(UNRESOLVED_ID_TYPE)

CHEBI_STRUCTURE = 'AWUCVROLDVIAJX-UHFFFAOYSA-N'
KEGG_STRUCTURE = 'WYOMHOATUARGQV-UHFFFAOYSA-N'


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
        SELECT canonical_identifier_type_id, canonical_identifier, status
        FROM entity_resolution
        WHERE entity_evidence_id = ?
        """,
        [ev_id],
    ).fetchone()


def test_unarbitrated_disagreement_stays_unresolved_and_is_recorded():
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
    # No identifier_authority row at all -- nothing can arbitrate.

    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id, status = _resolution(con, 'm1')
    assert status == 'unresolved'
    assert canonical_type == UNRESOLVED_TID

    conflicts = con.execute(
        """
        SELECT source, identifier_type, value_normalized,
               candidate_structure, candidate_source, skeleton
        FROM resolution_conflict
        WHERE entity_evidence_id = 'm1'
        ORDER BY candidate_structure
        """
    ).fetchall()
    assert len(conflicts) == 2

    skeletons = {row[5] for row in conflicts}
    assert skeletons == {CHEBI_STRUCTURE[:14], KEGG_STRUCTURE[:14]}

    structures = {row[3] for row in conflicts}
    assert structures == {CHEBI_STRUCTURE, KEGG_STRUCTURE}

    by_structure = {row[3]: row for row in conflicts}
    chebi_row = by_structure[CHEBI_STRUCTURE]
    assert chebi_row[0] == 'kegg'  # source of the mention
    assert chebi_row[1] == CHEBI_TYPE
    assert chebi_row[2] == 'CHEBI:16890'
    assert chebi_row[4] == 'kegg'  # candidate_source: who supplied this candidate

    kegg_row = by_structure[KEGG_STRUCTURE]
    assert kegg_row[1] == KEGG_COMPOUND_TYPE
    assert kegg_row[2] == 'C03189'
