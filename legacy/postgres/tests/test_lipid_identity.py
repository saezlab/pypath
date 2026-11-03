"""A lipid mention canonicalizes on its standardized nomenclature name, never
on structure (spec 011 WP7/US7, research R9, T116).

Needs the real utils Postgres attached (``OMNIPATH_BUILD_UTILS_PG_URL``) with
its ``lipid_name`` table populated (T113) -- an integration test, not a pure
unit test, since the resolution arm under test degrades to a no-op without a
live attach (verified separately by ``test_authority_arbitration.py``
continuing to pass with no utils DB configured at all).
"""

from __future__ import annotations

import os

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

UTILS_URL = os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL')
pytestmark = pytest.mark.skipif(
    not UTILS_URL,
    reason='set OMNIPATH_BUILD_UTILS_PG_URL to a built utils instance to run',
)

from omnipath_build.duckdb_load import (  # noqa: E402
    CHEBI_TYPE,
    CHEMICAL_ENTITY_TYPE,
    LIPID_NAME_TYPE,
    STANDARD_INCHI_KEY_TYPE,
    _canonicalize_loaded_duckdb,
    _create_duckdb_content_uuid_macro,
    _create_duckdb_evidence_tables,
)
from omnipath_build.resolver.identifier_types import (  # noqa: E402
    identifier_type_id,
)

LIPID_NAME_TID = identifier_type_id(LIPID_NAME_TYPE)
CHEBI_TID = identifier_type_id(CHEBI_TYPE)
INCHIKEY_TID = identifier_type_id(STANDARD_INCHI_KEY_TYPE)

# Two real, distinct raw names from the live utils lipid_name table that both
# canonicalize to the same Goslin rendering -- the format-variant convergence
# R9 acceptance scenario 2 describes, not a fixture invention.
VARIANT_A = 'PC 16:0/18:1'
VARIANT_B = 'PC(16:0/18:1)'
EXPECTED_CANONICAL = 'PC 16:0/18:1|sn_position|2|2'


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


def test_format_variants_of_one_lipid_converge_on_one_canonical_name():
    """Acceptance scenario 2: PC(16:0/18:1) and PC 16:0/18:1, written by two
    different (synthetic) resources, canonicalize identically."""
    con = _con()
    con.executemany(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        [
            ('sourceA', 'lipids', 1, 'm1', None, 'r', CHEMICAL_ENTITY_TYPE, None),
            ('sourceB', 'lipids', 1, 'm2', None, 'r', CHEMICAL_ENTITY_TYPE, None),
        ],
    )
    con.executemany(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        [
            ('sourceA', 'm1', 'i1', 'Name:OM:0202', VARIANT_A, VARIANT_A),
            ('sourceB', 'm2', 'i1', 'Name:OM:0202', VARIANT_B, VARIANT_B),
        ],
    )

    _canonicalize_loaded_duckdb(con)

    a = _resolution(con, 'm1')
    b = _resolution(con, 'm2')
    assert a is not None and b is not None, 'both mentions must resolve'
    assert a[0] == b[0] == LIPID_NAME_TID
    assert a[1] == b[1] == EXPECTED_CANONICAL
    assert a[2] == b[2] == 'resolved'
    assert a[3] == b[3] == 'lipid_name'


def test_a_lipid_with_a_structure_still_canonicalizes_on_the_name():
    """R9: "any InChIKey, ChEBI, ... identifier is attached as an ordinary
    identifier" -- a mention carrying BOTH a lipid-shaped name AND a real
    InChIKey still canonicalizes on the name, not the structure."""
    con = _con()
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('sourceC', 'lipids', 1, 'm3', None, 'r', CHEMICAL_ENTITY_TYPE, None),
    )
    con.executemany(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        [
            ('sourceC', 'm3', 'i1', 'Name:OM:0202', VARIANT_A, VARIANT_A),
            (
                'sourceC', 'm3', 'i2', 'Standard Inchi Key:MI:1101',
                'GHFVCXBIHHAKGP-BTJKTKAUSA-N', 'GHFVCXBIHHAKGP-BTJKTKAUSA-N',
            ),
        ],
    )

    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id, status, mechanism = _resolution(con, 'm3')
    assert status == 'resolved'
    assert mechanism == 'lipid_name'
    assert canonical_type == LIPID_NAME_TID
    assert canonical_id == EXPECTED_CANONICAL


def test_a_non_lipid_chemical_is_unaffected():
    """A chemical mention with no lipid-shaped name resolves exactly as
    before -- the lipid arm is a no-op for it."""
    con = _con()
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('sourceD', 'chemicals', 1, 'm4', None, 'r', CHEMICAL_ENTITY_TYPE, None),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        (
            'sourceD', 'm4', 'i1', 'Standard Inchi Key:MI:1101',
            'RYYVLZVUVIJVGH-UHFFFAOYSA-N', 'RYYVLZVUVIJVGH-UHFFFAOYSA-N',
        ),
    )

    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id, status, mechanism = _resolution(con, 'm4')
    assert status == 'resolved'
    assert mechanism == 'inchikey'
    assert canonical_type == INCHIKEY_TID
    assert canonical_id == 'RYYVLZVUVIJVGH-UHFFFAOYSA-N'
