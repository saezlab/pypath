"""Regression test for spec 011 WP1: the resolver-candidate join must key on
the normalized identifier, not the raw one.

T040/T041 taught ``entity_identifier_raw`` a chemical ``identifier_normalized``
column (e.g. Reactome cites ChEBI as a bare digit, ``30616``, normalized to the
CURIE form ``CHEBI:30616`` that ``resolver_chemical``/``resolver_lookup`` key
on). ``evidence_identifier_key`` (feeding the keyed fetch from utils Postgres)
already reads the coalesced normalized value. But ``entity_resolution_base``'s
``resolver_candidate`` CTE joined ``resolver_lookup`` back against the RAW
``ei.identifier`` -- so a mention citing the bare-digit form never matched its
own (correctly normalized) fetched resolver row, produced zero candidates, and
silently fell through to the chemical self-type fallback ('resolved', but
never promoted to the structure). This is the confirmed root cause of
Reactome's near-zero structure reach (SC-002/003) and the Reactome-HMDB
shared-metabolite regression (SC-004).
"""

from __future__ import annotations

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

from omnipath_build.duckdb_load import (  # noqa: E402
    CHEBI_TYPE,
    CHEMICAL_ENTITY_TYPE,
    STANDARD_INCHI_KEY_TYPE,
    _canonicalize_loaded_duckdb,
    _create_duckdb_content_uuid_macro,
    _create_duckdb_evidence_tables,
)
from omnipath_build.resolver.identifier_types import (  # noqa: E402
    identifier_type_id,
)


CHEBI_TID = identifier_type_id(CHEBI_TYPE)
INCHIKEY_TID = identifier_type_id(STANDARD_INCHI_KEY_TYPE)
INCHIKEY = 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J'


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
        SELECT canonical_identifier_type_id, canonical_identifier, status,
               resolution_mechanism
        FROM entity_resolution
        WHERE entity_evidence_id = ?
        """,
        [ev_id],
    ).fetchone()


def test_bare_digit_reactome_style_chebi_mention_resolves_via_normalized_key():
    con = _con()
    # Reactome's own emission: the raw evidence value is the bare ChEBI
    # digit ('30616'), while identifier_normalized carries the CURIE form
    # that resolver_lookup is actually keyed by ('CHEBI:30616') -- exactly
    # what omnipath_build.resolver.chemical_normalization._normalize_chebi
    # and resolver_chemical's own source_id convention produce.
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('reactome', 'reactions', 1, 'ev1', None, 'r', CHEMICAL_ENTITY_TYPE, None),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        ('reactome', 'ev1', 'i_ev1', CHEBI_TYPE, '30616', 'CHEBI:30616'),
    )
    con.execute(
        'INSERT INTO resolver_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (
            CHEMICAL_ENTITY_TYPE,
            CHEBI_TID,
            'CHEBI:30616',
            None,
            INCHIKEY_TID,
            INCHIKEY,
        ),
    )

    _canonicalize_loaded_duckdb(con)

    assert _resolution(con, 'ev1') == (
        INCHIKEY_TID,
        INCHIKEY,
        'resolved',
        'resolver',
    )
