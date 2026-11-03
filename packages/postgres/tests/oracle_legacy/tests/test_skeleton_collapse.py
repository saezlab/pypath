"""Candidates sharing a connectivity skeleton collapse to one representative
(spec 011 T048, WP2/User Story 2, acceptance scenario 1).

Fixture: the four ATP entities named in research.md R13 — one connectivity
skeleton (InChIKey block 1, ``ZKHQWZAMYRWXGA``) carrying four full InChIKeys
that differ only in stereo/protonation::

    ZKHQWZAMYRWXGA-KQYNXXCUSA-N   ATP           neutral
    ZKHQWZAMYRWXGA-KQYNXXCUSA-K   ATP(3-)       charged
    ZKHQWZAMYRWXGA-KQYNXXCUSA-J   ATP(4-)       charged
    ZKHQWZAMYRWXGA-UHFFFAOYSA-N   no stereo     neutral

A mention citing four different chemical identifier types that each resolve
(via the resolver) to one of these four InChIKeys currently sees
``candidate_count = 4`` and stays unresolved (an unarbitrated skeleton is
indistinguishable from a genuine cross-skeleton disagreement today). This test
asserts it instead resolves to one representative structure.
"""

from __future__ import annotations

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

from omnipath_build.duckdb_load import (  # noqa: E402
    CHEBI_TYPE,
    CHEMICAL_ENTITY_TYPE,
    HMDB_TYPE,
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

# The four ATP variants, one connectivity skeleton (ZKHQWZAMYRWXGA).
ATP_NEUTRAL_STEREO = 'ZKHQWZAMYRWXGA-KQYNXXCUSA-N'
ATP_CHARGED_3 = 'ZKHQWZAMYRWXGA-KQYNXXCUSA-K'
ATP_CHARGED_4 = 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J'
ATP_NEUTRAL_NO_STEREO = 'ZKHQWZAMYRWXGA-UHFFFAOYSA-N'
ATP_VARIANTS = (
    ATP_NEUTRAL_STEREO,
    ATP_CHARGED_3,
    ATP_CHARGED_4,
    ATP_NEUTRAL_NO_STEREO,
)


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


def test_four_candidates_sharing_a_skeleton_resolve_to_one_representative():
    con = _con()
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('reactome', 'reactions', 1, 'atp1', None, 'r', CHEMICAL_ENTITY_TYPE, None),
    )
    # Four different chemical identifier types on the same mention, each
    # resolving (via the resolver) to one of the four skeleton-sharing
    # InChIKeys -- exactly the "several candidates, one real molecule" shape
    # acceptance scenario 1 describes.
    con.executemany(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        [
            ('reactome', 'atp1', 'i1', CHEBI_TYPE, '15422', 'CHEBI:15422'),
            ('reactome', 'atp1', 'i2', HMDB_TYPE, 'HMDB0000538', 'HMDB0000538'),
            ('reactome', 'atp1', 'i3', KEGG_COMPOUND_TYPE, 'C00002', 'C00002'),
            ('reactome', 'atp1', 'i4', PUBCHEM_COMPOUND_TYPE, '5957', '5957'),
        ],
    )
    resolver_type = identifier_type_id(CHEBI_TYPE)
    con.executemany(
        'INSERT INTO resolver_lookup VALUES (?, ?, ?, ?, ?, ?)',
        [
            (
                CHEMICAL_ENTITY_TYPE,
                identifier_type_id(CHEBI_TYPE),
                'CHEBI:15422',
                None,
                INCHIKEY_TID,
                ATP_NEUTRAL_STEREO,
            ),
            (
                CHEMICAL_ENTITY_TYPE,
                identifier_type_id(HMDB_TYPE),
                'HMDB0000538',
                None,
                INCHIKEY_TID,
                ATP_CHARGED_3,
            ),
            (
                CHEMICAL_ENTITY_TYPE,
                identifier_type_id(KEGG_COMPOUND_TYPE),
                'C00002',
                None,
                INCHIKEY_TID,
                ATP_CHARGED_4,
            ),
            (
                CHEMICAL_ENTITY_TYPE,
                identifier_type_id(PUBCHEM_COMPOUND_TYPE),
                '5957',
                None,
                INCHIKEY_TID,
                ATP_NEUTRAL_NO_STEREO,
            ),
        ],
    )

    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id, status, mechanism = _resolution(con, 'atp1')
    assert status == 'resolved'
    assert canonical_type == INCHIKEY_TID
    assert canonical_id in ATP_VARIANTS
