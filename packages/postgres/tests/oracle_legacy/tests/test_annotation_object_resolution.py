"""An annotation-relation chemical object promotes through the resolver
instead of minting a permanent duplicate (spec 011, the cycle's "final
rebuild" bucket).

Previously ``annotation_object_entity`` self-typed the raw
``(object_id_type, object_id)`` verbatim -- a reaction annotated with a bare
ChEBI id (no ``entity_evidence`` mention of its own) minted its own entity
keyed by that ChEBI id, permanently separate from whatever InChIKey entity
the same real molecule resolves to via the normal mention pipeline.
Confirmed reproducibly with ``CHEBI:107644`` (piperonylic acid) across two
full rebuilds before this fix, identical orphaned UUID both times.

Fix: promote through ``annotation_object_resolution``, built from
``resolver_lookup`` the same way ``ontology_relation_endpoint_key`` already
promotes ontology relation endpoints (an established pattern for "a bare
external identifier that is not a full mention", not the heavier
per-mention candidate-arbitration machinery WP2 built for mentions with
possibly several identifiers to arbitrate between).
"""

from __future__ import annotations

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

from omnipath_build.duckdb_load import (  # noqa: E402
    CHEBI_TYPE,
    CHEMICAL_ENTITY_TYPE,
    STANDARD_INCHI_KEY_TYPE,
    _build_annotation_object_resolution,
    _canonicalize_loaded_duckdb,
    _create_duckdb_content_uuid_macro,
    _create_duckdb_evidence_tables,
)
from omnipath_build.resolver.identifier_types import (  # noqa: E402
    identifier_type_id,
)

INCHIKEY_TID = identifier_type_id(STANDARD_INCHI_KEY_TYPE)
CHEBI_TID = identifier_type_id(CHEBI_TYPE)
FIXTURE_INCHIKEY = 'RXWFNWIGFZWKRT-UHFFFAOYSA-N'  # synthetic, not the real structure -- this test only checks identity round-trips, never chemical accuracy

# Not read by the resolution logic under test, only needed to satisfy
# identifier_authority's foreign key (mirrors test_authority_arbitration.py).
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


def _annotation_resolution(con, object_id: str):
    return con.execute(
        """
        SELECT canonical_identifier_type_id, canonical_identifier
        FROM annotation_object_resolution
        WHERE entity_type = ? AND raw_identifier = ?
        """,
        [CHEMICAL_ENTITY_TYPE, object_id],
    ).fetchone()


def test_annotation_object_promotes_through_the_resolver_when_a_match_exists():
    """CHEBI:107644's real reproduction case, at unit scale: an annotation
    object citing a ChEBI id the resolver knows lands on that id's InChIKey,
    not the raw ChEBI id."""
    con = _con()
    con.execute(
        "INSERT INTO annotation_relation_evidence_raw VALUES "
        "('r1', 'reactome', 'reactions', 1, 'm1', 'has_participant', ?, ?, ?, 'reaction_metabolite')",
        [CHEMICAL_ENTITY_TYPE, CHEBI_TYPE, 'CHEBI:107644'],
    )
    con.execute(
        'INSERT INTO resolver_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (
            CHEMICAL_ENTITY_TYPE,
            CHEBI_TID,
            'CHEBI:107644',
            None,
            INCHIKEY_TID,
            FIXTURE_INCHIKEY,
        ),
    )

    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id = _annotation_resolution(con, 'CHEBI:107644')
    assert canonical_type == INCHIKEY_TID
    assert canonical_id == FIXTURE_INCHIKEY


def test_annotation_object_falls_through_to_the_raw_id_with_no_resolver_match():
    """No resolver coverage for this id at all -- falls through to the same
    raw self-typed behaviour as before the fix (never a crash, never NULL)."""
    con = _con()
    con.execute(
        "INSERT INTO annotation_relation_evidence_raw VALUES "
        "('r1', 'reactome', 'reactions', 1, 'm1', 'has_participant', ?, ?, ?, 'reaction_metabolite')",
        [CHEMICAL_ENTITY_TYPE, CHEBI_TYPE, 'CHEBI:99999999'],
    )

    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id = _annotation_resolution(con, 'CHEBI:99999999')
    assert canonical_type == CHEBI_TID
    assert canonical_id == 'CHEBI:99999999'


def test_a_mention_and_an_annotation_object_citing_the_same_id_land_on_one_entity():
    """The actual regression this fix closes: a real chemical mention and an
    annotation-relation object citing the SAME ChEBI id must resolve to the
    identical (canonical_identifier_type_id, canonical_identifier) pair --
    otherwise they mint two separate entities for one real molecule."""
    con = _con()
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('chebi', 'compounds', 1, 'm1', None, 'r', CHEMICAL_ENTITY_TYPE, None),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        ('chebi', 'm1', 'i1', CHEBI_TYPE, '107644', 'CHEBI:107644'),
    )
    con.execute(
        "INSERT INTO annotation_relation_evidence_raw VALUES "
        "('r1', 'reactome', 'reactions', 1, 'm2', 'has_participant', ?, ?, ?, 'reaction_metabolite')",
        [CHEMICAL_ENTITY_TYPE, CHEBI_TYPE, 'CHEBI:107644'],
    )
    con.execute(
        'INSERT INTO resolver_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (
            CHEMICAL_ENTITY_TYPE,
            CHEBI_TID,
            'CHEBI:107644',
            None,
            INCHIKEY_TID,
            FIXTURE_INCHIKEY,
        ),
    )

    _canonicalize_loaded_duckdb(con)

    mention_type, mention_id, _, _ = con.execute(
        """
        SELECT canonical_identifier_type_id, canonical_identifier, status,
               resolution_mechanism
        FROM entity_resolution
        WHERE entity_evidence_id = 'm1'
        """
    ).fetchone()
    annotation_type, annotation_id = _annotation_resolution(con, 'CHEBI:107644')

    assert mention_type == annotation_type == INCHIKEY_TID
    assert mention_id == annotation_id == FIXTURE_INCHIKEY


def test_annotation_object_resolution_is_a_pure_function_of_resolver_lookup():
    """Calling the builder directly (bypassing the heavier canonicalize pass)
    gives the same answer -- pins the table's own contract independent of
    where it's invoked from."""
    con = _con()
    con.execute(
        "INSERT INTO annotation_relation_evidence_raw VALUES "
        "('r1', 'reactome', 'reactions', 1, 'm1', 'has_participant', ?, ?, ?, 'reaction_metabolite')",
        [CHEMICAL_ENTITY_TYPE, CHEBI_TYPE, 'CHEBI:107644'],
    )
    con.execute(
        'INSERT INTO resolver_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (
            CHEMICAL_ENTITY_TYPE,
            CHEBI_TID,
            'CHEBI:107644',
            None,
            INCHIKEY_TID,
            FIXTURE_INCHIKEY,
        ),
    )
    from omnipath_build.duckdb_load import _create_duckdb_identifier_type_all_view

    _create_duckdb_identifier_type_all_view(con)

    _build_annotation_object_resolution(con)

    canonical_type, canonical_id = _annotation_resolution(con, 'CHEBI:107644')
    assert canonical_type == INCHIKEY_TID
    assert canonical_id == FIXTURE_INCHIKEY
