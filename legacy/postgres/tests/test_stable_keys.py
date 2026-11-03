"""Stable, deterministic entity keys: the base graph carries no hash-keyed
entity (an unresolvable entity is keyed by its best raw identifier), no
placeholder protein entity when a UniProt accession is available, and every
entity has a non-empty key."""

from __future__ import annotations

import pytest

duckdb = pytest.importorskip('duckdb')
pytest.importorskip('pkg_infra')

from omnipath_build.duckdb_load import (  # noqa: E402
    ENSEMBL_TYPE,
    GENE_ENTITY_TYPE,
    GENE_NAME_PRIMARY_TYPE,
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


def _insert(con, ev_id, entity_type, ids: list[tuple[str, str]], taxon='9606'):
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('s', 'd', 1, ev_id, None, 'r', entity_type, taxon),
    )
    for i, (id_type, value) in enumerate(ids):
        con.execute(
            'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
            ('s', ev_id, f'i_{ev_id}_{i}', id_type, value, None),
        )


def test_unresolved_key_is_best_identifier_not_md5():
    con = _con()
    # No resolver/fallback match: this record survives every path to the ultimate
    # unresolved branch. It carries an Ensembl and a gene symbol; the best-id
    # priority (UniProt > Ensembl > Entrez > symbol) picks the Ensembl id.
    _insert(
        con,
        'ultimate',
        PROTEIN_ENTITY_TYPE,
        [(GENE_NAME_PRIMARY_TYPE, 'FOO'), (ENSEMBL_TYPE, 'ENSP00000269305')],
    )
    _canonicalize_loaded_duckdb(con)

    canonical_type, canonical_id, status = con.execute(
        """
        SELECT canonical_identifier_type_id, canonical_identifier, status
        FROM entity_resolution WHERE entity_evidence_id = 'ultimate'
        """
    ).fetchone()
    assert status == 'unresolved'
    assert canonical_type == UNRESOLVED_TID
    # The key is the real Ensembl id, NOT an md5 hash of the identifier bag.
    assert canonical_id == 'ENSP00000269305'


def test_no_md5_hashed_entity_in_base_graph():
    con = _con()
    _insert(
        con,
        'ultimate',
        PROTEIN_ENTITY_TYPE,
        [(ENSEMBL_TYPE, 'ENSP00000269305')],
    )
    _canonicalize_loaded_duckdb(con)

    # No entity is keyed by a bare 32-char md5 hex digest.
    assert con.execute(
        r"""
        SELECT count(*) FROM canonical_entity
        WHERE regexp_matches(canonical_identifier, '^[0-9a-f]{32}$')
        """
    ).fetchone()[0] == 0


def test_every_entity_has_a_nonempty_key():
    con = _con()
    _insert(con, 'e1', PROTEIN_ENTITY_TYPE, [(ENSEMBL_TYPE, 'ENSP1')])
    _insert(con, 'e2', PROTEIN_ENTITY_TYPE, [(UNIPROT_TYPE, 'P00533')])
    con.execute(
        'INSERT INTO protein_uniprot_fallback_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (PROTEIN_ENTITY_TYPE, UNIPROT_TID, 'P00533', '9606', UNIPROT_TID, 'P00533'),
    )
    _canonicalize_loaded_duckdb(con)

    assert con.execute(
        """
        SELECT count(*) FROM canonical_entity
        WHERE canonical_identifier IS NULL OR canonical_identifier = ''
        """
    ).fetchone()[0] == 0


def test_protein_with_uniprot_leaves_no_protein_node():
    con = _con()
    _insert(con, 'known', PROTEIN_ENTITY_TYPE, [(UNIPROT_TYPE, 'P00533')])
    con.execute(
        'INSERT INTO protein_uniprot_fallback_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (PROTEIN_ENTITY_TYPE, UNIPROT_TID, 'P00533', '9606', UNIPROT_TID, 'P00533'),
    )
    _canonicalize_loaded_duckdb(con)

    assert con.execute(
        'SELECT count(*) FROM canonical_entity WHERE entity_type = ?',
        [PROTEIN_ENTITY_TYPE],
    ).fetchone()[0] == 0
    assert con.execute(
        'SELECT count(*) FROM canonical_entity WHERE entity_type = ?',
        [GENE_ENTITY_TYPE],
    ).fetchone()[0] == 1
