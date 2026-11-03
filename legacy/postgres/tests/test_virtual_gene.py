"""A gene product known only by a UniProt accession is typed as a gene keyed by
that accession, carrying the accession the source asserted as a protein state.
This holds whether the accession resolves through the identifier-resolution
database (primary or secondary) or is one the resolver does not carry at all (a
well-formed accession still names a gene product). No stray protein entity is left
behind, and the entity key is deterministic (re-runnable build)."""

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
    identifier_type_id,
)


UNIPROT_TID = identifier_type_id(UNIPROT_TYPE)


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


def _insert_protein(con, ev_id: str, accession: str, taxon: str = '9606') -> None:
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('s', 'd', 1, ev_id, None, 'r', PROTEIN_ENTITY_TYPE, taxon),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        ('s', ev_id, f'i_{ev_id}', UNIPROT_TYPE, accession, None),
    )


def _fallback(con, key: str, primary: str, taxon: str = '9606') -> None:
    con.execute(
        'INSERT INTO protein_uniprot_fallback_lookup VALUES (?, ?, ?, ?, ?, ?)',
        (PROTEIN_ENTITY_TYPE, UNIPROT_TID, key, taxon, UNIPROT_TID, primary),
    )


def _build_secondary_ac_case(con):
    # A source references the secondary UniProt accession Q15086; the resolver
    # maps it to the primary P04637. No gene candidate exists for this mention.
    _insert_protein(con, 'secondary_only', 'Q15086')
    _fallback(con, 'Q15086', 'P04637')
    _canonicalize_loaded_duckdb(con)


def test_secondary_uniprot_mints_gene_keyed_by_primary():
    con = _con()
    _build_secondary_ac_case(con)

    # Typed as a gene, keyed by the primary accession.
    assert con.execute(
        """
        SELECT entity_type, canonical_identifier_type_id, canonical_identifier,
               status, resolution_mechanism
        FROM entity_resolution WHERE entity_evidence_id = 'secondary_only'
        """
    ).fetchone() == (
        GENE_ENTITY_TYPE, UNIPROT_TID, 'P04637', 'resolved', 'unknown_gene',
    )


def test_no_stray_unresolved_protein_entity():
    con = _con()
    _build_secondary_ac_case(con)

    # The base graph is uniformly gene-typed: the mint produces a Gene entity and
    # NO Protein entity for the same mention.
    assert con.execute(
        """
        SELECT count(*) FROM canonical_entity
        WHERE entity_type = ? AND canonical_identifier = 'P04637'
        """,
        [GENE_ENTITY_TYPE],
    ).fetchone()[0] == 1
    assert con.execute(
        'SELECT count(*) FROM canonical_entity WHERE entity_type = ?',
        [PROTEIN_ENTITY_TYPE],
    ).fetchone()[0] == 0


def test_asserted_ac_captured_as_protein_state():
    con = _con()
    _build_secondary_ac_case(con)

    gene_entity_id = con.execute(
        'SELECT entity_id FROM canonical_entity WHERE canonical_identifier = ?',
        ['P04637'],
    ).fetchone()[0]
    # The exact accession the source asserted (the raw secondary Q15086) is
    # recorded as a protein state hanging off the gene.
    assert con.execute(
        """
        SELECT gene_entity_id, uniprot_ac, isoform
        FROM evidence_state_link WHERE entity_evidence_id = 'secondary_only'
        """
    ).fetchone() == (gene_entity_id, 'Q15086', None)


def test_representative_uniprot_arm_for_virtual_gene():
    con = _con()
    _build_secondary_ac_case(con)

    gene_entity_id = con.execute(
        'SELECT entity_id FROM canonical_entity WHERE canonical_identifier = ?',
        ['P04637'],
    ).fetchone()[0]
    # A UniProt-keyed gene is its own representative UniProt.
    assert con.execute(
        """
        SELECT representative_uniprot, uniprot_all
        FROM gene_protein_representative WHERE entity_id = ?
        """,
        [gene_entity_id],
    ).fetchone() == ('P04637', ['P04637'])


def test_resolver_unknown_uniprot_still_mints_gene():
    con = _con()
    # A well-formed UniProt accession the resolver carries no gene or fallback for
    # (e.g. a TrEMBL accession) still names a gene product: typed Gene keyed by the
    # accession, using the evidence taxon. No resolver/fallback rows are loaded.
    _insert_protein(con, 'trembl', 'A0A8M9QHZ6', taxon='7955')
    _canonicalize_loaded_duckdb(con)

    assert con.execute(
        """
        SELECT entity_type, taxonomy_id, canonical_identifier_type_id,
               canonical_identifier, status, resolution_mechanism
        FROM entity_resolution WHERE entity_evidence_id = 'trembl'
        """
    ).fetchone() == (
        GENE_ENTITY_TYPE, '7955', UNIPROT_TID, 'A0A8M9QHZ6',
        'resolved', 'unknown_gene',
    )
    assert con.execute(
        'SELECT count(*) FROM canonical_entity WHERE entity_type = ?',
        [PROTEIN_ENTITY_TYPE],
    ).fetchone()[0] == 0


def test_protein_uniprot_without_taxon_is_not_minted():
    con = _con()
    # Virtual genes must carry a real taxon. A UniProt-only protein with no
    # organism — and none recoverable from the accession (no accession->organism
    # lookup here) — cannot be safely typed as a gene, so it stays unresolved.
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('s', 'd', 1, 'notax', None, 'r', PROTEIN_ENTITY_TYPE, None),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        ('s', 'notax', 'i_notax', UNIPROT_TYPE, 'A0A8M9QHZ6', None),
    )
    _canonicalize_loaded_duckdb(con)

    assert con.execute(
        'SELECT entity_type, status FROM entity_resolution '
        "WHERE entity_evidence_id = 'notax'"
    ).fetchone() == (PROTEIN_ENTITY_TYPE, 'unresolved')


def test_no_taxon_protein_gains_organism_from_accession_then_mints():
    con = _con()
    # The same organism-less mention, but now its accession's organism is known
    # (the accession->organism lookup carries it): the organism is recovered onto
    # the mention, so it is typed as a gene keyed by the accession, carrying the
    # recovered organism — no bare protein node is left behind.
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('s', 'd', 1, 'recovered', None, 'r', PROTEIN_ENTITY_TYPE, None),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        ('s', 'recovered', 'i_recovered', UNIPROT_TYPE, 'A0A8M9QHZ6', None),
    )
    con.execute(
        'CREATE TABLE uniprot_taxon_lookup (uniprot VARCHAR, taxonomy_id VARCHAR)'
    )
    con.execute(
        "INSERT INTO uniprot_taxon_lookup VALUES ('A0A8M9QHZ6', '7955')"
    )
    _canonicalize_loaded_duckdb(con)

    assert con.execute(
        """
        SELECT entity_type, taxonomy_id, canonical_identifier, status,
               resolution_mechanism
        FROM entity_resolution WHERE entity_evidence_id = 'recovered'
        """
    ).fetchone() == (
        GENE_ENTITY_TYPE, '7955', 'A0A8M9QHZ6', 'resolved', 'unknown_gene',
    )
    assert con.execute(
        'SELECT count(*) FROM canonical_entity WHERE entity_type = ?',
        [PROTEIN_ENTITY_TYPE],
    ).fetchone()[0] == 0


def test_isoform_accession_recovers_organism_from_base_accession():
    con = _con()
    # A mention that cites an isoform accession (``-<n>`` suffix) is matched to the
    # lookup by its base accession, so the organism is still recovered.
    con.execute(
        'INSERT INTO entity_evidence_raw VALUES (?, ?, ?, ?, ?, ?, ?, ?)',
        ('s', 'd', 1, 'iso', None, 'r', PROTEIN_ENTITY_TYPE, None),
    )
    con.execute(
        'INSERT INTO entity_identifier_raw VALUES (?, ?, ?, ?, ?, ?)',
        ('s', 'iso', 'i_iso', UNIPROT_TYPE, 'A0A8M9QHZ6-2', None),
    )
    con.execute(
        'CREATE TABLE uniprot_taxon_lookup (uniprot VARCHAR, taxonomy_id VARCHAR)'
    )
    con.execute(
        "INSERT INTO uniprot_taxon_lookup VALUES ('A0A8M9QHZ6', '7955')"
    )
    _canonicalize_loaded_duckdb(con)

    assert con.execute(
        'SELECT taxonomy_id, status FROM entity_resolution '
        "WHERE entity_evidence_id = 'iso'"
    ).fetchone() == ('7955', 'resolved')


def test_malformed_uniprot_value_is_not_minted():
    con = _con()
    _insert_protein(con, 'bad', 'NOTANACCESSION', taxon='9606')
    _canonicalize_loaded_duckdb(con)

    assert con.execute(
        'SELECT entity_type, status FROM entity_resolution '
        "WHERE entity_evidence_id = 'bad'"
    ).fetchone() == (PROTEIN_ENTITY_TYPE, 'unresolved')


def test_virtual_gene_key_is_idempotent():
    # Deterministic keys: the same inputs, built independently, produce the same
    # entity_id (re-runnable build).
    con_a, con_b = _con(), _con()
    _build_secondary_ac_case(con_a)
    _build_secondary_ac_case(con_b)
    key_a = con_a.execute(
        'SELECT entity_id FROM canonical_entity WHERE canonical_identifier = ?',
        ['P04637'],
    ).fetchone()[0]
    key_b = con_b.execute(
        'SELECT entity_id FROM canonical_entity WHERE canonical_identifier = ?',
        ['P04637'],
    ).fetchone()[0]
    assert key_a == key_b
