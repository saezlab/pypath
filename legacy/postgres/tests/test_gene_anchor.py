"""Gene-anchored resolution invariants.

Asserts the post-build state of a gene-anchored OmniPath database:
- every gene-family entity is keyed by an **NCBI Gene (Entrez)** anchor;
- there is exactly **one gene entity per (gene, organism)** — no fragmentation by
  the id type a source happened to use (symbol / Entrez / Ensembl / UniProt /
  isoform all collapse to the same gene entity);
- the EGFR benchmark (human NCBI Gene 1956) is a single gene entity reachable by
  several id types;
- non-human organisms resolve too (orthologs are gene-anchored, not dropped).

Run against a built instance, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run pytest tests/test_gene_anchor.py -v

Skipped when DATABASE_URL is not set (no DB to test against).
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

GENE_ENTITY_TYPE = 'Gene:MI:0250'
ENTREZ_TYPE = 'Entrez:MI:0477'
EGFR_ENTREZ = '1956'  # human EGFR NCBI Gene id

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; gene-anchor test needs a built database',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _scalar(conn, query: str, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params or [])
        row = cur.fetchone()
        return row[0] if row else None


@pytest.fixture(scope='module')
def gene_type_id(conn):
    type_id = _scalar(
        conn,
        f'SELECT entity_type_id FROM {SCHEMA}.vocab_entity_type WHERE name = %s',
        [GENE_ENTITY_TYPE],
    )
    if type_id is None:
        pytest.skip('no Gene entity type — build predates the gene anchor')
    return type_id


def test_gene_entities_exist(conn, gene_type_id):
    assert (
        _scalar(
            conn,
            f'SELECT count(*) FROM {SCHEMA}.entity WHERE entity_type_id = %s',
            [gene_type_id],
        )
        > 0
    )


def test_every_gene_anchored_on_entrez(conn, gene_type_id):
    """No gene entity may be keyed by anything other than the Entrez anchor."""
    entrez_type_id = _scalar(
        conn,
        f'SELECT identifier_type_id FROM {SCHEMA}.vocab_identifier_type WHERE name = %s',
        [ENTREZ_TYPE],
    )
    assert entrez_type_id is not None, 'Entrez identifier type missing'
    non_entrez = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE entity_type_id = %s
          AND canonical_identifier_type_id IS DISTINCT FROM %s
        """,
        [gene_type_id, entrez_type_id],
    )
    assert non_entrez == 0, f'{non_entrez} gene entities not Entrez-anchored'


def test_one_gene_entity_per_gene_per_organism(conn, gene_type_id):
    """(canonical_identifier, taxonomy_id) is unique among gene entities."""
    dup = _scalar(
        conn,
        f"""
        SELECT count(*) FROM (
          SELECT canonical_identifier, taxonomy_id
          FROM {SCHEMA}.entity
          WHERE entity_type_id = %s
          GROUP BY canonical_identifier, taxonomy_id
          HAVING count(*) > 1
        ) d
        """,
        [gene_type_id],
    )
    assert dup == 0, f'{dup} (gene, organism) pairs map to >1 entity'


def test_egfr_is_one_gene_entity(conn, gene_type_id):
    """Human EGFR (NCBI Gene 1956) resolves to exactly one gene entity."""
    n = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE entity_type_id = %s
          AND canonical_identifier = %s
          AND taxonomy_id = 9606
        """,
        [gene_type_id, EGFR_ENTREZ],
    )
    assert n == 1, f'EGFR resolved to {n} gene entities (expected 1)'


def test_egfr_reachable_by_multiple_id_types(conn, gene_type_id):
    """The single EGFR gene entity carries several id types (cross-id collapse)."""
    egfr_id = _scalar(
        conn,
        f"""
        SELECT entity_id FROM {SCHEMA}.entity
        WHERE entity_type_id = %s AND canonical_identifier = %s AND taxonomy_id = 9606
        """,
        [gene_type_id, EGFR_ENTREZ],
    )
    if egfr_id is None:
        pytest.skip('EGFR gene entity absent in this (capped) build')
    id_types = _scalar(
        conn,
        f"""
        SELECT count(DISTINCT ie.identifier_type_id)
        FROM {SCHEMA}.entity_identifier ei
        JOIN {SCHEMA}.identifier_evidence ie ON ie.identifier_id = ei.identifier_id
        WHERE ei.entity_id = %s
        """,
        [egfr_id],
    )
    assert id_types >= 2, (
        f'EGFR entity has only {id_types} id type(s) — expected the gene to '
        f'collapse several (symbol / Entrez / UniProt / …)'
    )


def test_non_human_organisms_resolve(conn, gene_type_id):
    """Orthologs are gene-anchored: gene entities exist for non-human organisms."""
    non_human = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE entity_type_id = %s
          AND taxonomy_id IS NOT NULL
          AND taxonomy_id <> 9606
        """,
        [gene_type_id],
    )
    assert non_human > 0, 'no non-human gene entities — orthologs not resolving'
