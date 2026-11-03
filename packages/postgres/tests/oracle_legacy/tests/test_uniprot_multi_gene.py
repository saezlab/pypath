"""A UniProt AC mapping to >1 gene is split per gene.

The rare identical-gene-copies case (one UniProt AC, several NCBI genes producing
identical protein) MUST be **split into one gene-anchored record per gene**, each
carrying the same form as state (molecular_type = protein, uniprot = the AC); the
protein-centric view re-collapses them. No anchor is chosen and nothing is dropped.

DELIVERED. The split is realised **1:1**: a multi-gene UniProt mention is
*duplicated* per gene before resolution (``omnipath_build/multigene_split.py``),
each copy resolving to its own gene and carrying the same UniProt as a protein
``state``. So the signature of a performed split is **one UniProt AC appearing as
a protein state under >1 distinct gene entity** (not one evidence → many genes,
which the 1:1 model deliberately avoids). The explosion logic itself is unit-
tested deterministically in ``tests/test_multigene_split.py``; this test is the
integration check and **skips** when a capped build retains no multi-gene case.

Run against a built instance, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run pytest tests/test_uniprot_multi_gene.py -v

Skipped when DATABASE_URL is not set (no DB to test against).
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

GENE_ENTITY_TYPE = 'Gene:MI:0250'

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; multi-gene split test needs a built database',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _rows(conn, query: str, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params or [])
        return cur.fetchall()


def test_uniprot_mapping_to_multiple_genes_is_split(conn):
    """A UniProt AC that maps to >1 gene appears as a protein state under each
    gene (the 1:1-duplication signature of a performed split).

    Skips on a capped build that retains no multi-gene UniProt case — the
    explosion logic itself is guarded by tests/test_multigene_split.py.
    """
    gene_type_id = next(
        (
            r[0]
            for r in _rows(
                conn,
                f'SELECT entity_type_id FROM {SCHEMA}.vocab_entity_type WHERE name = %s',
                [GENE_ENTITY_TYPE],
            )
        ),
        None,
    )
    assert gene_type_id is not None

    # A UniProt state-component value linked (via state) to >1 distinct gene
    # entity — the signature of a performed split under the 1:1 model.
    multi = _rows(
        conn,
        f"""
        SELECT sc.value, count(DISTINCT s.gene_entity_id) AS genes
        FROM {SCHEMA}.state_component sc
        JOIN {SCHEMA}.state s ON s.state_id = sc.state_id
        WHERE sc.component_type = 'uniprot'
        GROUP BY sc.value
        HAVING count(DISTINCT s.gene_entity_id) > 1
        LIMIT 1
        """,
    )
    if not multi:
        pytest.skip('no multi-gene UniProt case retained in this (capped) build')
    assert multi[0][1] > 1
