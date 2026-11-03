"""Unit test for the multi-gene protein split explosion.

Builds a minimal synthetic set of the DuckDB raw tables — a protein mention
whose UniProt maps to two genes, participating in one relation (with a relation
annotation) — runs ``explode_multi_gene_protein_mentions``, and asserts the
1:1 duplication: the mention, its identifiers, the relation and its annotation
all fan out one-per-gene with fresh ids, and each copy gets a direct gene
resolution. A control mention (single-gene) and a control relation are left
untouched.

Needs ``duckdb`` only (no Postgres). Skipped when duckdb is absent.
"""

from __future__ import annotations

import hashlib

import pytest

duckdb = pytest.importorskip('duckdb')


def cu(payload: str) -> str:
    """Mirror of the build's content_uuid (md5 → 8-4-4-4-12 text)."""
    h = hashlib.md5(payload.encode()).hexdigest()
    return f'{h[0:8]}-{h[8:12]}-{h[12:16]}-{h[16:20]}-{h[20:32]}'

from omnipath_build.multigene_split import (  # noqa: E402
    ENTREZ_TYPE,
    UNIPROT_TYPE,
    explode_multi_gene_protein_mentions,
)

PROTEIN = 'Protein:MI:0326'
GENE = 'Gene:MI:0250'
TAX = '9606'
UNIPROT_TID = 10
ENTREZ_TID = 20


def _con():
    con = duckdb.connect(':memory:')
    con.execute(
        'CREATE TABLE identifier_type_all (identifier_type_id BIGINT, name VARCHAR)'
    )
    con.executemany(
        'INSERT INTO identifier_type_all VALUES (?, ?)',
        [(UNIPROT_TID, UNIPROT_TYPE), (ENTREZ_TID, ENTREZ_TYPE)],
    )
    con.execute(
        """
        CREATE TABLE needed_resolver_lookup (
          key_identifier_type_id BIGINT, key_value VARCHAR,
          evidence_entity_type VARCHAR, taxonomy_id VARCHAR,
          canonical_identifier_type_id BIGINT, canonical_identifier VARCHAR,
          taxonomy_optional_match BOOLEAN
        )
        """
    )
    con.executemany(
        'INSERT INTO needed_resolver_lookup VALUES (?, ?, ?, ?, ?, ?, ?)',
        [
            # P1 (multi-gene) → genes 100 and 200; P2 (single) → gene 300
            (UNIPROT_TID, 'P1', PROTEIN, TAX, ENTREZ_TID, '100', False),
            (UNIPROT_TID, 'P1', PROTEIN, TAX, ENTREZ_TID, '200', False),
            (UNIPROT_TID, 'P2', PROTEIN, TAX, ENTREZ_TID, '300', False),
        ],
    )
    con.execute(
        """
        CREATE TABLE entity_evidence_raw (
          source VARCHAR, dataset VARCHAR, row_id BIGINT,
          entity_evidence_id VARCHAR, parent_entity_evidence_id VARCHAR,
          entity_role VARCHAR, entity_type VARCHAR, taxonomy_id VARCHAR
        )
        """
    )
    con.executemany(
        'INSERT INTO entity_evidence_raw VALUES (?,?,?,?,?,?,?,?)',
        [
            ('s', 'd', 1, 'M', None, 'r', PROTEIN, TAX),   # multi-gene
            ('s', 'd', 2, 'P2ev', None, 'r', PROTEIN, TAX),  # single-gene
            ('s', 'd', 3, 'X', None, 'r', PROTEIN, TAX),   # partner (no uniprot)
        ],
    )
    con.execute(
        """
        CREATE TABLE entity_identifier_raw (
          source VARCHAR, entity_evidence_id VARCHAR, identifier_id VARCHAR,
          identifier_type VARCHAR, identifier VARCHAR, identifier_normalized VARCHAR
        )
        """
    )
    con.executemany(
        'INSERT INTO entity_identifier_raw VALUES (?,?,?,?,?,?)',
        [
            ('s', 'M', 'i1', UNIPROT_TYPE, 'P1', None),
            ('s', 'P2ev', 'i2', UNIPROT_TYPE, 'P2', None),
        ],
    )
    con.execute(
        """
        CREATE TABLE entity_annotation_raw (
          source VARCHAR, evidence_id VARCHAR, annotation_key VARCHAR,
          term VARCHAR, value VARCHAR, unit VARCHAR
        )
        """
    )
    con.executemany(
        'INSERT INTO entity_annotation_raw VALUES (?,?,?,?,?,?)',
        [('s', 'M', 'k', 't', 'v', None)],
    )
    con.execute(
        """
        CREATE TABLE relation_evidence_raw (
          source VARCHAR, dataset VARCHAR, row_id BIGINT,
          relation_evidence_id VARCHAR, subject_entity_evidence_id VARCHAR,
          predicate VARCHAR, object_entity_evidence_id VARCHAR,
          relation_category VARCHAR
        )
        """
    )
    con.executemany(
        'INSERT INTO relation_evidence_raw VALUES (?,?,?,?,?,?,?,?)',
        [
            ('s', 'd', 1, 'r1', 'X', 'interacts_with', 'M', 'c'),   # X–M(multi)
            ('s', 'd', 2, 'r2', 'X', 'interacts_with', 'P2ev', 'c'),  # control
        ],
    )
    con.execute(
        """
        CREATE TABLE relation_annotation_raw (
          source VARCHAR, evidence_id VARCHAR, annotation_key VARCHAR,
          annotation_scope VARCHAR, term VARCHAR, value VARCHAR, unit VARCHAR
        )
        """
    )
    con.executemany(
        'INSERT INTO relation_annotation_raw VALUES (?,?,?,?,?,?,?)',
        [('s', 'r1', 'k', 'sc', 't', 'v', None)],
    )
    con.execute(
        """
        CREATE TABLE annotation_relation_evidence_raw (
          relation_evidence_id VARCHAR, source VARCHAR, dataset VARCHAR,
          row_id BIGINT, subject_entity_evidence_id VARCHAR, predicate VARCHAR,
          object_entity_type VARCHAR, object_id_type VARCHAR, object_id VARCHAR,
          relation_category VARCHAR
        )
        """
    )
    con.executemany(
        'INSERT INTO annotation_relation_evidence_raw VALUES (?,?,?,?,?,?,?,?,?,?)',
        [('a1', 's', 'd', 1, 'M', 'located_in', 'Cv Term:OM:0012', 'go', 'GO:1', 'c')],
    )
    con.execute(
        """
        CREATE TABLE ontology_relation_raw (
          source VARCHAR, dataset VARCHAR, subject_entity_evidence_id VARCHAR,
          ontology_id VARCHAR, subject_entity_type VARCHAR,
          subject_identifier_type VARCHAR, subject_identifier VARCHAR,
          predicate VARCHAR, object_entity_type VARCHAR,
          object_identifier_type VARCHAR, object_identifier VARCHAR
        )
        """
    )
    return con


def _vals(con, sql):
    return con.execute(sql).fetchall()


def test_multi_gene_mention_splits_one_record_per_gene():
    con = _con()
    copies = explode_multi_gene_protein_mentions(con)
    assert copies == 2  # M → 2 genes
    m100, m200 = cu('M#mg=100'), cu('M#mg=200')

    # The mention M is replaced by one copy per gene (UUID ids); controls kept.
    ids = {r[0] for r in _vals(con, 'SELECT entity_evidence_id FROM entity_evidence_raw')}
    assert ids == {m100, m200, 'P2ev', 'X'}

    # Identifiers + entity annotation fan out, keeping the UniProt on each copy.
    uni = {tuple(r) for r in _vals(
        con,
        "SELECT entity_evidence_id, identifier FROM entity_identifier_raw "
        "WHERE identifier = 'P1'",
    )}
    assert uni == {(m100, 'P1'), (m200, 'P1')}
    ann = {r[0] for r in _vals(con, 'SELECT evidence_id FROM entity_annotation_raw')}
    assert ann == {m100, m200}


def test_relation_fans_out_with_fresh_ids():
    con = _con()
    explode_multi_gene_protein_mentions(con)

    rels = {tuple(r) for r in _vals(
        con,
        'SELECT relation_evidence_id, subject_entity_evidence_id, '
        'object_entity_evidence_id FROM relation_evidence_raw',
    )}
    m100, m200 = cu('M#mg=100'), cu('M#mg=200')
    r100, r200 = cu('r1#mgs=#mgo=100'), cu('r1#mgs=#mgo=200')
    # r1 (X–M) → two distinct relations to each gene copy; r2 control untouched.
    assert ('r2', 'X', 'P2ev') in rels
    assert (r100, 'X', m100) in rels
    assert (r200, 'X', m200) in rels
    # relation ids are unique (no collision after regeneration).
    all_ids = [r[0] for r in rels]
    assert len(all_ids) == len(set(all_ids))

    # The relation annotation on r1 follows both new relation ids.
    ra = {r[0] for r in _vals(con, 'SELECT evidence_id FROM relation_annotation_raw')}
    assert ra == {r100, r200}


def test_annotation_relation_subject_fans_out():
    con = _con()
    explode_multi_gene_protein_mentions(con)
    rows = {tuple(r) for r in _vals(
        con,
        'SELECT relation_evidence_id, subject_entity_evidence_id '
        'FROM annotation_relation_evidence_raw',
    )}
    assert rows == {
        (cu('a1#mgs=100'), cu('M#mg=100')),
        (cu('a1#mgs=200'), cu('M#mg=200')),
    }


def test_direct_gene_resolution_emitted():
    con = _con()
    explode_multi_gene_protein_mentions(con)
    res = {tuple(r) for r in _vals(
        con,
        'SELECT entity_evidence_id, entity_type, canonical_identifier '
        'FROM multigene_resolution',
    )}
    assert res == {
        (cu('M#mg=100'), GENE, '100'),
        (cu('M#mg=200'), GENE, '200'),
    }


def test_no_multi_gene_is_a_noop():
    con = _con()
    # Remove the second gene mapping for P1 so nothing is multi-gene.
    con.execute("DELETE FROM needed_resolver_lookup WHERE canonical_identifier = '200'")
    copies = explode_multi_gene_protein_mentions(con)
    assert copies == 0
    ids = {r[0] for r in _vals(con, 'SELECT entity_evidence_id FROM entity_evidence_raw')}
    assert ids == {'M', 'P2ev', 'X'}  # untouched


def test_explosion_preserves_identifier_normalized_column():
    # T040/T041 added entity_identifier_raw.identifier_normalized; the
    # multi-gene explosion rebuilds this table via an explicit column list
    # (_explode_one) and must carry it through unchanged, or every mention
    # caught up in a multi-gene split (any dataset with an isozyme-mapped
    # UniProt) silently loses the resolver's normalized join key downstream
    # in entity_resolution_base.
    con = _con()
    con.execute(
        "UPDATE entity_identifier_raw SET identifier_normalized = 'CHEBI:1' "
        "WHERE entity_evidence_id = 'M'"
    )
    copies = explode_multi_gene_protein_mentions(con)
    assert copies == 2

    cols = {
        r[0] for r in _vals(con, "DESCRIBE entity_identifier_raw")
    }
    assert 'identifier_normalized' in cols

    m100, m200 = cu('M#mg=100'), cu('M#mg=200')
    rows = {
        tuple(r) for r in _vals(
            con,
            'SELECT entity_evidence_id, identifier_normalized '
            'FROM entity_identifier_raw WHERE identifier = \'P1\'',
        )
    }
    assert rows == {(m100, 'CHEBI:1'), (m200, 'CHEBI:1')}
