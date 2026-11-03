"""Validate the resolution benchmark fixtures against the live DB.

Consumes :mod:`tests.fixtures.resolution_benchmarks` and asserts that the
gene-anchored entity build matches each pinned benchmark:

* **core** benchmarks must hold (EGFR human, alanine named, global invariants);
* **present-if-in-build** benchmarks (orthologs) are ``skip``-ped when the
  entity is absent in a capped ``MAX_RECORDS`` build.

Run against a built instance, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run pytest tests/test_resolution_benchmarks.py -v

Skipped when DATABASE_URL is not set (no DB to test against).
"""

from __future__ import annotations

import os

import pytest

from fixtures.resolution_benchmarks import (
    CHEMICAL_BENCHMARKS,
    GENE_BENCHMARKS,
    GENE_ENTITY_TYPE,
    HEX32_REGEX,
    ID_TYPE_NAMES,
    INCHIKEY_REGEX,
    ChemicalBenchmark,
    GeneBenchmark,
    looks_like_inchikey_or_hash,
)

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; benchmark test needs a built database',
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


def _gene_entity_id(conn, gene_type_id, entrez: str, taxonomy: int):
    return _scalar(
        conn,
        f"""
        SELECT entity_id FROM {SCHEMA}.entity
        WHERE entity_type_id = %s
          AND canonical_identifier = %s
          AND taxonomy_id = %s
        """,
        [gene_type_id, entrez, taxonomy],
    )


# ---------------------------------------------------------------------------
# Gene benchmarks
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'bm', GENE_BENCHMARKS, ids=lambda b: f'{b.expected_entrez}/{b.expected_taxonomy}'
)
def test_gene_benchmark(conn, gene_type_id, bm: GeneBenchmark):
    entity_id = _gene_entity_id(
        conn, gene_type_id, bm.expected_entrez, bm.expected_taxonomy,
    )
    if entity_id is None:
        if bm.core:
            pytest.fail(
                f'core gene benchmark absent: {bm.description} '
                f'(entrez={bm.expected_entrez}, taxon={bm.expected_taxonomy})'
            )
        pytest.skip(
            f'present-if-in-build gene absent in capped build: {bm.description}'
        )

    # label
    label = _scalar(
        conn,
        f'SELECT label FROM {SCHEMA}.entity WHERE entity_id = %s',
        [entity_id],
    )
    assert label == bm.expected_label, (
        f'{bm.description}: label {label!r} != {bm.expected_label!r}'
    )

    # representative UniProt (+ reviewed) via gene_output view
    row = None
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT uniprot, uniprot_is_reviewed, uniprot_all
            FROM {SCHEMA}.gene_output
            WHERE ncbi_gene_id = %s
            """,
            [bm.expected_entrez],
        )
        row = cur.fetchone()
    assert row is not None, f'{bm.description}: no gene_output row'
    uniprot, reviewed, uniprot_all = row
    assert uniprot == bm.expected_uniprot, (
        f'{bm.description}: representative UniProt {uniprot!r} != '
        f'{bm.expected_uniprot!r}'
    )
    assert reviewed == bm.expected_uniprot_reviewed, (
        f'{bm.description}: uniprot_is_reviewed {reviewed!r} != '
        f'{bm.expected_uniprot_reviewed!r}'
    )
    assert bm.expected_uniprot in (uniprot_all or []), (
        f'{bm.description}: {bm.expected_uniprot} not in uniprot_all {uniprot_all!r}'
    )

    # cross-id collapse: every supplied input id type must reach this same entity
    for alias, value in bm.inputs.items():
        id_type_name = ID_TYPE_NAMES[alias]
        reached = _scalar(
            conn,
            f"""
            SELECT count(DISTINCT ei.entity_id)
            FROM {SCHEMA}.identifier_evidence ie
            JOIN {SCHEMA}.entity_identifier ei
              ON ei.identifier_id = ie.identifier_id
            JOIN {SCHEMA}.entity e ON e.entity_id = ei.entity_id
            JOIN {SCHEMA}.vocab_identifier_type vit
              ON vit.identifier_type_id = ie.identifier_type_id
            WHERE e.entity_type_id = %s
              AND e.taxonomy_id = %s
              AND vit.name = %s
              AND ie.value = %s
            """,
            [gene_type_id, bm.expected_taxonomy, id_type_name, value],
        )
        assert reached == 1, (
            f'{bm.description}: input {alias}={value!r} reaches {reached} gene '
            f'entities for taxon {bm.expected_taxonomy} (expected exactly 1)'
        )
        # and that one entity is THIS entity
        is_this = _scalar(
            conn,
            f"""
            SELECT bool_or(ei.entity_id = %s)
            FROM {SCHEMA}.identifier_evidence ie
            JOIN {SCHEMA}.entity_identifier ei
              ON ei.identifier_id = ie.identifier_id
            JOIN {SCHEMA}.entity e ON e.entity_id = ei.entity_id
            JOIN {SCHEMA}.vocab_identifier_type vit
              ON vit.identifier_type_id = ie.identifier_type_id
            WHERE e.entity_type_id = %s
              AND e.taxonomy_id = %s
              AND vit.name = %s
              AND ie.value = %s
            """,
            [entity_id, gene_type_id, bm.expected_taxonomy, id_type_name, value],
        )
        assert is_this, (
            f'{bm.description}: input {alias}={value!r} does not reach the '
            f'expected gene entity {entity_id}'
        )


# ---------------------------------------------------------------------------
# Chemical benchmarks
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    'bm', CHEMICAL_BENCHMARKS, ids=lambda b: b.expected_labels[0]
)
def test_chemical_benchmark(conn, bm: ChemicalBenchmark):
    chem_type_id = _scalar(
        conn,
        f'SELECT entity_type_id FROM {SCHEMA}.vocab_entity_type WHERE name = %s',
        [bm.expected_entity_type],
    )
    if chem_type_id is None:
        pytest.skip(f'no {bm.expected_entity_type} entity type')

    # At least one of the expected labels must exist with an acceptable rule.
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT label, label_rule
            FROM {SCHEMA}.entity
            WHERE entity_type_id = %s
              AND label = ANY(%s)
            """,
            [chem_type_id, list(bm.expected_labels)],
        )
        rows = cur.fetchall()
    if not rows:
        if bm.core:
            pytest.fail(
                f'core chemical benchmark absent: {bm.description} '
                f'(none of {bm.expected_labels})'
            )
        pytest.skip(f'present-if-in-build chemical absent: {bm.description}')

    found_label_rule = False
    for label, label_rule in rows:
        # invariant: never an InChIKey / opaque hash
        assert not looks_like_inchikey_or_hash(label), (
            f'{bm.description}: label {label!r} looks like an InChIKey/hash'
        )
        if label_rule in bm.expected_label_rules:
            found_label_rule = True
    assert found_label_rule, (
        f'{bm.description}: no entity among {bm.expected_labels} uses an '
        f'expected label_rule {bm.expected_label_rules} (got '
        f'{sorted({r for _, r in rows})})'
    )


# ---------------------------------------------------------------------------
# Global invariants
# ---------------------------------------------------------------------------


def test_genes_anchored_on_entrez(conn, gene_type_id):
    """Every gene entity is keyed by the Entrez anchor."""
    entrez_type_id = _scalar(
        conn,
        f'SELECT identifier_type_id FROM {SCHEMA}.vocab_identifier_type WHERE name = %s',
        [ID_TYPE_NAMES['entrez']],
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


def test_no_chemical_labelled_by_inchikey_or_hash(conn):
    """No chemical is left labelled by a raw InChIKey or 32-hex hash."""
    chem_type_id = _scalar(
        conn,
        "SELECT entity_type_id FROM {0}.vocab_entity_type WHERE name = %s".format(
            SCHEMA,
        ),
        [CHEMICAL_BENCHMARKS[0].expected_entity_type],
    )
    if chem_type_id is None:
        pytest.skip('no Chemical entity type')
    raw_like = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE entity_type_id = %s
          AND (label ~ %s OR label ~ %s)
        """,
        [chem_type_id, INCHIKEY_REGEX, HEX32_REGEX],
    )
    assert raw_like == 0, (
        f'{raw_like} chemicals labelled by a raw InChIKey or hash (orphaned)'
    )
