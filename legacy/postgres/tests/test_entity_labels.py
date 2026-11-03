"""Every entity has a stored human-readable label.

Asserts the post-``derive`` label state:
- **100%** of entities have a non-empty stored ``label``;
- gene entities are labelled by their **gene symbol** (``label_rule = gene_symbol``),
  and the EGFR benchmark shows ``EGFR``;
- the producing rule is recorded in ``label_rule``.

The chemical brevity-first cascade — a chemical shows a short name, not a raw
InChIKey or id, where a name exists — is a separate piece of work.

Run against a built instance, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55432/omnipath \
        uv run pytest tests/test_entity_labels.py -v

Skipped when DATABASE_URL is not set (no DB to test against).
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

GENE_ENTITY_TYPE = 'Gene:MI:0250'
CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; entity-label test needs a built database',
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


def _type_id(conn, name: str):
    return _scalar(
        conn,
        f'SELECT entity_type_id FROM {SCHEMA}.vocab_entity_type WHERE name = %s',
        [name],
    )


def test_entity_label_column_exists(conn):
    present = _scalar(
        conn,
        """
        SELECT count(*) FROM information_schema.columns
        WHERE table_schema = %s AND table_name = 'entity'
          AND column_name IN ('label', 'label_rule')
        """,
        [SCHEMA],
    )
    if present != 2:
        pytest.skip('entity.label/label_rule absent — build predates labels')


def test_all_entities_have_a_label(conn):
    """100% of entities carry a non-empty stored label."""
    missing = _scalar(
        conn,
        f"SELECT count(*) FROM {SCHEMA}.entity WHERE label IS NULL OR label = ''",
    )
    assert missing == 0, f'{missing} entities without a non-empty label'


def test_label_rule_recorded(conn):
    """Every labelled entity records the rule that produced its label."""
    missing_rule = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE label IS NOT NULL AND label <> ''
          AND (label_rule IS NULL OR label_rule = '')
        """,
    )
    assert missing_rule == 0, f'{missing_rule} labels without a label_rule'


def test_gene_entities_labelled_by_symbol(conn):
    """Most gene entities are labelled by gene symbol; the rest fall back to id."""
    gene_type_id = _type_id(conn, GENE_ENTITY_TYPE)
    if gene_type_id is None:
        pytest.skip('no Gene entity type')
    total = _scalar(
        conn,
        f'SELECT count(*) FROM {SCHEMA}.entity WHERE entity_type_id = %s',
        [gene_type_id],
    )
    by_symbol = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE entity_type_id = %s AND label_rule = 'gene_symbol'
        """,
        [gene_type_id],
    )
    assert total > 0
    # The vast majority resolve to a primary symbol; only symbol-less genes fall
    # back to the NCBI Gene id.
    assert by_symbol >= 0.9 * total, (
        f'only {by_symbol}/{total} gene entities labelled by symbol'
    )


def test_egfr_labelled_egfr(conn):
    gene_type_id = _type_id(conn, GENE_ENTITY_TYPE)
    if gene_type_id is None:
        pytest.skip('no Gene entity type')
    label = _scalar(
        conn,
        f"""
        SELECT label FROM {SCHEMA}.entity
        WHERE entity_type_id = %s AND canonical_identifier = '1956' AND taxonomy_id = 9606
        """,
        [gene_type_id],
    )
    if label is None:
        pytest.skip('EGFR gene entity absent in this (capped) build')
    assert label == 'EGFR', f'EGFR labelled {label!r} (expected EGFR)'


def test_chemical_has_a_name_not_a_raw_identifier(conn):
    """No chemical is labelled by a raw InChIKey or an opaque hash."""
    chemical_type_id = _type_id(conn, CHEMICAL_ENTITY_TYPE)
    if chemical_type_id is None:
        pytest.skip('no Chemical entity type')
    # An InChIKey-shaped label (e.g. ``XXXXXXXXXXXXXX-YYYYYYYYYY-Z``) or the
    # 32-hex ``unresolved_entity_key`` hash means no human-readable name/id was
    # selected — the label cascade must avoid both.
    raw_like = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE entity_type_id = %s
          AND (label ~ '^[A-Z]{{14}}-[A-Z]{{10}}-[A-Z]$'
               OR label ~ '^[0-9a-f]{{32}}$')
        """,
        [chemical_type_id],
    )
    assert raw_like == 0, (
        f'{raw_like} chemicals labelled by a raw InChIKey or hash'
    )


def test_chemicals_labelled_by_cascade_rules(conn):
    """Every chemical leaves the universal fallback for a cascade rule.

    Robust to the metabo Goslin lipid layer, which may overwrite some
    chemical labels with ``goslin_lipid`` after the build's derive — the
    invariant is that **no** chemical is left on the universal
    ``identifier_fallback`` (the opaque hash/InChIKey last resort).
    """
    chemical_type_id = _type_id(conn, CHEMICAL_ENTITY_TYPE)
    if chemical_type_id is None:
        pytest.skip('no Chemical entity type')
    total = _scalar(
        conn,
        f'SELECT count(*) FROM {SCHEMA}.entity WHERE entity_type_id = %s',
        [chemical_type_id],
    )
    if not total:
        pytest.skip('no chemical entities in this build')
    on_universal_fallback = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE entity_type_id = %s AND label_rule = 'identifier_fallback'
        """,
        [chemical_type_id],
    )
    assert on_universal_fallback == 0, (
        f'{on_universal_fallback}/{total} chemicals still on the universal '
        'identifier fallback (cascade did not relabel them)'
    )
    by_name = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity
        WHERE entity_type_id = %s
          AND label_rule IN
            ('chemical_name', 'chemical_iupac_name', 'goslin_lipid')
        """,
        [chemical_type_id],
    )
    # The name cascade produces real human-readable names (sanity floor). The
    # named *fraction* is data-composition-dependent — many builds carry large
    # id-only sets (FooDB/PTFI stubs) — so we don't assert a ratio here; the
    # hard invariants are "no universal fallback" (above) and "no InChIKey/hash
    # label" (test_chemical_has_a_name_not_a_raw_identifier).
    assert by_name > 0, 'the chemical name cascade produced no name-based labels'
