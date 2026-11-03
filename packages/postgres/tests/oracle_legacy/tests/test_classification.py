"""Integration tests for chemical-class classification.

Run against a built instance after `derive`, e.g. on beauty::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_classification.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')
CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; classification test needs a built database',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _scalar(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchone()[0]


def test_vocab_chemical_class_seeded(conn):
    """The chemical-class CV is seeded (6 classes with precedence)."""
    assert _scalar(conn, f'SELECT count(*) FROM {SCHEMA}.vocab_chemical_class') >= 5
    names = {
        row
        for (row,) in _rows(conn, f'SELECT name FROM {SCHEMA}.vocab_chemical_class')
    }
    assert {'metabolite', 'lipid', 'drug', 'food', 'other'} <= names


def test_every_chemical_has_a_class(conn):
    """Every Chemical:OM:0037 entity carries a chemical_class_id."""
    unclassified = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity e
        JOIN {SCHEMA}.vocab_entity_type v ON v.entity_type_id = e.entity_type_id
        WHERE v.name = %s AND e.chemical_class_id IS NULL
        """,
        [CHEMICAL_ENTITY_TYPE],
    )
    assert unclassified == 0


def test_multiple_chemical_classes_present(conn):
    """Chemicals resolve into several classes, not one undifferentiated bucket."""
    n_classes = _scalar(
        conn,
        f"""
        SELECT count(DISTINCT e.chemical_class_id) FROM {SCHEMA}.entity e
        JOIN {SCHEMA}.vocab_entity_type v ON v.entity_type_id = e.entity_type_id
        WHERE v.name = %s
        """,
        [CHEMICAL_ENTITY_TYPE],
    )
    assert n_classes >= 2


def test_chemical_class_facet_present(conn):
    """The chemical_class facet is built into facet_entity_bitmap."""
    assert (
        _scalar(
            conn,
            f"SELECT count(*) FROM {SCHEMA}.facet_entity_bitmap "
            f"WHERE facet_name = 'chemical_class'",
        )
        > 0
    )


# --- metabolic-domain + interaction-class ---------------------------------


def test_vocab_metabolic_domain_seeded(conn):
    """The metabolic-domain CV is seeded with the coarse buckets."""
    names = {
        row
        for (row,) in _rows(
            conn, f'SELECT name FROM {SCHEMA}.vocab_metabolic_domain'
        )
    }
    assert {
        'lipid',
        'nucleotide',
        'amino_acid',
        'carbohydrate',
        'cofactor_vitamin',
        'other',
    } <= names


def test_every_chemical_has_a_metabolic_domain(conn):
    """Every Chemical:OM:0037 entity carries a metabolic_domain_id."""
    unclassified = _scalar(
        conn,
        f"""
        SELECT count(*) FROM {SCHEMA}.entity e
        JOIN {SCHEMA}.vocab_entity_type v ON v.entity_type_id = e.entity_type_id
        WHERE v.name = %s AND e.metabolic_domain_id IS NULL
        """,
        [CHEMICAL_ENTITY_TYPE],
    )
    assert unclassified == 0


def test_multiple_metabolic_domains_present(conn):
    """Chemicals resolve into several metabolic domains, not just `other`."""
    n_domains = _scalar(
        conn,
        f"""
        SELECT count(DISTINCT e.metabolic_domain_id) FROM {SCHEMA}.entity e
        JOIN {SCHEMA}.vocab_entity_type v ON v.entity_type_id = e.entity_type_id
        WHERE v.name = %s
        """,
        [CHEMICAL_ENTITY_TYPE],
    )
    assert n_domains >= 2


def test_metabolic_domain_facet_present(conn):
    """The metabolic_domain facet is built into facet_entity_bitmap."""
    assert (
        _scalar(
            conn,
            f"SELECT count(*) FROM {SCHEMA}.facet_entity_bitmap "
            f"WHERE facet_name = 'metabolic_domain'",
        )
        > 0
    )


def test_vocab_interaction_class_seeded(conn):
    """The interaction-class CV is seeded, by snake_case slug."""
    names = {
        row
        for (row,) in _rows(
            conn, f'SELECT name FROM {SCHEMA}.vocab_interaction_class'
        )
    }
    assert {'signaling', 'transport', 'other'} <= names


def test_every_predicate_has_an_interaction_class(conn):
    """No vocab_relation_predicate row is left with a NULL interaction_class_id."""
    null_predicates = _scalar(
        conn,
        f'SELECT count(*) FROM {SCHEMA}.vocab_relation_predicate '
        f'WHERE interaction_class_id IS NULL',
    )
    assert null_predicates == 0


def _rows(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()
