"""Integration tests for the interaction-type registry.

``vocab_interaction_class`` is the single canonical interaction-type
vocabulary: a snake_case ``name`` slug for storage and filtering, a capitalised
``label`` for display, and a nullable ``controlled_vocabulary_mapping`` slot for
an ontology CURIE. The legacy OmniPath dataset names are preset scopes, never
rows here.

Run against a built instance after ``derive``, e.g.::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_interaction_types.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os
import re

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

SNAKE_CASE = re.compile(r'^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$')

# The canonical eight classes of `vocab_interaction_class`.
EXPECTED_SLUGS = {
    'signaling',
    'tf_target',
    'allosteric',
    'orthosteric',
    'transport',
    'ligand_receptor',
    'maturation',
    'other',
}

# Legacy OmniPath dataset names: preset identities, not type rows.
LEGACY_DATASET_NAMES = {
    'post_translational',
    'transcriptional',
    'post_transcriptional',
    'mirna_transcriptional',
    'lncrna_post_transcriptional',
}

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; interaction-type test needs a built database',
)


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    # Read-only tests: autocommit keeps one failing query from aborting the
    # transaction and masking the rest with InFailedSqlTransaction.
    connection.autocommit = True
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture(scope='module')
def classes(conn):
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT name, label, controlled_vocabulary_mapping '
            f'FROM {SCHEMA}.vocab_interaction_class ORDER BY name'
        )
        return cur.fetchall()


def test_registry_holds_the_canonical_classes(classes):
    """The eight canonical interaction classes are present, by slug."""
    assert {name for name, _, _ in classes} == EXPECTED_SLUGS


def test_every_name_is_a_snake_case_slug(classes):
    """`name` is the filterable slug, so it never carries a display spelling."""
    offenders = [name for name, _, _ in classes if not SNAKE_CASE.match(name)]
    assert offenders == []


def test_every_class_has_a_distinct_capitalised_label(classes):
    """Every class carries a display label, capitalised and unique."""
    labels = [label for _, label, _ in classes]
    assert all(label for label in labels), 'a class has no display label'
    assert all(label[0].isupper() for label in labels), (
        f'labels not capitalised: {[la for la in labels if not la[0].isupper()]}'
    )
    assert len(set(labels)) == len(labels), 'display labels are not distinct'


def test_slug_and_label_are_separate_fields(classes):
    """No class is filterable only by its display form."""
    assert all(name != label for name, label, _ in classes)


def test_controlled_vocabulary_mapping_slot_exists(conn):
    """The CV/ontology mapping slot exists and is nullable."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT is_nullable FROM information_schema.columns
            WHERE table_schema = %s
              AND table_name = 'vocab_interaction_class'
              AND column_name = 'controlled_vocabulary_mapping'
            """,
            [SCHEMA],
        )
        row = cur.fetchone()
    assert row is not None, 'controlled_vocabulary_mapping column is missing'
    assert row[0] == 'YES'


def test_no_legacy_dataset_name_is_a_type_row(classes):
    """Legacy dataset names are preset scopes, not type-vocabulary rows."""
    present = {name for name, _, _ in classes} & LEGACY_DATASET_NAMES
    assert present == set(), f'legacy dataset names used as types: {present}'
    labels = {(label or '').lower().replace('-', '_') for _, label, _ in classes}
    assert labels & LEGACY_DATASET_NAMES == set()
