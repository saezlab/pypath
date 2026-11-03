"""The network-view framework, after both datasets became presets.

This file used to check fifteen materialized views and the combined contract
they unioned into. There is no such contract left to check: a dataset is a
``network_registry`` row now, and the query that answers it is assembled from
that row by the serving side. So what is left here is the framework's own
promise — that a dataset is metadata, that a dataset with a recipe stores the
recipe, and that nothing is materialized on any dataset's behalf.

**The old views are still on disk, and this file deliberately does not test
them.** Converting a dataset stops the framework managing its views; it does
not drop them, because the only ``DROP`` for each relation sits inside the SQL
file the derive has stopped executing. They stand frozen at their last refresh,
serving data the build has since corrected. Asserting their contents here would
claim they are maintained, which is exactly what they are not. Their retirement
is a step of its own, and until it runs they are a rollback path.

What replaced the content assertions: the metabolite gate, the inline
annotations and the curation now belong to the served dataset, and the
api-service checks them against the record the dataset is drawn from.

Run against a built instance, e.g. on dev4::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_network_views.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os

import pytest

from omnipath_build.network_views import NETWORKS

DATABASE_URL = os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; network-view test needs a built database',
)

REGISTERED_NETWORKS = ('metalinksdb', 'liana')

# The operations a stored recipe may name. Two of their orders are binding and
# the registry has to be able to state them.
COMPOSITION_OPERATIONS = ('union', 'collapse', 'exclude', 'annotate')


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


def _rows(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()


def _scalar(conn, query, params=None):
    return _rows(conn, query, params)[0][0]


def _registry_rows(conn, names=REGISTERED_NETWORKS):
    """The registry rows for the named datasets, in the order they were asked for."""
    by_name = {
        row[0]: row
        for row in _rows(
            conn,
            'SELECT name, kind, schema_name, combined_relation, included_sources, '
            'interaction_class_scope, default_attributes, composition '
            'FROM public.network_registry WHERE name = ANY(%s)',
            [list(names)],
        )
    }
    missing = [name for name in names if name not in by_name]
    assert not missing, f'not registered: {", ".join(missing)}'
    return [by_name[name] for name in names]


def test_registry_lists_both_networks(conn):
    names = {row for (row,) in _rows(conn, 'SELECT name FROM public.network_registry')}
    assert set(REGISTERED_NETWORKS) <= names


def test_every_registered_dataset_is_a_preset(conn):
    """No dataset owns a relation any more.

    A row that still names a schema and a combined relation is a dataset the
    serving side would read from a view instead of from the record — a second
    copy of the same data, refreshed on a different schedule, and only one of
    the two carries the build's corrections.
    """
    owning = [
        (name, schema, combined)
        for name, _kind, schema, combined, *_rest in _registry_rows(conn)
        if schema is not None or combined is not None
    ]
    assert not owning, (
        f'{owning} still name a relation of their own; the dataset would be '
        f'served from a copy rather than from the record'
    )


def test_registry_metadata_well_formed(conn):
    """Every row carries what a caller needs to resolve the dataset.

    An empty interaction-class scope is a legitimate answer — it means the
    dataset spans every class, which `metalinksdb` does — so it is not asserted
    to be non-empty. What must be there is the resource list and the attributes
    the dataset returns when the caller asks for none.
    """
    for name, kind, _schema, _combined, sources, _classes, attributes, _recipe in (
        _registry_rows(conn)
    ):
        assert kind, f'{name} has no kind'
        assert isinstance(sources, list) and len(sources) >= 1, (
            f'{name} contributes from no source'
        )
        assert attributes, f'{name} returns nothing when the caller asks for nothing'


def test_a_dataset_with_a_recipe_stores_it_whole(conn):
    """A composition is components plus ordered steps, or it is not one."""
    composed = [
        (name, recipe)
        for name, *_rest, recipe in _registry_rows(conn)
        if recipe is not None
    ]
    if not composed:
        pytest.skip('no registered dataset is a composition')
    for name, recipe in composed:
        assert recipe.get('operation') in COMPOSITION_OPERATIONS, (
            f'{name} names the operation {recipe.get("operation")!r}'
        )
        assert recipe.get('components'), f'{name} composes no component'
        steps = [step.get('operation') for step in recipe.get('steps') or []]
        if 'exclude' in steps and 'collapse' in steps:
            assert steps.index('exclude') < steps.index('collapse'), (
                f'{name} excludes after it folds, which leaves the excluded '
                f'resource inside source_count, the references and the sign flags'
            )


def test_the_framework_materializes_nothing_for_a_dataset(conn):
    """No definition carries SQL, so no build step creates a view for one."""
    carrying = [
        definition.name for definition in NETWORKS
        if definition.sql_files or definition.matviews
    ]
    assert not carrying, (
        f'{carrying} still carry curated SQL; a dataset that materialises '
        f'something is the mechanism this cycle retires, and the next dataset '
        f'onboarded beside it would copy the pattern'
    )


def test_absent_view_is_detectable(conn):
    """A missing relation resolves to NULL via to_regclass (API → 503)."""
    assert _scalar(conn, "SELECT to_regclass('custom_views.no_such_network')") is None
