"""License filtering over the interaction record.

A license is not a name. `pypath/internals/license.py` models it as three
**ordinal** levels — ``purpose``, ``sharing``, ``attrib`` — and defines
``enables(other)`` as ``self >= other``, so "usable commercially" is a
comparison over three smallints and never a match on ``license_name``. Two
resources under different names can permit the same use, and one name under two
versions can not. `data_source_license` stores the levels, and a license
question resolves to a **set of resources** before anything touches the
interaction tables.

Two properties this file exists to pin:

**Unknown terms are exclusions, never defaults.** A resource with ``is_known =
false`` MUST NOT appear in a license-filtered result. This is the one failure
mode of the table that cannot be detected downstream — an admitted unknown
looks exactly like an admitted known — so it gets its own test, and the fixture
is adversarial: the unknown resource is stored with the *most permissive levels
the vocabulary has*. If the resolution consults the levels rather than
``is_known``, it admits it, and the test says so.

**A license filter restricts the resource set, so the scope rule applies.** The
surviving summaries must be recomputed by collapsing
``interaction_fact_resource`` over the surviving resources. Reusing a fold
computed over a **wider** set returns the right interactions carrying numbers
that describe resources the license excluded. The last test holds that
shortcut to be a defect by showing the two answers differ.

**Amended 2026-08-21**: the wider fold used to be a table —
``interaction_fact_combined``, the all-resources scope materialised — and it is
removed, so the last test folds the record over every resource instead of
reading one. The rule survives the table because it was never about the table:
it forbids reusing a fold computed over a wider set, whether that fold is
stored or computed a moment ago.

Run::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_interactions_license.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os

import pytest

from tests.fixtures.interaction_graph import (
    SOURCE_A,
    SOURCE_B,
    SOURCE_C,
    SOURCE_LR,
    SOURCE_NAMES,
    build_interaction_fixture,
)

DATABASE_URL = os.environ.get('DATABASE_URL')
SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_INTERACTIONS_LICENSE',
    'interactions_license_test',
)

# The ordinal vocabularies of `pypath/internals/license.py`, as
# `data_source_license` stores them. Ascending: a higher level enables
# everything below it.
PURPOSE = {
    'ignore': 0,
    'academic': 5,
    'nonprofit': 10,
    'commercial': 15,
    'free': 20,
    'composite': 25,
}
SHARING = {
    'ignore': 0,
    'noshare': 5,
    'noderiv': 10,
    'alike': 15,
    'share': 20,
    'free': 25,
}
ATTRIB = {'ignore': 0, 'attrib': 5, 'free': 10}

# The license terms the four fixture resources hold. Deliberately arranged so
# that the requested level (`commercial`) matches **no** stored `license_name`:
# the two resources that qualify sit *above* it under two unrelated names, and
# the one that does not sits below it. A name match cannot produce this answer.
FIXTURE_LICENSES = {
    SOURCE_A: {
        'license_name': 'CC-BY-4.0',
        'license_full_name': 'Creative Commons Attribution 4.0 International',
        'license_url': 'https://creativecommons.org/licenses/by/4.0/',
        'purpose_level': PURPOSE['free'],
        'sharing_level': SHARING['free'],
        'attrib_level': ATTRIB['attrib'],
        'is_known': True,
    },
    SOURCE_B: {
        'license_name': 'Academic-only',
        'license_full_name': 'Free for academic use',
        'license_url': None,
        'purpose_level': PURPOSE['academic'],
        'sharing_level': SHARING['noshare'],
        'attrib_level': ATTRIB['attrib'],
        'is_known': True,
    },
    SOURCE_C: {
        'license_name': 'Composite',
        'license_full_name': 'Composite of the terms of its contributors',
        'license_url': None,
        'purpose_level': PURPOSE['composite'],
        'sharing_level': SHARING['free'],
        'attrib_level': ATTRIB['free'],
        'is_known': True,
    },
    # No license record maps to this resource. The levels stored beside the
    # flag are the *most permissive the vocabulary offers*, on purpose: the
    # exclusion must follow from `is_known`, not from the numbers. A resolution
    # that reads the levels first admits this resource and fails the test.
    SOURCE_LR: {
        'license_name': None,
        'license_full_name': None,
        'license_url': None,
        'purpose_level': PURPOSE['composite'],
        'sharing_level': SHARING['free'],
        'attrib_level': ATTRIB['free'],
        'is_known': False,
    },
}

# The filter every test below asks: "usable commercially".
COMMERCIAL = {'purpose_level': PURPOSE['commercial']}

NAME_A = SOURCE_NAMES[SOURCE_A]
NAME_B = SOURCE_NAMES[SOURCE_B]
NAME_C = SOURCE_NAMES[SOURCE_C]
NAME_LR = SOURCE_NAMES[SOURCE_LR]

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the license-filter test needs a Postgres',
)


@pytest.fixture(scope='module')
def built():
    """The fixture graph, projected by the derive step, in a scratch schema."""
    import psycopg2

    from omnipath_build.db import schema as build_schema
    from omnipath_build.db.derived_tables import rebuild_interaction_tables

    connection = psycopg2.connect(DATABASE_URL)
    try:
        build_schema.ensure_schema(
            connection,
            schema=SCRATCH,
            drop_existing=True,
        )
        connection.commit()
        build_interaction_fixture(connection, SCRATCH)
        rebuild_interaction_tables(connection, schema=SCRATCH)
        # From here the tests only read, and write the license fixture. In
        # autocommit a failing statement aborts nothing, so a missing table is
        # reported once per test rather than masked by InFailedSqlTransaction
        # on every test after the first.
        connection.autocommit = True
        yield connection
    finally:
        with connection.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {SCRATCH} CASCADE')
        connection.close()


def _load_fixture_licenses(conn) -> None:
    """Record the four fixture resources' terms in `data_source_license`.

    Called from each test rather than from the module fixture, so that a
    missing table is reported as a failing test rather than as a collection
    error — the point of the task is that these tests fail on the table that
    does not exist yet. The connection is in autocommit, so a statement that
    fails leaves no aborted transaction to mask the next test's reason.
    """
    for source_id, terms in FIXTURE_LICENSES.items():
        with conn.cursor() as cur:
            cur.execute(
                f"""
                INSERT INTO {SCRATCH}.data_source_license
                  (source_id, license_name, license_full_name, license_url,
                   purpose_level, sharing_level, attrib_level, is_known)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (source_id) DO UPDATE SET
                  license_name = EXCLUDED.license_name,
                  license_full_name = EXCLUDED.license_full_name,
                  license_url = EXCLUDED.license_url,
                  purpose_level = EXCLUDED.purpose_level,
                  sharing_level = EXCLUDED.sharing_level,
                  attrib_level = EXCLUDED.attrib_level,
                  is_known = EXCLUDED.is_known
                """,
                (
                    source_id,
                    terms['license_name'],
                    terms['license_full_name'],
                    terms['license_url'],
                    terms['purpose_level'],
                    terms['sharing_level'],
                    terms['attrib_level'],
                    terms['is_known'],
                ),
            )


def _resolve(conn, **minimum_levels: int) -> dict[str, object]:
    """Resolve a license question to a resource set, and say what it dropped.

    The reference semantics the implementation owes, written out so the
    assertions below are about behaviour rather than about an import:

    * a resource whose license is **unknown** is excluded first, whatever its
      stored levels say;
    * a resource whose license is known is admitted when **every** requested
      level is met under ``enables(other) == self >= other``.

    Returns ``{'admitted': [names], 'excluded': {name: reason}}``, so the
    caller can name the resources the filter removed.
    """
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT s.name, l.is_known, l.purpose_level, l.sharing_level,
                   l.attrib_level
            FROM {SCRATCH}.data_source AS s
            LEFT JOIN {SCRATCH}.data_source_license AS l
              ON l.source_id = s.source_id
            ORDER BY s.name
            """
        )
        rows = cur.fetchall()

    admitted: list[str] = []
    excluded: dict[str, str] = {}
    for name, is_known, purpose, sharing, attrib in rows:
        if not is_known:
            # Includes the source with no row at all: absence of terms is not
            # permission.
            excluded[name] = 'license_unknown'
            continue
        stored = {
            'purpose_level': purpose,
            'sharing_level': sharing,
            'attrib_level': attrib,
        }
        below = [
            level
            for level, minimum in minimum_levels.items()
            if stored[level] is None or stored[level] < minimum
        ]
        if below:
            excluded[name] = f'below_{below[0]}'
        else:
            admitted.append(name)
    return {'admitted': admitted, 'excluded': excluded}


def _collapse(conn, sources: list[str] | None) -> dict[tuple[str, str], dict]:
    """Collapse the record over ``sources``; ``None`` means every resource.

    This is the collapse the scope rule requires: the summaries are recomputed
    from ``interaction_fact_resource`` over the rows that survived the
    restriction. ``bool_or`` gives the three-valued answer, because no record
    row ever carries an asserted ``false``.
    """
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT subject.canonical_identifier,
                   object.canonical_identifier,
                   array_agg(DISTINCT s.name ORDER BY s.name) AS sources,
                   count(DISTINCT r.source_id) AS source_count,
                   bool_or(r.is_directed) AS is_directed,
                   bool_or(r.is_stimulation) AS is_stimulation,
                   bool_or(r.is_inhibition) AS is_inhibition,
                   count(DISTINCT r.source_id) FILTER (
                     WHERE r.is_stimulation IS NOT NULL
                        OR r.is_inhibition IS NOT NULL
                   ) AS sign_source_count,
                   count(DISTINCT r.source_id) FILTER (
                     WHERE r.is_directed IS NOT NULL
                   ) AS direction_source_count
            FROM {SCRATCH}.interaction_fact_resource AS r
            JOIN {SCRATCH}.data_source AS s ON s.source_id = r.source_id
            JOIN {SCRATCH}.entity AS subject
              ON subject.entity_id = r.subject_entity_id
            JOIN {SCRATCH}.entity AS object
              ON object.entity_id = r.object_entity_id
            WHERE %(sources)s::text[] IS NULL
               OR s.name = ANY(%(sources)s::text[])
            GROUP BY 1, 2, r.interaction_class_id
            """,
            {'sources': sources},
        )
        rows = cur.fetchall()

    keys = (
        'sources',
        'source_count',
        'is_directed',
        'is_stimulation',
        'is_inhibition',
        'sign_source_count',
        'direction_source_count',
    )
    return {
        (subject.removeprefix('FIXTURE_'), object_.removeprefix('FIXTURE_')): dict(
            zip(keys, rest)
        )
        for subject, object_, *rest in rows
    }


def test_the_license_table_stores_ordinal_levels_not_only_a_name(built):
    """`data_source_license` holds three levels and `is_known` per resource."""
    conn = built
    _load_fixture_licenses(conn)
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT column_name, data_type
            FROM information_schema.columns
            WHERE table_schema = %s AND table_name = 'data_source_license'
            """,
            (SCRATCH,),
        )
        columns = dict(cur.fetchall())

    assert set(columns) >= {
        'source_id',
        'license_name',
        'purpose_level',
        'sharing_level',
        'attrib_level',
        'is_known',
    }, f'data_source_license is missing columns: {sorted(columns)}'
    for level in ('purpose_level', 'sharing_level', 'attrib_level'):
        assert columns[level] == 'smallint', (
            f'{level} is {columns[level]}, not an ordinal smallint — '
            'a license question is a comparison, not a name match'
        )
    assert columns['is_known'] == 'boolean'


def test_the_filter_is_a_comparison_over_levels_not_a_name_match(built):
    """`enables(other)` is `self >= other`, so `free` and `composite` qualify.

    The requested level, `commercial`, is the name of no stored license. The
    two resources that qualify sit above it under two unrelated names, and the
    one that does not sits below it under a third. Only an ordinal comparison
    produces this partition.
    """
    conn = built
    _load_fixture_licenses(conn)

    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT s.name, l.license_name, l.purpose_level
            FROM {SCRATCH}.data_source AS s
            JOIN {SCRATCH}.data_source_license AS l
              ON l.source_id = s.source_id
            WHERE l.is_known
            """
        )
        stored = {name: (license_name, level) for name, license_name, level in cur}

    threshold = PURPOSE['commercial']
    # Above the bar, admitted: `free` (20) and `composite` (25).
    assert stored[NAME_A][1] >= threshold
    assert stored[NAME_C][1] >= threshold
    # Below the bar, excluded: `academic` (5).
    assert stored[NAME_B][1] < threshold
    # And no name match could have said so.
    assert 'commercial' not in {value[0] for value in stored.values()}
    assert stored[NAME_A][0] != stored[NAME_C][0], (
        'the two admitted resources must carry different license names, '
        'or the test cannot tell a comparison from a name match'
    )

    resolved = _resolve(conn, **COMMERCIAL)
    assert set(resolved['admitted']) == {NAME_A, NAME_C}
    assert NAME_B not in resolved['admitted']


def test_the_filter_names_the_resources_it_excluded(built):
    """The excluded resources are named, with the reason each fell."""
    conn = built
    _load_fixture_licenses(conn)

    resolved = _resolve(conn, **COMMERCIAL)
    assert resolved['excluded'] == {
        NAME_B: 'below_purpose_level',
        NAME_LR: 'license_unknown',
    }, f'the filter did not name what it dropped: {resolved["excluded"]}'


def test_an_unknown_license_is_never_admitted_by_a_permissive_default(built):
    """`is_known = false` excludes, whatever the stored levels say.

    The fixture stores this resource with the highest level in every
    vocabulary, so a resolution that consults the numbers admits it under any
    request at all. Its own test, because an admitted unknown is
    indistinguishable from an admitted known once it reaches a caller.
    """
    conn = built
    _load_fixture_licenses(conn)

    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT is_known, purpose_level, sharing_level, attrib_level
            FROM {SCRATCH}.data_source_license
            WHERE source_id = %s
            """,
            (SOURCE_LR,),
        )
        is_known, purpose, sharing, attrib = cur.fetchone()
    assert is_known is False
    assert (purpose, sharing, attrib) == (
        PURPOSE['composite'],
        SHARING['free'],
        ATTRIB['free'],
    ), 'the adversarial fixture stopped being adversarial'

    # The most permissive request the vocabulary can express: every level at
    # `ignore`. It still must not admit a resource with unknown terms.
    permissive = _resolve(
        conn,
        purpose_level=PURPOSE['ignore'],
        sharing_level=SHARING['ignore'],
        attrib_level=ATTRIB['ignore'],
    )
    assert NAME_LR not in permissive['admitted'], (
        'a resource with no mapped license was admitted under a permissive '
        'default — the one failure mode of this table that cannot be '
        'detected downstream'
    )
    assert permissive['excluded'][NAME_LR] == 'license_unknown'

    # And under the commercial filter, likewise.
    assert NAME_LR not in _resolve(conn, **COMMERCIAL)['admitted']


# --- the filtered result is a strict subset, and it is recomputed ------------


def test_a_license_filtered_result_is_a_strict_subset(built):
    """A license-filtered result is strictly smaller, and every survivor was
    there before."""
    conn = built
    _load_fixture_licenses(conn)

    unfiltered = _collapse(conn, None)
    filtered = _collapse(conn, sorted(_resolve(conn, **COMMERCIAL)['admitted']))

    assert set(filtered) < set(unfiltered), (
        'the license-filtered result is not a proper subset of the '
        'unfiltered one'
    )
    assert len(filtered) < len(unfiltered), (
        f'the filter removed nothing: {len(filtered)} rows either way'
    )


def test_the_interactions_the_filter_removed_are_named(built):
    """The ligand/receptor pair goes, because its only resource is unknown."""
    conn = built
    _load_fixture_licenses(conn)

    unfiltered = _collapse(conn, None)
    filtered = _collapse(conn, sorted(_resolve(conn, **COMMERCIAL)['admitted']))
    removed = set(unfiltered) - set(filtered)

    assert ('a', 'b') in removed, (
        'a->b is contributed only by the unknown-license resource and must '
        'not survive a license filter'
    )
    assert unfiltered[('a', 'b')]['sources'] == [NAME_LR]
    # Nothing survives on the strength of a resource the filter excluded.
    for key, row in filtered.items():
        assert NAME_LR not in row['sources'], (
            f'{key} still credits the unknown-license resource'
        )
        assert NAME_B not in row['sources'], (
            f'{key} still credits an academic-only resource under a '
            'commercial filter'
        )


def test_the_surviving_summaries_are_recomputed_over_the_survivors(built):
    """The scope rule: c->d loses a resource, and its numbers change.

    Three resources report c->d: one asserts stimulation, one inhibition, one
    neither. The inhibiting resource is academic-only, so a commercial filter
    drops it — and with it the inhibition. A result that still says
    `is_inhibition = true` is describing a resource the caller excluded.
    """
    conn = built
    _load_fixture_licenses(conn)

    filtered = _collapse(conn, sorted(_resolve(conn, **COMMERCIAL)['admitted']))
    row = filtered[('c', 'd')]

    assert sorted(row['sources']) == [NAME_A, NAME_C]
    assert row['source_count'] == 2
    assert row['is_stimulation'] is True
    assert row['is_inhibition'] is None, (
        'the inhibition survived the removal of the only resource asserting '
        'it — the summary was read off a wider scope'
    )
    assert row['sign_source_count'] == 1


def test_reusing_a_wider_fold_reports_the_wrong_numbers(built):
    """Why the shortcut is a defect, not an optimisation.

    The all-resources fold of c->d says three resources and an inhibition.
    Filtering *that* row down to the license's resource set — with `sources &&
    ARRAY['fixture_res_a','fixture_res_c']`, or by any other predicate over an
    already-folded row — returns the right interaction carrying a source count
    and a sign that describe three resources, one of which the license
    excluded. Only re-folding over the surviving resources answers the
    question that was asked.

    The test now computes the wider fold rather than reading it from
    `interaction_fact_combined`, which is removed. That changes where the wrong
    numbers come from and not whether they are wrong.
    """
    conn = built
    _load_fixture_licenses(conn)

    wider = _collapse(conn, None)[('c', 'd')]
    assert wider['source_count'] == 3
    assert wider['is_inhibition'] is True

    scoped = _collapse(
        conn,
        sorted(_resolve(conn, **COMMERCIAL)['admitted']),
    )[('c', 'd')]

    assert scoped['source_count'] != wider['source_count']
    assert scoped['is_inhibition'] != wider['is_inhibition'], (
        'the scoped collapse agreed with the all-resources fold, so this '
        'fixture no longer proves the scope rule bites'
    )
