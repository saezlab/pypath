"""What ConnectomeDB2025 puts into the record, and under what name.

The reduced ligand-receptor drop needs no ingestion and no download: the
release is already loaded, and the derive projects it into
``interaction_fact_resource`` like any other resource. So the build-side
question is not "did it arrive" but "did it arrive intact, and named for the
release it is".

Both halves of that matter.

**The name.** A slug that reads ``connectomedb`` names a resource *family*
while denoting one particular release. Callers filter on the slug, so the
ambiguity is not cosmetic — a query written today against ``connectomedb``
would silently change meaning the day a second ConnectomeDB release is
onboarded, and the two would collide in the source table. The build therefore
carries ``connectomedb2025`` and nothing shorter.

**The fields.** Whatever the record says about these rows is what the resource
said, because there is no second resource on them to blur the picture:

- **Evidence type.** The release marks each pair as directly observed or
  inferred from another pair. That is the resource's main quality signal and it
  lands in ``curation_flags``.
- **Taxon.** The loaded release is the all-species drop, not the human subset.
  Human is a minority of it, so a projection that lost the taxon or narrowed to
  human would throw most of the data away invisibly.
- **Direction.** A ligand acts on a receptor. The class names its endpoints
  asymmetrically, so the ordered pair carries meaning and the rows are directed
  however coarse the predicate the ingest layer recorded.
- **No sign.** The release publishes no stimulation or inhibition at all. The
  sign columns are therefore NULL — never ``false``, which would be a positive
  claim that the pair is known not to stimulate. Nobody made that claim.
- **Partner roles.** Which protein is the ligand and which the receptor is
  recorded, so the serving layer has something to project.

Run against a built instance after ``derive``::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_liana_slice.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason = 'DATABASE_URL not set; the ligand-receptor slice needs a build',
)

#: The slug the release is loaded under.
RESOURCE = 'connectomedb2025'

#: The family name it must never be loaded under.
AMBIGUOUS_SLUG = 'connectomedb'

#: The one interaction class the release contributes.
INTERACTION_CLASS = 'ligand_receptor'

#: How the release says a pair was established.
EVIDENCE_VALUES = {'Direct', 'Inferred'}

#: NCBI taxon id of human — a small part of this release.
HUMAN = 9606

#: ``interaction_party.role_flag``: which end of the pair a participant is.
LIGAND_FLAG = 1
RECEPTOR_FLAG = 2

#: The per-entity annotation terms the release's own loader emits.
LIGAND_TERM = 'Ligand:OM:7777'
RECEPTOR_TERM = 'Receptor:OM:7778'


@pytest.fixture(scope = 'module')
def conn():
    """An open connection to the built database."""

    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)

    try:

        yield connection

    finally:

        connection.close()


def _rows(conn, query, parameters = None):
    """Run one statement and return every row of it.

    Args:
        conn: An open connection.
        query: The statement.
        parameters: Its bound parameters, if any.

    Returns:
        The rows, as tuples.
    """

    with conn.cursor() as cur:

        cur.execute(query, parameters)

        return cur.fetchall()


def _one(conn, query, parameters = None):
    """Run one statement expected to return a single row.

    Args:
        conn: An open connection.
        query: The statement.
        parameters: Its bound parameters, if any.

    Returns:
        The single row, as a tuple.
    """

    rows = _rows(conn, query, parameters)

    assert len(rows) == 1, f'expected one row, got {len(rows)}'

    return rows[0]


@pytest.fixture(scope = 'module')
def record_rows(conn) -> int:
    """How many record rows the release contributes.

    Every assertion below is over this set, so a build in which the release
    contributed nothing would pass them all vacuously. This fixture is where
    that is caught once, rather than in each test.

    Args:
        conn: An open connection.

    Returns:
        The row count.
    """

    count = _one(
        conn,
        f'SELECT count(*) FROM {SCHEMA}.interaction_fact_resource f '
        f'JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id '
        f'WHERE ds.name = %s',
        (RESOURCE,),
    )[0]

    if not count:

        pytest.fail(
            f'{RESOURCE} contributes no row to '
            f'{SCHEMA}.interaction_fact_resource. Either the derive did not '
            f'project it or the slug changed; every assertion in this file '
            f'would otherwise pass over an empty set'
        )

    return count


def test_the_release_reaches_the_record(record_rows):
    """The drop's data is projected, so the serving layer has rows to fold."""

    assert record_rows > 0


def test_every_row_of_the_release_is_a_ligand_receptor_pair(conn, record_rows):
    """The release publishes ligand-receptor pairs and nothing else.

    The class is not read off the predicate — the ingest layer records a coarse
    verb that would put most of this resource in `other`. It comes from the
    participant roles, which is why the whole release lands in one class.
    """

    classes = sorted(
        name for name, in _rows(
            conn,
            f'SELECT DISTINCT c.name '
            f'FROM {SCHEMA}.interaction_fact_resource f '
            f'JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id '
            f'JOIN {SCHEMA}.vocab_interaction_class c '
            f'  ON c.interaction_class_id = f.interaction_class_id '
            f'WHERE ds.name = %s',
            (RESOURCE,),
        )
    )

    assert classes == [INTERACTION_CLASS], (
        f'the release landed in classes {classes}; a row outside '
        f'{INTERACTION_CLASS!r} means the class was derived from the predicate '
        f'rather than from the participant roles'
    )


def test_only_the_release_named_slug_exists_in_the_source_table(conn):
    """No version-ambiguous ConnectomeDB name survives in the source table.

    This is the assertion the rename exists for. A bare family slug alongside
    the release-named one is worse than either alone: callers cannot tell which
    they are filtering on, and rows would be split between the two.
    """

    names = sorted(
        name for name, in _rows(
            conn,
            f'SELECT name FROM {SCHEMA}.data_source WHERE name ILIKE %s',
            (f'{AMBIGUOUS_SLUG}%',),
        )
    )

    assert names == [RESOURCE], (
        f'{SCHEMA}.data_source carries {names}; only the release-named slug '
        f'{RESOURCE!r} may exist, because {AMBIGUOUS_SLUG!r} names a family '
        f'and would collide with any other release of it'
    )


def test_the_resource_registries_carry_no_bare_family_slug(conn):
    """The rename reaches the registries callers and presets resolve through.

    Renaming the source row alone leaves the old name in the places that
    *reference* it, where it is just as load-bearing: a preset scoped to the
    bare slug resolves to nothing, and a resource inventory keyed on it
    describes a resource the source table no longer has.

    The frozen legacy ligand-receptor matview is deliberately outside this
    sweep. It is unmanaged, unrefreshed and retired separately; rewriting its
    rows here would only hide that it is stale in every other respect too.
    """

    stale = _rows(
        conn,
        f"""
        SELECT 'resources.resource_id', resource_id
        FROM {SCHEMA}.resources
        WHERE resource_id ILIKE %(family)s AND resource_id <> %(release)s
        UNION ALL
        SELECT 'network_registry.included_sources', name
        FROM {SCHEMA}.network_registry
        WHERE %(bare)s = ANY(included_sources)
        """,
        {
            'family': f'{AMBIGUOUS_SLUG}%',
            'release': RESOURCE,
            'bare': AMBIGUOUS_SLUG,
        },
    )

    assert stale == [], (
        f'the bare slug {AMBIGUOUS_SLUG!r} still appears in {stale}; a '
        f'reference to it resolves to no source at all'
    )


def test_the_evidence_type_lands_on_every_row(conn, record_rows):
    """Direct observation and inference are distinct claims, and both are kept.

    The release states, for every pair, whether it was seen directly or
    inferred from another pair. Losing that would present a measured
    interaction and an inferred one as the same statement.
    """

    empty = _one(
        conn,
        f'SELECT count(*) FROM {SCHEMA}.interaction_fact_resource f '
        f'JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id '
        f'WHERE ds.name = %s '
        f'  AND (f.curation_flags IS NULL OR cardinality(f.curation_flags) = 0)',
        (RESOURCE,),
    )[0]

    assert empty == 0, (
        f'{empty} of {record_rows} rows carry no evidence type in '
        f'`curation_flags`; the release states one for every pair'
    )


def test_the_evidence_vocabulary_is_the_releases_own(conn):
    """Both values appear, and nothing the release does not publish."""

    seen = {
        value for value, in _rows(
            conn,
            f'SELECT DISTINCT value '
            f'FROM {SCHEMA}.interaction_fact_resource f '
            f'JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id, '
            f'     unnest(f.curation_flags) AS value '
            f'WHERE ds.name = %s',
            (RESOURCE,),
        )
    }

    assert seen == EVIDENCE_VALUES, (
        f'`curation_flags` carries {sorted(seen)} for this release, against '
        f'the {sorted(EVIDENCE_VALUES)} it publishes; a missing value means '
        f'the distinction was collapsed, an extra one that something else '
        f'was written into the same column'
    )


def test_every_row_names_the_taxon_on_both_endpoints(conn, record_rows):
    """A pair without a species is not a usable ligand-receptor statement."""

    unnamed = _one(
        conn,
        f'SELECT count(*) FROM {SCHEMA}.interaction_fact_resource f '
        f'JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id '
        f'WHERE ds.name = %s '
        f'  AND (f.subject_organism IS NULL OR f.object_organism IS NULL)',
        (RESOURCE,),
    )[0]

    assert unnamed == 0, (
        f'{unnamed} of {record_rows} rows leave an endpoint without a taxon'
    )


def test_the_release_is_the_all_species_drop(conn, record_rows):
    """Many species, with human a minority — not a human-only resource.

    The loaded release covers every species the source publishes. A projection
    that kept human only would drop the majority of the rows, and would do it
    without any signal that it had happened.
    """

    taxa, human = _one(
        conn,
        f'SELECT count(DISTINCT f.subject_organism), '
        f'       count(*) FILTER (WHERE f.subject_organism = %s) '
        f'FROM {SCHEMA}.interaction_fact_resource f '
        f'JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id '
        f'WHERE ds.name = %s',
        (HUMAN, RESOURCE),
    )

    assert taxa > 1, (
        f'the release landed under {taxa} taxon; it is the all-species drop '
        f'and narrowing it to one is a silent loss of most of its rows'
    )
    assert 0 < human < record_rows, (
        f'human contributes {human} of {record_rows} rows, which is either '
        f'nothing or everything; it should be a minority of a multi-species '
        f'release'
    )


def test_every_row_is_directed(conn, record_rows):
    """A ligand acts on a receptor, so the ordered pair carries meaning.

    Direction here is a property of the class rather than something the
    resource asserts per row: the class names its two ends asymmetrically, and
    the reverse pair would be a different statement.
    """

    undirected = _one(
        conn,
        f'SELECT count(*) FROM {SCHEMA}.interaction_fact_resource f '
        f'JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id '
        f'WHERE ds.name = %s AND f.is_directed IS DISTINCT FROM true',
        (RESOURCE,),
    )[0]

    assert undirected == 0, (
        f'{undirected} of {record_rows} rows are not directed; the ordered '
        f'roles of the class fix the direction whatever predicate the ingest '
        f'layer recorded'
    )


def test_no_row_carries_a_sign_it_was_never_given(conn, record_rows):
    """The release publishes no sign, so the columns stay NULL — not false.

    This is the distinction the whole three-valued design exists for. NULL says
    no resource in scope asserts anything about the sign. ``false`` says the
    pair is known not to stimulate, or known not to inhibit. The release makes
    no such claim about any pair, so a column defaulted from the first to the
    second turns silence into a finding across every row of the drop, and
    nothing downstream can tell the two apart afterwards.
    """

    negative, asserted = _one(
        conn,
        f'SELECT count(*) FILTER (WHERE f.is_stimulation IS false '
        f'                           OR f.is_inhibition IS false), '
        f'       count(*) FILTER (WHERE f.is_stimulation IS NOT NULL '
        f'                           OR f.is_inhibition IS NOT NULL) '
        f'FROM {SCHEMA}.interaction_fact_resource f '
        f'JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id '
        f'WHERE ds.name = %s',
        (RESOURCE,),
    )

    assert negative == 0, (
        f'{negative} of {record_rows} rows carry an asserted `false` sign; '
        f'the release publishes no sign at all, so `false` here is a claim '
        f'nobody made'
    )
    assert asserted == 0, (
        f'{asserted} of {record_rows} rows carry a sign the release does not '
        f'publish'
    )


def test_each_partner_carries_its_role_annotation(conn, record_rows):
    """Ligand and receptor are annotated on the entities of every pair.

    The release's own loader emits these terms alongside the interactions, so
    the serving layer has something to project when a caller asks which end of
    the pair to look for in the sending cell. The bridge from the annotation to
    the entity runs through the evidence resolution, because the annotation is
    attached to the raw evidence rather than to the resolved entity.
    """

    subject_ligand, object_receptor = _one(
        conn,
        f"""
        WITH src AS (
          SELECT source_id FROM {SCHEMA}.data_source WHERE name = %(resource)s
        ),
        annotated AS (
          SELECT DISTINCT res.entity_id, a.term
          FROM {SCHEMA}.entity_evidence_annotation ea
          JOIN src ON src.source_id = ea.source_id
          JOIN {SCHEMA}.annotation a ON a.annotation_key = ea.annotation_key
          JOIN {SCHEMA}.entity_evidence_resolution res
            ON res.source_id = ea.source_id
           AND res.entity_evidence_id = ea.entity_evidence_id
          WHERE a.term IN (%(ligand)s, %(receptor)s)
        )
        SELECT
          count(*) FILTER (
            WHERE EXISTS (
              SELECT 1 FROM annotated
              WHERE annotated.entity_id = f.subject_entity_id
                AND annotated.term = %(ligand)s
            )
          ),
          count(*) FILTER (
            WHERE EXISTS (
              SELECT 1 FROM annotated
              WHERE annotated.entity_id = f.object_entity_id
                AND annotated.term = %(receptor)s
            )
          )
        FROM {SCHEMA}.interaction_fact_resource f
        JOIN src ON src.source_id = f.source_id
        """,
        {
            'resource': RESOURCE,
            'ligand': LIGAND_TERM,
            'receptor': RECEPTOR_TERM,
        },
    )

    assert subject_ligand == record_rows, (
        f'{record_rows - subject_ligand} of {record_rows} rows have a subject '
        f'the release does not annotate as a ligand'
    )
    assert object_receptor == record_rows, (
        f'{record_rows - object_receptor} of {record_rows} rows have an object '
        f'the release does not annotate as a receptor'
    )


def test_the_role_within_the_interaction_is_flagged(conn):
    """The role a partner plays *in one pair*, distinct from what it is overall.

    The per-entity annotation is a property of the protein across every pair it
    appears in, so a protein that is a ligand in one pair and a receptor in
    another carries both terms and cannot say which it is here. The party rows
    can: they record the role for this interaction. Both roles must be
    reachable, and no third value may appear in a two-role vocabulary.
    """

    flags = {
        flag: count for flag, count in _rows(
            conn,
            f"""
            SELECT p.role_flag, count(*)
            FROM {SCHEMA}.interaction_fact_resource f
            JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id
            JOIN {SCHEMA}.interaction_party p
              ON p.interaction_id = f.interaction_id
            WHERE ds.name = %s
            GROUP BY 1
            """,
            (RESOURCE,),
        )
    }

    assert set(flags) == {LIGAND_FLAG, RECEPTOR_FLAG}, (
        f'the party rows of this release carry role flags {sorted(flags)}; '
        f'the class has exactly two roles and both must be reachable'
    )


@pytest.mark.xfail(
    strict = True,
    reason = (
        'known open defect: the pairwise role flag marks ligand on object-side '
        'party rows and receptor on no subject-side row at all, leaving 3,844 '
        'pairs of this release with no receptor end. The per-entity annotation '
        'route disagrees and is complete, so nothing serving is blocked, but '
        'the flag must not be projected as a role until this is settled - '
        'either these are genuine ligand-ligand pairs the class flattens, or '
        'the flag is wrong. When it is fixed this test passes and the strict '
        'marker fails the run, which is the reminder to delete the marker.'
    ),
)
def test_every_pair_names_both_of_its_roles(conn, record_rows):
    """Each pair records a ligand end and a receptor end, not just one.

    A pair with only one flagged participant leaves the other end's role to be
    guessed from column order, which is exactly what the flag exists to avoid.
    """

    incomplete = _one(
        conn,
        f"""
        SELECT count(*) FROM (
          SELECT f.interaction_id
          FROM {SCHEMA}.interaction_fact_resource f
          JOIN {SCHEMA}.data_source ds ON ds.source_id = f.source_id
          JOIN {SCHEMA}.interaction_party p
            ON p.interaction_id = f.interaction_id
          WHERE ds.name = %s
          GROUP BY f.interaction_id
          HAVING NOT (bool_or(p.role_flag = %s) AND bool_or(p.role_flag = %s))
        ) incomplete
        """,
        (RESOURCE, LIGAND_FLAG, RECEPTOR_FLAG),
    )[0]

    assert incomplete == 0, (
        f'{incomplete} interactions of this release flag only one of their two '
        f'roles, so the other end of the pair has no recorded role at all '
        f'(over {record_rows} record rows)'
    )
