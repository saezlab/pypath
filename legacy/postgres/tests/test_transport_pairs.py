"""A transported metabolite and its transporter are one binary interaction.

The build has identified transport since the loaders did: ``Transport:OM:0035``
is an entity type of its own, and Recon3D and Human-GEM each publish a
transport dataset. None of that reached a caller. A transport is stored as a
**star** — one relation per member, whose subject is the event — so selecting
it yields metabolite-to-event rows, which are not metabolite-protein
interactions whatever class they carry. Measured on the build before this
projection existed: recon3d and metatlas contributed no ``transport`` relation
at all, and rhea contributed 12.

This file asserts the pair that owes. A metabolite stated as a reactant in one
compartment and as a product in another, inside a reaction some protein
catalyses, is that protein transporting that metabolite across that boundary,
and it projects to **one** ordinary binary record naming both ends and both
compartments.

Two things the projection has to get right beyond producing a row:

* **The two compartments have to differ.** They are what a transport *is*, and
  an earlier projection folded a member's compartment across its roles before
  it split them, so both sides of every transported metabolite read the same
  value. A row asserting a move from the cytosol to the cytosol is worse than
  no row.
* **Both roles is not enough.** A cofactor consumed and regenerated in one
  compartment holds the same pair of roles and crosses nothing.

Run::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_transport_pairs.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os

import pytest

from tests.fixtures.interaction_graph import ENTITY, build_interaction_fixture

DATABASE_URL = os.environ.get('DATABASE_URL')
SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_TRANSPORT_PAIRS',
    'transport_pairs_test',
)

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the transport projection needs a Postgres',
)


@pytest.fixture(scope='module')
def conn():
    """A connection to the database the scratch schema is built in."""
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture(scope='module')
def scratch(conn):
    """The fixture graph, projected, in a throwaway namespace."""
    from omnipath_build.db import schema as build_schema
    from omnipath_build.db.derived_tables import rebuild_interaction_tables

    build_schema.ensure_schema(conn, schema=SCRATCH, drop_existing=True)
    conn.commit()
    build_interaction_fixture(conn, SCRATCH)
    rebuild_interaction_tables(conn, schema=SCRATCH)
    conn.commit()
    try:
        yield SCRATCH
    finally:
        with conn.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {SCRATCH} CASCADE')
        conn.commit()


def _rows(conn, statement: str, parameters: list | None = None):
    with conn.cursor() as cur:
        cur.execute(statement, parameters or [])
        return cur.fetchall()


def _records_between(conn, schema: str, subject: str, obj: str) -> list[tuple]:
    """Every record row for one ordered pair, whatever class it carries."""
    return _rows(
        conn,
        f"""
        SELECT vic.name, ds.name, f.is_directed, f.attributes, f.interaction_id
        FROM {schema}.interaction_fact_resource f
        JOIN {schema}.vocab_interaction_class vic
          ON vic.interaction_class_id = f.interaction_class_id
        JOIN {schema}.data_source ds ON ds.source_id = f.source_id
        WHERE f.subject_entity_id = %s AND f.object_entity_id = %s
        ORDER BY vic.name, ds.name
        """,
        [ENTITY[subject], ENTITY[obj]],
    )


class TestTheTransportBecomesAPair:
    """The fixture's transport, read as the pair it is.

    Two of one cargo leave the cytosol, one arrives outside, and ``enz_tb`` is
    the protein that moves it.
    """

    @pytest.fixture(scope='class')
    def record(self, conn, scratch):
        found = _records_between(conn, scratch, 'enz_tb', 'met_cargo')
        assert len(found) == 1, (
            'the transport projects to one record row per resource; '
            f'found {len(found)}'
        )
        return found[0]

    def test_the_pair_exists_at_all(self, record):
        """Before this projection the transporter reached no metabolite.

        The star joins the cargo to the **event**, so every pair reading of it
        named an abstract node on one end.
        """
        assert record is not None

    def test_the_class_is_transport(self, record):
        """The class is stated here, not inherited from the star.

        The star's own class resolves to ``other``, which is a gap in the
        classification map. The compartment change is the transport statement,
        so the pair says so rather than inheriting the gap.
        """
        interaction_class, _source, _directed, _attributes, _header = record
        assert interaction_class == 'transport'

    def test_the_resource_is_the_one_that_stated_the_movement(self, record):
        _class, source, _directed, _attributes, _header = record
        assert source == 'fixture_res_a'

    def test_the_pair_is_directed(self, record):
        """A transporter moving a cargo is asymmetric.

        Every transport row the build already holds carries the flag.
        """
        _class, _source, is_directed, _attributes, _header = record
        assert is_directed is True

    def test_both_compartments_are_recorded(self, record):
        _class, _source, _directed, attributes, _header = record
        assert attributes == {'transport': [{'from': 'c', 'to': 'e'}]}

    def test_the_compartments_differ(self, record):
        """The movement is the whole content of a transport.

        A projection that resolved the compartment per member rather than per
        role reads `min('c', 'e')` on both sides and states that the cargo
        went from the cytosol to the cytosol.
        """
        _class, _source, _directed, attributes, _header = record
        movement = attributes['transport'][0]
        assert movement['from'] != movement['to']

    def test_the_record_reaches_a_header_of_its_own(
        self,
        conn,
        scratch,
        record,
    ):
        """An ordinary pair reaches the ordinary header.

        Two participants, and the event entity is not one of them.
        """
        _class, _source, _directed, _attributes, header = record
        assert header is not None
        arity, sources = _rows(
            conn,
            f'SELECT arity, sources FROM {scratch}.interaction '
            'WHERE interaction_id = %s',
            [header],
        )[0]
        assert arity == 2
        assert sources == ['fixture_res_a']
        participants = _rows(
            conn,
            f'SELECT entity_id::text FROM {scratch}.interaction_party '
            'WHERE interaction_id = %s ORDER BY entity_id',
            [header],
        )
        assert sorted(row[0] for row in participants) == sorted(
            [ENTITY['enz_tb'], ENTITY['met_cargo']]
        )

    def test_the_star_rows_are_untouched(self, conn, scratch):
        """The pair is added beside the reaction reading, not instead of it.

        The record keeps one row per `(event, member)` and the catalyst keeps
        its `controls` edge, so a caller reading the transport as a reaction
        still finds every participant of it.
        """
        spokes = _rows(
            conn,
            f"""
            SELECT count(*)
            FROM {scratch}.interaction_fact_resource
            WHERE subject_entity_id = %s
            """,
            [ENTITY['trn_b']],
        )[0][0]
        assert spokes == 2
        catalysis = _records_between(conn, scratch, 'enz_tb', 'trn_b')
        assert len(catalysis) == 1


class TestWhatIsNotATransport:
    """Three shapes that hold some of a transport's parts and not the move."""

    def test_a_member_with_one_role_moves_nowhere(self, conn, scratch):
        """The fuel of the same transport is consumed and never regenerated.

        It sits on one side of the arrow, so there is no second compartment to
        compare and nothing was carried across anything.
        """
        assert _records_between(conn, scratch, 'enz_tb', 'met_fuel') == []

    def test_two_roles_in_one_compartment_move_nothing(self, conn, scratch):
        """A cofactor turned over in the cytosol holds both roles.

        This is the shape that makes the compartment change load-bearing
        rather than decorative: read the roles alone and this reaction states
        a transporter that transports nothing.
        """
        assert _records_between(conn, scratch, 'enz_still', 'met_still') == []

    def test_two_sides_that_are_two_metabolites_are_not_a_move(
        self,
        conn,
        scratch,
    ):
        """The other transport in the fixture names one metabolite per side.

        That is how Rhea publishes a transport — one membership for the inside
        and another for the outside — and the two are different entities, so
        nothing here says the same molecule appeared in both places. The
        reaction reading keeps it; the pair reading cannot invent it.
        """
        assert _records_between(conn, scratch, 'enz_t', 'met_in') == []
        assert _records_between(conn, scratch, 'enz_t', 'met_out') == []


class TestTheRestOfTheRecordIsUnchanged:
    """The projection adds rows. It must not annotate the ones already there."""

    def test_only_the_transport_pairs_carry_attributes(self, conn, scratch):
        annotated = _rows(
            conn,
            f"""
            SELECT vic.name, count(*)
            FROM {scratch}.interaction_fact_resource f
            JOIN {scratch}.vocab_interaction_class vic
              ON vic.interaction_class_id = f.interaction_class_id
            WHERE f.attributes IS NOT NULL
            GROUP BY 1
            """,
        )
        assert annotated == [('transport', 1)]

    def test_the_directly_stated_transport_keeps_its_shape(
        self,
        conn,
        scratch,
    ):
        """A resource publishing a ``transports`` verb already made a pair.

        That pair gains no compartments the resource never stated.
        """
        found = _records_between(conn, scratch, 'o', 'p')
        assert len(found) == 1
        interaction_class, _source, is_directed, attributes, _header = found[0]
        assert interaction_class == 'transport'
        assert is_directed is True
        assert attributes is None
