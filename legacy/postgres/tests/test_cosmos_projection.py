"""A reaction becomes substrate and product edges around a per-reaction node.

The causal-reasoning prior knowledge network is a **binary** graph: every edge
is one metabolite and one gene, and the reaction itself survives only as the
node the two sides meet at. Projecting the build's N-ary reaction headers into
it is therefore a rewrite, not a flattening, and this file asserts the four
things the rewrite has to get right.

* **The enzyme node is per reaction.** The same protein catalysing five
  reactions is five nodes, ``Gene1__P00533`` through ``Gene5__P00533``, because
  a reasoning method treats each of them as an independently switchable step.
  The index is a modelling artifact, which is why it lives in its own table and
  never in ``entity``, and it has to be **stable across rebuilds** — a label
  that moved between runs would make every downstream comparison meaningless.
* **Only the two sides of the arrow project.** A substrate reaches the enzyme
  node and the enzyme node reaches a product. Nothing joins a substrate to a
  product, nothing joins the enzyme to itself, and a party row holding a role
  that names no side of the arrow — a cofactor, a regulator — is left alone
  and counted, because a genome-scale model has no such role and reading one as
  a substrate would state something the resource did not.
* **A reaction nobody names an enzyme for keeps its edges.** It is the majority
  shape in a metabolic model. The reaction stands in for the unknown catalyst,
  the row says so in ``orphan``, and no reaction is dropped for want of a
  catalyst.
* **The labels reach the namespaces the network reads, or say that they did
  not.** COSMOS is written against UniProt on the gene side and ChEBI on the
  metabolite side, and the build hands the projection whatever identifier its
  resources agreed on. A metabolite typed ChEBI gains the prefix the build
  does not store and the network does; a catalyst a mapping answers for
  becomes an accession; a catalyst nothing answers for keeps its own
  identifier, and both endpoints of every edge record the namespace their
  label actually carries, so a fallen-back node cannot be mistaken for a
  translated one. A run with no mapping database at all falls every label
  back and still emits the whole network.
* **A reversible reaction runs both ways, and only a reversible one does.** The
  resource states the direction on the reaction event, under two spellings per
  value. A reversible conversion gains the mirrored pair — every product into
  the enzyme, every reactant out of it — on a ``_rev`` copy of the gene node
  alone, never of the metabolite. A one-way conversion gains nothing. A
  conversion nobody published a direction for gains nothing either, and says
  so with a null rather than by looking one-way, because a source that is
  silent and a source that says left to right are different claims.

Run::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_cosmos_projection.py -v

Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os
import re

import pytest

from tests.fixtures.interaction_graph import (
    ENTITY,
    SOURCE_NAMES,
    build_interaction_fixture,
)

# The projection serves a named set of metabolic resources. The fixture
# publishes under its own names, so the tests hand it those rather than
# renaming the fixture after real resources and reading as though the build
# held them.
FIXTURE_SOURCES = tuple(SOURCE_NAMES.values())

DATABASE_URL = os.environ.get('DATABASE_URL')
SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_COSMOS',
    'cosmos_projection_test',
)
# The translated projection gets a namespace of its own rather than a second
# run over the first: the two differ in every label the mappings touch, and one
# table holding whichever of them ran last would make every assertion in this
# file depend on the order pytest happened to build its fixtures in.
TRANSLATION_SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_COSMOS_TRANSLATION',
    'cosmos_translation_test',
)
# The mapping database, staged. The real one lives on another server, holds 618
# million rows and changes when the utils build runs, so a test that asked it
# would be asserting today's coverage rather than the rule. This one is four
# mappings and one resolver row, in the shape the utils schema publishes, and
# the step reads it through exactly the queries it uses against the real thing.
MAPPINGS_SCRATCH = os.environ.get(
    'OMNIPATH_TEST_SCRATCH_SCHEMA_COSMOS_MAPPINGS',
    'cosmos_mappings_test',
)

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the binary projection needs a Postgres',
)

# The node strings the fixture's ordinary conversion has to produce. They are
# written out rather than composed, because a helper that built them from the
# same rule the projection uses would agree with it however wrong the rule is.
MET_A = 'Metab__FIXTURE_cos_met_a_c'
MET_B = 'Metab__FIXTURE_cos_met_b'
MET_C = 'Metab__FIXTURE_cos_met_c_m'
ENZYME_ID = 'FIXTURE_cos_enz'
ORPHAN_REACTION_ID = 'FIXTURE_cos_rxn_orphan'

# The node strings of the reversible conversion, both sides of the mirror.
BACK_IN_A = 'Metab__FIXTURE_cos_back_in_a_c'
BACK_IN_B = 'Metab__FIXTURE_cos_back_in_b_c'
BACK_OUT = 'Metab__FIXTURE_cos_back_out_c'
BACK_ENZYME_ID = 'FIXTURE_cos_enz_back'
ORPHAN_BACK_REACTION_ID = 'FIXTURE_cos_rxn_orphan_back'

GENE_NODE = re.compile(r'^Gene(\d+)__(.+)$')
REVERSE_SUFFIX = '_rev'
ORPHAN_NODE_PREFIX = 'orphanReac'

COLUMNS = (
    'source_label',
    'target_label',
    'source_entity_id',
    'target_entity_id',
    'source_type',
    'target_type',
    'source_compartment',
    'target_compartment',
    'mor',
    'interaction_type',
    'reaction_entity_id',
    'interaction_id',
    'reaction_index',
    'orphan',
    'reverse',
    'direction',
    'source_id_type',
    'target_id_type',
    'sources',
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
        # A failed projection leaves the transaction aborted, and the drop
        # would then leak the schema into the next run.
        conn.rollback()
        with conn.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {SCRATCH} CASCADE')
        conn.commit()


@pytest.fixture(scope='module')
def projection(conn, scratch):
    """The binary edges, built once, and what the step reported building."""
    from omnipath_build.cosmos import (
        build_cosmos_projection,
        ensure_cosmos_edge_table,
    )

    ensure_cosmos_edge_table(conn, schema=scratch)
    conn.commit()
    try:
        stats = build_cosmos_projection(
            conn, schema=scratch, reaction_sources=FIXTURE_SOURCES,
        )
    except Exception:
        conn.rollback()
        raise
    conn.commit()
    return stats


def _rows(conn, statement: str, parameters: list | None = None):
    with conn.cursor() as cur:
        cur.execute(statement, parameters or [])
        return cur.fetchall()


def _edges(
    conn,
    schema: str,
    where: str = 'true',
    parameters: list | None = None,
) -> list[dict]:
    """The projected edges matching ``where``, as column-named rows."""
    from omnipath_build.cosmos import TABLE

    rows = _rows(
        conn,
        f"""
        SELECT {', '.join(COLUMNS)}
        FROM {schema}.{TABLE}
        WHERE {where}
        """,
        parameters,
    )
    return [dict(zip(COLUMNS, row, strict=True)) for row in rows]


def _reaction_edges(conn, schema: str, name: str) -> list[dict]:
    """The edges the named fixture reaction produced, connectors aside."""
    return _edges(
        conn,
        schema,
        "reaction_entity_id = %s AND interaction_type <> 'connector'",
        [ENTITY[name]],
    )


def _labels(edges: list[dict]) -> set[tuple[str, str]]:
    return {(edge['source_label'], edge['target_label']) for edge in edges}


def _gene_nodes(edges: list[dict]) -> set[str]:
    """Every enzyme node the given edges touch, on either end."""
    return {
        label
        for edge in edges
        for label in (edge['source_label'], edge['target_label'])
        if GENE_NODE.match(label)
    }


def _one_gene_node(edges: list[dict]) -> str:
    nodes = _gene_nodes(edges)
    assert len(nodes) == 1, (
        f'these edges should meet at one enzyme node; found {sorted(nodes)}'
    )
    return nodes.pop()


class TestTheEdgesAreTheSubstrateAndProductPairs:
    """The ordinary conversion: two substrates, one product, one enzyme."""

    @pytest.fixture(scope='class')
    def edges(self, conn, scratch, projection):
        found = _reaction_edges(conn, scratch, 'cos_rxn')
        assert found, 'the catalysed conversion produced no edges at all'
        return found

    @pytest.fixture(scope='class')
    def gene_label(self, edges):
        return _one_gene_node(edges)

    def test_one_edge_per_substrate_and_one_per_product(self, edges):
        """Two substrates in, one product out, and nothing besides."""
        assert len(edges) == 3

    def test_every_substrate_points_at_the_enzyme(self, edges, gene_label):
        incoming = {
            edge['source_label']
            for edge in edges
            if edge['target_label'] == gene_label
        }
        assert incoming == {MET_A, MET_B}

    def test_the_enzyme_points_at_the_product(self, edges, gene_label):
        outgoing = {
            edge['target_label']
            for edge in edges
            if edge['source_label'] == gene_label
        }
        assert outgoing == {MET_C}

    def test_no_substrate_reaches_a_product_directly(self, edges):
        """The enzyme node is the only path across the arrow.

        A substrate-to-product edge would say the conversion happens whether or
        not the gene is expressed, which is the one thing the network exists to
        deny.
        """
        assert (MET_A, MET_C) not in _labels(edges)
        assert (MET_B, MET_C) not in _labels(edges)

    def test_the_enzyme_does_not_point_at_itself(self, edges, gene_label):
        assert (gene_label, gene_label) not in _labels(edges)

    def test_the_edges_carry_no_sign(self, edges):
        """The projection knows the direction of the arrow, not its sign."""
        assert {edge['mor'] for edge in edges} == {1}

    def test_the_edges_are_catalysis(self, edges):
        assert {edge['interaction_type'] for edge in edges} == {'catalysis'}

    def test_the_two_sides_keep_their_entity_and_their_type(
        self,
        edges,
        gene_label,
    ):
        """A metabolite side resolves to an entity, and says it is one."""
        metabolite_sides = [
            (edge['source_entity_id'], edge['source_type'])
            for edge in edges
            if edge['source_label'] != gene_label
        ]
        assert metabolite_sides
        assert {kind for _entity_id, kind in metabolite_sides} == {'metabolite'}
        assert all(entity_id for entity_id, _kind in metabolite_sides)

    def test_the_edges_credit_the_resource_that_published_them(self, edges):
        for edge in edges:
            assert 'fixture_res_a' in edge['sources']


class TestTheNodeStringsAreTheOnesTheNetworkReads:
    """A consumer matches these strings literally. They are the contract."""

    @pytest.fixture(scope='class')
    def edges(self, conn, scratch, projection):
        return _reaction_edges(conn, scratch, 'cos_rxn')

    def test_a_metabolite_states_the_compartment_it_sits_in(self, edges):
        assert MET_A in {edge['source_label'] for edge in edges}
        assert MET_C in {edge['target_label'] for edge in edges}

    def test_a_metabolite_with_no_compartment_carries_no_separator(
        self,
        edges,
    ):
        """``Metab__x_`` is a different node from ``Metab__x``.

        The compartment is joined on with an underscore, so a metabolite whose
        location nobody stated has to end at the identifier rather than at the
        separator that would have introduced one.
        """
        labels = {edge['source_label'] for edge in edges}
        assert MET_B in labels
        assert 'Metab__FIXTURE_cos_met_b_' not in labels

    def test_the_metabolite_prefix_is_five_letters(
        self,
        conn,
        scratch,
        projection,
    ):
        """``Metab__``, not ``Metabo__``.

        The prose contract says the longer one. The pinned consumer tests say
        the shorter, and the consumer is what the strings are for.
        """
        from omnipath_build.cosmos import TABLE

        wrong = _rows(
            conn,
            f"""
            SELECT count(*)
            FROM {scratch}.{TABLE}
            WHERE left(source_label, 8) = 'Metabo__'
               OR left(target_label, 8) = 'Metabo__'
            """,
        )
        assert wrong[0][0] == 0

        right = _rows(
            conn,
            f"""
            SELECT count(*)
            FROM {scratch}.{TABLE}
            WHERE left(source_label, 7) = 'Metab__'
               OR left(target_label, 7) = 'Metab__'
            """,
        )
        assert right[0][0] > 0

    def test_the_enzyme_node_names_the_reaction_it_belongs_to(self, edges):
        """``Gene{N}__<id>``, with the N the row reports in its own column."""
        for edge in edges:
            for side in ('source', 'target'):
                match = GENE_NODE.match(edge[f'{side}_label'])
                if match is None:
                    continue
                assert edge['reaction_index'] is not None
                assert edge[f'{side}_label'] == (
                    f'Gene{edge["reaction_index"]}__{ENZYME_ID}'
                )

    def test_the_suffix_and_the_column_say_the_same_thing(
        self,
        conn,
        scratch,
        projection,
    ):
        """A ``_rev`` label and a ``reverse`` row are one claim, not two.

        A consumer reads the direction off the label and the build reads it off
        the column. A row where the two disagree would be a network whose
        answer depends on which of them the reader trusted.
        """
        from omnipath_build.cosmos import TABLE

        disagreeing = _rows(
            conn,
            f"""
            SELECT source_label, target_label, reverse
            FROM {scratch}.{TABLE}
            WHERE (
                right(source_label, 4) = '_rev'
                OR right(target_label, 4) = '_rev'
              ) <> reverse
              AND interaction_type <> 'connector'
            """,
        )
        assert disagreeing == []

    def test_only_the_enzyme_side_of_a_reverse_row_is_suffixed(
        self,
        conn,
        scratch,
        projection,
    ):
        """The mirrored metabolite is the same metabolite.

        Suffixing it would mint a second node for one molecule, and the
        forward and reverse halves of a reversible reaction would stop meeting
        anywhere in the graph.
        """
        from omnipath_build.cosmos import TABLE

        suffixed = _rows(
            conn,
            f"""
            SELECT source_label, target_label
            FROM {scratch}.{TABLE}
            WHERE (
                left(source_label, 7) = 'Metab__'
                AND right(source_label, 4) = '_rev'
              )
               OR (
                left(target_label, 7) = 'Metab__'
                AND right(target_label, 4) = '_rev'
              )
            """,
        )
        assert suffixed == []


class TestTheCompartmentIsAColumnOfItsOwn:
    """Reading it back out of the label is parsing, not querying."""

    @pytest.fixture(scope='class')
    def edges(self, conn, scratch, projection):
        return _reaction_edges(conn, scratch, 'cos_rxn')

    def test_the_substrate_compartment_is_readable_without_the_label(
        self,
        edges,
    ):
        compartments = {
            edge['source_label']: edge['source_compartment']
            for edge in edges
        }
        assert compartments[MET_A] == 'c'

    def test_the_product_compartment_is_readable_without_the_label(
        self,
        edges,
    ):
        compartments = {
            edge['target_label']: edge['target_compartment']
            for edge in edges
        }
        assert compartments[MET_C] == 'm'

    def test_an_unstated_compartment_stays_unstated(self, edges):
        """NULL says nobody said. An empty string says somebody said nothing."""
        compartments = {
            edge['source_label']: edge['source_compartment']
            for edge in edges
        }
        assert compartments[MET_B] is None

    def test_the_enzyme_node_claims_no_compartment(self, edges):
        for edge in edges:
            for side in ('source', 'target'):
                if GENE_NODE.match(edge[f'{side}_label']):
                    assert edge[f'{side}_compartment'] is None


class TestAReactionWithNoNamedEnzymeKeepsItsEdges:
    """Most reactions in a metabolic model have no catalyst on record."""

    @pytest.fixture(scope='class')
    def edges(self, conn, scratch, projection):
        found = _reaction_edges(conn, scratch, 'cos_rxn_orphan')
        assert found, (
            'the conversion with no named enzyme was dropped; it is the '
            'majority case, not an edge case'
        )
        return found

    def test_both_sides_of_the_arrow_survive(self, edges):
        assert len(edges) == 2

    def test_the_reaction_stands_in_for_the_unknown_catalyst(self, edges):
        labels = {
            label
            for edge in edges
            for label in (edge['source_label'], edge['target_label'])
            if GENE_NODE.match(label)
        }
        assert len(labels) == 1
        label = labels.pop()
        index = {edge['reaction_index'] for edge in edges}
        assert len(index) == 1
        assert label == (
            f'Gene{index.pop()}__orphanReac{ORPHAN_REACTION_ID}'
        )

    def test_the_missing_enzyme_is_stated_in_the_data(self, edges):
        """A consumer has to be able to tell a real gene from a placeholder."""
        assert {edge['orphan'] for edge in edges} == {True}

    def test_the_placeholder_resolves_to_no_entity(self, edges):
        """It is not a protein, so it holds no protein's id."""
        for edge in edges:
            for side in ('source', 'target'):
                if GENE_NODE.match(edge[f'{side}_label']):
                    assert edge[f'{side}_entity_id'] is None

    def test_a_catalysed_reaction_is_not_marked_orphan(
        self,
        conn,
        scratch,
        projection,
    ):
        catalysed = _reaction_edges(conn, scratch, 'cos_rxn')
        assert {edge['orphan'] for edge in catalysed} == {False}


class TestAPartyOnNeitherSideOfTheArrowIsNotProjected:
    """A cofactor and a regulator hold roles a binarised model has no word for.

    Every species in a genome-scale model is a reactant or a product by the
    sign of its coefficient. Reading a regulator as a substrate would put a
    statement into the network that no resource made, so the step leaves those
    party rows alone and reports how many it left.
    """

    @pytest.fixture(scope='class')
    def edges(self, conn, scratch, projection):
        return _reaction_edges(conn, scratch, 'cos_rxn_helper')

    def test_only_the_substrate_and_the_product_project(self, edges):
        assert len(edges) == 2

    def test_the_cofactor_is_not_read_as_a_substrate(
        self,
        conn,
        scratch,
        projection,
    ):
        from omnipath_build.cosmos import TABLE

        rows = _rows(
            conn,
            f"""
            SELECT count(*)
            FROM {scratch}.{TABLE}
            WHERE source_label LIKE %s OR target_label LIKE %s
            """,
            ['%FIXTURE_cos_cofactor%'] * 2,
        )
        assert rows[0][0] == 0

    def test_the_regulator_is_not_read_as_an_enzyme(
        self,
        conn,
        scratch,
        projection,
    ):
        from omnipath_build.cosmos import TABLE

        rows = _rows(
            conn,
            f"""
            SELECT count(*)
            FROM {scratch}.{TABLE}
            WHERE source_label LIKE %s OR target_label LIKE %s
            """,
            ['%FIXTURE_cos_regulator%'] * 2,
        )
        assert rows[0][0] == 0

    def test_the_step_reports_what_it_left_out(self, projection):
        """A silent skip is indistinguishable from a projection that failed."""
        assert projection.skipped_parties >= 2


class TestAReversibleReactionRunsBothWays:
    """A reversible conversion is two reactions sharing one enzyme.

    The resource states the direction on the reaction event. Where it says
    reversible, the network needs the mirrored pair as well, on a node of its
    own: a reasoning method has to be able to switch the two directions
    independently, and one node carrying both would tie them together.
    """

    @pytest.fixture(scope='class')
    def edges(self, conn, scratch, projection):
        found = _reaction_edges(conn, scratch, 'cos_rxn_back')
        assert found, 'the reversible conversion produced no edges at all'
        return found

    @pytest.fixture(scope='class')
    def forward(self, edges):
        return [edge for edge in edges if not edge['reverse']]

    @pytest.fixture(scope='class')
    def backward(self, edges):
        found = [edge for edge in edges if edge['reverse']]
        assert found, (
            'the conversion is published as reversible and gained no '
            'mirrored edges'
        )
        return found

    def test_the_forward_pair_is_what_it_was(self, forward):
        gene = _one_gene_node(forward)
        assert _labels(forward) == {
            (BACK_IN_A, gene),
            (BACK_IN_B, gene),
            (gene, BACK_OUT),
        }

    def test_the_product_gains_an_edge_into_the_enzyme(self, backward):
        gene = _one_gene_node(backward)
        assert (BACK_OUT, gene) in _labels(backward)

    def test_every_reactant_gains_an_edge_out_of_the_enzyme(self, backward):
        gene = _one_gene_node(backward)
        assert (gene, BACK_IN_A) in _labels(backward)
        assert (gene, BACK_IN_B) in _labels(backward)

    def test_the_mirror_holds_nothing_else(self, backward):
        """Run backwards, a product is a substrate and a substrate a product.

        The mirrored side is the forward side with the arrow turned round, not
        a second chance to join every metabolite to every other.
        """
        assert len(backward) == 3

    def test_the_mirrored_node_is_the_forward_node_suffixed(
        self,
        forward,
        backward,
    ):
        assert _one_gene_node(backward) == (
            _one_gene_node(forward) + REVERSE_SUFFIX
        )

    def test_the_mirrored_node_still_names_its_enzyme(self, backward):
        edge = backward[0]
        assert _one_gene_node(backward) == (
            f'Gene{edge["reaction_index"]}__{BACK_ENZYME_ID}{REVERSE_SUFFIX}'
        )

    def test_the_metabolites_are_the_same_molecules(self, backward):
        """Literally the same strings, so the two halves meet in the graph."""
        labels = {
            label
            for edge in backward
            for label in (edge['source_label'], edge['target_label'])
            if not GENE_NODE.match(label)
        }
        assert labels == {BACK_IN_A, BACK_IN_B, BACK_OUT}

    def test_the_forward_rows_are_not_marked_mirrored(self, forward):
        assert {edge['reverse'] for edge in forward} == {False}

    def test_the_direction_is_on_every_row(self, edges):
        assert {edge['direction'] for edge in edges} == {'reversible'}

    def test_both_halves_belong_to_one_reaction(self, edges):
        """The index is what says the two directions are one conversion."""
        assert len({edge['reaction_index'] for edge in edges}) == 1


class TestAOneWayReactionGainsNothing:
    """A resource that says left to right has said which way it runs."""

    @pytest.fixture(scope='class')
    def edges(self, conn, scratch, projection):
        found = _reaction_edges(conn, scratch, 'cos_rxn_fwd')
        assert found, 'the one-way conversion produced no edges at all'
        return found

    def test_it_keeps_only_the_forward_pair(self, edges):
        assert len(edges) == 2
        assert {edge['reverse'] for edge in edges} == {False}

    def test_nothing_it_holds_is_suffixed(self, edges):
        for edge in edges:
            assert not edge['source_label'].endswith(REVERSE_SUFFIX)
            assert not edge['target_label'].endswith(REVERSE_SUFFIX)

    def test_the_direction_says_which_way(self, edges):
        assert {edge['direction'] for edge in edges} == {'left_to_right'}


class TestAnUnpublishedDirectionIsNotAOneWayReaction:
    """Silence and a stated direction are different claims about a reaction.

    One resource states no direction on any of its reactions. Recording that as
    left to right would turn a gap in what the source says into a statement the
    source never made, and nothing downstream could tell the two apart again.
    """

    @pytest.fixture(scope='class')
    def edges(self, conn, scratch, projection):
        return _reaction_edges(conn, scratch, 'cos_rxn')

    def test_nothing_is_mirrored(self, edges):
        assert {edge['reverse'] for edge in edges} == {False}

    def test_the_direction_is_null(self, edges):
        assert {edge['direction'] for edge in edges} == {None}

    def test_a_one_way_reaction_is_still_told_apart_from_a_silent_one(
        self,
        conn,
        scratch,
        edges,
    ):
        one_way = _reaction_edges(conn, scratch, 'cos_rxn_fwd')
        assert {edge['direction'] for edge in one_way} == {'left_to_right'}
        assert {edge['direction'] for edge in edges} == {None}


class TestTheTwoSpellingsOfADirectionAgree:
    """Each direction reaches the build in upper and in lower case.

    The upper-case forms carry most of the build and the lower-case forms the
    rest, so reading one spelling and not the other would silently drop a
    resource rather than fail.
    """

    @pytest.fixture(scope='class')
    def by_reaction(self, conn, scratch, projection):
        return {
            name: _reaction_edges(conn, scratch, name)
            for name in (
                'cos_rxn_back',
                'cos_rxn_back_lower',
                'cos_rxn_fwd',
                'cos_rxn_fwd_lower',
            )
        }

    def test_both_spellings_of_reversible_read_the_same(self, by_reaction):
        for name in ('cos_rxn_back', 'cos_rxn_back_lower'):
            edges = by_reaction[name]
            assert edges, f'{name} produced no edges'
            assert {edge['direction'] for edge in edges} == {'reversible'}
            assert any(edge['reverse'] for edge in edges), (
                f'{name} is published as reversible and gained no mirror'
            )

    def test_both_spellings_of_one_way_read_the_same(self, by_reaction):
        for name in ('cos_rxn_fwd', 'cos_rxn_fwd_lower'):
            edges = by_reaction[name]
            assert edges, f'{name} produced no edges'
            assert {edge['direction'] for edge in edges} == {'left_to_right'}
            assert not any(edge['reverse'] for edge in edges)


class TestAReversibleReactionWithNoNamedEnzyme:
    """The two situations meet: a placeholder node that also mirrors."""

    @pytest.fixture(scope='class')
    def edges(self, conn, scratch, projection):
        found = _reaction_edges(conn, scratch, 'cos_rxn_orphan_back')
        assert found, (
            'a reversible conversion with no named enzyme produced no edges'
        )
        return found

    def test_the_mirror_is_there_as_well(self, edges):
        assert len(edges) == 4
        assert len([edge for edge in edges if edge['reverse']]) == 2

    def test_the_placeholder_takes_the_suffix(self, edges):
        backward = [edge for edge in edges if edge['reverse']]
        index = {edge['reaction_index'] for edge in backward}
        assert len(index) == 1
        assert _one_gene_node(backward) == (
            f'Gene{index.pop()}__orphanReac{ORPHAN_BACK_REACTION_ID}'
            f'{REVERSE_SUFFIX}'
        )

    def test_it_is_still_an_orphan_on_both_sides(self, edges):
        assert {edge['orphan'] for edge in edges} == {True}

    def test_the_direction_is_recorded(self, edges):
        assert {edge['direction'] for edge in edges} == {'reversible'}


class TestEveryEnzymeNodeIsReachableFromItsBareIdentifier:
    """The per-reaction index is the network's, not the caller's.

    A caller holding a UniProt accession cannot know which N the projection
    gave it, so each enzyme node carries a connector from the bare identifier
    to the node. Without it the gene side of the network is unaddressable.
    """

    @pytest.fixture(scope='class')
    def connectors(self, conn, scratch, projection):
        return _edges(conn, scratch, "interaction_type = 'connector'")

    @pytest.fixture(scope='class')
    def gene_nodes(self, conn, scratch, projection):
        return _gene_nodes(
            _edges(conn, scratch, "interaction_type <> 'connector'"),
        )

    @pytest.fixture(scope='class')
    def forward_gene_nodes(self, conn, scratch, projection):
        return _gene_nodes(
            _edges(
                conn,
                scratch,
                "interaction_type <> 'connector' AND NOT reverse",
            ),
        )

    def test_the_connector_joins_the_bare_identifier_to_the_node(
        self,
        connectors,
    ):
        """The identifier a caller holds carries no index and no suffix."""
        for edge in connectors:
            match = GENE_NODE.match(edge['target_label'])
            assert match is not None, (
                f'a connector points at {edge["target_label"]}, which is not '
                'an enzyme node'
            )
            bare = match.group(2)
            if bare.endswith(REVERSE_SUFFIX):
                bare = bare[: -len(REVERSE_SUFFIX)]
            # A reaction standing in for an unknown enzyme is reachable by the
            # identifier the resource published, not by the node string the
            # projection built around it.
            if bare.startswith(ORPHAN_NODE_PREFIX):
                bare = bare[len(ORPHAN_NODE_PREFIX):]
            assert edge['source_label'] == bare

    def test_every_forward_enzyme_node_has_one(
        self,
        connectors,
        forward_gene_nodes,
        gene_nodes,
    ):
        """Whether the mirrored node gets its own connector is the step's call.

        What is not optional is that every enzyme node the forward graph holds
        is reachable, and that a connector never points at something that is
        not an enzyme node at all.
        """
        assert forward_gene_nodes
        targets = {edge['target_label'] for edge in connectors}
        assert forward_gene_nodes <= targets
        assert targets <= gene_nodes

    def test_the_named_enzyme_reaches_its_node(self, connectors):
        joined = {
            (edge['source_label'], edge['target_label'])
            for edge in connectors
        }
        assert any(
            source == ENZYME_ID and target.endswith(f'__{ENZYME_ID}')
            for source, target in joined
        )


class TestASecondRunLandsOnTheSameGraph:
    """The index is embedded in every enzyme label, so it has to be stable.

    Ranking the reactions by their entity id gives the same N on every rebuild.
    Numbering them by the order the rows happened to arrive in does not, and a
    network whose node names move between builds cannot be compared with
    anything, including itself.
    """

    @pytest.fixture(scope='class')
    def rebuilt(self, conn, scratch, projection):
        from omnipath_build.cosmos import build_cosmos_projection

        before = _edges(conn, scratch)
        build_cosmos_projection(
            conn, schema=scratch, reaction_sources=FIXTURE_SOURCES,
        )
        conn.commit()
        return before, _edges(conn, scratch)

    def test_the_rows_are_not_written_twice(self, rebuilt):
        before, after = rebuilt
        assert len(after) == len(before)

    def test_every_label_comes_back_the_same(self, rebuilt):
        before, after = rebuilt
        assert _labels(after) == _labels(before)

    def test_the_reaction_index_does_not_move(self, rebuilt):
        before, after = rebuilt
        assert sorted(
            (edge['source_label'], edge['target_label'], edge['reaction_index'])
            for edge in after
        ) == sorted(
            (edge['source_label'], edge['target_label'], edge['reaction_index'])
            for edge in before
        )


def test_the_edges_live_in_their_own_table():
    """The enzyme nodes are modelling artifacts and stay out of the entities.

    One protein in five hundred reactions is five hundred nodes here. Writing
    those into the canonical entity table would move every count the rest of
    the build reports, so the projection owns a table instead.
    """
    from omnipath_build.cosmos import TABLE

    assert TABLE == 'cosmos_edge'


# The identifier types the staged mapping database knows, numbered as the utils
# build numbers them. The step reads the numbers back by name rather than
# assuming them, and these are here so that it has something to read.
STAGED_ID_TYPES = (
    (1, 'uniprot'),
    (18, 'entrez'),
    (44, 'chebi'),
    (52, 'inchikey'),
)

# The mappings themselves, in the utils table's own shape: a source type, a
# target type, the taxonomy the pair holds under, and the two identifiers.
#
# Two of them are deliberately ambiguous, and the two ambiguities resolve in
# opposite directions. One Entrez gene answering two accessions takes the lower
# one; one InChIKey two ChEBI entries claim keeps the InChIKey.
#
# The ChEBI side is stored **prefixed**, which is how the utils build stores it
# and is not how the OmniPath build stores it. A lookup that compared the two
# spellings directly would match nothing at all and report a mapping database
# that answers everything as one that answers nothing.
STAGED_MAPPINGS = (
    (18, 1, 9606, '7157', 'Q00000'),
    (18, 1, 9606, '7158', 'Q99999'),
    (18, 1, 9606, '7158', 'P11111'),
    # Published ChEBI to InChIKey, which is the only direction the real mapping
    # table holds it in, so reading it backwards is the whole of the chemical
    # lookup rather than an optimisation of it.
    (44, 52, 0, 'CHEBI:30000', 'FIXTUREINCHIKA-FIXTUREKEY-N'),
    (44, 52, 0, 'CHEBI:30001', 'FIXTUREINCHIKA-FIXTUREKEY-N'),
)

# The protein resolver, which is the gene side's first oracle: it covers 32,838
# organisms against the mapping table's two, and where both answer it returns
# the reviewed accession while the mapping table lists every accession the gene
# ever reached.
#
# `7157` is here **and** in the mapping table, under two different accessions,
# because the order the two are asked in is a claim the projection makes rather
# than an implementation detail. `7158` is in the mapping table alone, so the
# table still answers for what the resolver does not hold.
STAGED_PROTEINS = (
    (9606, 'entrez', '7157', 'P04637'),
)

# Which bare numerals the mapping database knows as ChEBI identifiers. This is
# the oracle an untyped numeral is checked against, and the reason a numeral it
# does not answer for stays a numeral.
STAGED_CHEMICALS = (
    ('chebi', 'CHEBI:17234', 'FIXTUREKNOWNKA-FIXTUREKEY-N'),
)

# What the translation is expected to produce, per fixture entity: the
# identifier the label carries and the namespace the columns report. Written
# out rather than derived, because a table computed from the same rules the
# step applies would agree with it however wrong the rules are.
TRANSLATED = {
    # Typed ChEBI. Only the prefix changes, and it has to.
    'cos_chebi_in': ('CHEBI:15422', 'chebi'),
    'cos_chebi_out': ('CHEBI:16810', 'chebi'),
    'cos_shaped_out': ('CHEBI:100', 'chebi'),
    'cos_ambiguous_out': ('CHEBI:200', 'chebi'),
    # An untyped numeral the mapping database knows, and one it does not.
    'cos_numeral_out': ('CHEBI:17234', 'chebi'),
    'cos_numeral_in': ('90000001', 'unresolved'),
    # An InChIKey two ChEBI entries claim.
    'cos_ambiguous_in': ('FIXTUREINCHIKA-FIXTUREKEY-N', 'inchikey'),
    # An untyped string that is no identifier of any namespace.
    'cos_shaped_in': ('FIXTURE_cos_shaped_in', 'unresolved'),
    # The gene side: one mapped, one ambiguous and resolved low, one with no
    # mapping at all, one never typed but plainly an accession.
    'cos_entrez_enz': ('P04637', 'uniprot'),
    'cos_entrez_many_enz': ('P11111', 'uniprot'),
    'cos_entrez_missing_enz': ('90000002', 'entrez'),
    'cos_shaped_enz': ('P12345', 'uniprot'),
}

# The shapes a namespace column claims about the string beside it. A column
# saying `uniprot` over a bare numeral would be worse than one saying `entrez`,
# because a consumer filtering on the column would take the row.
NAMESPACE_SHAPE = {
    'chebi': re.compile(r'^CHEBI:[0-9]+$'),
    'entrez': re.compile(r'^[0-9]+$'),
    'uniprot': re.compile(
        r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]'
        r'|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$'
    ),
    'inchikey': re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'),
}

METABOLITE_PREFIX = 'Metab__'


def _stage_mappings(conn, schema: str) -> None:
    """Build a mapping database in the shape the utils schema publishes.

    Three tables and ten rows, named and keyed exactly as the real ones are, so
    the projection reaches them through the queries it runs in production
    rather than through a seam opened for the test.
    """
    from psycopg2 import sql

    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL('DROP SCHEMA IF EXISTS {} CASCADE').format(schema_id)
        )
        cur.execute(sql.SQL('CREATE SCHEMA {}').format(schema_id))
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE {}.id_type (
                  id smallint PRIMARY KEY,
                  name text NOT NULL
                )
                """
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE {}.id_mapping (
                  source_type_id smallint NOT NULL,
                  target_type_id smallint NOT NULL,
                  ncbi_tax_id integer NOT NULL,
                  source_id varchar(64) NOT NULL,
                  target_id varchar(64) NOT NULL
                )
                """
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE {}.resolver_protein (
                  ncbi_tax_id integer,
                  source_type varchar(64),
                  source_id varchar(64),
                  uniprot varchar(64)
                )
                """
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE {}.resolver_chemical (
                  source_type varchar(64),
                  source_id varchar(64),
                  inchikey varchar(64)
                )
                """
            ).format(schema_id)
        )
        cur.executemany(
            sql.SQL('INSERT INTO {}.id_type (id, name) VALUES (%s, %s)')
            .format(schema_id)
            .as_string(conn),
            STAGED_ID_TYPES,
        )
        cur.executemany(
            sql.SQL(
                'INSERT INTO {}.id_mapping (source_type_id, target_type_id, '
                'ncbi_tax_id, source_id, target_id) VALUES (%s, %s, %s, %s, %s)'
            )
            .format(schema_id)
            .as_string(conn),
            STAGED_MAPPINGS,
        )
        cur.executemany(
            sql.SQL(
                'INSERT INTO {}.resolver_protein '
                '(ncbi_tax_id, source_type, source_id, uniprot) '
                'VALUES (%s, %s, %s, %s)'
            )
            .format(schema_id)
            .as_string(conn),
            STAGED_PROTEINS,
        )
        cur.executemany(
            sql.SQL(
                'INSERT INTO {}.resolver_chemical '
                '(source_type, source_id, inchikey) VALUES (%s, %s, %s)'
            )
            .format(schema_id)
            .as_string(conn),
            STAGED_CHEMICALS,
        )
    conn.commit()


@pytest.fixture(scope='module')
def mappings(conn):
    """The staged mapping database, for the length of the module."""
    from psycopg2 import sql

    _stage_mappings(conn, MAPPINGS_SCRATCH)
    try:
        yield MAPPINGS_SCRATCH
    finally:
        conn.rollback()
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL('DROP SCHEMA IF EXISTS {} CASCADE').format(
                    sql.Identifier(MAPPINGS_SCRATCH)
                )
            )
        conn.commit()


@pytest.fixture(scope='module')
def translated_scratch(conn):
    """A second copy of the fixture graph, for the translated projection."""
    from omnipath_build.db import schema as build_schema
    from omnipath_build.db.derived_tables import rebuild_interaction_tables

    build_schema.ensure_schema(
        conn, schema=TRANSLATION_SCRATCH, drop_existing=True
    )
    conn.commit()
    build_interaction_fixture(conn, TRANSLATION_SCRATCH)
    rebuild_interaction_tables(conn, schema=TRANSLATION_SCRATCH)
    conn.commit()
    try:
        yield TRANSLATION_SCRATCH
    finally:
        conn.rollback()
        with conn.cursor() as cur:
            cur.execute(
                f'DROP SCHEMA IF EXISTS {TRANSLATION_SCRATCH} CASCADE'
            )
        conn.commit()


@pytest.fixture(scope='module')
def translated(conn, translated_scratch, mappings):
    """The projection with the mappings in reach, and what it reported.

    The mapping database is another schema of the same Postgres rather than
    another server, which is a difference the step cannot see: it opens a
    second connection either way and reads the same three tables through it.
    """
    from omnipath_build.cosmos import (
        build_cosmos_projection,
        ensure_cosmos_edge_table,
    )

    ensure_cosmos_edge_table(conn, schema=translated_scratch)
    conn.commit()
    try:
        stats = build_cosmos_projection(
            conn,
            schema=translated_scratch,
            utils_db_url=DATABASE_URL,
            utils_schema=mappings,
            reaction_sources=FIXTURE_SOURCES,
        )
    except Exception:
        conn.rollback()
        raise
    conn.commit()
    return stats


def _translated_edges(conn, translated_scratch, name: str) -> list[dict]:
    """The reaction edges of one named fixture reaction, after translation."""
    return _reaction_edges(conn, translated_scratch, name)


def _endpoint_identifier(label: str, compartment: str | None) -> str:
    """The bare identifier inside a node label.

    A metabolite label is the prefix, the identifier and the compartment; a
    gene label is the prefix with the reaction index, the identifier and, on
    the mirrored half, the suffix. Everything else — the bare side of a
    connector — is the identifier already.
    """
    if label.startswith(METABOLITE_PREFIX):
        body = label[len(METABOLITE_PREFIX):]
        if compartment and body.endswith(f'_{compartment}'):
            body = body[: -len(compartment) - 1]
        return body
    gene = GENE_NODE.match(label)
    if gene:
        body = gene.group(2)
        if body.endswith(REVERSE_SUFFIX):
            body = body[: -len(REVERSE_SUFFIX)]
        return body
    return label


def _endpoints(edges: list[dict]) -> list[tuple[str, str | None]]:
    """Every (identifier, namespace) pair the given edges claim."""
    return [
        (
            _endpoint_identifier(edge[f'{side}_label'],
                                 edge[f'{side}_compartment']),
            edge[f'{side}_id_type'],
        )
        for edge in edges
        for side in ('source', 'target')
    ]


class TestAChebiMetaboliteCarriesThePrefixTheNetworkReads:
    """The one change that is not optional.

    The build stores a ChEBI identifier as a bare numeral, `15422`. Every other
    database that holds one, the delivered network included, writes
    `CHEBI:15422`. Without the prefix every metabolite node in the projection
    differs from its counterpart in the network for a reason that has nothing
    to do with the chemistry, and a consumer joining the two matches nothing.
    """

    @pytest.fixture(scope='class')
    def edges(self, conn, translated_scratch, translated):
        found = _translated_edges(conn, translated_scratch, 'cos_rxn_named')
        assert found, 'the translated conversion produced no edges at all'
        return found

    def test_the_substrate_label_is_the_prefixed_identifier(self, edges):
        assert 'Metab__CHEBI:15422_c' in {
            edge['source_label'] for edge in edges
        }

    def test_the_product_label_is_the_prefixed_identifier(self, edges):
        assert 'Metab__CHEBI:16810_m' in {
            edge['target_label'] for edge in edges
        }

    def test_the_bare_numeral_is_nowhere_in_the_labels(self, edges):
        """A label the build could have emitted before, and must not now."""
        labels = {
            label
            for edge in edges
            for label in (edge['source_label'], edge['target_label'])
        }
        assert 'Metab__15422_c' not in labels
        assert 'Metab__16810_m' not in labels

    def test_the_prefix_is_applied_once(self, edges):
        """Prefixing an already prefixed identifier is its own bug."""
        for edge in edges:
            for label in (edge['source_label'], edge['target_label']):
                assert label.count('CHEBI:') <= 1

    def test_the_namespace_says_chebi(self, edges):
        namespaces = {
            namespace
            for identifier, namespace in _endpoints(edges)
            if identifier.startswith('CHEBI:')
        }
        assert namespaces == {'chebi'}


class TestACatalystTranslatesWhereAMappingAnswers:
    """An Entrez gene becomes the accession COSMOS is written against."""

    @pytest.fixture(scope='class')
    def edges(self, conn, translated_scratch, translated):
        return _translated_edges(conn, translated_scratch, 'cos_rxn_named')

    @pytest.fixture(scope='class')
    def gene_label(self, edges):
        return _one_gene_node(edges)

    def test_the_node_carries_the_accession(self, gene_label):
        assert gene_label.endswith('__P04637')

    def test_the_gene_identifier_is_gone_from_the_node(self, gene_label):
        assert '7157' not in gene_label

    def test_the_namespace_says_uniprot(self, edges, gene_label):
        namespaces = {
            edge[f'{side}_id_type']
            for edge in edges
            for side in ('source', 'target')
            if edge[f'{side}_label'] == gene_label
        }
        assert namespaces == {'uniprot'}

    def test_the_connector_starts_from_the_translated_identifier(
        self,
        conn,
        translated_scratch,
        gene_label,
    ):
        """The attachment point has to be the identifier the node carries.

        A connector running from `7157` to `Gene3__P04637` would leave the node
        reachable only by a caller holding the identifier the build happened to
        canonicalise to, which is the one thing the connector exists to stop.
        """
        connectors = _edges(
            conn,
            translated_scratch,
            "interaction_type = 'connector' AND target_label = %s",
            [gene_label],
        )
        assert [edge['source_label'] for edge in connectors] == ['P04637']

    def test_the_resolver_answers_before_the_mapping_table(self, gene_label):
        """Both oracles hold this gene, under different accessions.

        The order is a claim, not an accident. The mapping table lists every
        accession an Entrez gene ever reached, so choosing among them by sort
        order returns an unreviewed isoform; the resolver returns the one
        accession the gene is known by. Where both answer, that one stands.
        """
        assert not gene_label.endswith('__Q00000')

    def test_one_gene_reaches_one_node_even_where_two_accessions_answer(
        self,
        conn,
        translated_scratch,
        translated,
    ):
        """Ambiguity resolves low, and it resolves to a single node.

        Two accessions for one Entrez gene are two names for the enzyme, not
        two enzymes. Emitting both would say the reaction runs twice, and
        picking whichever the mapping table returned first would move the label
        between rebuilds; the lowest by sort order does neither.
        """
        edges = _translated_edges(
            conn, translated_scratch, 'cos_rxn_ambiguous'
        )
        assert _one_gene_node(edges).endswith('__P11111')


class TestACatalystWithNoMappingKeepsWhatTheBuildHolds:
    """Entrez is an identifier. An unanswered lookup is not a reason to lie."""

    @pytest.fixture(scope='class')
    def edges(self, conn, translated_scratch, translated):
        return _translated_edges(conn, translated_scratch, 'cos_rxn_unmapped')

    @pytest.fixture(scope='class')
    def gene_label(self, edges):
        return _one_gene_node(edges)

    def test_the_node_keeps_the_gene_identifier(self, gene_label):
        assert gene_label.endswith('__90000002')

    def test_the_namespace_says_entrez(self, edges, gene_label):
        namespaces = {
            edge[f'{side}_id_type']
            for edge in edges
            for side in ('source', 'target')
            if edge[f'{side}_label'] == gene_label
        }
        assert namespaces == {'entrez'}

    def test_the_reaction_keeps_both_sides_of_its_arrow(self, edges):
        """A catalyst nothing could translate loses no chemistry."""
        assert len(edges) == 2


class TestAnUntypedAccessionIsReportedAsUniprot:
    """The largest population on the gene side is this one.

    140,312 of the build's 140,345 untyped catalysts are UniProt accessions
    that merely failed to be typed as such. Reading the declared type and
    calling them unresolved would report the gene side as almost entirely
    untranslated when almost all of it is already in the namespace COSMOS
    wants. The column follows the string, not the type.
    """

    @pytest.fixture(scope='class')
    def edges(self, conn, translated_scratch, translated):
        return _translated_edges(conn, translated_scratch, 'cos_rxn_shaped')

    @pytest.fixture(scope='class')
    def gene_label(self, edges):
        return _one_gene_node(edges)

    def test_the_label_is_unchanged(self, gene_label):
        assert gene_label.endswith('__P12345')

    def test_the_namespace_says_uniprot(self, edges, gene_label):
        namespaces = {
            edge[f'{side}_id_type']
            for edge in edges
            for side in ('source', 'target')
            if edge[f'{side}_label'] == gene_label
        }
        assert namespaces == {'uniprot'}

    def test_a_string_of_no_namespace_stays_unresolved(self, edges):
        """The same untyped population, where the string is not an accession."""
        found = {
            namespace
            for identifier, namespace in _endpoints(edges)
            if identifier == 'FIXTURE_cos_shaped_in'
        }
        assert found == {'unresolved'}


class TestAnUntypedNumeralIsChebiOnlyWhereTheMappingSaysSo:
    """Shape is not evidence.

    About 28 per cent of the build's untyped bare numerals are ChEBI
    identifiers and the rest are model-local accessions that merely look like
    them. Prefixing all of them would put a wrong identifier into three
    published labels out of four, which is worse than an honest fallback.
    """

    @pytest.fixture(scope='class')
    def endpoints(self, conn, translated_scratch, translated):
        return dict(
            _endpoints(
                _translated_edges(
                    conn, translated_scratch, 'cos_rxn_unmapped'
                )
            )
        )

    def test_a_confirmed_numeral_becomes_a_chebi_identifier(self, endpoints):
        assert endpoints['CHEBI:17234'] == 'chebi'

    def test_an_unconfirmed_numeral_is_left_alone(self, endpoints):
        assert endpoints['90000001'] == 'unresolved'

    def test_the_unconfirmed_numeral_never_gains_the_prefix(self, endpoints):
        assert 'CHEBI:90000001' not in endpoints


class TestAnAmbiguousChemicalMappingKeepsTheBuildsIdentifier:
    """The chemical side resolves an ambiguity the other way from the gene one.

    Two ChEBI entries claiming one InChIKey is a disagreement about
    stereochemistry, and choosing between them would state a structure neither
    database did. An InChIKey is an identifier a consumer can work with, so the
    honest answer keeps it and says so.
    """

    @pytest.fixture(scope='class')
    def endpoints(self, conn, translated_scratch, translated):
        return dict(
            _endpoints(
                _translated_edges(
                    conn, translated_scratch, 'cos_rxn_ambiguous'
                )
            )
        )

    def test_the_inchikey_survives(self, endpoints):
        assert (
            endpoints['FIXTUREINCHIKA-FIXTUREKEY-N'] == 'inchikey'
        )

    def test_neither_candidate_was_picked(self, endpoints):
        assert 'CHEBI:30000' not in endpoints
        assert 'CHEBI:30001' not in endpoints


class TestEveryTranslatedLabelIsTheOneTheTableWasToldToExpect:
    """One assertion per fixture entity, against a written-out answer."""

    @pytest.fixture(scope='class')
    def endpoints(self, conn, translated_scratch, translated):
        return dict(_endpoints(_edges(conn, translated_scratch)))

    @pytest.mark.parametrize('name', sorted(TRANSLATED))
    def test_the_label_and_the_namespace_are_both_what_was_expected(
        self,
        endpoints,
        name,
    ):
        identifier, namespace = TRANSLATED[name]
        assert identifier in endpoints, (
            f'{name} should be labelled {identifier}; '
            f'the projection holds {sorted(endpoints)}'
        )
        assert endpoints[identifier] == namespace


class TestTheNamespaceColumnAgreesWithTheLabel:
    """The column is the whole point, so it has to be true of every row.

    A consumer that needs UniProt takes the endpoints whose namespace says
    `uniprot` and leaves the rest, without parsing one label. That only works
    if a column never claims a namespace the string beside it does not belong
    to, which is asserted here over every endpoint of every edge rather than
    over the ones the other classes name.
    """

    @pytest.fixture(scope='class')
    def endpoints(self, conn, translated_scratch, translated):
        return _endpoints(_edges(conn, translated_scratch))

    def test_every_endpoint_names_a_namespace(self, endpoints):
        assert all(namespace for _identifier, namespace in endpoints)

    def test_no_column_holds_a_vocabulary_name(self, endpoints):
        """`chebi`, not `Chebi:MI:0474`.

        The build names an identifier type after the controlled-vocabulary term
        it came from. That is the right name inside the build and the wrong one
        in a published network, where the column is read by a consumer who has
        never heard of the vocabulary.
        """
        for _identifier, namespace in endpoints:
            assert ':' not in namespace
            assert namespace == namespace.lower()

    def test_a_claimed_namespace_matches_the_string_beside_it(
        self,
        endpoints,
    ):
        for identifier, namespace in endpoints:
            shape = NAMESPACE_SHAPE.get(namespace)
            if shape is None:
                continue
            assert shape.match(identifier), (
                f'{identifier} is labelled {namespace} and does not look '
                f'like one'
            )

    def test_the_step_counts_what_it_reached(self, translated):
        """Each side reports what it reached and what it fell back to."""
        assert translated.identifier_translation is True
        assert translated.gene_labels_mapped >= 2
        assert translated.gene_labels_translated >= (
            translated.gene_labels_mapped
        )
        assert translated.chemical_labels_translated >= 5
        assert translated.chemical_labels_fallback >= 1


class TestWithoutAMappingDatabaseTheProjectionIsStillWhole:
    """A missing mapping database narrows the output and never fails it.

    This is the path every run of this file outside these classes takes, and it
    is the path a build machine off the lab network takes. Every label keeps
    the identifier the build canonicalised it to, every namespace column says
    which one that is, and the edge set is the same edge set.
    """

    @pytest.fixture(scope='class')
    def endpoints(self, conn, scratch, projection):
        return dict(_endpoints(_edges(conn, scratch)))

    def test_the_step_reports_that_nothing_was_reachable(self, projection):
        assert projection.identifier_translation is False

    def test_nothing_was_mapped(self, projection):
        assert projection.gene_labels_mapped == 0
        assert projection.chemical_labels_mapped == 0

    def test_the_edge_set_is_the_same_one(self, conn, scratch, projection,
                                          translated_scratch, translated):
        """Same reactions, same arrows, same count. Only the labels differ."""
        assert projection.edges == translated.edges
        assert projection.connectors == translated.connectors
        assert projection.reactions == translated.reactions

    def test_a_catalyst_keeps_its_gene_identifier(self, endpoints):
        assert endpoints['7157'] == 'entrez'
        assert 'P04637' not in endpoints

    def test_an_untyped_numeral_is_not_guessed_at(self, endpoints):
        """With nothing to ask, a numeral cannot be confirmed as a ChEBI."""
        assert endpoints['17234'] == 'unresolved'
        assert 'CHEBI:17234' not in endpoints

    def test_the_chebi_prefix_is_still_applied(self, endpoints):
        """The prefix is a spelling of an identifier the build already holds.

        Nothing has to be looked up to know that the build's `15422` and the
        network's `CHEBI:15422` are the same identifier, so the one change the
        network cannot do without survives having no mapping database at all.
        """
        assert endpoints['CHEBI:15422'] == 'chebi'

    def test_an_untyped_accession_is_still_read_as_one(self, endpoints):
        """Reading the shape of a string needs no database either."""
        assert endpoints['P12345'] == 'uniprot'


def test_a_mapping_database_that_answers_with_an_error_falls_back():
    """A reachable mapping database is not the same as a usable one.

    An older utils build with one of these tables missing, or a connection
    without the rights to read them, answers every query with an error. The
    step cannot tell that from an absent database and must not treat it as a
    reason to fail a projection over content the build already wrote, so both
    land on the same fallback and the step records that nothing answered.

    Pointed at a schema that holds no mapping tables at all, on a connection of
    its own: the failing statements abort the transaction they run in, and
    lending them the connection the rest of this file builds on would take the
    other fixtures down with them.
    """
    import psycopg2

    from omnipath_build.cosmos.translate import (
        IdentifierMappings,
        LabelEntity,
        translate_identifiers,
    )

    entity_id = ENTITY['cos_entrez_enz']
    entities = [
        LabelEntity(entity_id, '7157', 'Entrez:MI:0477', 9606, 'gene'),
    ]
    connection = psycopg2.connect(DATABASE_URL)
    connection.set_session(readonly=True, autocommit=True)
    try:
        mappings = IdentifierMappings(connection, schema='pg_catalog')
        labels = translate_identifiers(entities, mappings)
    finally:
        connection.close()

    assert mappings.failed is True
    assert labels[entity_id].identifier == '7157'
    assert labels[entity_id].namespace == 'entrez'
    assert labels[entity_id].mapped is False


class TestOnlyTheNamedResourcesReachTheNetwork:
    """The network is a metabolic one, and the resource list is what says so.

    A general reaction database states its events in proteins and complexes
    rather than small molecules. One reaching the output would be labelled a
    metabolite by the role it plays rather than by what it is, so the gate sits
    on the resource that published the reaction.
    """

    def test_a_reaction_outside_the_named_resources_is_absent(
        self,
        conn,
        scratch,
    ):
        """Narrowing the list to one resource leaves the others' reactions out."""
        from omnipath_build.cosmos import (
            build_cosmos_projection,
            ensure_cosmos_edge_table,
        )

        ensure_cosmos_edge_table(conn, schema=scratch)
        narrowed = build_cosmos_projection(
            conn,
            schema=scratch,
            reaction_sources=(),
        )
        conn.commit()
        assert narrowed.edges == 0, (
            'no resource is in scope, so the projection has nothing to publish'
        )
        assert narrowed.reactions == 0

        restored = build_cosmos_projection(
            conn,
            schema=scratch,
            reaction_sources=FIXTURE_SOURCES,
        )
        conn.commit()
        assert restored.edges > 0, (
            'and naming the resources again brings the whole network back'
        )

    def test_the_gate_reads_the_resources_that_published_the_reaction(
        self,
        conn,
        scratch,
    ):
        """A reaction two resources report survives on either one's word.

        The overlap rather than a containment test is what keeps a reaction
        whose second contributor is out of scope.
        """
        rows = _rows(
            conn,
            f'SELECT sources FROM {scratch}.interaction '
            'WHERE arity > 2 AND sources IS NOT NULL LIMIT 5',
        )
        assert rows, 'the fixture holds reaction headers carrying resources'
        for (sources,) in rows:
            assert set(sources) & set(FIXTURE_SOURCES), (
                'every fixture reaction names a resource the tests scope to'
            )
