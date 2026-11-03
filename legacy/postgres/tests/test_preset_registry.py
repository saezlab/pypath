"""The preset registry: ``network_registry`` carries a full preset spec.

A dataset is a preset — a ``network_registry`` row that filters the interaction
fact table — not a bespoke matview. These tests assert the registry round-trips
every field of that spec (class scope, evidence scope, default and mandatory
attributes, labels, curation, attribute sources) and that the live table carries
the columns to hold them.

**Amended 2026-08-20** with the two ``network_registry`` columns the grain
amendment adds:

``collapse_mode``
    How the preset collapses the per-resource record over **its own** resource
    scope: ``none`` (one row per resource assertion), ``assertion`` (fold the
    resources that agree on sign and direction) or ``endpoints`` (fold to the
    collapsed key). ``endpoints`` is the default, reproducing the legacy
    one-row-per-interaction contract. A single-resource preset collapses
    nothing whatever the mode says, and that gets a test of its own.
``license_scope``
    The minimum ``purpose`` / ``sharing`` / ``attrib`` levels a resource must
    meet to contribute. No ``license_scope`` means unrestricted; a scope that
    is set excludes a resource whose license is unknown, however permissive its
    recorded levels look.

``grain``
    What the preset's groups are keyed on: ``interaction`` (the default) groups
    on the binary triple of subject entity, object entity and interaction
    class, and ``participant`` groups on the ``interaction_id`` header alone,
    so one group is one interaction whatever its arity. The grain is read
    before ``collapse_mode`` and can silence it — each collapse mode folds an
    ordered endpoint pair, and a reaction is not one — so a preset at
    ``participant`` grain records ``none`` as its mode, which states no fold
    rather than the wrong one. Nothing enforces that across the two columns;
    it is how the pair is read, and the tests here say so.

**Amended 2026-08-21** with ``composition``, the column a preset carries when
it is not one query: the ordered component list and the operation that joins
them (``union``, ``collapse``, ``exclude``, ``annotate``). A component is
either a parameter set or the **name of another preset**, which is what gives
``nichenet`` its per-component override. NULL means the preset is a single
parameter set — the common case, and the shape every other test in this
module registers.

Two orders in that value are binding rather than stylistic, and the column has
to be able to state them: the ``collapse`` runs after the ``union`` and over the
union's own resolved scope, and the ``exclude`` runs **before** the ``collapse``
— dropping a resource after the fold leaves its contribution inside
``source_count``, ``references`` and the sign flags. The algebra itself lives
in the api-service; what is asserted here is that a registry round-trip
preserves the order it was written in.

A third column, ``materialize_collapse``, was proposed and **withdrawn on
2026-08-20**: both interaction tables are built unconditionally, so no preset
carries a materialisation flag. Nothing here tests for it.

The license assertions here are about the *preset spec* — that the registry
stores enough to resolve a scope, and that resolving it excludes unknown terms.
The query path carries its own version of that rule.

Run against a build database, e.g. on dev4::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:5404/omnipath \
        uv run --with pytest pytest tests/test_preset_registry.py -v

The round-trip runs in a throwaway schema, so it never touches the presets the
build registered in ``public``. Skipped when DATABASE_URL is not set.
"""

from __future__ import annotations

import os

import pytest

from omnipath_build.network_views import (
    NetworkDefinition,
    ensure_network_registry,
    register_network,
)

DATABASE_URL = os.environ.get('DATABASE_URL')
SCHEMA = os.environ.get('OMNIPATH_PG_SCHEMA', 'public')
TEST_SCHEMA = 'preset_registry_roundtrip_test'

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; preset-registry test needs a database',
)


# A preset spec with every field populated — the shape `network_registry` has
# to carry. `interaction_class_scope` holds `vocab_interaction_class` slugs,
# never legacy dataset names; the class derivation itself is another task.
# Built by a factory rather than at import, so a missing field fails the test
# that needs it instead of the whole module at collection.
def _full_preset(**overrides: object) -> NetworkDefinition:
    spec: dict[str, object] = {
        'name': '_roundtrip_preset',
        'kind': 'ligand_receptor',
        'included_sources': ('connectomedb2025', 'cellphonedb'),
        'interaction_class_scope': ('ligand_receptor', 'transport'),
        'evidence_scope': {
            'evidence_type': ['experimental', 'curated'],
            'predicate': ['binds'],
            'min_confidence': 2,
        },
        'default_attributes': ('endpoints', 'references'),
        'mandatory_attributes': ('label', 'evidence'),
        'labels': {
            'preset': 'Round-trip preset',
            'columns': {'label': 'Interaction label'},
        },
        'curation': {
            'moa_only': True,
            'affinity_cutoff': 6.0,
            'metabolite_class_gate': True,
        },
        'attribute_sources': {
            'protein_localization': {'stage': 'interim', 'source': 'uniprot'},
        },
        'collapse_mode': 'assertion',
        'license_scope': {'purpose': 15, 'sharing': 0, 'attrib': 0},
    }
    spec.update(overrides)
    return NetworkDefinition(**spec)


PRESET_COLUMNS = {
    'interaction_class_scope': '_text',
    'evidence_scope': 'jsonb',
    'default_attributes': '_text',
    'mandatory_attributes': '_text',
    'labels': 'jsonb',
    'curation': 'jsonb',
    'attribute_sources': 'jsonb',
    # The grain amendment.
    'collapse_mode': 'text',
    'license_scope': 'jsonb',
    # The composition algebra.
    'composition': 'jsonb',
    # The key the preset's groups are folded on.
    'grain': 'text',
}

COLLAPSE_MODES = ('none', 'assertion', 'endpoints')

GRAINS = ('interaction', 'participant')

# The group key each collapse mode folds to. The derive and the query path
# share one collapse routine; this mirrors its keys so the registry test can
# state what a mode *means* without importing it.
COLLAPSE_KEYS = {
    'none': (
        'subject_entity_id, object_entity_id, interaction_class_id, source_id, '
        'is_directed, is_stimulation, is_inhibition'
    ),
    'assertion': (
        'subject_entity_id, object_entity_id, interaction_class_id, '
        'is_directed, is_stimulation, is_inhibition'
    ),
    'endpoints': 'subject_entity_id, object_entity_id, interaction_class_id',
}


# The composition of `metalinksdb`: a union of three parameter sets, then the
# exclusion, then the fold, then the annotation layer. The two binding orders
# are visible in the value itself — the `collapse` follows the `union` and the
# `exclude` precedes the `collapse`.
METALINKSDB_COMPOSITION = {
    'operation': 'union',
    'components': [
        {
            'parameters': {
                'resources': ['chembl'],
                'curation_flags': ['chembl_mechanism'],
            }
        },
        {'parameters': {'resources': ['cellinker', 'guidetopharma', 'stitch']}},
        {
            'parameters': {
                'resources': ['recon3d', 'humangem'],
                'interaction_classes': ['transport'],
            }
        },
    ],
    'steps': [
        {'operation': 'exclude', 'resources': ['bindingdb']},
        {'operation': 'collapse', 'collapse_mode': 'endpoints'},
        {'operation': 'annotate', 'layers': ['entity_intercell_annotation']},
    ],
}

# The composition of `nichenet`: the components are *named presets*, not
# parameter sets. That is where the per-component override comes from —
# an override replaces one named component and leaves the recipe alone.
NICHENET_COMPOSITION = {
    'operation': 'union',
    'components': [
        {'preset': 'curated_ligand_receptor'},
        {'preset': 'omnipath'},
        {'preset': 'collectri'},
    ],
    'steps': [{'operation': 'collapse', 'collapse_mode': 'endpoints'}],
}


def _component_names(composition):
    """Each entry's named preset, or None where the component is a parameter set."""
    return [component.get('preset') for component in composition['components']]


def _step_operations(composition):
    return [step['operation'] for step in composition.get('steps', ())]


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture(scope='module')
def registry(conn):
    """A throwaway registry schema, so the build's own presets stay untouched."""
    with conn.cursor() as cur:
        cur.execute(f'DROP SCHEMA IF EXISTS {TEST_SCHEMA} CASCADE')
        cur.execute(f'CREATE SCHEMA {TEST_SCHEMA}')
    conn.commit()
    ensure_network_registry(conn, registry_schema=TEST_SCHEMA)
    yield TEST_SCHEMA
    with conn.cursor() as cur:
        cur.execute(f'DROP SCHEMA IF EXISTS {TEST_SCHEMA} CASCADE')
    conn.commit()


@pytest.fixture(autouse=True)
def _clean_transaction(conn):
    """Leave the shared connection usable after a test that hit a SQL error."""
    yield
    conn.rollback()


@pytest.fixture(scope='module')
def licenses(conn, registry):
    """A throwaway license catalogue, shaped like ``data_source_license``.

    ``mystery_db`` is the trap: its recorded levels are maximal, but nothing
    maps a license to it, so a license-filtered scope must drop it.
    """
    with conn.cursor() as cur:
        cur.execute(
            f"""
            CREATE TABLE {registry}.data_source (
              source_id bigint PRIMARY KEY,
              name text NOT NULL UNIQUE
            );
            CREATE TABLE {registry}.data_source_license (
              source_id bigint PRIMARY KEY
                REFERENCES {registry}.data_source (source_id),
              license_name text,
              purpose_level smallint NOT NULL,
              sharing_level smallint NOT NULL,
              attrib_level smallint NOT NULL,
              is_known boolean NOT NULL
            );
            INSERT INTO {registry}.data_source (source_id, name) VALUES
              (1, 'signor'), (2, 'cellphonedb'), (3, 'mystery_db');
            INSERT INTO {registry}.data_source_license VALUES
              (1, 'CC-BY', 20, 25, 5, true),
              (2, 'academic-only', 5, 5, 5, true),
              (3, NULL, 25, 25, 10, false);
            """
        )
    conn.commit()
    return registry


@pytest.fixture(scope='module')
def records(conn, registry, licenses):
    """A throwaway record table, shaped like ``interaction_fact_resource``.

    Two triples. ``signor`` asserts each of them once; ``cellphonedb`` asserts
    the first one too, with the same signature. So a two-resource scope folds
    under ``endpoints`` and a single-resource scope has nothing to fold.
    """
    with conn.cursor() as cur:
        cur.execute(
            f"""
            CREATE TABLE {registry}.interaction_fact_resource (
              subject_entity_id bigint NOT NULL,
              object_entity_id bigint NOT NULL,
              interaction_class_id int NOT NULL,
              source_id bigint NOT NULL,
              is_directed boolean,
              is_stimulation boolean,
              is_inhibition boolean
            );
            INSERT INTO {registry}.interaction_fact_resource VALUES
              (10, 20, 1, 1, true, true, NULL),
              (30, 40, 1, 1, true, NULL, NULL),
              (10, 20, 1, 2, true, true, NULL);
            """
        )
    conn.commit()
    return registry


def _preset_row(conn, schema, name):
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT kind, included_sources, interaction_class_scope,
                   evidence_scope, default_attributes, mandatory_attributes,
                   labels, curation, attribute_sources,
                   collapse_mode, license_scope, composition
            FROM {schema}.network_registry WHERE name = %s
            """,
            [name],
        )
        return cur.fetchone()


def _collapsed_row_count(conn, schema, mode, sources):
    """Rows a scope yields once the preset's collapse mode has folded it."""
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT count(*) FROM (
              SELECT {COLLAPSE_KEYS[mode]}
              FROM {schema}.interaction_fact_resource f
              JOIN {schema}.data_source d USING (source_id)
              WHERE d.name = ANY(%s)
              GROUP BY {COLLAPSE_KEYS[mode]}
            ) collapsed
            """,
            [list(sources)],
        )
        return cur.fetchone()[0]


def _contributing_sources(conn, schema, preset_name):
    """The resources a registered preset admits, after its ``license_scope``.

    The scope is a comparison over three ordinal levels, and a resource whose
    license is unknown fails it however permissive its levels read.
    """
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT included_sources, license_scope
            FROM {schema}.network_registry WHERE name = %s
            """,
            [preset_name],
        )
        included_sources, license_scope = cur.fetchone()
        if license_scope is None:
            return set(included_sources)
        cur.execute(
            f"""
            SELECT d.name
            FROM {schema}.data_source d
            JOIN {schema}.data_source_license l USING (source_id)
            WHERE d.name = ANY(%s)
              AND l.is_known
              AND l.purpose_level >= %s
              AND l.sharing_level >= %s
              AND l.attrib_level >= %s
            """,
            [
                list(included_sources),
                license_scope.get('purpose', 0),
                license_scope.get('sharing', 0),
                license_scope.get('attrib', 0),
            ],
        )
        return {name for (name,) in cur.fetchall()}


def test_registry_round_trips_the_full_preset_spec(conn, registry):
    """Every field of the preset spec survives registration unchanged."""
    preset = _full_preset()
    register_network(conn, preset, registry_schema=registry)
    row = _preset_row(conn, registry, preset.name)
    assert row is not None, 'the preset was not registered'
    (
        kind,
        included_sources,
        class_scope,
        evidence_scope,
        default_attributes,
        mandatory_attributes,
        labels,
        curation,
        attribute_sources,
        collapse_mode,
        license_scope,
        composition,
    ) = row
    assert kind == preset.kind
    assert tuple(included_sources) == preset.included_sources
    assert tuple(class_scope) == preset.interaction_class_scope
    assert evidence_scope == preset.evidence_scope
    assert tuple(default_attributes) == preset.default_attributes
    assert tuple(mandatory_attributes) == preset.mandatory_attributes
    assert labels == preset.labels
    assert curation == preset.curation
    assert attribute_sources == preset.attribute_sources
    assert collapse_mode == preset.collapse_mode
    assert license_scope == preset.license_scope
    # The full spec is still one parameter set: a composition is the exception,
    # and its absence is NULL rather than an empty recipe.
    assert composition is None


def test_re_registration_upserts_the_preset(conn, registry):
    """``register_network`` keeps its upsert: one row per preset, spec replaced."""
    preset = _full_preset()
    register_network(conn, preset, registry_schema=registry)
    revised = _full_preset(
        included_sources=('connectomedb2025',),
        interaction_class_scope=('ligand_receptor',),
        evidence_scope={'predicate': ['binds']},
        default_attributes=('endpoints',),
        mandatory_attributes=('label',),
        labels={'preset': 'Revised'},
        curation={'moa_only': False},
        attribute_sources={'protein_localization': {'stage': 'intercell'}},
        collapse_mode='none',
        license_scope=None,
    )
    register_network(conn, revised, registry_schema=registry)
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT count(*) FROM {registry}.network_registry WHERE name = %s',
            [preset.name],
        )
        assert cur.fetchone()[0] == 1
    row = _preset_row(conn, registry, preset.name)
    assert tuple(row[2]) == revised.interaction_class_scope
    assert row[3] == revised.evidence_scope
    assert row[6] == revised.labels
    # The amendment's columns are replaced by the upsert like any other field:
    # a revised preset may drop a license restriction it once carried.
    assert row[9] == revised.collapse_mode
    assert row[10] is None


def test_preset_without_matview_registers(conn, registry):
    """A preset is metadata over the fact table — it needs no schema/matview."""
    bare = NetworkDefinition(
        name='_roundtrip_preset_bare',
        kind='signaling',
        included_sources=('signor',),
        interaction_class_scope=('signaling',),
    )
    register_network(conn, bare, registry_schema=registry)
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT schema_name, combined_relation
            FROM {registry}.network_registry WHERE name = %s
            """,
            [bare.name],
        )
        schema_name, combined_relation = cur.fetchone()
    # The matview-era columns stay empty for a fact-table preset, and go once
    # the last bespoke matview retires.
    assert schema_name is None
    assert combined_relation is None


def test_collapse_mode_defaults_to_endpoints(conn, registry):
    """A preset that names no mode collapses to the endpoints key.

    That is the legacy one-row-per-interaction contract, so a preset written
    before the amendment keeps behaving as it did.
    """
    legacy = NetworkDefinition(
        name='_roundtrip_preset_default_mode',
        kind='signaling',
        included_sources=('signor',),
    )
    register_network(conn, legacy, registry_schema=registry)
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT collapse_mode FROM {registry}.network_registry WHERE name = %s',
            [legacy.name],
        )
        assert cur.fetchone()[0] == 'endpoints'


@pytest.mark.parametrize('mode', COLLAPSE_MODES)
def test_registry_round_trips_every_collapse_mode(conn, registry, mode):
    """Every ``network_registry.collapse_mode`` value survives registration."""
    preset = _full_preset(
        name=f'_roundtrip_preset_mode_{mode}',
        collapse_mode=mode,
    )
    register_network(conn, preset, registry_schema=registry)
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT collapse_mode FROM {registry}.network_registry WHERE name = %s',
            [preset.name],
        )
        assert cur.fetchone()[0] == mode


def test_single_resource_preset_collapses_nothing_whatever_the_mode(
    conn, registry, records
):
    """One resource in scope means every collapse group holds one row.

    The mode is a real choice only where the scope holds resources that can
    disagree — the two-resource control below shows the same data folding.
    """
    single = _full_preset(
        name='_roundtrip_preset_single_resource',
        included_sources=('signor',),
        license_scope=None,
    )
    register_network(conn, single, registry_schema=registry)
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT count(*) FROM {registry}.interaction_fact_resource f
            JOIN {registry}.data_source d USING (source_id)
            WHERE d.name = 'signor'
            """
        )
        record_rows = cur.fetchone()[0]
    counts = {
        mode: _collapsed_row_count(conn, registry, mode, single.included_sources)
        for mode in COLLAPSE_MODES
    }
    assert set(counts.values()) == {record_rows}, (
        f'a single-resource preset collapsed something: {counts} '
        f'over {record_rows} record rows'
    )
    # Control: with a second resource in scope the mode does change the answer,
    # so the assertion above is about the scope and not about empty data.
    two = ('signor', 'cellphonedb')
    assert _collapsed_row_count(conn, registry, 'endpoints', two) < (
        _collapsed_row_count(conn, registry, 'none', two)
    )


def test_grain_defaults_to_interaction(conn, registry):
    """A preset that names no grain groups on the binary triple.

    That is the pair every dataset was written around, so a definition from
    before the grain existed is registered unchanged rather than acquiring a
    key nobody chose for it.
    """
    legacy = NetworkDefinition(
        name='_roundtrip_preset_default_grain',
        kind='signaling',
        included_sources=('signor',),
    )
    register_network(conn, legacy, registry_schema=registry)
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT grain FROM {registry}.network_registry WHERE name = %s',
            [legacy.name],
        )
        assert cur.fetchone()[0] == 'interaction'


@pytest.mark.parametrize('grain', GRAINS)
def test_registry_round_trips_every_grain(conn, registry, grain):
    """Every ``network_registry.grain`` value survives registration."""
    preset = _full_preset(
        name=f'_roundtrip_preset_grain_{grain}',
        grain=grain,
    )
    register_network(conn, preset, registry_schema=registry)
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT grain FROM {registry}.network_registry WHERE name = %s',
            [preset.name],
        )
        assert cur.fetchone()[0] == grain


def test_participant_grain_round_trips_with_no_fold_claimed(conn, registry):
    """A participant-grain preset stores its grain and claims no fold.

    One group is one interaction whatever its arity, so there is no ordered
    endpoint pair for a collapse mode to describe. ``none`` is the value that
    says so; the pair has to come back as it was written, because a reader who
    trusted ``collapse_mode`` here would be reading a fold that never runs.
    """
    reactions = _full_preset(
        name='_roundtrip_preset_participant_grain',
        grain='participant',
        collapse_mode='none',
    )
    register_network(conn, reactions, registry_schema=registry)
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT grain, collapse_mode
            FROM {registry}.network_registry WHERE name = %s
            """,
            [reactions.name],
        )
        grain, collapse_mode = cur.fetchone()
    assert grain == 'participant'
    assert collapse_mode == 'none'


def test_grain_refuses_a_value_outside_the_two(conn, registry):
    """The registry takes the two grains and nothing else.

    A grain nobody implemented is not a preset that returns something unusual;
    it is a preset whose result shape no query path knows how to produce. The
    CHECK is what keeps that out of the table instead of leaving it for the
    api-service to discover.
    """
    import psycopg2

    invented = _full_preset(
        name='_roundtrip_preset_invented_grain',
        grain='reaction',
    )
    with pytest.raises(psycopg2.errors.CheckViolation):
        register_network(conn, invented, registry_schema=registry)


def test_preset_without_license_scope_is_unrestricted(conn, registry, licenses):
    """No ``license_scope`` restricts nothing — every included resource contributes."""
    unrestricted = _full_preset(
        name='_roundtrip_preset_no_license',
        included_sources=('signor', 'cellphonedb', 'mystery_db'),
        license_scope=None,
    )
    register_network(conn, unrestricted, registry_schema=registry)
    row = _preset_row(conn, registry, unrestricted.name)
    assert row[10] is None, 'an absent license scope must be stored as NULL'
    assert _contributing_sources(conn, registry, unrestricted.name) == {
        'signor',
        'cellphonedb',
        'mystery_db',
    }


def test_license_scope_excludes_an_unknown_license_resource(conn, registry, licenses):
    """Unknown terms are an exclusion, never a permissive default.

    ``mystery_db`` records the most permissive levels in the catalogue and is
    still dropped, because ``is_known`` is false. ``cellphonedb`` is dropped on
    the comparison itself — academic-only fails a commercial purpose level.
    """
    restricted = _full_preset(
        name='_roundtrip_preset_licensed',
        included_sources=('signor', 'cellphonedb', 'mystery_db'),
        license_scope={'purpose': 15, 'sharing': 0, 'attrib': 0},
    )
    register_network(conn, restricted, registry_schema=registry)
    row = _preset_row(conn, registry, restricted.name)
    assert row[10] == restricted.license_scope
    contributing = _contributing_sources(conn, registry, restricted.name)
    assert 'mystery_db' not in contributing, (
        'a resource with an unknown license was admitted to a license-scoped '
        'preset'
    )
    assert contributing == {'signor'}


def test_registry_round_trips_a_composition(conn, registry):
    """A preset that is not one query round-trips its whole recipe.

    Value equality is not enough on its own here: the recipe is *ordered*, so
    the component list and the step list are compared as sequences.
    """
    composed = _full_preset(
        name='_roundtrip_preset_composition',
        composition=METALINKSDB_COMPOSITION,
    )
    register_network(conn, composed, registry_schema=registry)
    stored = _preset_row(conn, registry, composed.name)[11]
    assert stored == METALINKSDB_COMPOSITION
    assert stored['operation'] == 'union'
    assert [
        component['parameters']['resources']
        for component in stored['components']
    ] == [
        component['parameters']['resources']
        for component in METALINKSDB_COMPOSITION['components']
    ], 'the component order did not survive the round trip'


def test_composition_keeps_the_two_binding_orders(conn, registry):
    """The stored recipe states both order rules, not only the operations.

    The ``collapse`` follows the ``union`` and works over the union's own
    resolved scope, and the ``exclude`` precedes the ``collapse`` — an exclusion
    applied after the fold leaves the dropped resource inside ``source_count``,
    ``references`` and the sign flags, which is the scoped-collapse rule under
    another name. The api-service executes this; the column has to say it.
    """
    composed = _full_preset(
        name='_roundtrip_preset_composition_order',
        composition=METALINKSDB_COMPOSITION,
    )
    register_network(conn, composed, registry_schema=registry)
    stored = _preset_row(conn, registry, composed.name)[11]
    operations = _step_operations(stored)
    assert stored['operation'] == 'union', 'the components are joined by the union'
    assert 'collapse' in operations, 'a multi-resource composition folds'
    # The union is the top-level operation and the collapse is a step, so the
    # fold cannot precede the union by construction; what has to be checked is
    # the exclusion against the fold.
    assert operations.index('exclude') < operations.index('collapse'), (
        f'the exclude runs after the collapse: {operations}'
    )


def test_composition_component_can_name_another_preset(conn, registry):
    """A component is a parameter set **or** the name of another preset.

    That is where the per-component override comes from: replacing one named
    component leaves the other two exactly as the recipe wrote them.
    """
    nichenet = _full_preset(
        name='_roundtrip_preset_named_components',
        composition=NICHENET_COMPOSITION,
    )
    register_network(conn, nichenet, registry_schema=registry)
    stored = _preset_row(conn, registry, nichenet.name)[11]
    assert _component_names(stored) == [
        'curated_ligand_receptor',
        'omnipath',
        'collectri',
    ]
    # The override: one named component swapped, the rest of the recipe alone.
    overridden = {
        **NICHENET_COMPOSITION,
        'components': [
            NICHENET_COMPOSITION['components'][0],
            NICHENET_COMPOSITION['components'][1],
            {'preset': 'dorothea'},
        ],
    }
    register_network(
        conn,
        _full_preset(name=nichenet.name, composition=overridden),
        registry_schema=registry,
    )
    stored = _preset_row(conn, registry, nichenet.name)[11]
    assert _component_names(stored) == [
        'curated_ligand_receptor',
        'omnipath',
        'dorothea',
    ]
    assert _step_operations(stored) == ['collapse'], (
        'the override changed the rest of the recipe'
    )


def test_preset_without_composition_stores_null(conn, registry):
    """One parameter set is the common case, and it stores NULL.

    NULL and not an empty recipe: a reader must be able to tell "no composition"
    from "a composition with no components" without guessing.
    """
    single = _full_preset(name='_roundtrip_preset_no_composition')
    register_network(conn, single, registry_schema=registry)
    assert _preset_row(conn, registry, single.name)[11] is None


def test_live_registry_carries_the_preset_columns(conn):
    """The build database's own ``network_registry`` holds the preset columns."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT column_name, udt_name, is_nullable
            FROM information_schema.columns
            WHERE table_schema = %s AND table_name = 'network_registry'
            """,
            [SCHEMA],
        )
        columns = {name: (udt, nullable) for name, udt, nullable in cur.fetchall()}
    for column, udt in PRESET_COLUMNS.items():
        assert column in columns, f'{SCHEMA}.network_registry lacks {column}'
        assert columns[column][0] == udt, f'{column} is {columns[column][0]}, want {udt}'
    # A license restriction is optional; its absence is NULL, not a level of 0.
    assert columns['license_scope'][1] == 'YES'
    # So is a composition: NULL is the single-parameter-set preset, which
    # is what a database written before the amendment holds in every row.
    assert columns['composition'][1] == 'YES'
    # Matview-era columns: nullable through the transition, dropped once the
    # last bespoke matview retires.
    assert columns['schema_name'][1] == 'YES'
    assert columns['combined_relation'][1] == 'YES'


# ── adding a dataset is a declarative change ────────────────────────────────
#
# The claim: a new dataset costs a definition and nothing else. No view is
# created for it, no column is added to the record for it, and it is queryable
# the moment it is registered. Every part of that is easy to believe and easy
# to lose — the first bespoke view added "just for this one dataset" is how the
# fifteen MetaLinksDB views started. So the invariants are measured across a
# registration rather than trusted.


def _matviews(conn):
    """Every materialized view in the database, as `schema.name` strings."""
    with conn.cursor() as cur:
        cur.execute('SELECT schemaname, matviewname FROM pg_matviews')
        return {f'{schema}.{name}' for schema, name in cur.fetchall()}


def _record_columns(conn, schema):
    """The column names of the interaction record, in a set."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT column_name FROM information_schema.columns
            WHERE table_schema = %s AND table_name = 'interaction_fact_resource'
            """,
            [schema],
        )
        return {name for (name,) in cur.fetchall()}


def test_registering_a_dataset_creates_no_view(conn, registry):
    """A dataset is a row. Nothing is materialized on its behalf."""
    before = _matviews(conn)
    register_network(
        conn,
        _full_preset(name='_declarative_preset'),
        registry_schema=registry,
    )
    after = _matviews(conn)
    assert after == before, (
        f'registering a preset created {sorted(after - before)}; a dataset '
        f'that needs a view of its own is the mechanism this cycle retires'
    )


def test_registering_a_dataset_changes_no_record_column(conn, registry):
    """The record's shape is a property of the build, not of any dataset."""
    before = _record_columns(conn, SCHEMA)
    if not before:
        pytest.skip('the build database holds no interaction record to compare')
    register_network(
        conn,
        _full_preset(name='_declarative_preset_columns'),
        registry_schema=registry,
    )
    after = _record_columns(conn, SCHEMA)
    assert after == before, (
        f'registering a preset changed the record columns by '
        f'{sorted(after ^ before)}; a dataset that widens the fact table makes '
        f'every other dataset pay for it'
    )


def test_a_registered_dataset_is_enumerable_and_scoped(conn, registry):
    """Registration alone makes a dataset discoverable with its full spec."""
    definition = _full_preset(name='_declarative_preset_enumerable')
    register_network(conn, definition, registry_schema=registry)
    with conn.cursor() as cur:
        cur.execute(
            f'SELECT name, kind, included_sources, interaction_class_scope, '
            f'default_attributes, mandatory_attributes '
            f'FROM {registry}.network_registry WHERE name = %s',
            [definition.name],
        )
        row = cur.fetchone()
    assert row is not None, 'the registered dataset is not enumerable'
    name, kind, sources, classes, defaults, mandatory = row
    assert kind == definition.kind
    assert tuple(sources) == definition.included_sources
    assert tuple(classes) == definition.interaction_class_scope
    assert tuple(defaults) == definition.default_attributes
    assert tuple(mandatory) == definition.mandatory_attributes


def test_a_registered_dataset_resolves_to_record_rows(conn, records):
    """Queryable, not only listed: its scope selects rows of the record."""
    definition = _full_preset(
        name='_declarative_preset_queryable',
        included_sources=('signor', 'cellphonedb'),
        license_scope=None,
    )
    register_network(conn, definition, registry_schema=records)
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT count(*)
            FROM {records}.interaction_fact_resource f
            JOIN {records}.data_source d USING (source_id)
            WHERE d.name = ANY(
              SELECT unnest(included_sources)
              FROM {records}.network_registry WHERE name = %s
            )
            """,
            [definition.name],
        )
        selected = cur.fetchone()[0]
    assert selected, (
        'the registered dataset selects no record row through its own scope; '
        'enumerable and queryable are one claim, and this is the second half'
    )
