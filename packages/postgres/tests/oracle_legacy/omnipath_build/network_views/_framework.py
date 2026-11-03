"""Declarative network-view framework, carrying presets since cycle 008.

A dataset is described by a :class:`NetworkDefinition` and registered as one
``network_registry`` row. Two generations of definition share that table:

* the **preset** (cycle 008) — metadata over the interaction fact table: which
  sources contribute, which interaction classes and evidence it scopes to, which
  attributes it returns, how it is labelled and curated. It materialises nothing;
  registering it is the whole build step.
* the **matview network** — metadata plus the curated SQL that
  materialises per-source and combined matviews under its own schema.

The framework manages both uniformly — create the schema and apply the SQL where
there is any, refresh the matviews, and upsert the registry row — so a uniform
API discovers and serves either. Adding a dataset is a definition, never bespoke
framework or API code.

**Transitional columns.** ``schema_name`` and ``combined_relation`` describe a
matview and mean nothing for a preset; they stay nullable while both generations
coexist and are dropped when the last bespoke matview retires. They are
recorded here so they are not left behind as silent dead columns.

**Grain amendment (2026-08-20).** A preset also says how it collapses the
per-resource record over its own resource scope (``collapse_mode``) and which
license terms a resource must meet to contribute (``license_scope``). Neither
buys a preset a table of its own: a third column, ``materialize_collapse``, was
proposed and withdrawn the same day, because interaction queries are served from
one precomputed record table rather than from N per-dataset ones, and a derived
table built only when a flag says so is the silently skipped phase constitution
Principle V rules out. There are two interaction tables and both are built
unconditionally — a scoped preset collapses the record at query time.

**Composition amendment (2026-08-21).** A preset is a parameter set or a
composition of them, and never a third thing. A composition names its components
in order and the operation that joins them (``composition``); a component is
either a parameter set or the name of another preset, which is where a
per-component override comes from. Most presets are one parameter set and leave
the column NULL. The algebra runs in the api-service, not here — the
build stores the recipe and nothing executes it.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

from psycopg2 import sql
from psycopg2.extras import Json
import psycopg2.extensions

_SQL_DIR = Path(__file__).with_name('sql')

logger = logging.getLogger(__name__)

# The matview-era registry columns. A preset leaves them NULL; they are dropped
# once the last bespoke matview retires — see the module docstring.
MATVIEW_ERA_COLUMNS = ('schema_name', 'combined_relation')

# How a preset folds the per-resource record over its own resource scope, as
# the `network_registry.collapse_mode` column records it. `endpoints` is the
# default because it reproduces the legacy one-row-per-interaction contract,
# so a preset written before the amendment keeps behaving as it did.
COLLAPSE_MODES = ('none', 'assertion', 'endpoints')
DEFAULT_COLLAPSE_MODE = 'endpoints'

# What a preset's groups are keyed on, as the `network_registry.grain` column
# records it. `interaction` — the default — groups on the binary triple of
# subject entity, object entity and interaction class, so a group is an ordered
# pair. `participant` groups on the `interaction_id` header alone, so one group
# is one interaction whatever its arity and its members come back as a
# participant list rather than as a first and a second endpoint.
#
# **The grain is read before the collapse mode, and it can silence it.** Each
# of the three collapse modes describes a fold over an ordered endpoint pair,
# and a reaction is not one, so at `participant` grain nothing reads
# `collapse_mode` at all and whatever it holds makes no claim. A preset at that
# grain records `none` there — the closest available value, which states no
# fold rather than the wrong one. This is a reading rule and not a constraint:
# no CHECK enforces it across the two columns.
GRAINS = ('interaction', 'participant')
DEFAULT_GRAIN = 'interaction'

# The operations a composition joins its components with. `union` is the
# operation of the composition itself; the rest are the steps that follow it.
COMPOSITION_OPERATIONS = ('union', 'collapse', 'exclude', 'annotate')

AMENDMENT_COLUMN_COMMENTS = {
    'collapse_mode': (
        'How this preset folds the per-resource interaction record over its '
        "own resource scope (cycle 008): 'none' keeps one row per "
        "resource assertion, 'assertion' folds resources agreeing on sign and "
        "direction, 'endpoints' folds to the collapsed key. Default "
        "'endpoints', the legacy one-row-per-interaction contract. No preset "
        'gets a materialisation of its own: the scoped collapse happens at '
        'query time.'
    ),
    'grain': (
        'What this preset groups on (cycle 008): '
        "'interaction' — the default — groups on the binary triple of subject "
        "entity, object entity and interaction class, and 'participant' groups "
        'on the interaction_id header alone, so one group is one interaction '
        'whatever its arity and its members are returned as a participant list '
        'rather than as two endpoints. Read this column before collapse_mode: '
        'each collapse mode folds an ordered endpoint pair and a reaction is '
        "not one, so at 'participant' grain nothing reads collapse_mode and "
        "its value makes no claim. Such a preset records 'none' there, which "
        'states no fold rather than the wrong one. A reading rule, not a '
        'constraint: no CHECK ties the two columns together.'
    ),
    'license_scope': (
        'Minimum purpose/sharing/attrib levels a resource must meet to '
        'contribute to this preset (cycle 008). NULL is no license '
        'restriction. A resource whose license is unknown is excluded, never '
        'admitted under a permissive default.'
    ),
    'composition': (
        'For a preset that is not one query (cycle 008): the ordered '
        "component list and the operation joining them ('union'), followed by "
        "the ordered steps ('exclude', 'collapse', 'annotate'). A component is "
        'a parameter set or the name of another preset, which is what gives a '
        'per-component override. NULL means one parameter set, the common '
        'case. Two orders are binding: the collapse runs after the '
        "union and over the union's own resolved scope, and the exclude runs "
        'before the collapse, or the dropped resource stays inside '
        'source_count, references and the sign flags.'
    ),
}


@dataclass(frozen=True)
class NetworkDefinition:
    """A dataset: a preset over the fact table, or a matview-backed network.

    The preset fields — the ``network_registry`` columns — describe *what a
    query for this dataset selects and returns*; they carry no SQL:

    ``interaction_class_scope``
        Interaction-class slugs from ``vocab_interaction_class``
        (``signaling``, ``tf_target``, ``ligand_receptor``, …) the preset
        restricts to; empty means all classes.
        A class is derived from the resource annotations, never from a legacy
        dataset name — a legacy name is preset identity, not a class.
    ``evidence_scope``
        The evidence-type / predicate / confidence filter separating the
        legacy datasets that share a class.
    ``default_attributes`` / ``mandatory_attributes``
        Returned when the caller asks for nothing, and returned always even
        when unrequested.
    ``labels``
        Display labels for the preset and its columns.
    ``curation``
        Configurable thresholds and flags — MoA-only, affinity cut-off,
        metabolite-class gate — as config, never inline SQL.
    ``attribute_sources``
        Which source supplies each mandatory attribute, carrying the
        interim-vs-Intercell provenance.
    ``grain``
        What the preset's groups are keyed on: ``interaction`` — the default —
        groups on the binary triple of subject entity, object entity and
        interaction class, so a group is an ordered pair; ``participant``
        groups on the ``interaction_id`` header alone, so one group is one
        interaction whatever its arity and its members come back as a
        participant list instead of a first and a second endpoint. It is read
        before ``collapse_mode`` and can silence it — see that field.
    ``collapse_mode``
        How the preset folds the per-resource record over *its own* resource
        scope: ``none`` keeps one row per resource assertion, ``assertion``
        folds the resources that agree on sign and direction, and ``endpoints``
        — the default — folds to the collapsed key, the legacy
        one-row-per-interaction contract. A single-resource preset collapses
        nothing whatever the mode says, because every group holds one row. At
        ``participant`` grain nothing reads this field at all, because each of
        the three modes folds an ordered endpoint pair and a reaction is not
        one; a preset at that grain records ``none``, the closest available
        value, which states no fold rather than the wrong one.
    ``license_scope``
        The minimum ``purpose`` / ``sharing`` / ``attrib`` levels a resource
        must meet to contribute. ``None`` is no restriction at all; a scope that
        is set resolves to a resource set before the query runs, and a resource
        whose license is unknown fails it however permissive its recorded
        levels read — never admitted by default.
    ``composition``
        The recipe of a preset that is not one query: the ordered
        ``components`` and the ``operation`` joining them, then the ordered
        ``steps`` that follow. A component is a parameter set or the name of
        another preset — the second form is what lets an override replace one
        component and leave the rest of the recipe alone. ``None`` — the common
        case — means the preset is a single parameter set. The order is part of
        the value: the collapse follows the union and folds over the union's own
        resolved scope, and the exclude precedes the collapse, because a
        resource dropped after the fold still counts towards ``source_count``,
        the references and the sign flags. The api-service executes the algebra;
        the registry only stores it.

    ``schema``, ``combined_relation``, ``matviews`` and ``sql_files`` are the
    matview-era fields: a preset leaves them empty.
    """

    name: str
    kind: str
    included_sources: tuple[str, ...] = ()
    interaction_class_scope: tuple[str, ...] = ()
    evidence_scope: Mapping[str, Any] | None = None
    default_attributes: tuple[str, ...] = ()
    mandatory_attributes: tuple[str, ...] = ()
    labels: Mapping[str, Any] | None = None
    curation: Mapping[str, Any] | None = None
    attribute_sources: Mapping[str, Any] | None = None
    # The grain amendment. Both default, so no existing definition
    # changes: no mode named is the legacy collapse, no scope named is no
    # license restriction.
    collapse_mode: str = DEFAULT_COLLAPSE_MODE
    license_scope: Mapping[str, Any] | None = None
    # The key the groups are folded on, defaulted to the binary triple, so a
    # definition written before the grain existed keeps the pair it had.
    grain: str = DEFAULT_GRAIN
    # The composition amendment. Defaulted too: a preset that is one
    # parameter set names no composition at all.
    composition: Mapping[str, Any] | None = None
    # Matview-era fields — retire with the columns above.
    schema: str | None = None
    combined_relation: str | None = None
    matviews: tuple[str, ...] = ()  # refresh order: per-source → combined → annotations
    sql_files: tuple[str, ...] = ()  # applied in order, relative to sql/

    @property
    def is_preset(self) -> bool:
        """True when the dataset is metadata over the fact table (no matview)."""
        return not self.sql_files and not self.matviews

    def sql_text(self) -> str:
        return '\n'.join(
            (_SQL_DIR / name).read_text(encoding='utf-8') for name in self.sql_files
        )


def ensure_network_registry(
    conn: psycopg2.extensions.connection,
    *,
    registry_schema: str = 'public',
) -> None:
    """Create (or extend) the discovery table the network API reads.

    Fresh databases get the full preset table. Databases carrying the older
    matview-only descriptor are migrated in place: the preset columns are added, and
    ``schema_name`` / ``combined_relation`` lose their NOT NULL so a preset can
    leave them empty until they are dropped.

    The amendment's two columns migrate the same way. ``collapse_mode`` arrives
    NOT NULL with the ``endpoints`` default, so rows registered before it keep
    the legacy collapse rather than acquiring an undefined grain, and
    ``license_scope`` arrives nullable, because no license restriction is NULL
    and not a level of zero. So does ``composition``: a database written
    before the amendment holds one parameter set per row, which is exactly what
    a NULL composition says, so the migration rewrites no existing row.

    ``grain`` follows ``collapse_mode``: NOT NULL with the ``interaction``
    default, so every row written before it keeps the binary triple it was
    already grouped on rather than acquiring an undefined key.
    """
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.network_registry (
                  name text PRIMARY KEY,
                  kind text NOT NULL,
                  schema_name text,
                  combined_relation text,
                  included_sources text[] NOT NULL,
                  interaction_class_scope text[],
                  evidence_scope jsonb,
                  default_attributes text[],
                  mandatory_attributes text[],
                  labels jsonb,
                  curation jsonb,
                  attribute_sources jsonb,
                  collapse_mode text NOT NULL DEFAULT 'endpoints',
                  license_scope jsonb,
                  composition jsonb,
                  grain text NOT NULL DEFAULT 'interaction',
                  built_at timestamptz NOT NULL DEFAULT now()
                )
                """
            ).format(sql.Identifier(registry_schema))
        )
        cur.execute(
            sql.SQL(
                """
                ALTER TABLE {}.network_registry
                  ADD COLUMN IF NOT EXISTS interaction_class_scope text[],
                  ADD COLUMN IF NOT EXISTS evidence_scope jsonb,
                  ADD COLUMN IF NOT EXISTS default_attributes text[],
                  ADD COLUMN IF NOT EXISTS mandatory_attributes text[],
                  ADD COLUMN IF NOT EXISTS labels jsonb,
                  ADD COLUMN IF NOT EXISTS curation jsonb,
                  ADD COLUMN IF NOT EXISTS attribute_sources jsonb,
                  ADD COLUMN IF NOT EXISTS collapse_mode text NOT NULL
                    DEFAULT 'endpoints',
                  ADD COLUMN IF NOT EXISTS license_scope jsonb,
                  ADD COLUMN IF NOT EXISTS composition jsonb,
                  ADD COLUMN IF NOT EXISTS grain text NOT NULL
                    DEFAULT 'interaction',
                  ALTER COLUMN schema_name DROP NOT NULL,
                  ALTER COLUMN combined_relation DROP NOT NULL
                """
            ).format(sql.Identifier(registry_schema))
        )
        # `COLLAPSE_MODES` names three modes and no more. Dropping the
        # constraint before adding it keeps the statement idempotent and lets
        # the set of modes change with the definition rather than only on a
        # fresh database.
        cur.execute(
            sql.SQL(
                """
                ALTER TABLE {}.network_registry
                  DROP CONSTRAINT IF EXISTS network_registry_collapse_mode_check,
                  ADD CONSTRAINT network_registry_collapse_mode_check
                    CHECK (collapse_mode = ANY({}))
                """
            ).format(
                sql.Identifier(registry_schema),
                sql.Literal(list(COLLAPSE_MODES)),
            )
        )
        # `GRAINS` names two keys and no more, and the same drop-then-add shape
        # keeps the statement idempotent while letting the set of grains change
        # with the definition. Nothing here ties the grain to `collapse_mode`:
        # that a `participant` preset records `none` is how the pair is read,
        # not something the table refuses to store otherwise.
        cur.execute(
            sql.SQL(
                """
                ALTER TABLE {}.network_registry
                  DROP CONSTRAINT IF EXISTS network_registry_grain_check,
                  ADD CONSTRAINT network_registry_grain_check
                    CHECK (grain = ANY({}))
                """
            ).format(
                sql.Identifier(registry_schema),
                sql.Literal(list(GRAINS)),
            )
        )
        # A composition is a named operation over an ordered component list, or
        # it is nothing at all. The same drop-then-add shape, so the set
        # of operations can change with the definition. Only the shape is
        # checked here: what the components mean is the api-service's business.
        cur.execute(
            sql.SQL(
                """
                ALTER TABLE {}.network_registry
                  DROP CONSTRAINT IF EXISTS network_registry_composition_check,
                  ADD CONSTRAINT network_registry_composition_check
                    CHECK (
                      composition IS NULL
                      OR (
                        composition->>'operation' = ANY({})
                        AND jsonb_typeof(composition->'components') = 'array'
                      )
                    )
                """
            ).format(
                sql.Identifier(registry_schema),
                sql.Literal(list(COMPOSITION_OPERATIONS)),
            )
        )
        # The transition is recorded in the database too, so nobody meets these
        # two columns without learning they are on their way out.
        for column in MATVIEW_ERA_COLUMNS:
            cur.execute(
                sql.SQL(
                    'COMMENT ON COLUMN {}.network_registry.{} IS %s'
                ).format(
                    sql.Identifier(registry_schema), sql.Identifier(column)
                ),
                [
                    'Matview-era column: NULL for a preset over the interaction '
                    'fact table. Dropped when the last bespoke matview retires '
                    '(cycle 008).'
                ],
            )
        # The amendment's two columns say the same thing to a reader with psql
        # and no specification open.
        for column, comment in AMENDMENT_COLUMN_COMMENTS.items():
            cur.execute(
                sql.SQL(
                    'COMMENT ON COLUMN {}.network_registry.{} IS %s'
                ).format(
                    sql.Identifier(registry_schema), sql.Identifier(column)
                ),
                [comment],
            )
    conn.commit()


def apply_network(
    conn: psycopg2.extensions.connection,
    definition: NetworkDefinition,
) -> None:
    """Create the network's schema + (re)materialise its views from curated SQL.

    A preset materialises nothing — registration is its whole build step.
    """
    if definition.is_preset:
        logger.debug(
            'network-views: %s is a preset over the fact table; nothing to apply',
            definition.name,
        )
        return
    schema_id = sql.Identifier(definition.schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL('CREATE SCHEMA IF NOT EXISTS {}').format(schema_id)
        )
        # Unqualified objects in the curated SQL land in the network schema;
        # public.* (the canonical graph) still resolves.
        cur.execute(
            sql.SQL('SET search_path = {}, public').format(schema_id)
        )
        cur.execute(definition.sql_text())
        cur.execute('RESET search_path')
    conn.commit()


def refresh_network(
    conn: psycopg2.extensions.connection,
    definition: NetworkDefinition,
) -> None:
    """Refresh the network's matviews in dependency order (per-source → combined)."""
    if definition.is_preset:
        logger.debug(
            'network-views: %s is a preset over the fact table; nothing to refresh',
            definition.name,
        )
        return
    schema_id = sql.Identifier(definition.schema)
    with conn.cursor() as cur:
        for matview in definition.matviews:
            cur.execute(
                sql.SQL('REFRESH MATERIALIZED VIEW {}.{}').format(
                    schema_id, sql.Identifier(matview)
                )
            )
    conn.commit()


def register_network(
    conn: psycopg2.extensions.connection,
    definition: NetworkDefinition,
    *,
    registry_schema: str = 'public',
) -> None:
    """Upsert the dataset's row in ``network_registry`` (stamps ``built_at``).

    The whole preset spec travels with the row, so the API resolves a dataset
    from the registry alone.

    A resource the definition names and the build has not loaded is normal — a
    dataset is registered once and its resources arrive over several cycles,
    and a capped development build loads a fraction of them. It is also
    invisible: the dataset registers, the resource contributes nothing, and an
    empty contribution looks exactly like a resource with nothing to say. So
    the missing names are logged, and logged individually, before the row goes
    in. The scope keeps them: they are what the dataset is waiting for, and
    dropping them would turn a known gap into a forgotten one.
    """
    _warn_unloaded(conn, definition, registry_schema=registry_schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {}.network_registry
                  (name, kind, schema_name, combined_relation, included_sources,
                   interaction_class_scope, evidence_scope, default_attributes,
                   mandatory_attributes, labels, curation, attribute_sources,
                   collapse_mode, license_scope, composition, grain, built_at)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                        %s, %s, now())
                ON CONFLICT (name) DO UPDATE SET
                  kind = EXCLUDED.kind,
                  schema_name = EXCLUDED.schema_name,
                  combined_relation = EXCLUDED.combined_relation,
                  included_sources = EXCLUDED.included_sources,
                  interaction_class_scope = EXCLUDED.interaction_class_scope,
                  evidence_scope = EXCLUDED.evidence_scope,
                  default_attributes = EXCLUDED.default_attributes,
                  mandatory_attributes = EXCLUDED.mandatory_attributes,
                  labels = EXCLUDED.labels,
                  curation = EXCLUDED.curation,
                  attribute_sources = EXCLUDED.attribute_sources,
                  collapse_mode = EXCLUDED.collapse_mode,
                  license_scope = EXCLUDED.license_scope,
                  composition = EXCLUDED.composition,
                  grain = EXCLUDED.grain,
                  built_at = now()
                """
            ).format(sql.Identifier(registry_schema)),
            [
                definition.name,
                definition.kind,
                definition.schema,
                definition.combined_relation,
                list(definition.included_sources),
                list(definition.interaction_class_scope) or None,
                Json(definition.evidence_scope)
                if definition.evidence_scope is not None
                else None,
                list(definition.default_attributes),
                list(definition.mandatory_attributes),
                Json(definition.labels) if definition.labels is not None else None,
                Json(definition.curation) if definition.curation is not None else None,
                Json(definition.attribute_sources)
                if definition.attribute_sources is not None
                else None,
                definition.collapse_mode,
                Json(definition.license_scope)
                if definition.license_scope is not None
                else None,
                Json(definition.composition)
                if definition.composition is not None
                else None,
                definition.grain,
            ],
        )
    conn.commit()


def _warn_unloaded(
    conn: psycopg2.extensions.connection,
    definition: NetworkDefinition,
    *,
    registry_schema: str = 'public',
) -> tuple[str, ...]:
    """Log every resource this dataset names that the build has not loaded.

    Args:
        conn: An open connection.
        definition: The dataset being registered.
        registry_schema: The schema holding ``data_source``.

    Returns:
        The unloaded resource names, sorted. Empty when every one is present,
        and empty too when the schema carries no ``data_source`` at all — a
        database without one cannot say what is loaded, and guessing there
        would produce a warning for every resource of every dataset.
    """
    named = list(definition.included_sources)
    if not named:
        return ()
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL('SELECT to_regclass({})').format(
                sql.Literal(f'{registry_schema}.data_source')
            )
        )
        if cur.fetchone()[0] is None:
            return ()
        cur.execute(
            sql.SQL(
                'SELECT name FROM {}.data_source WHERE name = ANY(%s)'
            ).format(sql.Identifier(registry_schema)),
            [named],
        )
        loaded = {name for (name,) in cur.fetchall()}
    missing = tuple(sorted(set(named) - loaded))
    for name in missing:
        logger.warning(
            'the dataset %r names the resource %r, which this build has not '
            'loaded; it contributes no row and the dataset is short of '
            'whatever it would have contributed',
            definition.name,
            name,
        )
    return missing


@dataclass(frozen=True)
class NetworkViewStats:
    applied: tuple[str, ...] = field(default_factory=tuple)


def apply_all(
    conn: psycopg2.extensions.connection,
    definitions: list[NetworkDefinition],
    *,
    registry_schema: str = 'public',
    log=lambda *_: None,
) -> NetworkViewStats:
    """Apply + register every dataset (the build hook). Idempotent."""
    ensure_network_registry(conn, registry_schema=registry_schema)
    applied: list[str] = []
    for definition in definitions:
        log(f'[network-views] apply {definition.name}')
        apply_network(conn, definition)
        register_network(conn, definition, registry_schema=registry_schema)
        applied.append(definition.name)
    return NetworkViewStats(applied=tuple(applied))


def refresh_all(
    conn: psycopg2.extensions.connection,
    definitions: list[NetworkDefinition],
    *,
    registry_schema: str = 'public',
    log=lambda *_: None,
) -> NetworkViewStats:
    """Refresh + re-register every dataset (fast path; views already exist)."""
    refreshed: list[str] = []
    for definition in definitions:
        log(f'[network-views] refresh {definition.name}')
        refresh_network(conn, definition)
        register_network(conn, definition, registry_schema=registry_schema)
        refreshed.append(definition.name)
    return NetworkViewStats(applied=tuple(refreshed))
