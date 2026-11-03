"""Maintain the API-facing resource summary table.

The resource table joins static pypath resource configuration with the content
that is actually present in PostgreSQL. Counts can come either from direct
evidence/graph tables or from bitmap tables when they have been refreshed.
"""

from __future__ import annotations

from typing import Any, Mapping, NamedTuple, Sequence
from pathlib import Path
from functools import lru_cache
import json
import hashlib
import logging
import subprocess
from dataclasses import dataclass
import importlib.util

from psycopg2 import sql
from psycopg2.extras import Json
import psycopg2.extensions

from omnipath_build.cv_terms import CV_TERM_ENTITY_TYPE

logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class ResourceTableStats:
    """Summary counts from resource metadata sync."""

    resources: int = 0


@dataclass(frozen=True)
class BuildManifestStats:
    """Identity of the manifest written at the end of a build."""

    build_id: str = ''
    partial_build: bool = False
    resources: int = 0


def emit_build_manifest(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    inputs_package: str = 'pypath.inputs_v2',
    partial_build: bool = False,
    derive_cost: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None = None,
    scope_cost: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None = None,
    deferral_cost: Mapping[str, Any] | None = None,
    utils_db_url: str | None = None,
) -> BuildManifestStats:
    """Write the single self-describing ``build_manifest`` row.

    ``build_id`` is the SHA-256 of the canonical sorted-key JSON of
    ``{package_commits, resources}`` (12 hex chars) — reproducible for identical
    content, independent of ``built_at``. One current row per build (TRUNCATE +
    INSERT). Projects per-resource provenance from the ``resources`` table, so it
    runs in ``derive`` right after ``sync_resources_table``.

    ``derive_cost`` is what the interaction derive steps cost this run — seconds
    and rows per step (record fact, assay, party, reaction projection,
    intercell), either as ``{step: {seconds, rows}}`` or as a sequence of such
    records. The **interaction table is named apart** under
    ``interaction_tables`` as ``record``, so the cost the derive-time budget is
    argued against is readable without knowing which step wrote the table.

    ``scope_cost`` is what materialising each query scope cost — seconds, rows
    and whether it was materialised at all, per scope — either as
    ``{scope: {...}}`` or as a sequence of ``{scope, ...}`` records. It is
    recorded for **every** scope measured, so that declining to materialise a
    scope stays an available answer to a cost overrun.

    ``deferral_cost`` is what deferring the foreign keys and the secondary
    indexes over the load bought: the seconds it saved, and the seconds the
    drop, the restore and the revalidation cost — see
    :func:`_interactions_deferral_cost`. All three are optional: a build that
    recorded none of them still emits a manifest.

    The cost is written to ``interactions_derive_cost``, the deferral to
    ``interactions_deferral_cost``, and both, with the ``network_presets``
    inventory read from ``network_registry``, stay **outside** the hashed
    payload: all of it is volatile — timings vary run to run and the preset list
    changes without the built content changing — and none of it may move a
    build's identity (the cycle-007 lesson). The hash covers
    ``{package_commits, resources}`` and nothing else.
    """

    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.build_manifest (
                  build_id text PRIMARY KEY,
                  built_at timestamptz NOT NULL DEFAULT now(),
                  package_commits jsonb NOT NULL,
                  resources jsonb NOT NULL,
                  partial_build boolean NOT NULL
                )
                """
            ).format(schema_id)
        )
        # Descriptive columns, every one of them deliberately outside the
        # build_id content hash: `translation_tables` names the
        # identifier-resolution relations the build consumed,
        # `canonicalization_coverage` counts how many records resolved vs stayed
        # unresolved (and why), `interactions_derive_cost` holds this
        # run's per-step seconds and rows, the interaction table named apart and
        # the per-scope materialisation cost, `interactions_deferral_cost` what
        # deferring the constraints over the load saved, and `network_presets`
        # the registered preset inventory — volatile numbers and metadata that must not change a
        # build's identity. They are added as columns rather than folded into
        # the hashed payload, which is what keeps them out of the hash by
        # construction: `payload` below is built from two names only.
        cur.execute(
            sql.SQL(
                """
                ALTER TABLE {}.build_manifest
                  ADD COLUMN IF NOT EXISTS translation_tables jsonb,
                  ADD COLUMN IF NOT EXISTS canonicalization_coverage jsonb,
                  ADD COLUMN IF NOT EXISTS interactions_derive_cost jsonb,
                  ADD COLUMN IF NOT EXISTS interactions_deferral_cost jsonb,
                  ADD COLUMN IF NOT EXISTS network_presets jsonb,
                  ADD COLUMN IF NOT EXISTS capabilities jsonb,
                  ADD COLUMN IF NOT EXISTS chemical_resolution_coverage jsonb
                """
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                SELECT
                  resource_id,
                  (entity_count + interaction_count + association_count
                   + identifier_count + ontology_term_count) AS record_count,
                  input_module_commit,
                  input_module_dirty
                FROM {}.resources
                ORDER BY resource_id
                """
            ).format(schema_id)
        )
        resources = [
            {
                'name': name,
                'record_count': int(record_count),
                'version': None,
                'input_module_commit': commit,
                'input_module_dirty': bool(dirty),
            }
            for name, record_count, commit, dirty in cur.fetchall()
        ]
        package_commits = {
            'omnipath_build': _package_commit('omnipath_build'),
            'omnipath_resources': _package_commit(inputs_package.split('.')[0]),
        }
        payload = {'package_commits': package_commits, 'resources': resources}
        build_id = hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()
        ).hexdigest()[:12]
        translation_tables = _translation_tables_used(cur, schema)
        coverage = _canonicalization_coverage(cur, schema)
        interactions_derive_cost = _interactions_derive_cost(derive_cost, scope_cost)
        # The sign-conflict measurement — how often a collapsed row hides
        # resources disagreeing on sign — rides next to the cost, read from
        # what the derive step recorded rather than passed in, so it reaches
        # the manifest without the orchestration having to forward it. It is a
        # measurement of the build, never part of its identity, so it lands in
        # the same non-hashed column.
        sign_conflict = _interaction_sign_conflict(cur, schema)
        if sign_conflict is not None:
            interactions_derive_cost = {
                **(interactions_derive_cost or {}),
                'sign_conflict': sign_conflict,
            }
        interactions_deferral_cost = _interactions_deferral_cost(deferral_cost)
        network_presets = _network_preset_inventory(cur, schema)
        capabilities = _build_capabilities(cur, schema, utils_db_url=utils_db_url)
        _populate_build_capability_table(cur, schema, capabilities)
        chemical_resolution_coverage = _chemical_resolution_coverage_manifest(
            cur, schema,
        )
        cur.execute(sql.SQL('TRUNCATE {}.build_manifest').format(schema_id))
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {}.build_manifest
                  (build_id, package_commits, resources, partial_build,
                   translation_tables, canonicalization_coverage,
                   interactions_derive_cost, interactions_deferral_cost,
                   network_presets, capabilities, chemical_resolution_coverage)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                """
            ).format(schema_id),
            [
                build_id,
                Json(package_commits),
                Json(resources),
                partial_build,
                Json(translation_tables),
                Json(coverage),
                Json(interactions_derive_cost)
                if interactions_derive_cost is not None
                else None,
                Json(interactions_deferral_cost)
                if interactions_deferral_cost is not None
                else None,
                Json(network_presets) if network_presets is not None else None,
                Json(capabilities),
                Json(chemical_resolution_coverage)
                if chemical_resolution_coverage is not None
                else None,
            ],
        )
    conn.commit()
    return BuildManifestStats(
        build_id=build_id,
        partial_build=partial_build,
        resources=len(resources),
    )


# The identifier-resolution relations the build reads to canonicalise evidence.
# Kept here (rather than imported from the build pipeline) so the manifest layer
# stays free of the pipeline's heavy imports.
_TRANSLATION_TABLES = (
    'resolver_gene',
    'resolver_gene_protein_global',
    'resolver_protein',
    'resolver_chemical',
)


def _translation_tables_used(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> list[dict[str, Any]]:
    """Identifier-resolution relations consulted while canonicalising evidence.

    These relations live in the separate identifier-resolution database the build
    reads from, so only their names are recorded here; whether each was present
    and usable at build time is reported by the build's resolver pre-flight (and
    its warnings) in the build log.
    """

    return [{'name': name} for name in _TRANSLATION_TABLES]


def _canonicalization_coverage(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> dict[str, Any] | None:
    """How completely evidence was canonicalised.

    Counts how many evidence records resolved to an entity versus stayed
    unresolved, and — for the unresolved ones — why. Returns ``None`` if the
    resolution table is not present. A single grouped scan run once per build;
    the counts are descriptive and excluded from the build's identity hash.
    """

    schema_id = sql.Identifier(schema)
    # A savepoint, not a bare rollback: the manifest's own DDL was issued in this
    # same transaction, and a missing resolution table must not undo it.
    cur.execute('SAVEPOINT canonicalization_coverage_probe')
    try:
        cur.execute(
            sql.SQL(
                """
                SELECT rs.name AS status, rr.name AS reason, count(*) AS n
                FROM {0}.entity_evidence_resolution eer
                JOIN {0}.vocab_resolution_status rs
                  ON rs.resolution_status_id = eer.status_id
                LEFT JOIN {0}.vocab_resolution_reason rr
                  ON rr.resolution_reason_id = eer.reason_id
                GROUP BY rs.name, rr.name
                """
            ).format(schema_id)
        )
    except psycopg2.Error:
        cur.execute('ROLLBACK TO SAVEPOINT canonicalization_coverage_probe')
        logger.debug('build manifest: no resolution table in %s; coverage skipped', schema)
        return None
    # Read the rows out before releasing the savepoint: executing the RELEASE
    # replaces the cursor's result set, and fetching afterwards raises "no
    # results to fetch" on every schema that actually has a resolution table.
    rows = cur.fetchall()
    cur.execute('RELEASE SAVEPOINT canonicalization_coverage_probe')

    by_status: dict[str, int] = {}
    by_reason: dict[str, int] = {}
    total = 0
    for status, reason, n in rows:
        n = int(n)
        total += n
        by_status[status] = by_status.get(status, 0) + n
        if reason is not None:
            by_reason[reason] = by_reason.get(reason, 0) + n
    return {'total': total, 'by_status': by_status, 'by_reason': by_reason}


def _chemical_resolution_coverage_manifest(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> list[dict[str, Any]] | None:
    """spec 011 T128: the per (resource, namespace, role) coverage
    :func:`derived_tables._populate_chemical_resolution_coverage` already
    wrote, read into the manifest as-is -- so two builds compare from
    their manifests alone (data-model.md section 9), without re-running
    analysis. Returns ``None`` if the coverage table is not present (an
    older build, or one where derive has not run yet).
    """

    schema_id = sql.Identifier(schema)
    cur.execute('SAVEPOINT chemical_resolution_coverage_probe')
    try:
        cur.execute(
            sql.SQL(
                """
                SELECT
                  crc.source_id, ds.name AS source, crc.identifier_type_id,
                  vit.name AS namespace, crc.role, crc.mentions, crc.entities,
                  crc.reached_structure, crc.reached_name, crc.unresolved,
                  crc.conflicted
                FROM {0}.chemical_resolution_coverage crc
                JOIN {0}.data_source ds ON ds.source_id = crc.source_id
                JOIN {0}.vocab_identifier_type vit
                  ON vit.identifier_type_id = crc.identifier_type_id
                ORDER BY ds.name, vit.name, crc.role
                """
            ).format(schema_id)
        )
    except psycopg2.Error:
        cur.execute('ROLLBACK TO SAVEPOINT chemical_resolution_coverage_probe')
        logger.debug(
            'build manifest: no chemical_resolution_coverage table in %s; skipped',
            schema,
        )
        return None
    rows = cur.fetchall()
    cur.execute('RELEASE SAVEPOINT chemical_resolution_coverage_probe')

    return [
        {
            'source_id': source_id,
            'source': source,
            'identifier_type_id': identifier_type_id,
            'namespace': namespace,
            'role': role,
            'mentions': int(mentions),
            'entities': int(entities),
            'reached_structure': int(reached_structure),
            'reached_name': int(reached_name),
            'unresolved': int(unresolved),
            'conflicted': int(conflicted),
        }
        for (
            source_id, source, identifier_type_id, namespace, role,
            mentions, entities, reached_structure, reached_name,
            unresolved, conflicted,
        ) in rows
    ]


def _build_capabilities(
    cur: psycopg2.extensions.cursor,
    schema: str,
    *,
    utils_db_url: str | None,
) -> list[dict[str, Any]]:
    """What the build could do this run, and why not when it could not
    (spec 011 data-model section 8, research R14).

    Absence never changes the identity of an entity that resolves. It
    only reduces coverage or switches off a feature -- and says so,
    rather than failing the build or quietly returning less.
    """

    return [
        _detect_structure_substrate(cur, schema),
        _detect_in_database_structure_key(cur),
        _detect_lipid_nomenclature(cur, schema),
        _detect_structure_key_computation(utils_db_url),
    ]


def _detect_structure_substrate(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> dict[str, Any]:
    """The metabo phase's RDKit structure substrate.

    Checked by row count, not just table existence. Schema setup creates
    the table up front. Only the metabo phase populates it, once it
    actually runs against this build.
    """

    from omnipath_build.db.schema import _table_exists

    provider = 'metabo phase'
    if not _table_exists(cur, schema, 'metabo_entity_structure'):
        return {
            'capability': 'structure_substrate', 'available': False,
            'provider': provider,
            'reason': 'metabo_entity_structure does not exist',
        }
    cur.execute(
        sql.SQL('SELECT count(*) FROM {}.metabo_entity_structure')
        .format(sql.Identifier(schema))
    )
    populated = int(cur.fetchone()[0]) > 0
    return {
        'capability': 'structure_substrate', 'available': populated,
        'provider': provider,
        'reason': None if populated else (
            'metabo_entity_structure is empty. The metabo phase has not '
            'run against this build'
        ),
    }


def _detect_in_database_structure_key(
    cur: psycopg2.extensions.cursor,
) -> dict[str, Any]:
    """The rdkit Postgres cartridge, for ad-hoc in-database structure
    queries. A convenience, not a prerequisite (plan.md). Structure
    computation itself runs through the library binding -- which the
    metabo phase already uses -- whether or not this extension is present.
    """

    cur.execute("SELECT EXISTS (SELECT 1 FROM pg_extension WHERE extname = 'rdkit')")
    available = bool(cur.fetchone()[0])
    return {
        'capability': 'in_database_structure_key', 'available': available,
        'provider': 'the container image',
        'reason': None if available else (
            'the rdkit extension is not installed in this Postgres image'
        ),
    }


def _detect_lipid_nomenclature(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> dict[str, Any]:
    """The nomenclature grammar (an omnipath-utils optional extra) that
    assigns lipid identity. When absent, ``lipid_name_node`` and
    ``lipid_name_edge`` exist (schema setup creates them) but stay empty,
    and lipids fall back to a normalized raw name as identity.
    """

    from omnipath_build.db.schema import _table_exists

    provider = 'the utils build, nomenclature extra'
    if not _table_exists(cur, schema, 'lipid_name_node'):
        return {
            'capability': 'lipid_nomenclature', 'available': False,
            'provider': provider,
            'reason': 'lipid_name_node does not exist',
        }
    cur.execute(
        sql.SQL('SELECT count(*) FROM {}.lipid_name_node')
        .format(sql.Identifier(schema))
    )
    populated = int(cur.fetchone()[0]) > 0
    return {
        'capability': 'lipid_nomenclature', 'available': populated,
        'provider': provider,
        'reason': None if populated else (
            'lipid_name_node is empty. The nomenclature grammar has not '
            'assigned lipid identity for this build'
        ),
    }


def _detect_structure_key_computation(
    utils_db_url: str | None,
) -> dict[str, Any]:
    """The chemistry toolkit (an omnipath-utils optional extra) that
    computes a structure key -- InChIKey, canonical SMILES -- from a raw
    structure. Structure-to-anything translation routes report themselves
    unavailable without it. Identifier-to-identifier translation is
    unaffected.

    The real signal (spec 011 T082) is the utils build's own record of
    whether it ran with the toolkit importable
    (``DatabaseBuilder.record_structure_key_capability``,
    ``build_info`` row ``table_name='capability'``) -- this build never
    imports rdkit itself, only reads what the *utils* build already
    determined about its own environment.
    """

    provider = 'the utils build, chemistry extra'
    if not utils_db_url:
        return {
            'capability': 'structure_key_computation', 'available': False,
            'provider': provider,
            'reason': 'OMNIPATH_BUILD_UTILS_PG_URL is not set',
        }

    import psycopg2

    try:
        conn = psycopg2.connect(utils_db_url)
        try:
            with conn.cursor() as cur:
                cur.execute(
                    "SELECT status FROM omnipath_utils.build_info"
                    " WHERE table_name = 'capability'"
                    " AND source_type = 'structure_key_computation'"
                    " ORDER BY built_at DESC LIMIT 1"
                )
                row = cur.fetchone()
        finally:
            conn.close()
    except Exception as exc:  # noqa: BLE001 -- report unavailable, don't crash
        return {
            'capability': 'structure_key_computation', 'available': False,
            'provider': provider,
            'reason': f'could not reach the utils database: {exc}',
        }

    if row is None:
        return {
            'capability': 'structure_key_computation', 'available': False,
            'provider': provider,
            'reason': (
                'the utils build has never run '
                'DatabaseBuilder.record_structure_key_capability'
            ),
        }

    available = row[0] == 'available'
    return {
        'capability': 'structure_key_computation', 'available': available,
        'provider': provider,
        'reason': (
            None if available
            else 'the utils build ran without the chemistry extra installed'
        ),
    }


def _populate_build_capability_table(
    cur: psycopg2.extensions.cursor,
    schema: str,
    capabilities: list[dict[str, Any]],
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(sql.SQL('TRUNCATE {}.build_capability').format(schema_id))
    cur.executemany(
        sql.SQL(
            """
            INSERT INTO {}.build_capability
              (capability, available, provider, reason)
            VALUES (%s, %s, %s, %s)
            """
        ).format(schema_id).as_string(cur.connection),
        [
            (c['capability'], c['available'], c['provider'], c['reason'])
            for c in capabilities
        ],
    )


# The one interaction table the derive writes: one row per (subject, object,
# class, resource) plus the assertion signature. It is named apart from the
# step list because it is the table the derive-time budget is argued against,
# and a reader should not have to know which step wrote it. The collapse
# `interaction_fact_combined` was removed, so there is no second half to name
# and no all-resources scope to reserve a name for — the fold is what a query
# does now, and it materialises nothing.
INTERACTION_RECORD_STEP = 'interaction_fact_resource'

# The interaction derive steps whose cost the manifest reports, in build order.
_DERIVE_STEPS = (
    'interaction_header',
    'interaction_party',
    INTERACTION_RECORD_STEP,
    'interaction_assay',
    'reaction_projection',
    'intercell',
)


def _interaction_sign_conflict(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> dict[str, Any] | None:
    """The sign-conflict summary the interaction derive step recorded.

    ``None`` when no interaction derive step has run in this schema, so a build
    without one says nothing rather than claiming zero conflicts.
    """
    cur.execute(
        'SELECT to_regclass(%s)',
        [f'{schema}.interaction_sign_conflict'],
    )
    if cur.fetchone()[0] is None:
        return None
    cur.execute(
        sql.SQL(
            """
            SELECT fact_rows, signed_rows, both_flags_rows,
                   both_flags_percent, single_resource_rows,
                   cross_resource_rows
            FROM {}.interaction_sign_conflict
            ORDER BY measured_at DESC
            LIMIT 1
            """
        ).format(sql.Identifier(schema))
    )
    row = cur.fetchone()
    if row is None:
        return None
    names = (
        'fact_rows',
        'signed_rows',
        'both_flags_rows',
        'both_flags_percent',
        'single_resource_rows',
        'cross_resource_rows',
    )
    return dict(zip(names, row))


def _interactions_derive_cost(
    derive_cost: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
    scope_cost: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None = None,
) -> dict[str, Any] | None:
    """Normalise what the interaction derive cost into one manifest record.

    Returns, as available:

    ``steps``
        ``[{step, seconds, rows}]`` ordered by the build order of the steps it
        knows, unknown steps last. Fed from ``{step: {seconds, rows}}`` or a
        sequence of ``{step, seconds, rows}`` records.
    ``interaction_tables``
        the interaction table named apart from the step list — see
        :func:`_interaction_table_cost`.
    ``scopes``
        what each measured query scope cost to materialise — see
        :func:`_scope_materialisation_cost`.

    Returns ``None`` when nothing was reported at all, so a build that ran no
    interaction derive step says so rather than claiming zeros.

    These are timings and row counts — the most volatile thing a build produces.
    They are recorded next to the identity hash and never inside it.
    """

    scopes = _scope_materialisation_cost(scope_cost)
    if not derive_cost:
        return {'scopes': scopes} if scopes else None
    if isinstance(derive_cost, Mapping):
        records = [
            {'step': step, **dict(cost)} for step, cost in derive_cost.items()
        ]
    else:
        records = [dict(cost) for cost in derive_cost]

    steps: list[dict[str, Any]] = []
    for record in records:
        step = record.get('step')
        if step is None:
            logger.warning('build manifest: derive-cost record without a step: %r', record)
            continue
        seconds = record.get('seconds')
        rows = record.get('rows')
        steps.append(
            {
                'step': str(step),
                'seconds': float(seconds) if seconds is not None else None,
                'rows': int(rows) if rows is not None else None,
            }
        )
    if not steps:
        return {'scopes': scopes} if scopes else None

    order = {step: n for n, step in enumerate(_DERIVE_STEPS)}
    steps.sort(key=lambda entry: (order.get(entry['step'], len(order)), entry['step']))
    cost: dict[str, Any] = {'steps': steps}
    tables = _interaction_table_cost(steps)
    if tables is not None:
        cost['interaction_tables'] = tables
    if scopes:
        cost['scopes'] = scopes
    return cost


def _interaction_table_cost(
    steps: Sequence[Mapping[str, Any]],
) -> dict[str, Any] | None:
    """The interaction table's cost, named apart from the step list.

    ``{'record': {table, seconds, rows}}`` — ``record`` is
    ``interaction_fact_resource``, the one table the projection writes. A
    seconds or rows the step did not report stays ``None`` rather than becoming
    a zero: not measured and measured as free are different facts, and the
    ceiling this field will be argued against is decided by which.

    Returns ``None`` when no step wrote the table, so a build that projected no
    interactions says nothing rather than naming a table it does not hold.

    **The collapse half is gone**, together with the ``total_seconds`` that
    summed the two halves, because the table they described was removed. A
    ``collapse: null`` would read as a table that ran unmeasured, which is a
    claim about a build that has no such table.
    """

    entry = next(
        (
            entry
            for entry in steps
            if entry['step'] == INTERACTION_RECORD_STEP
        ),
        None,
    )
    if entry is None:
        return None
    return {
        'record': {
            'table': INTERACTION_RECORD_STEP,
            'seconds': entry.get('seconds'),
            'rows': entry.get('rows'),
        }
    }


def _scope_materialisation_cost(
    scope_cost: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
) -> list[dict[str, Any]] | None:
    """What each query scope cost to materialise.

    Accepts ``{scope: {...}}`` or a sequence of ``{scope, ...}`` records and
    normalises each to ``{scope, table, materialised, seconds, rows, sources}``:

    ``scope``
        the scope's name.
    ``table``
        the relation it was materialised into, ``None`` when it was measured
        without being materialised.
    ``materialised``
        whether this build actually stores it. Inferred from ``table`` when the
        record does not say. A measured-but-not-materialised scope is the whole
        point of the field: declining to materialise has to stay an available
        response to a cost overrun, and that argument needs the number for the
        scope that was declined.
    ``seconds`` / ``rows``
        derive time and row count, ``None`` when not measured.
    ``sources``
        the resource set the scope resolves to, so a later reader can tell which
        scope this was without re-resolving the preset or license filter.

    Ordered materialised scopes first, then measured ones, each by name.
    Returns ``None`` when nothing was reported, so a build that measured no
    scope says so rather than claiming it materialised none.
    """

    if not scope_cost:
        return None
    if isinstance(scope_cost, Mapping):
        records = [
            {'scope': scope, **dict(cost)} for scope, cost in scope_cost.items()
        ]
    else:
        records = [dict(cost) for cost in scope_cost]

    scopes: list[dict[str, Any]] = []
    for record in records:
        scope = record.get('scope')
        if scope is None:
            logger.warning('build manifest: scope-cost record without a scope: %r', record)
            continue
        table = record.get('table')
        materialised = record.get('materialised')
        seconds = record.get('seconds')
        rows = record.get('rows')
        sources = record.get('sources')
        scopes.append(
            {
                'scope': str(scope),
                'table': str(table) if table is not None else None,
                # Absent means "as the table says": a scope written somewhere
                # was materialised, one measured into nothing was not.
                'materialised': bool(materialised)
                if materialised is not None
                else table is not None,
                'seconds': float(seconds) if seconds is not None else None,
                'rows': int(rows) if rows is not None else None,
                'sources': [str(source) for source in sources]
                if sources is not None
                else None,
            }
        )
    if not scopes:
        return None

    scopes.sort(key=lambda entry: (not entry['materialised'], entry['scope']))
    return scopes


# What a deferred load reports, read back by kind: seconds as floats, object
# counts as ints, the two claims as booleans. The names are the contract the
# derive's deferred load fills; a field it did not measure is simply absent
# from what it hands over.
_DEFERRAL_SECONDS = (
    'seconds_saved',
    'drop_seconds',
    'restore_seconds',
    'revalidate_seconds',
    # ``load_seconds`` is this build's own load window, and it is recorded so a
    # later build can subtract it. ``seconds_saved`` is a difference against a
    # build that ran without the deferral, which is not a thing any single build
    # can measure about itself: the first deferred build reports it as null and
    # the next one reads this number off the manifest. Dropping it would leave
    # the baseline recoverable only by parsing the step list.
    'load_seconds',
)
_DEFERRAL_COUNTS = ('constraints_deferred', 'indexes_deferred')
_DEFERRAL_FLAGS = ('deferred', 'catalogue_unchanged')


def _interactions_deferral_cost(
    deferral_cost: Mapping[str, Any] | None,
) -> dict[str, Any] | None:
    """What deferring the constraints over the interaction load bought.

    The deferral is taken: the load runs with the 13 foreign keys and the 18
    secondary indexes off and puts them back validated, at 709.7 s against
    1,814.7 s, with revalidation set-based at 44.6 s against 726.3 s of per-row
    triggers. Recording it per build is what makes a regression in the deferral
    show up in the manifest rather than in somebody's stopwatch.

    Normalises one record to::

        {deferred, seconds_saved, drop_seconds, restore_seconds,
         revalidate_seconds, constraints_deferred, indexes_deferred,
         catalogue_unchanged}

    ``seconds_saved``
        the load without the deferral less the load with it. The producer
        computes it; the manifest records it, because the baseline is not a
        thing this build ran.
    ``drop_seconds`` / ``restore_seconds`` / ``revalidate_seconds``
        the three halves of the mechanism's own cost, kept apart from the load
        in ``build_manifest`` — the revalidation is the one that could quietly
        be skipped, and a deferral that returns a `NOT VALID` key has saved
        nothing.
    ``catalogue_unchanged``
        whether the catalogue on the far side matched the near side.

    Every field is ``None`` when the producer did not report it, and the whole
    record is ``None`` when it reported nothing readable — **a build that ran
    without the deferral says nothing rather than claiming zero**. Zero seconds
    saved is a measurement, and it is a regression to chase; no deferral at all
    is not, and the two must not read alike.

    A field that cannot be read as its kind is warned about and dropped, never
    raised: a malformed number from the derive must not cost the build its
    manifest. ``deferred`` is inferred from the presence of a measurement when
    the record does not say, the way ``materialised`` is inferred from
    ``table``.
    """

    if not deferral_cost:
        return None

    record: dict[str, Any] = {}

    def read(field: str, cast: Any) -> Any:
        value = deferral_cost.get(field)
        if value is None:
            return None
        try:
            return cast(value)
        except (TypeError, ValueError):
            logger.warning(
                'build manifest: unreadable deferral-cost field %s: %r',
                field,
                value,
            )
            return None

    for field in _DEFERRAL_SECONDS:
        record[field] = read(field, float)
    for field in _DEFERRAL_COUNTS:
        record[field] = read(field, int)
    for field in _DEFERRAL_FLAGS:
        record[field] = read(field, bool)

    measured = [
        value
        for field, value in record.items()
        if field != 'deferred' and value is not None
    ]
    if not measured and record['deferred'] is None:
        logger.warning(
            'build manifest: deferral-cost record with nothing readable in it: %r',
            deferral_cost,
        )
        return None
    if record['deferred'] is None:
        # Something was measured, so something was deferred.
        record['deferred'] = True
    return {
        'deferred': record['deferred'],
        **{field: record[field] for field in _DEFERRAL_SECONDS},
        **{field: record[field] for field in _DEFERRAL_COUNTS},
        'catalogue_unchanged': record['catalogue_unchanged'],
    }


def _network_preset_inventory(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> list[dict[str, Any]] | None:
    """The datasets registered in ``network_registry`` at the end of this build.

    One record per preset — its scope, the attributes it returns, its curation
    and where its attributes come from — so the manifest says which datasets this
    build serves without a second query. Returns ``None`` when the registry is
    absent (a build that registered nothing). Descriptive only: a preset is
    metadata over the fact table, so the inventory changes without the built
    content changing, and it stays out of the identity hash.
    """

    schema_id = sql.Identifier(schema)
    # Same savepoint discipline as the coverage probe: a missing registry must
    # not undo the manifest DDL issued earlier in this transaction.
    cur.execute('SAVEPOINT network_preset_inventory_probe')
    try:
        cur.execute(
            sql.SQL(
                """
                SELECT name, kind, included_sources, interaction_class_scope,
                       evidence_scope, default_attributes, mandatory_attributes,
                       curation, attribute_sources
                FROM {}.network_registry
                ORDER BY name
                """
            ).format(schema_id)
        )
    except psycopg2.Error:
        cur.execute('ROLLBACK TO SAVEPOINT network_preset_inventory_probe')
        logger.debug('build manifest: no network_registry in %s; presets skipped', schema)
        return None
    rows = cur.fetchall()
    cur.execute('RELEASE SAVEPOINT network_preset_inventory_probe')

    return [
        {
            'name': name,
            'kind': kind,
            'included_sources': list(included_sources or ()),
            'interaction_class_scope': list(class_scope or ()),
            'evidence_scope': evidence_scope,
            'default_attributes': list(default_attributes or ()),
            'mandatory_attributes': list(mandatory_attributes or ()),
            'curation': curation,
            'attribute_sources': attribute_sources,
        }
        for (
            name,
            kind,
            included_sources,
            class_scope,
            evidence_scope,
            default_attributes,
            mandatory_attributes,
            curation,
            attribute_sources,
        ) in rows
    ]


@lru_cache(maxsize=32)
def _package_commit(module: str) -> dict[str, Any]:
    """``{commit, dirty}`` for a package's git checkout (contract-shaped, never empty).

    Unlike :func:`_module_git_metadata` (whose dirty check is scoped to a single
    input-module file), this reports the HEAD and the dirty state of the whole
    package subtree — the right granularity for a per-package manifest entry.
    """
    spec = importlib.util.find_spec(module)
    origin = getattr(spec, 'origin', None)
    if not origin or origin == 'built-in':
        return {'commit': None, 'dirty': False}
    package_dir = Path(origin).parent
    try:
        commit = subprocess.check_output(
            ['git', '-C', str(package_dir), 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        status = subprocess.check_output(
            ['git', '-C', str(package_dir), 'status', '--porcelain', '--', '.'],
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return {'commit': None, 'dirty': False}
    return {'commit': commit or None, 'dirty': bool(status.strip())}


def sync_resources_table(
    conn: psycopg2.extensions.connection,
    discovered: dict[str, list[object]],
    *,
    schema: str = 'public',
    prefer_bitmaps: bool = False,
) -> ResourceTableStats:
    """Upsert discovered resource metadata and omnipath_build source-level counts."""

    with conn.cursor() as cur:
        cur.execute('SET LOCAL max_parallel_workers_per_gather = 0')
        _ensure_resources_metadata_columns(cur, schema)
        present_sources = _present_sources(
            cur,
            schema=schema,
            prefer_bitmaps=prefer_bitmaps,
        )
        rows = [
            _resource_row(
                cur,
                schema=schema,
                source=source,
                functions=functions,
                present_sources=present_sources,
                prefer_bitmaps=prefer_bitmaps,
            )
            for source, functions in sorted(discovered.items())
        ]
        _validate_resource_names(rows)
        cur.executemany(
            sql.SQL(
                """
                INSERT INTO {}.resources (
                  resource_id,
                  resource_name,
                  resource_short,
                  resource_full,
                  synonyms,
                  description,
                  homepage_url,
                  license,
                  license_label,
                  pubmed_id,
                  resource_kind,
                  input_module,
                  input_module_commit,
                  input_module_dirty,
                  categories,
                  annotation_ontologies,
                  entity_count,
                  interaction_count,
                  association_count,
                  identifier_count,
                  ontology_term_count,
                  total_size_bytes,
                  last_downloaded_at,
                  last_built_at,
                  build_status
                )
                VALUES (
                  %(resource_id)s,
                  %(resource_name)s,
                  %(resource_short)s,
                  %(resource_full)s,
                  %(synonyms)s,
                  %(description)s,
                  %(homepage_url)s,
                  %(license)s,
                  %(license_label)s,
                  %(pubmed_id)s,
                  %(resource_kind)s,
                  %(input_module)s,
                  %(input_module_commit)s,
                  %(input_module_dirty)s,
                  %(categories)s,
                  %(annotation_ontologies)s,
                  %(entity_count)s,
                  %(interaction_count)s,
                  %(association_count)s,
                  %(identifier_count)s,
                  %(ontology_term_count)s,
                  %(total_size_bytes)s,
                  %(last_downloaded_at)s,
                  %(last_built_at)s,
                  %(build_status)s
                )
                ON CONFLICT (resource_id) DO UPDATE SET
                  resource_name = EXCLUDED.resource_name,
                  resource_short = EXCLUDED.resource_short,
                  resource_full = EXCLUDED.resource_full,
                  synonyms = EXCLUDED.synonyms,
                  description = EXCLUDED.description,
                  homepage_url = EXCLUDED.homepage_url,
                  license = EXCLUDED.license,
                  license_label = EXCLUDED.license_label,
                  pubmed_id = EXCLUDED.pubmed_id,
                  resource_kind = EXCLUDED.resource_kind,
                  input_module = EXCLUDED.input_module,
                  input_module_commit = EXCLUDED.input_module_commit,
                  input_module_dirty = EXCLUDED.input_module_dirty,
                  categories = EXCLUDED.categories,
                  annotation_ontologies = EXCLUDED.annotation_ontologies,
                  entity_count = EXCLUDED.entity_count,
                  interaction_count = EXCLUDED.interaction_count,
                  association_count = EXCLUDED.association_count,
                  identifier_count = EXCLUDED.identifier_count,
                  ontology_term_count = EXCLUDED.ontology_term_count,
                  total_size_bytes = EXCLUDED.total_size_bytes,
                  last_downloaded_at = EXCLUDED.last_downloaded_at,
                  last_built_at = EXCLUDED.last_built_at,
                  build_status = EXCLUDED.build_status
                """
            )
            .format(sql.Identifier(schema))
            .as_string(cur.connection),
            rows,
        )
    conn.commit()
    return ResourceTableStats(resources=len(rows))


def _ensure_resources_metadata_columns(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.resources
            ADD COLUMN IF NOT EXISTS input_module text,
            ADD COLUMN IF NOT EXISTS input_module_commit text,
            ADD COLUMN IF NOT EXISTS input_module_dirty boolean
              NOT NULL DEFAULT false,
            ADD COLUMN IF NOT EXISTS license_label text,
            ADD COLUMN IF NOT EXISTS resource_short text,
            ADD COLUMN IF NOT EXISTS resource_full text,
            ADD COLUMN IF NOT EXISTS synonyms text[] NOT NULL
              DEFAULT ARRAY[]::text[]
            """
        ).format(sql.Identifier(schema))
    )


def _resource_row(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source: str,
    functions: list[object],
    present_sources: set[str],
    prefer_bitmaps: bool,
) -> dict[str, Any]:
    config = _resource_config(functions)
    module_metadata = _input_module_metadata(functions)
    snapshot_metadata = _resource_snapshot_metadata(source, functions)
    counts = (
        _source_counts(
            cur,
            schema=schema,
            source=source,
            prefer_bitmaps=prefer_bitmaps,
        )
        if source in present_sources
        else _empty_source_counts()
    )
    last_built_at = (
        counts['last_built_at'] or snapshot_metadata.get('last_built_at')
        if counts['has_rows']
        else None
    )
    categories = _resource_categories(
        primary_category=getattr(config, 'primary_category', None),
        interaction_count=counts['interaction_count'],
        association_count=counts['association_count'],
        ontology_term_count=counts['ontology_term_count'],
    )
    names = _resource_names(config, source)
    return {
        'resource_id': source,
        'resource_name': getattr(config, 'name', source),
        'resource_short': names.short,
        'resource_full': names.full,
        'synonyms': list(names.synonyms),
        'description': getattr(config, 'description', None),
        'homepage_url': getattr(config, 'url', None),
        'license': _text_or_none(getattr(config, 'license', None)),
        'license_label': _cv_label(getattr(config, 'license', None)),
        'pubmed_id': getattr(config, 'pubmed', None),
        'resource_kind': getattr(config, 'resource_kind', 'data_resource'),
        'input_module': module_metadata['module'],
        'input_module_commit': module_metadata['commit'],
        'input_module_dirty': module_metadata['dirty'],
        'categories': Json(categories),
        'annotation_ontologies': Json(_ontology_labels(config)),
        'entity_count': counts['entity_count'],
        'interaction_count': counts['interaction_count'],
        'association_count': counts['association_count'],
        'identifier_count': counts['identifier_count'],
        'ontology_term_count': counts['ontology_term_count'],
        'total_size_bytes': snapshot_metadata.get('total_size_bytes', 0),
        'last_downloaded_at': snapshot_metadata.get('last_downloaded_at'),
        'last_built_at': last_built_at,
        'build_status': 'success' if counts['has_rows'] else 'not_built',
    }


def _validate_resource_names(rows: list[dict[str, Any]]) -> None:
    """Flag resource-name rule violations at build time.

    The ``_`` character is reserved for ``primary_secondary`` composite labels
    (``PhosphoSite_SIGNOR``), so no single resource may carry one in any of its
    three names; the slug must also be lowercase alphanumeric and space-free.
    """
    try:
        from pypath.inputs_v2.resource_names import validate_resource_name
    except ImportError:
        return  # pypath without the 3-name model (older pin); skip validation

    violations: list[tuple[str, list[str]]] = []
    for row in rows:
        errors = validate_resource_name(
            row['resource_id'],
            row.get('resource_short') or '',
            row.get('resource_full') or '',
        )
        if errors:
            violations.append((row['resource_id'], errors))
    if violations:
        print(
            f'[resources] resource-name rule violations: {len(violations)}',
            flush=True,
        )
        for resource_id, errors in violations[:30]:
            print(f'[resources]   {resource_id}: {"; ".join(errors)}', flush=True)


class _FallbackNames(NamedTuple):
    short: str
    full: str
    synonyms: tuple[str, ...]


def _resource_names(config: object, source: str):
    """Resolve the resource 3-name model.

    Prefers the resource's own ``ResourceConfig.names()`` — which resolves via the
    canonical ``ResourceCv`` slug against the authoritative ``resources.json`` —
    and otherwise resolves the build's ``resource_id`` (``source``) through the
    same registry / filter index.
    """
    try:
        from pypath.inputs_v2.resource_names import resolve_filter, resolve_names
    except ImportError:
        # pypath without the 3-name model (older pin): short/full fall back to
        # the resource's existing name, no synonyms.
        name = getattr(config, 'name', None) or source
        return _FallbackNames(short=name, full=name, synonyms=())

    if config is not None and hasattr(config, 'names'):
        try:
            return config.names()
        except Exception:  # pragma: no cover - defensive
            pass

    # No usable config: resolve the source key via the registry/filter index.
    canonical = resolve_filter(source)
    return resolve_names(name=source, slug=canonical)


def _present_sources(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    prefer_bitmaps: bool = False,
) -> set[str]:
    if prefer_bitmaps:
        cur.execute(
            sql.SQL(
                """
                SELECT facet_value
                FROM {}.facet_entity_bitmap
                WHERE facet_name = 'source'
                UNION
                SELECT facet_value
                FROM {}.facet_relation_bitmap
                WHERE facet_name = 'source'
                """
            ).format(
                sql.Identifier(schema),
                sql.Identifier(schema),
            )
        )
        return {row[0] for row in cur.fetchall()}

    cur.execute(
        sql.SQL(
            """
            SELECT ds.name
            FROM {}.entity_evidence ee
            JOIN {}.data_source ds
              ON ds.source_id = ee.source_id
            UNION
            SELECT ds.name
            FROM {}.relation_evidence re
            JOIN {}.data_source ds
              ON ds.source_id = re.source_id
            UNION
            SELECT ds.name
            FROM {}.ontology_terms ot
            JOIN {}.data_source ds
              ON ds.source_id = ot.source_id
            """
        ).format(
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
        )
    )
    return {row[0] for row in cur.fetchall()}


def _empty_source_counts() -> dict[str, Any]:
    return {
        'entity_count': 0,
        'identifier_count': 0,
        'interaction_count': 0,
        'association_count': 0,
        'ontology_term_count': 0,
        'last_built_at': None,
        'has_rows': False,
    }


def _resource_config(functions: list[object]) -> object | None:
    for fn in functions:
        if getattr(fn, 'function_name', None) != 'resource':
            continue
        config = getattr(getattr(fn, 'call', None), 'config', None)
        if config is not None:
            return config
    return None


def _input_module_metadata(functions: list[object]) -> dict[str, Any]:
    fn = next(
        (
            fn
            for fn in functions
            if getattr(fn, 'function_name', None) == 'resource'
        ),
        functions[0] if functions else None,
    )
    module = getattr(fn, 'qualified_module', None) if fn is not None else None
    metadata = _module_git_metadata(module) if module else {}
    return {
        'module': module,
        'commit': metadata.get('commit'),
        'dirty': bool(metadata.get('dirty', False)),
    }


@lru_cache(maxsize=512)
def _module_git_metadata(module: str) -> dict[str, Any]:
    spec = importlib.util.find_spec(module)
    origin = getattr(spec, 'origin', None)
    if not origin or origin == 'built-in':
        return {}
    path = Path(origin)
    try:
        commit = subprocess.check_output(
            ['git', '-C', str(path.parent), 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        status = subprocess.check_output(
            [
                'git',
                '-C',
                str(path.parent),
                'status',
                '--porcelain',
                '--',
                str(path),
            ],
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return {}
    return {
        'commit': commit or None,
        'dirty': bool(status.strip()),
    }


def _resource_snapshot_metadata(
    source: str,
    functions: list[object],
) -> dict[str, Any]:
    del source, functions
    return {}


def _source_counts(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source: str,
    prefer_bitmaps: bool = False,
) -> dict[str, Any]:
    if prefer_bitmaps:
        bitmap_counts = _bitmap_source_counts(
            cur,
            schema=schema,
            source=source,
        )
        if bitmap_counts is not None:
            return bitmap_counts

    cur.execute(
        sql.SQL(
            """
            WITH selected_source AS (
              SELECT source_id
              FROM {}.data_source
              WHERE name = %s
            ),
            source_relations AS (
              SELECT DISTINCT rer.relation_id
              FROM {}.relation_evidence_relation rer
              JOIN {}.relation_evidence re
                ON re.source_id = rer.source_id
               AND re.relation_evidence_id = rer.relation_evidence_id
              WHERE re.source_id = (SELECT source_id FROM selected_source)
              UNION
              SELECT DISTINCT ear.relation_id
              FROM {}.entity_annotation_relation ear
              WHERE ear.source_id = (SELECT source_id FROM selected_source)
            )
            SELECT
              (
                SELECT COUNT(DISTINCT r.entity_id)::bigint
                FROM {}.entity_evidence_resolution r
                JOIN {}.entity_evidence ee
                  ON ee.source_id = r.source_id
                 AND ee.entity_evidence_id = r.entity_evidence_id
                WHERE ee.source_id = (SELECT source_id FROM selected_source)
                  AND r.entity_id IS NOT NULL
              ) AS entity_count,
              (
                SELECT 0::bigint
              ) AS identifier_count,
              (
                SELECT COUNT(DISTINCT rel.relation_id)::bigint
                FROM source_relations sr
                JOIN {}.relation rel
                  ON rel.relation_id = sr.relation_id
                JOIN {}.vocab_relation_category rc
                  ON rc.relation_category_id = rel.relation_category_id
                WHERE rc.name = 'interaction'
              ) AS interaction_count,
              (
                SELECT COUNT(DISTINCT rel.relation_id)::bigint
                FROM source_relations sr
                JOIN {}.relation rel
                  ON rel.relation_id = sr.relation_id
                JOIN {}.vocab_relation_category rc
                  ON rc.relation_category_id = rel.relation_category_id
                WHERE rc.name = 'association'
              ) AS association_count,
              (
                SELECT COUNT(*)::bigint
                FROM {}.ontology_terms ot
                WHERE ot.source_id = (SELECT source_id FROM selected_source)
              ) AS ontology_term_count,
              (
                SELECT NULL::timestamptz
              ) AS last_built_at,
              (
                SELECT EXISTS (
                  SELECT 1
                  FROM {}.entity_evidence
                  WHERE source_id = (
                    SELECT source_id FROM selected_source
                  )
                  UNION ALL
                  SELECT 1
                  FROM {}.ontology_terms
                  WHERE source_id = (
                    SELECT source_id FROM selected_source
                  )
                )
              ) AS has_rows
            """
        ).format(
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
        ),
        [source],
    )
    row = cur.fetchone()
    return {
        'entity_count': int(row[0] or 0),
        'identifier_count': int(row[1] or 0),
        'interaction_count': int(row[2] or 0),
        'association_count': int(row[3] or 0),
        'ontology_term_count': int(row[4] or 0),
        'last_built_at': row[5],
        'has_rows': bool(row[6]),
    }


def _bitmap_source_counts(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source: str,
) -> dict[str, Any] | None:
    cur.execute(
        sql.SQL(
            """
            WITH
            source_entity AS (
              SELECT entity_bitmap, entity_count
              FROM {}.facet_entity_bitmap
              WHERE facet_name = 'source'
                AND facet_value = %s
            ),
            ontology_entity AS (
              SELECT entity_bitmap
              FROM {}.facet_entity_bitmap
              WHERE facet_name = 'entity_type'
                AND facet_value = %s
            ),
            source_relation AS (
              SELECT relation_bitmap, relation_count
              FROM {}.facet_relation_bitmap
              WHERE facet_name = 'source'
                AND facet_value = %s
            ),
            interaction_relation AS (
              SELECT rb_or_agg(relation_bitmap) AS relation_bitmap
              FROM {}.facet_relation_bitmap
              WHERE facet_name = 'predicate'
                AND facet_category = 'interaction'
            ),
            association_relation AS (
              SELECT rb_or_agg(relation_bitmap) AS relation_bitmap
              FROM {}.facet_relation_bitmap
              WHERE facet_name = 'predicate'
                AND facet_category = 'association'
            )
            SELECT
              COALESCE(
                (SELECT entity_count FROM source_entity),
                0
              )::bigint AS entity_count,
              (
                SELECT 0::bigint
              ) AS identifier_count,
              COALESCE(
                (
                  SELECT rb_and_cardinality(
                    source_relation.relation_bitmap,
                    interaction_relation.relation_bitmap
                  )
                  FROM source_relation
                  CROSS JOIN interaction_relation
                ),
                0
              )::bigint AS interaction_count,
              COALESCE(
                (
                  SELECT rb_and_cardinality(
                    source_relation.relation_bitmap,
                    association_relation.relation_bitmap
                  )
                  FROM source_relation
                  CROSS JOIN association_relation
                ),
                0
              )::bigint AS association_count,
              COALESCE(
                (
                  SELECT rb_and_cardinality(
                    source_entity.entity_bitmap,
                    ontology_entity.entity_bitmap
                  )
                  FROM source_entity
                  CROSS JOIN ontology_entity
                ),
                0
              )::bigint AS ontology_term_count,
              NULL::timestamptz AS last_built_at,
              EXISTS (SELECT 1 FROM source_entity)
                OR EXISTS (SELECT 1 FROM source_relation) AS has_rows
            """
        ).format(
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
            sql.Identifier(schema),
        ),
        [
            source,
            CV_TERM_ENTITY_TYPE,
            source,
        ],
    )
    row = cur.fetchone()
    if row is None or not row[6]:
        return None
    return {
        'entity_count': int(row[0] or 0),
        'identifier_count': int(row[1] or 0),
        'interaction_count': int(row[2] or 0),
        'association_count': int(row[3] or 0),
        'ontology_term_count': int(row[4] or 0),
        'last_built_at': row[5],
        'has_rows': bool(row[6]),
    }


def _resource_categories(
    *,
    primary_category: object,
    interaction_count: int,
    association_count: int,
    ontology_term_count: int,
) -> list[str]:
    categories = []
    primary = _text_or_none(primary_category)
    if primary:
        categories.append(primary)
    if interaction_count:
        categories.append('interaction')
    if association_count:
        categories.append('association')
    if ontology_term_count:
        categories.append('ontology')
    return sorted(set(categories))


def _ontology_labels(config: object | None) -> list[str]:
    values = []
    for ontology in getattr(config, 'annotation_ontologies', ()) or ():
        label = getattr(ontology, 'definition', None) or str(ontology)
        values.append(str(label))
    return values


def _text_or_none(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _cv_label(value: object) -> str | None:
    label = getattr(value, 'definition', None) or getattr(value, 'label', None)
    if label:
        return str(label)
    return _text_or_none(value)
