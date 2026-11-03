"""Build the MetSigDB membership substrate from the build database.

The step reads what the core build already wrote. It parses no upstream
resource and downloads nothing: `relation_evidence`,
`entity_evidence_resolution` and `entity_ontology_relation` already hold every
membership, resolved against the canonical entity layer.

Each resource runs the same three phases. An extraction query stages one row
per published pair, a shared projection reads the metabolite side off the
canonical layer, and a shared upsert writes the rows under their deterministic
identity. Only the extraction differs per resource.
"""

from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path

import psycopg2.errors
import psycopg2.extensions
from psycopg2 import sql

from omnipath_build.metsigdb.mapping import (
    CHEMICAL_ENTITY_TYPE,
    KEGG_OVERVIEW_MAPS,
    RESOURCES,
    ResourceRule,
)

_SQL_DIR = Path(__file__).with_name('sql')

logger = logging.getLogger(__name__)

TABLE = 'metsigdb_membership'

# Temp-table headroom for the extractions. A temp table lives in the backend's
# local buffers, and the default 8 MB is far too small for intermediates of a
# few million rows: Postgres reports "no empty local buffer available" rather
# than spilling. The setting only takes effect before a session touches its
# first temp table, so the build sets it once per connection and lets a later
# call fail quietly.
TEMP_BUFFERS = os.environ.get('METSIGDB_TEMP_BUFFERS', '1GB')


@dataclass(frozen=True)
class ResourceLoadStats:
    """What one resource cost and what it published."""

    resource: str
    rows: int
    sets: int
    metabolites: int
    removed: int
    seconds: float


def _sql_text(name: str) -> str:
    return (_SQL_DIR / name).read_text(encoding='utf-8')


def ensure_membership_table(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> None:
    """Create the membership table and its filter indexes if they are absent.

    The DDL is idempotent. Applying it to a populated table leaves both the
    schema and the rows untouched, so a rebuild can call it first without a
    guard.
    """
    with conn.cursor() as cur:
        # Unqualified names in the DDL land in the target schema. public.*
        # still resolves, so the table can reference the canonical layer.
        cur.execute(
            sql.SQL('SET search_path = {}, public').format(sql.Identifier(schema))
        )
        cur.execute(_sql_text('membership_table.sql'))
        cur.execute('RESET search_path')
    conn.commit()
    logger.info('metsigdb: %s.%s is present', schema, TABLE)


def _widen_temp_buffers(cur) -> None:
    """Give the backend room for the extraction intermediates.

    Postgres refuses the change once the session has touched a temp table, so a
    connection that already loaded one resource keeps the value it got then.
    That is the wanted behavior, and the refusal is not an error.

    The refusal arrives as `invalid_parameter_value`, and it aborts the
    transaction it happened in, so the attempt runs inside a savepoint and the
    caller keeps a usable connection.

    Tolerating the refusal is right **within one build**, where the session
    already carries the widened value from the first resource. It is not a
    licence to run on any connection: a session that reaches here still on the
    8 MB default cannot stage KEGG or ClassyFire, and fails with "no empty
    local buffer available". Callers must pass a connection that has touched
    no temp table — which is why the derive opens one of its own rather than
    lending the step the session it has been building on.
    """
    savepoint = not cur.connection.autocommit
    if savepoint:
        cur.execute('SAVEPOINT metsigdb_temp_buffers')
    try:
        cur.execute(f"SET temp_buffers = '{TEMP_BUFFERS}'")
    except (
        psycopg2.errors.InvalidParameterValue,
        psycopg2.errors.ActiveSqlTransaction,
    ):
        if savepoint:
            cur.execute('ROLLBACK TO SAVEPOINT metsigdb_temp_buffers')
        logger.debug('metsigdb: temp_buffers already fixed for this session')
    else:
        if savepoint:
            cur.execute('RELEASE SAVEPOINT metsigdb_temp_buffers')


def _scalar(cur, query: str, params=None):
    cur.execute(query, params)
    row = cur.fetchone()
    return None if row is None else row[0]


def _source_id(cur, name: str) -> int:
    source_id = _scalar(cur, 'SELECT source_id FROM data_source WHERE name = %s', [name])
    if source_id is None:
        raise LookupError(f'metsigdb: source {name!r} is not loaded in this build')
    return source_id


def _entity_type_id(cur, name: str) -> int:
    type_id = _scalar(
        cur, 'SELECT entity_type_id FROM vocab_entity_type WHERE name = %s', [name]
    )
    if type_id is None:
        raise LookupError(f'metsigdb: entity type {name!r} is absent')
    return type_id


def build_id(conn: psycopg2.extensions.connection) -> str:
    """The build the substrate is stamped with, read from the manifest.

    A membership row records which build produced it, so a consumer can tell a
    stale row from a current one. The manifest holds exactly one row.
    """
    with conn.cursor() as cur:
        stamp = _scalar(cur, 'SELECT build_id FROM build_manifest')
    if stamp is None:
        raise LookupError('metsigdb: build_manifest is empty; run the derive first')
    return stamp


def _provenance_source(cur, source_name: str) -> str:
    """The input module and commit that parsed the resource.

    The MetSigDB step downloads nothing, so it cites the provenance the core
    build recorded rather than restating a URL of its own.
    """
    cur.execute(
        'SELECT input_module, input_module_commit FROM resources WHERE resource_id = %s',
        [source_name],
    )
    row = cur.fetchone()
    if row is None or not row[0]:
        return source_name
    module, commit = row
    return f'{module}@{commit[:12]}' if commit else module


def load_resource(
    conn: psycopg2.extensions.connection,
    rule: ResourceRule,
    *,
    stamp: str,
    max_records: int | None = None,
) -> ResourceLoadStats:
    """Extract, project and publish one resource into the substrate.

    Rows this build did not write are removed for that resource alone, so a
    source that shrinks upstream shrinks the substrate instead of leaving
    memberships behind.
    """
    started = time.monotonic()
    with conn.cursor() as cur:
        _widen_temp_buffers(cur)
        params = {
            'source_id': _source_id(cur, rule.source_name),
            'set_entity_type_id': _entity_type_id(cur, rule.set_entity_type),
            'chemical_entity_type_id': _entity_type_id(cur, CHEMICAL_ENTITY_TYPE),
            'max_records': max_records,
        }
        if rule.name == 'KEGG':
            params['overview_maps'] = list(KEGG_OVERVIEW_MAPS)
        if rule.hierarchy_source_name:
            params['hierarchy_source_id'] = _source_id(cur, rule.hierarchy_source_name)

        cur.execute(_sql_text(rule.extraction), params)
        cur.execute('CREATE INDEX ON metsigdb_stage (metabolite_entity_id)')
        cur.execute('CREATE INDEX ON metsigdb_stage (set_source_id)')
        cur.execute(
            'CREATE INDEX ON metsigdb_stage (set_source_id, metabolite_entity_id)'
        )
        cur.execute('ANALYZE metsigdb_stage')
        staged = _scalar(cur, 'SELECT count(*) FROM metsigdb_stage')

        cur.execute('CREATE INDEX ON metsigdb_stage (set_entity_id)')
        cur.execute(_sql_text('publish_membership.sql'), params)

        cur.execute(
            _sql_text('upsert_membership.sql'),
            {
                'resource': rule.name,
                'set_type': rule.set_type,
                'provenance_source': _provenance_source(cur, rule.source_name),
                'build_id': stamp,
            },
        )

        # After the upsert, the resource's rows are whatever it staged and
        # nothing else. An anti-join against the stage is exact for every case
        # the build has: a membership that vanished upstream, a metabolite that
        # re-anchored on a merged entity, and a capped run that must hold only
        # what it loaded. Matching on the build stamp would miss all three
        # whenever the stamp did not change.
        cur.execute(
            f"""
            DELETE FROM {TABLE} m
            WHERE m.resource = %(resource)s
              AND NOT EXISTS (
                SELECT 1 FROM metsigdb_stage s
                WHERE s.set_source_id = m.set_source_id
                  AND s.metabolite_entity_id = m.metabolite_entity_id
              )
            """,
            {'resource': rule.name},
        )
        removed = cur.rowcount

        cur.execute(
            f"""
            SELECT count(*), count(DISTINCT set_source_id),
                   count(DISTINCT metabolite_entity_id)
            FROM {TABLE} WHERE resource = %s
            """,
            [rule.name],
        )
        rows, sets, metabolites = cur.fetchone()

        cur.execute('DROP TABLE IF EXISTS metsigdb_stage')
        cur.execute('DROP TABLE IF EXISTS metsigdb_projection')
    conn.commit()

    stats = ResourceLoadStats(
        resource=rule.name,
        rows=rows,
        sets=sets,
        metabolites=metabolites,
        removed=removed,
        seconds=round(time.monotonic() - started, 1),
    )
    logger.info(
        'metsigdb: %s staged=%s rows=%s sets=%s metabolites=%s removed=%s seconds=%s',
        stats.resource, staged, stats.rows, stats.sets,
        stats.metabolites, stats.removed, stats.seconds,
    )
    return stats


@dataclass(frozen=True)
class BuildStats:
    """What a whole MetSigDB build published."""

    build_id: str
    partial: bool
    rows: int
    seconds: float
    resources: tuple[ResourceLoadStats, ...]


def build_metsigdb(
    conn: psycopg2.extensions.connection,
    *,
    database_url: str,
    resources: tuple[ResourceRule, ...] = RESOURCES,
    max_records: int | None = None,
    schema: str = 'public',
) -> BuildStats:
    """Publish the whole MetSigDB substrate.

    Applies the DDL, loads each resource in turn, and stamps every row with the
    current build. A non-zero ``max_records`` caps each resource and marks the
    build manifest partial, because a capped substrate is not authoritative.
    """
    started = time.monotonic()
    ensure_membership_table(conn, schema=schema)
    stamp = build_id(conn)

    loaded = tuple(
        load_resource(conn, rule, stamp=stamp, max_records=max_records)
        for rule in resources
    )

    partial = bool(max_records)
    if partial:
        _mark_partial_build(conn)

    _vacuum(conn, database_url=database_url, schema=schema)

    with conn.cursor() as cur:
        rows = _scalar(cur, f'SELECT count(*) FROM {TABLE}')

    stats = BuildStats(
        build_id=stamp,
        partial=partial,
        rows=rows,
        seconds=round(time.monotonic() - started, 1),
        resources=loaded,
    )
    logger.info(
        'metsigdb: build=%s rows=%s partial=%s seconds=%s',
        stats.build_id, stats.rows, stats.partial, stats.seconds,
    )
    return stats


def _vacuum(
    conn: psycopg2.extensions.connection,
    *,
    database_url: str,
    schema: str = 'public',
) -> None:
    """Reclaim the rebuild's dead rows and refresh the planner's statistics.

    A rebuild updates every row it already had, so Postgres leaves one dead
    version per row: a second full build takes the table from 2.2 GB to 4.3 GB
    on disk. Plain VACUUM does not shrink the file, it makes that space
    reusable, which is what keeps the table at a steady size instead of growing
    with every build. Reclaiming the file itself needs VACUUM FULL and an
    exclusive lock, which is an operator's decision, not the build's.

    ANALYZE is not optional here. The build replaces the whole table, and the
    serving layer plans its filters against these statistics.

    The fix first tried toggling ``conn.autocommit`` on the load connection.
    That did not reliably close an already-open transaction on this
    driver/server combination. VACUUM still refused with "cannot run inside a
    transaction block", even right after an explicit ``conn.commit()``. A
    dedicated connection sidesteps the ambiguity: it starts in autocommit
    mode, so Postgres never sees an open transaction to refuse.
    """
    conn.commit()  # make the load visible to the dedicated connection below.
    vacuum_conn = psycopg2.connect(database_url)
    vacuum_conn.autocommit = True  # set before any statement runs on it.
    try:
        with vacuum_conn.cursor() as cur:
            cur.execute(
                sql.SQL('VACUUM (ANALYZE) {}.{}').format(
                    sql.Identifier(schema), sql.Identifier(TABLE)
                )
            )
    finally:
        vacuum_conn.close()
    logger.info('metsigdb: vacuumed and analyzed %s.%s', schema, TABLE)


def _mark_partial_build(conn: psycopg2.extensions.connection) -> None:
    """Flag the manifest, because a capped run did not load everything.

    The flag is never cleared here. Only a full build may claim to be complete,
    and that claim belongs to the derive that writes the manifest.
    """
    with conn.cursor() as cur:
        cur.execute('UPDATE build_manifest SET partial_build = true')
    conn.commit()
    logger.info('metsigdb: capped run; build_manifest.partial_build set')
