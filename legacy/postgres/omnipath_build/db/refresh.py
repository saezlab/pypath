"""Source-scoped content detection and deletion helpers.

Plain source loads can skip already-loaded content. Reloads remove evidence rows
and evidence annotations before a source is streamed again. Graph entities and
relations that were only reachable through that source are then
garbage-collected, while shared canonical graph rows backed by other sources are
preserved.
"""

from __future__ import annotations

from dataclasses import dataclass

from psycopg2 import sql
import psycopg2.extensions

from omnipath_build.db.schema import (
    SOURCE_PARTITIONED_TABLES,
    SOURCE_PARTITION_DROP_ORDER,
    _source_partition_suffix,
)

@dataclass(frozen=True)
class SourceContentDropStats:
    """Summary of a source refresh/delete cleanup."""

    source: str
    source_id: int | None = None
    strategy: str = 'missing_source'
    affected_relations: int = 0
    affected_entities: int = 0
    affected_identifiers: int = 0
    partitions_dropped: int = 0
    deleted_relations: int = 0
    deleted_entities: int = 0
    deleted_identifiers: int = 0
    deleted_annotations: int = 0
    refreshed_relation_counts: int = 0


def delete_source_content(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    source: str,
    drop_partitions: bool = True,
) -> SourceContentDropStats:
    """Delete source-scoped evidence and orphaned resolved graph rows.

    Source-specific evidence normally lives in list partitions. When all source
    partitions are present, dropping those child tables is much faster than
    deleting row by row from each partitioned parent. Older databases may have
    source rows in default partitions, so this keeps a row-delete fallback.
    """

    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                """
                SELECT source_id
                FROM {}.data_source
                WHERE name = %s
                """
            ).format(schema_id),
            [source],
        )
        row = cur.fetchone()
        if row is None:
            conn.commit()
            return SourceContentDropStats(source=source)
        source_id = int(row[0])
        affected = _create_source_cleanup_scope(
            cur,
            schema=schema,
            source_id=source_id,
        )
        use_partition_drop = (
            drop_partitions
            and _source_partitions_exist(
                cur,
                schema=schema,
                source=source,
            )
        )
        partitions_dropped = (
            _drop_source_partitions(cur, schema=schema, source=source)
            if use_partition_drop
            else 0
        )
        if not use_partition_drop:
            _delete_source_evidence_rows(
                cur,
                schema=schema,
                source_id=source_id,
            )

        deleted = _garbage_collect_source_cleanup(
            cur,
            schema=schema,
            source=source,
        )
    conn.commit()
    return SourceContentDropStats(
        source=source,
        source_id=source_id,
        strategy='partition_drop' if use_partition_drop else 'row_delete',
        affected_relations=affected['relations'],
        affected_entities=affected['entities'],
        affected_identifiers=affected['identifiers'],
        partitions_dropped=partitions_dropped,
        deleted_relations=deleted['relations'],
        deleted_entities=deleted['entities'],
        deleted_identifiers=deleted['identifiers'],
        deleted_annotations=deleted['annotations'],
        refreshed_relation_counts=deleted['relation_counts'],
    )


def source_has_content(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    source: str,
) -> bool:
    """Return whether a source already has source-scoped content loaded."""

    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                """
                SELECT source_id
                FROM {}.data_source
                WHERE name = %s
                """
            ).format(schema_id),
            [source],
        )
        row = cur.fetchone()
        if row is None:
            return False
        source_id = int(row[0])
        for table in SOURCE_PARTITIONED_TABLES:
            cur.execute(
                sql.SQL(
                    """
                    SELECT 1
                    FROM {}.{}
                    WHERE source_id = %s
                    LIMIT 1
                    """
                ).format(schema_id, sql.Identifier(table)),
                [source_id],
            )
            if cur.fetchone() is not None:
                return True
    return False


def _create_source_cleanup_scope(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source_id: int,
) -> dict[str, int]:
    """Snapshot graph and dimension rows that may become orphaned."""

    schema_id = sql.Identifier(schema)
    cur.execute('DROP TABLE IF EXISTS _refresh_removed_relation')
    cur.execute(
        sql.SQL(
            """
            CREATE TEMP TABLE _refresh_removed_relation ON COMMIT DROP AS
            SELECT DISTINCT rer.relation_id
            FROM {}.relation_evidence_relation rer
            WHERE rer.source_id = %s
            UNION
            SELECT DISTINCT ear.relation_id
            FROM {}.entity_annotation_relation ear
            WHERE ear.source_id = %s
            """
        ).format(schema_id, schema_id),
        [source_id, source_id],
    )
    cur.execute(
        """
        CREATE UNIQUE INDEX IF NOT EXISTS _refresh_removed_relation_idx
        ON _refresh_removed_relation (relation_id)
        """
    )
    cur.execute('DROP TABLE IF EXISTS _refresh_removed_entity')
    cur.execute(
        sql.SQL(
            """
            CREATE TEMP TABLE _refresh_removed_entity ON COMMIT DROP AS
            SELECT DISTINCT ot.term_entity_id AS entity_id
            FROM {}.ontology_terms ot
            WHERE ot.source_id = %s
            UNION
            SELECT DISTINCT eor.subject_entity_id AS entity_id
            FROM {}.entity_ontology_relation eor
            WHERE eor.source_id = %s
            UNION
            SELECT DISTINCT eor.object_entity_id AS entity_id
            FROM {}.entity_ontology_relation eor
            WHERE eor.source_id = %s
            UNION
            SELECT DISTINCT r.entity_id
            FROM {}.entity_evidence_resolution r
            WHERE r.source_id = %s
              AND r.entity_id IS NOT NULL
            UNION
            SELECT DISTINCT re.subject_entity_id AS entity_id
            FROM {}.relation_evidence re
            WHERE re.source_id = %s
              AND re.subject_entity_id IS NOT NULL
            UNION
            SELECT DISTINCT re.object_entity_id AS entity_id
            FROM {}.relation_evidence re
            WHERE re.source_id = %s
              AND re.object_entity_id IS NOT NULL
            """
        ).format(
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
        ),
        [source_id, source_id, source_id, source_id, source_id, source_id],
    )
    cur.execute(
        """
        CREATE UNIQUE INDEX IF NOT EXISTS _refresh_removed_entity_idx
        ON _refresh_removed_entity (entity_id)
        """
    )
    cur.execute('DROP TABLE IF EXISTS _refresh_removed_identifier')
    cur.execute(
        sql.SQL(
            """
            CREATE TEMP TABLE _refresh_removed_identifier ON COMMIT DROP AS
            SELECT DISTINCT eei.identifier_id
            FROM {}.entity_evidence_identifier eei
            WHERE eei.source_id = %s
            UNION
            SELECT DISTINCT ei.identifier_id
            FROM {}.entity_identifier ei
            WHERE ei.source_id = %s
            """
        ).format(schema_id, schema_id),
        [source_id, source_id],
    )
    cur.execute(
        """
        CREATE UNIQUE INDEX IF NOT EXISTS _refresh_removed_identifier_idx
        ON _refresh_removed_identifier (identifier_id)
        """
    )
    cur.execute('DROP TABLE IF EXISTS _refresh_dirty_entity')
    cur.execute(
        sql.SQL(
            """
            CREATE TEMP TABLE _refresh_dirty_entity ON COMMIT DROP AS
            SELECT entity_id
            FROM _refresh_removed_entity
            UNION
            SELECT r.subject_entity_id AS entity_id
            FROM {}.relation r
            JOIN _refresh_removed_relation rr
              ON rr.relation_id = r.relation_id
            UNION
            SELECT r.object_entity_id AS entity_id
            FROM {}.relation r
            JOIN _refresh_removed_relation rr
              ON rr.relation_id = r.relation_id
            """
        ).format(schema_id, schema_id)
    )
    cur.execute(
        """
        CREATE UNIQUE INDEX IF NOT EXISTS _refresh_dirty_entity_idx
        ON _refresh_dirty_entity (entity_id)
        """
    )
    return {
        'relations': _count_temp_table(cur, '_refresh_removed_relation'),
        'entities': _count_temp_table(cur, '_refresh_removed_entity'),
        'identifiers': _count_temp_table(cur, '_refresh_removed_identifier'),
    }


def _source_partitions_exist(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source: str,
) -> bool:
    suffix = _source_partition_suffix(source)
    for table in SOURCE_PARTITIONED_TABLES:
        cur.execute(
            'SELECT to_regclass(%s) IS NOT NULL',
            [f'{schema}.{table}_{suffix}'],
        )
        if not bool(cur.fetchone()[0]):
            return False
    return True


def _drop_source_partitions(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source: str,
) -> int:
    suffix = _source_partition_suffix(source)
    dropped = 0
    for table in SOURCE_PARTITION_DROP_ORDER:
        partition = f'{table}_{suffix}'
        cur.execute('SELECT to_regclass(%s) IS NOT NULL', [f'{schema}.{partition}'])
        if not bool(cur.fetchone()[0]):
            continue
        if table == 'entity_evidence':
            cur.execute(
                sql.SQL('DELETE FROM {}.{}').format(
                    sql.Identifier(schema),
                    sql.Identifier(partition),
                )
            )
        cur.execute(
            sql.SQL('ALTER TABLE {}.{} DETACH PARTITION {}.{}').format(
                sql.Identifier(schema),
                sql.Identifier(table),
                sql.Identifier(schema),
                sql.Identifier(partition),
            )
        )
        cur.execute(
            sql.SQL(
                'DROP TABLE {}.{}'
            ).format(sql.Identifier(schema), sql.Identifier(partition))
        )
        dropped += 1
    return dropped


def _delete_source_evidence_rows(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source_id: int,
) -> None:
    schema_id = sql.Identifier(schema)
    for table in (
        'ontology_terms',
        'entity_ontology_relation',
        'relation_evidence_annotation',
        'relation_evidence_relation',
        'entity_annotation_relation',
        'relation_evidence',
        'entity_evidence_annotation',
        'entity_evidence_resolution',
        # evidence_state (the molecular state each record asserted) is
        # source-partitioned like the rows above; it must be deleted here too,
        # or a reload of the same source leaves the old per-source rows behind
        # and the re-staged rows collide on the evidence_state_source_<id>
        # primary key at derive time. Delete it before entity_evidence (it
        # references entity_evidence). `state` / `state_component` are
        # content-hashed and shared across sources, so they are intentionally
        # not source-scoped here.
        'evidence_state',
        'entity_identifier',
        'entity_evidence_identifier',
        'entity_evidence',
    ):
        cur.execute(
            sql.SQL('DELETE FROM {}.{} WHERE source_id = %s').format(
                schema_id,
                sql.Identifier(table),
            ),
            [source_id],
        )


def _garbage_collect_source_cleanup(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source: str,
) -> dict[str, int]:
    schema_id = sql.Identifier(schema)
    deleted: dict[str, int] = {}
    cur.execute(
        sql.SQL(
            """
            DELETE FROM {}.identifier_evidence i
            USING _refresh_removed_identifier ri
            WHERE i.identifier_id = ri.identifier_id
              AND NOT EXISTS (
                SELECT 1
                FROM {}.entity_evidence_identifier eei
                WHERE eei.identifier_id = i.identifier_id
              )
              AND NOT EXISTS (
                SELECT 1
                FROM {}.entity_identifier ei
                WHERE ei.identifier_id = i.identifier_id
              )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    deleted['identifiers'] = int(cur.rowcount)
    cur.execute(
        sql.SQL(
            """
            DELETE FROM {}.relation r
            USING _refresh_removed_relation rr
            WHERE r.relation_id = rr.relation_id
              AND NOT EXISTS (
                SELECT 1
                FROM {}.relation_evidence_relation rer
                WHERE rer.relation_id = r.relation_id
              )
              AND NOT EXISTS (
                SELECT 1
                FROM {}.entity_annotation_relation ear
                WHERE ear.relation_id = r.relation_id
              )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    deleted['relations'] = int(cur.rowcount)
    deleted['relation_counts'] = _refresh_dirty_relation_counts(cur, schema)
    cur.execute(
        sql.SQL(
            """
            DELETE FROM {}.entity e
            USING _refresh_removed_entity re
            WHERE e.entity_id = re.entity_id
              AND NOT EXISTS (
                SELECT 1
                FROM {}.entity_evidence_resolution r
                WHERE r.entity_id = e.entity_id
              )
              AND NOT EXISTS (
                SELECT 1
                FROM {}.relation rel
                WHERE rel.subject_entity_id = e.entity_id
              )
              AND NOT EXISTS (
                SELECT 1
                FROM {}.relation rel
                WHERE rel.object_entity_id = e.entity_id
              )
            """
        ).format(schema_id, schema_id, schema_id, schema_id)
    )
    deleted['entities'] = int(cur.rowcount)
    cur.execute(
        sql.SQL(
            """
            DELETE FROM {}.annotation a
            WHERE NOT EXISTS (
                SELECT 1
                FROM {}.entity_evidence_annotation eea
                WHERE eea.annotation_key = a.annotation_key
              )
              AND NOT EXISTS (
                SELECT 1
                FROM {}.relation_evidence_annotation rea
                WHERE rea.annotation_key = a.annotation_key
              )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    deleted['annotations'] = int(cur.rowcount)
    _clear_source_derived_facets(cur, schema=schema, source=source)
    _mark_resource_not_built(cur, schema=schema, source=source)
    return deleted


def _refresh_dirty_relation_counts(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> int:
    if not _table_exists(cur, schema, 'entity_relation_counts'):
        return 0
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            DELETE FROM {}.entity_relation_counts c
            USING _refresh_dirty_entity d
            WHERE c.entity_id = d.entity_id
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            INSERT INTO {}.entity_relation_counts (
              entity_id,
              relation_count,
              ontology_annotated_entity_count,
              ontology_annotated_relation_count,
              search_count
            )
            SELECT
              d.entity_id,
              COUNT(DISTINCT endpoints.relation_id)::bigint,
              0::bigint,
              0::bigint,
              COUNT(DISTINCT endpoints.relation_id)::bigint
            FROM _refresh_dirty_entity d
            JOIN (
              SELECT subject_entity_id AS entity_id, relation_id
              FROM {}.relation
              UNION ALL
              SELECT object_entity_id AS entity_id, relation_id
              FROM {}.relation
            ) endpoints
              ON endpoints.entity_id = d.entity_id
            GROUP BY d.entity_id
            """
        ).format(schema_id, schema_id, schema_id)
    )
    return int(cur.rowcount)


def _clear_source_derived_facets(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source: str,
) -> None:
    schema_id = sql.Identifier(schema)
    for table in ('facet_entity_bitmap', 'facet_relation_bitmap'):
        if not _table_exists(cur, schema, table):
            continue
        cur.execute(
            sql.SQL(
                """
                DELETE FROM {}.{}
                WHERE facet_name = 'source'
                  AND facet_value = %s
                """
            ).format(schema_id, sql.Identifier(table)),
            [source],
        )


def _mark_resource_not_built(
    cur: psycopg2.extensions.cursor,
    *,
    schema: str,
    source: str,
) -> None:
    if not _table_exists(cur, schema, 'resources'):
        return
    cur.execute(
        sql.SQL(
            """
            UPDATE {}.resources
            SET entity_count = 0,
                interaction_count = 0,
                association_count = 0,
                identifier_count = 0,
                ontology_term_count = 0,
                total_size_bytes = 0,
                last_built_at = NULL,
                build_status = 'not_built'
            WHERE resource_id = %s
            """
        ).format(sql.Identifier(schema)),
        [source],
    )


def _table_exists(
    cur: psycopg2.extensions.cursor,
    schema: str,
    table: str,
) -> bool:
    cur.execute('SELECT to_regclass(%s) IS NOT NULL', [f'{schema}.{table}'])
    return bool(cur.fetchone()[0])


def _count_temp_table(
    cur: psycopg2.extensions.cursor,
    table: str,
) -> int:
    cur.execute(sql.SQL('SELECT COUNT(*) FROM {}').format(sql.Identifier(table)))
    return int(cur.fetchone()[0])
