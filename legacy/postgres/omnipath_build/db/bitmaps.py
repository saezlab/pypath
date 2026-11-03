"""Roaring bitmap indexes for fast entity and relation filtering.

Bitmap tables encode ontology annotations, source facets, entity types, and
relation categories as compressed sets of canonical entity or relation IDs.
They are derived from canonical graph tables and refreshed in full by the
derive phase.
"""

from __future__ import annotations

from dataclasses import dataclass
import time

import psycopg2.extensions
from psycopg2 import sql


@dataclass(frozen=True)
class BitmapStats:
    """Summary counts from bitmap table population."""

    annotation_term_entities: int = 0
    annotation_term_direct_relations: int = 0
    entity_relations: int = 0
    entity_facets: int = 0
    relation_facets: int = 0


def rebuild_bitmap_tables(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    progress: bool = False,
) -> BitmapStats:
    """Create and fully rebuild bitmap index tables."""

    started = time.perf_counter()
    with conn.cursor() as cur:
        _log(progress, 'create_tables', 'start', schema=schema)
        step_started = time.perf_counter()
        _create_bitmap_tables(cur, schema)
        _log(
            progress,
            'create_tables',
            'done',
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )

        _log(progress, 'bitmap_ids', 'start')
        step_started = time.perf_counter()
        _refresh_bitmap_id_tables(cur, schema)
        _log(
            progress,
            'bitmap_ids',
            'done',
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )

        _log(progress, 'annotation_term_entity_bitmap', 'start')
        step_started = time.perf_counter()
        annotation_term_entities = _populate_annotation_term_entity_bitmap(
            cur, schema
        )
        _log(
            progress,
            'annotation_term_entity_bitmap',
            'done',
            rows=annotation_term_entities,
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )

        _log(progress, 'entity_relation_bitmap', 'start')
        step_started = time.perf_counter()
        entity_relations = _populate_entity_relation_bitmap(cur, schema)
        _log(
            progress,
            'entity_relation_bitmap',
            'done',
            rows=entity_relations,
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )

        _log(progress, 'annotation_term_direct_relation_bitmap', 'start')
        step_started = time.perf_counter()
        annotation_term_direct_relations = (
            _populate_annotation_term_direct_relation_bitmap(cur, schema)
        )
        _log(
            progress,
            'annotation_term_direct_relation_bitmap',
            'done',
            rows=annotation_term_direct_relations,
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )

        _log(progress, 'entity_search_counts', 'start')
        step_started = time.perf_counter()
        _populate_entity_search_counts(cur, schema)
        _log(
            progress,
            'entity_search_counts',
            'done',
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )

        _log(progress, 'facet_entity_bitmap', 'start')
        step_started = time.perf_counter()
        entity_facets = _populate_facet_entity_bitmap(cur, schema)
        _log(
            progress,
            'facet_entity_bitmap',
            'done',
            rows=entity_facets,
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )

        _log(progress, 'facet_relation_bitmap', 'start')
        step_started = time.perf_counter()
        relation_facets = _populate_facet_relation_bitmap(cur, schema)
        _log(
            progress,
            'facet_relation_bitmap',
            'done',
            rows=relation_facets,
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )

        _log(progress, 'indexes', 'start')
        step_started = time.perf_counter()
        _create_bitmap_indexes(cur, schema)
        _log(
            progress,
            'indexes',
            'done',
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )
    conn.commit()
    _log(progress, 'all', 'done', seconds=f'{time.perf_counter() - started:.3f}')
    return BitmapStats(
        annotation_term_entities=annotation_term_entities,
        annotation_term_direct_relations=annotation_term_direct_relations,
        entity_relations=entity_relations,
        entity_facets=entity_facets,
        relation_facets=relation_facets,
    )


def _log(progress: bool, step: str, event: str, **fields: object) -> None:
    if not progress:
        return
    details = ' '.join(f'{key}={value}' for key, value in fields.items())
    print(
        f'[derive-bitmaps] step={step} event={event}'
        + (f' {details}' if details else ''),
        flush=True,
    )


def _create_bitmap_tables(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute('CREATE EXTENSION IF NOT EXISTS roaringbitmap')
    cur.execute(
        sql.SQL('DROP TABLE IF EXISTS {}.annotation_term_relation_bitmap').format(
            schema_id
        )
    )
    for table in (
        'annotation_term_entity_bitmap',
        'annotation_term_direct_relation_bitmap',
        'entity_relation_bitmap',
    ):
        cur.execute(
            """
            SELECT data_type
            FROM information_schema.columns
            WHERE table_schema = %s
              AND table_name = %s
              AND column_name = 'term_entity_id'
            """,
            [schema, table],
        )
        row = cur.fetchone()
        if row is not None and row[0] != 'uuid':
            cur.execute(
                sql.SQL('DROP TABLE {}.{}').format(
                    schema_id,
                    sql.Identifier(table),
                )
            )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_bitmap_id (
              entity_id uuid PRIMARY KEY,
              bitmap_id integer NOT NULL
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.relation_bitmap_id (
              relation_id uuid PRIMARY KEY,
              bitmap_id integer NOT NULL
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.annotation_term_entity_bitmap (
              term_entity_id uuid PRIMARY KEY,
              entity_bitmap roaringbitmap NOT NULL,
              global_count integer NOT NULL
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.annotation_term_direct_relation_bitmap (
              term_entity_id uuid PRIMARY KEY,
              relation_bitmap roaringbitmap NOT NULL,
              global_count integer NOT NULL
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_relation_bitmap (
              entity_id uuid PRIMARY KEY,
              relation_bitmap roaringbitmap NOT NULL,
              global_count integer NOT NULL
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_relation_counts (
              entity_id uuid PRIMARY KEY,
              relation_count bigint NOT NULL,
              ontology_annotated_entity_count bigint NOT NULL,
              ontology_annotated_relation_count bigint NOT NULL,
              search_count bigint NOT NULL
            )
            """
        ).format(schema_id)
    )
    _drop_legacy_bitmap_constraints(cur, schema)
    for column_name in (
        'ontology_annotated_entity_count',
        'ontology_annotated_relation_count',
        'search_count',
    ):
        cur.execute(
            sql.SQL(
                'ALTER TABLE {}.entity_relation_counts '
                'ADD COLUMN IF NOT EXISTS {} bigint NOT NULL DEFAULT 0'
            ).format(schema_id, sql.Identifier(column_name))
        )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.facet_entity_bitmap (
              facet_name text NOT NULL,
              facet_value text NOT NULL,
              entity_bitmap roaringbitmap NOT NULL,
              entity_count integer NOT NULL,
              PRIMARY KEY (facet_name, facet_value)
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.facet_relation_bitmap (
              facet_name text NOT NULL,
              facet_value text NOT NULL,
              facet_category text NOT NULL DEFAULT '',
              relation_bitmap roaringbitmap NOT NULL,
              relation_count integer NOT NULL,
              PRIMARY KEY (facet_name, facet_value, facet_category)
            )
            """
        ).format(schema_id)
    )


def _drop_legacy_bitmap_constraints(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """Remove old defensive constraints that slow full bitmap rebuilds."""

    schema_id = sql.Identifier(schema)
    for table, constraint in (
        ('entity_bitmap_id', 'entity_bitmap_id_entity_id_fkey'),
        ('entity_bitmap_id', 'entity_bitmap_id_bitmap_id_key'),
        ('relation_bitmap_id', 'relation_bitmap_id_relation_id_fkey'),
        ('relation_bitmap_id', 'relation_bitmap_id_bitmap_id_key'),
        (
            'annotation_term_entity_bitmap',
            'annotation_term_entity_bitmap_term_entity_id_fkey',
        ),
        (
            'annotation_term_direct_relation_bitmap',
            'annotation_term_direct_relation_bitmap_term_entity_id_fkey',
        ),
        ('entity_relation_bitmap', 'entity_relation_bitmap_entity_id_fkey'),
        ('entity_relation_counts', 'entity_relation_counts_entity_id_fkey'),
    ):
        cur.execute(
            sql.SQL('ALTER TABLE {}.{} DROP CONSTRAINT IF EXISTS {}').format(
                schema_id,
                sql.Identifier(table),
                sql.Identifier(constraint),
            )
        )


def _refresh_bitmap_id_tables(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(sql.SQL('TRUNCATE {}.entity_bitmap_id').format(schema_id))
    cur.execute(
        sql.SQL(
            """
            INSERT INTO {}.entity_bitmap_id (entity_id, bitmap_id)
            SELECT
              entity_id,
              row_number() OVER (ORDER BY entity_id)::integer
            FROM {}.entity
            """
        ).format(schema_id, schema_id)
    )
    cur.execute(sql.SQL('TRUNCATE {}.relation_bitmap_id').format(schema_id))
    cur.execute(
        sql.SQL(
            """
            INSERT INTO {}.relation_bitmap_id (relation_id, bitmap_id)
            SELECT
              relation_id,
              row_number() OVER (ORDER BY relation_id)::integer
            FROM {}.relation
            """
        ).format(schema_id, schema_id)
    )


def _populate_annotation_term_entity_bitmap(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> int:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL('TRUNCATE {}.annotation_term_entity_bitmap').format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            INSERT INTO {}.annotation_term_entity_bitmap (
              term_entity_id,
              entity_bitmap,
              global_count
            )
            WITH ontology_terms AS MATERIALIZED (
              SELECT DISTINCT term_entity_id
              FROM {}.entity_ontology_term
            ),
            term_entities AS MATERIALIZED (
              SELECT
                r.subject_entity_id AS term_entity_id,
                r.object_entity_id AS entity_id
              FROM {}.relation r
              JOIN ontology_terms terms
                ON terms.term_entity_id = r.subject_entity_id
              JOIN {}.vocab_relation_category rc
                ON rc.relation_category_id = r.relation_category_id
              WHERE rc.name = 'association'
              UNION
              SELECT
                term_entity_id,
                term_entity_id AS entity_id
              FROM ontology_terms
              UNION
              SELECT
                eor.object_entity_id AS term_entity_id,
                eor.subject_entity_id AS entity_id
              FROM {}.entity_ontology_relation eor
            )
            SELECT
              term_entities.term_entity_id,
              rb_build_agg(DISTINCT object_bitmap.bitmap_id),
              COUNT(DISTINCT term_entities.entity_id)::integer
            FROM term_entities
            JOIN {}.entity_bitmap_id object_bitmap
              ON object_bitmap.entity_id = term_entities.entity_id
            GROUP BY term_entities.term_entity_id
            """
        ).format(
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
        )
    )
    return int(cur.rowcount)


def _populate_entity_relation_bitmap(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> int:
    schema_id = sql.Identifier(schema)
    cur.execute(sql.SQL('TRUNCATE {}.entity_relation_bitmap').format(schema_id))
    cur.execute(
        sql.SQL(
            """
            WITH subject_relation_bitmap AS (
              SELECT
                r.subject_entity_id AS entity_id,
                rb_build_agg(relation_bitmap.bitmap_id) AS relation_bitmap
              FROM {}.relation r
              JOIN {}.relation_bitmap_id relation_bitmap
                ON relation_bitmap.relation_id = r.relation_id
              GROUP BY r.subject_entity_id
            ),
            object_relation_bitmap AS (
              SELECT
                r.object_entity_id AS entity_id,
                rb_build_agg(relation_bitmap.bitmap_id) AS relation_bitmap
              FROM {}.relation r
              JOIN {}.relation_bitmap_id relation_bitmap
                ON relation_bitmap.relation_id = r.relation_id
              GROUP BY r.object_entity_id
            ),
            endpoint_relation_bitmap AS (
              SELECT
                entity_id,
                rb_or_agg(relation_bitmap) AS relation_bitmap
              FROM (
                SELECT entity_id, relation_bitmap
                FROM subject_relation_bitmap
                UNION ALL
                SELECT entity_id, relation_bitmap
                FROM object_relation_bitmap
              ) endpoint
              GROUP BY entity_id
            )
            INSERT INTO {}.entity_relation_bitmap (
              entity_id,
              relation_bitmap,
              global_count
            )
            SELECT
              entity_id,
              relation_bitmap,
              rb_cardinality(relation_bitmap)::integer
            FROM endpoint_relation_bitmap
            """
        ).format(schema_id, schema_id, schema_id, schema_id, schema_id)
    )
    return int(cur.rowcount)


def _populate_annotation_term_direct_relation_bitmap(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> int:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL('TRUNCATE {}.annotation_term_direct_relation_bitmap').format(
            schema_id
        )
    )
    cur.execute(
        sql.SQL(
            """
            WITH term_lookup AS MATERIALIZED (
              SELECT DISTINCT term_entity_id, term_id AS term
              FROM {}.entity_ontology_term
              WHERE term_id <> ''
              UNION
              SELECT DISTINCT term_entity_id, alias.value AS term
              FROM {}.entity_ontology_term terms
              CROSS JOIN LATERAL unnest(terms.term_aliases) AS alias(value)
              WHERE alias.value <> ''
            ),
            direct_relation_terms AS MATERIALIZED (
              SELECT DISTINCT
                relation_link.relation_id,
                annotation.term
              FROM {}.relation_evidence_annotation annotation_link
              JOIN {}.annotation annotation
                ON annotation.annotation_key = annotation_link.annotation_key
              JOIN {}.relation_evidence_relation relation_link
                ON relation_link.source_id = annotation_link.source_id
               AND relation_link.relation_evidence_id =
                   annotation_link.relation_evidence_id
              WHERE annotation.term <> ''
            ),
            direct_relation_annotations AS (
              SELECT
                term_lookup.term_entity_id,
                rb_build_agg(DISTINCT relation_bitmap.bitmap_id) AS relation_bitmap
              FROM direct_relation_terms
              JOIN term_lookup
                ON term_lookup.term = direct_relation_terms.term
              JOIN {}.relation_bitmap_id relation_bitmap
                ON relation_bitmap.relation_id =
                   direct_relation_terms.relation_id
              GROUP BY term_lookup.term_entity_id
            )
            INSERT INTO {}.annotation_term_direct_relation_bitmap (
              term_entity_id,
              relation_bitmap,
              global_count
            )
            SELECT
              term_entity_id,
              relation_bitmap,
              rb_cardinality(relation_bitmap)::integer
            FROM direct_relation_annotations
            """
        ).format(
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
        )
    )
    return int(cur.rowcount)


def _populate_entity_search_counts(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> int:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            UPDATE {}.entity_relation_counts
            SET
              ontology_annotated_entity_count = 0,
              ontology_annotated_relation_count = 0,
              search_count = relation_count
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
              e.entity_id,
              COALESCE(rc.relation_count, 0)::bigint AS relation_count,
              COALESCE(entity_bitmap.global_count, 0)::bigint
                AS ontology_annotated_entity_count,
              COALESCE(relation_bitmap.global_count, 0)::bigint
                AS ontology_annotated_relation_count,
              (
                COALESCE(rc.relation_count, 0)
                + COALESCE(entity_bitmap.global_count, 0)
                + COALESCE(relation_bitmap.global_count, 0)
              )::bigint AS search_count
            FROM {}.entity e
            LEFT JOIN {}.entity_relation_counts rc
              ON rc.entity_id = e.entity_id
            LEFT JOIN {}.annotation_term_entity_bitmap entity_bitmap
              ON entity_bitmap.term_entity_id = e.entity_id
            LEFT JOIN {}.annotation_term_direct_relation_bitmap relation_bitmap
              ON relation_bitmap.term_entity_id = e.entity_id
            WHERE entity_bitmap.term_entity_id IS NOT NULL
               OR relation_bitmap.term_entity_id IS NOT NULL
            ON CONFLICT (entity_id) DO UPDATE
            SET
              ontology_annotated_entity_count =
                EXCLUDED.ontology_annotated_entity_count,
              ontology_annotated_relation_count =
                EXCLUDED.ontology_annotated_relation_count,
              search_count =
                {}.entity_relation_counts.relation_count
                + EXCLUDED.ontology_annotated_entity_count
                + EXCLUDED.ontology_annotated_relation_count
            """
        ).format(
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
        )
    )
    return int(cur.rowcount)


def _populate_facet_entity_bitmap(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> int:
    schema_id = sql.Identifier(schema)
    cur.execute(sql.SQL('TRUNCATE {}.facet_entity_bitmap').format(schema_id))
    total = 0
    statements = [
        sql.SQL(
            """
            INSERT INTO {}.facet_entity_bitmap (
              facet_name,
              facet_value,
              entity_bitmap,
              entity_count
            )
            SELECT
              'entity_type',
              et.name,
              rb_build_agg(bitmap.bitmap_id),
              COUNT(*)::integer
            FROM {}.entity e
            JOIN {}.entity_bitmap_id bitmap
              ON bitmap.entity_id = e.entity_id
            JOIN {}.vocab_entity_type et
              ON et.entity_type_id = e.entity_type_id
            GROUP BY et.name
            """
        ).format(schema_id, schema_id, schema_id, schema_id),
        sql.SQL(
            """
            INSERT INTO {}.facet_entity_bitmap (
              facet_name,
              facet_value,
              entity_bitmap,
              entity_count
            )
            SELECT
              'chemical_class',
              vcc.name,
              rb_build_agg(bitmap.bitmap_id),
              COUNT(*)::integer
            FROM {}.entity e
            JOIN {}.entity_bitmap_id bitmap
              ON bitmap.entity_id = e.entity_id
            JOIN {}.vocab_chemical_class vcc
              ON vcc.chemical_class_id = e.chemical_class_id
            GROUP BY vcc.name
            """
        ).format(schema_id, schema_id, schema_id, schema_id),
        sql.SQL(
            """
            INSERT INTO {}.facet_entity_bitmap (
              facet_name,
              facet_value,
              entity_bitmap,
              entity_count
            )
            SELECT
              'metabolic_domain',
              vmd.name,
              rb_build_agg(bitmap.bitmap_id),
              COUNT(*)::integer
            FROM {}.entity e
            JOIN {}.entity_bitmap_id bitmap
              ON bitmap.entity_id = e.entity_id
            JOIN {}.vocab_metabolic_domain vmd
              ON vmd.metabolic_domain_id = e.metabolic_domain_id
            GROUP BY vmd.name
            """
        ).format(schema_id, schema_id, schema_id, schema_id),
        sql.SQL(
            """
            INSERT INTO {}.facet_entity_bitmap (
              facet_name,
              facet_value,
              entity_bitmap,
              entity_count
            )
            SELECT
              'taxonomy_id',
              taxonomy_id,
              rb_build_agg(bitmap.bitmap_id),
              COUNT(*)::integer
            FROM {}.entity
            JOIN {}.entity_bitmap_id bitmap
              ON bitmap.entity_id = entity.entity_id
            WHERE taxonomy_id IS NOT NULL
            GROUP BY taxonomy_id
            """
        ).format(schema_id, schema_id, schema_id),
        sql.SQL(
            """
            INSERT INTO {}.facet_entity_bitmap (
              facet_name,
              facet_value,
              entity_bitmap,
              entity_count
            )
            WITH entity_sources AS (
              SELECT DISTINCT r.entity_id, ds.name AS source
              FROM {}.entity_evidence_resolution r
              JOIN {}.entity_evidence ee
                ON ee.source_id = r.source_id
               AND ee.entity_evidence_id = r.entity_evidence_id
              JOIN {}.data_source ds
                ON ds.source_id = ee.source_id
              WHERE r.entity_id IS NOT NULL
              UNION
              SELECT DISTINCT re.subject_entity_id AS entity_id, ds.name AS source
              FROM {}.relation_evidence re
              JOIN {}.data_source ds
                ON ds.source_id = re.source_id
              WHERE re.subject_entity_id IS NOT NULL
              UNION
              SELECT DISTINCT re.object_entity_id AS entity_id, ds.name AS source
              FROM {}.relation_evidence re
              JOIN {}.data_source ds
                ON ds.source_id = re.source_id
              WHERE re.object_entity_id IS NOT NULL
              UNION
              SELECT DISTINCT ot.term_entity_id AS entity_id, source.value AS source
              FROM {}.entity_ontology_term ot
              CROSS JOIN LATERAL unnest(ot.sources) AS source(value)
            )
            SELECT
              'source',
              source,
              rb_build_agg(bitmap.bitmap_id),
              COUNT(*)::integer
            FROM entity_sources
            JOIN {}.entity_bitmap_id bitmap
              ON bitmap.entity_id = entity_sources.entity_id
            WHERE source IS NOT NULL
              AND source <> ''
            GROUP BY source
            """
        ).format(
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
        ),
        sql.SQL(
            """
            INSERT INTO {}.facet_entity_bitmap (
              facet_name,
              facet_value,
              entity_bitmap,
              entity_count
            )
            SELECT
              'ontology_id',
              ot.ontology_id,
              rb_build_agg(DISTINCT bitmap.bitmap_id),
              COUNT(DISTINCT ot.term_entity_id)::integer
            FROM {}.entity_ontology_term ot
            JOIN {}.entity_bitmap_id bitmap
              ON bitmap.entity_id = ot.term_entity_id
            WHERE COALESCE(ot.ontology_id, '') <> ''
            GROUP BY ot.ontology_id
            """
        ).format(
            schema_id,
            schema_id,
            schema_id,
        ),
    ]
    for statement in statements:
        cur.execute(statement)
        total += int(cur.rowcount)
    return total


def _populate_facet_relation_bitmap(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> int:
    schema_id = sql.Identifier(schema)
    cur.execute(sql.SQL('TRUNCATE {}.facet_relation_bitmap').format(schema_id))
    total = 0
    statements = [
        sql.SQL(
            """
            INSERT INTO {}.facet_relation_bitmap (
              facet_name,
              facet_value,
              facet_category,
              relation_bitmap,
              relation_count
            )
            SELECT
              'predicate',
              rp.name,
              COALESCE(rc.name, ''),
              rb_build_agg(bitmap.bitmap_id),
              COUNT(*)::integer
            FROM {}.relation r
            JOIN {}.relation_bitmap_id bitmap
              ON bitmap.relation_id = r.relation_id
            JOIN {}.vocab_relation_predicate rp
              ON rp.relation_predicate_id = r.predicate_id
            LEFT JOIN {}.vocab_relation_category rc
              ON rc.relation_category_id = r.relation_category_id
            GROUP BY rp.name, rc.name
            """
        ).format(schema_id, schema_id, schema_id, schema_id, schema_id),
        sql.SQL(
            """
            INSERT INTO {}.facet_relation_bitmap (
              facet_name,
              facet_value,
              relation_bitmap,
              relation_count
            )
            WITH relation_sources AS (
              SELECT DISTINCT rer.relation_id, ds.name AS source
              FROM {}.relation_evidence_relation rer
              JOIN {}.data_source ds
                ON ds.source_id = rer.source_id
              UNION
              SELECT DISTINCT ear.relation_id, ds.name AS source
              FROM {}.entity_annotation_relation ear
              JOIN {}.data_source ds
                ON ds.source_id = ear.source_id
            )
            SELECT
              'source',
              source,
              rb_build_agg(bitmap.bitmap_id),
              COUNT(*)::integer
            FROM relation_sources
            JOIN {}.relation_bitmap_id bitmap
              ON bitmap.relation_id = relation_sources.relation_id
            WHERE source IS NOT NULL
              AND source <> ''
            GROUP BY source
            """
        ).format(
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
        ),
        sql.SQL(
            """
            INSERT INTO {}.facet_relation_bitmap (
              facet_name,
              facet_value,
              relation_bitmap,
              relation_count
            )
            WITH participant_types AS (
              SELECT r.relation_id, subject_type.name AS vocab_entity_type
              FROM {}.relation r
              JOIN {}.entity subject
                ON subject.entity_id = r.subject_entity_id
              JOIN {}.vocab_entity_type subject_type
                ON subject_type.entity_type_id = subject.entity_type_id
              UNION
              SELECT r.relation_id, object_type.name AS vocab_entity_type
              FROM {}.relation r
              JOIN {}.entity object
                ON object.entity_id = r.object_entity_id
              JOIN {}.vocab_entity_type object_type
                ON object_type.entity_type_id = object.entity_type_id
            )
            SELECT
              'participant_type',
              vocab_entity_type,
              rb_build_agg(bitmap.bitmap_id),
              COUNT(DISTINCT participant_types.relation_id)::integer
            FROM participant_types
            JOIN {}.relation_bitmap_id bitmap
              ON bitmap.relation_id = participant_types.relation_id
            WHERE vocab_entity_type IS NOT NULL
            GROUP BY vocab_entity_type
            """
        ).format(
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
        ),
        sql.SQL(
            """
            INSERT INTO {}.facet_relation_bitmap (
              facet_name,
              facet_value,
              relation_bitmap,
              relation_count
            )
            WITH relation_taxonomy AS (
              SELECT r.relation_id, subject.taxonomy_id
              FROM {}.relation r
              JOIN {}.entity subject
                ON subject.entity_id = r.subject_entity_id
              WHERE subject.taxonomy_id IS NOT NULL
              UNION
              SELECT r.relation_id, object.taxonomy_id
              FROM {}.relation r
              JOIN {}.entity object
                ON object.entity_id = r.object_entity_id
              WHERE object.taxonomy_id IS NOT NULL
            )
            SELECT
              'taxonomy_id',
              taxonomy_id,
              rb_build_agg(bitmap.bitmap_id),
              COUNT(*)::integer
            FROM relation_taxonomy
            JOIN {}.relation_bitmap_id bitmap
              ON bitmap.relation_id = relation_taxonomy.relation_id
            GROUP BY taxonomy_id
            """
        ).format(
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
            schema_id,
        ),
    ]
    for statement in statements:
        cur.execute(statement)
        total += int(cur.rowcount)
    return total


def _create_bitmap_indexes(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    statements = [
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS annotation_term_entity_count_idx
            ON {}.annotation_term_entity_bitmap (global_count DESC)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS annotation_term_direct_relation_count_idx
            ON {}.annotation_term_direct_relation_bitmap (global_count DESC)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS entity_relation_count_idx
            ON {}.entity_relation_bitmap (global_count DESC)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS entity_relation_counts_search_count_idx
            ON {}.entity_relation_counts (search_count DESC, entity_id ASC)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS facet_entity_count_idx
            ON {}.facet_entity_bitmap (
              facet_name,
              entity_count DESC,
              facet_value
            )
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS facet_relation_count_idx
            ON {}.facet_relation_bitmap (
              facet_name,
              relation_count DESC,
              facet_value
            )
            """
        ).format(schema_id),
    ]
    for statement in statements:
        cur.execute(statement)
