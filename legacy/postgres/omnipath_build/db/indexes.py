"""Secondary PostgreSQL indexes for query-ready pipeline outputs."""

from __future__ import annotations

from psycopg2 import sql
import psycopg2.extensions

def create_secondary_indexes(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> None:
    """Create indexes for derived tables, annotations, and graph lookups."""

    schema_id = sql.Identifier(schema)
    statements = [
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS entity_evidence_source_dataset_row_idx
            ON {}.entity_evidence (source_id, dataset_id, row_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS entity_evidence_type_taxonomy_idx
            ON {}.entity_evidence (entity_type_id, taxonomy_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS entity_type_taxonomy_idx
            ON {}.entity (entity_type_id, taxonomy_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS canonical_relation_subject_entity_idx
            ON {}.relation (subject_entity_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS canonical_relation_object_entity_idx
            ON {}.relation (object_entity_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS relation_predicate_category_idx
            ON {}.relation (predicate_id, relation_category_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS relation_category_subject_idx
            ON {}.relation (relation_category_id, subject_entity_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS relation_category_object_idx
            ON {}.relation (relation_category_id, object_entity_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS relation_evidence_predicate_category_idx
            ON {}.relation_evidence (predicate_id, relation_category_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS relation_evidence_source_dataset_row_idx
            ON {}.relation_evidence (source_id, dataset_id, row_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS relation_evidence_relation_relation_id_idx
            ON {}.relation_evidence_relation (relation_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS entity_annotation_relation_relation_id_idx
            ON {}.entity_annotation_relation (relation_id)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS relation_evidence_annotation_annotation_key_idx
            ON {}.relation_evidence_annotation (annotation_key)
            """
        ).format(schema_id),
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS entity_evidence_annotation_annotation_key_idx
            ON {}.entity_evidence_annotation (annotation_key)
            """
        ).format(schema_id),
    ]
    with conn.cursor() as cur:
        for statement in statements:
            cur.execute(statement)
    conn.commit()
