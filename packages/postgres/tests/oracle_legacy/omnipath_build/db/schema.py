"""PostgreSQL DDL for source evidence, resolution, and graph tables.

The schema separates refresh-loaded source evidence from canonical graph
materialization. Evidence tables are source-partitioned and keep source row
identity, participant occurrence identity, identifiers, annotations, and
relation endpoints. Resolution tables map those evidence occurrences to
canonical entities, while relation mapping tables connect source relation
evidence to deduplicated graph relations.

Several high-volume secondary indexes are treated as deferred content indexes.
Scratch builds can create the schema without them, stream source evidence
quickly, and then create the indexes once before canonicalization or derivation.
"""

from __future__ import annotations

import logging

from psycopg2 import sql
import psycopg2.extensions

# Named for the module rather than taken from ``__name__``: this module is only
# ever imported, so the two agree, but naming it keeps the logger inside the
# ``omnipath_build`` subtree that ``_session`` configures.
_logger = logging.getLogger('omnipath_build.db.schema')

CONTENT_TABLES: tuple[str, ...] = (
    # Published datasets over the core content come first: they are keyed on
    # entity ids the truncate below invalidates, and they carry no foreign key
    # that would make a stale row announce itself.
    'metsigdb_membership',
    'annotation_term_entity_bitmap',
    'annotation_term_direct_relation_bitmap',
    'entity_relation_bitmap',
    'facet_entity_bitmap',
    'facet_relation_bitmap',
    'entity_bitmap_id',
    'relation_bitmap_id',
    'entity_identifier_lookup',
    'entity_identifier',
    'entity_relation_counts',
    'entity_ontology_term',
    'entity_ontology_relation',
    'ontology_terms',
    'entity_name',
    'relation_annotation',
    'relation_evidence_annotation',
    'entity_annotation_relation',
    'entity_annotation',
    'entity_evidence_annotation',
    'relation_evidence_relation',
    'entity_evidence_resolution',
    'resolution_conflict',
    # Gene-anchored state model (spec 011's predecessor, "spec-003 Phase 6"):
    # not previously listed here, so reset-content never truncated these --
    # stale rows survived every wipe and eventually collided with a later
    # PRIMARY KEY add on evidence_state. Found live on 2026-09-09.
    'evidence_state',
    'state_component',
    'state',
    'gene_protein_representative',
    'annotation',
    'interaction_fact_resource',
    'interaction_party',
    'interaction',
    'relation',
    'relation_evidence',
    'entity_evidence_identifier',
    'entity_evidence',
    'identifier_evidence',
    'entity',
    'resources',
    'dataset',
    'data_source_license',
    'entity_annotation_relation_default',
    'identifier_authority',
    'identifier_role',
    'chemical_resolution_coverage',
    'lipid_name_edge',
    'lipid_name_node',
    'resource_overlap_summary',
    # omnipath-metabo's WP6 QC layer (structure_consistency_finding
    # references data_source; found live during the cycle's final
    # reset-content -- same "not previously listed here" bug class as
    # evidence_state/state_component above, T098's tables just never got
    # added when they landed in the other repo).
    'structure_consistency_finding',
    'structure_consistency_summary',
    'data_source',
)

SOURCE_PARTITIONED_TABLES: tuple[str, ...] = (
    'ontology_terms',
    'entity_ontology_relation',
    'entity_name',
    'entity_evidence',
    'entity_evidence_identifier',
    'entity_identifier',
    'relation_evidence',
    'entity_evidence_annotation',
    'relation_evidence_annotation',
    'entity_evidence_resolution',
    'resolution_conflict',
    'relation_evidence_relation',
    'entity_annotation_relation',
    'evidence_state',
)

SOURCE_PARTITION_DROP_ORDER: tuple[str, ...] = (
    'ontology_terms',
    'entity_ontology_relation',
    'entity_name',
    'relation_evidence_annotation',
    'relation_evidence_relation',
    'entity_annotation_relation',
    'relation_evidence',
    'entity_evidence_annotation',
    'entity_evidence_resolution',
    'resolution_conflict',
    'entity_identifier',
    'entity_evidence_identifier',
    'evidence_state',
    'entity_evidence',
)

CONTENT_PRIMARY_KEYS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ('identifier_evidence', ('identifier_id',)),
    ('annotation', ('annotation_key',)),
    ('entity', ('entity_id',)),
    ('entity_evidence', ('source_id', 'entity_evidence_id')),
    (
        'entity_evidence_identifier',
        ('source_id', 'entity_evidence_id', 'identifier_id'),
    ),
    ('entity_identifier', ('source_id', 'entity_id', 'identifier_id')),
    (
        'entity_evidence_annotation',
        ('source_id', 'entity_evidence_id', 'annotation_key'),
    ),
    ('relation_evidence', ('source_id', 'relation_evidence_id')),
    (
        'relation_evidence_annotation',
        (
            'source_id',
            'relation_evidence_id',
            'annotation_key',
            'annotation_scope_id',
        ),
    ),
    ('entity_evidence_resolution', ('source_id', 'entity_evidence_id')),
    (
        'resolution_conflict',
        ('source_id', 'entity_evidence_id', 'candidate_structure'),
    ),
    ('relation', ('relation_id',)),
    (
        'relation_evidence_relation',
        ('source_id', 'relation_evidence_id'),
    ),
    ('ontology_terms', ('source_id', 'term_entity_id')),
    (
        'entity_ontology_relation',
        (
            'source_id',
            'subject_entity_id',
            'predicate_id',
            'object_entity_id',
            'ontology_id',
        ),
    ),
    (
        'entity_name',
        ('source_id', 'entity_id', 'name_normalized'),
    ),
    ('gene_protein_representative', ('entity_id',)),
    ('state', ('state_id',)),
    ('state_component', ('state_id', 'component_type', 'value')),
    ('evidence_state', ('source_id', 'entity_evidence_id', 'state_id')),
)


def ensure_schema(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    drop_existing: bool = False,
    progress: bool = False,
    indexes: bool = True,
) -> None:
    """Create or refresh the omnipath_build evidence and resolution schema."""

    def log_step(message: str) -> None:
        if progress:
            print(f'[schema] {message}', flush=True)

    started = None
    if progress:
        import time

        started = time.perf_counter()
        log_step(
            f'ensure start schema={schema} drop_existing={drop_existing} '
            f'indexes={indexes}'
        )

    with conn.cursor() as cur:
        if drop_existing:
            log_step('drop existing schema')
            cur.execute(
                sql.SQL('DROP SCHEMA IF EXISTS {} CASCADE').format(
                    sql.Identifier(schema)
                )
            )

        log_step('create schema')
        cur.execute(
            sql.SQL('CREATE SCHEMA IF NOT EXISTS {}').format(
                sql.Identifier(schema)
            )
        )

        log_step('rename legacy vocabulary tables')
        _rename_legacy_vocab_tables(cur, schema)
        log_step('create identifier type table')
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.vocab_identifier_type (
                  identifier_type_id bigint PRIMARY KEY,
                  name text NOT NULL UNIQUE
                )
                """
            ).format(sql.Identifier(schema))
        )
        _ensure_static_identifier_types(cur, schema)
        log_step('create identifier evidence table')
        cur.execute(
            sql.SQL(
                """
                DO $$
                BEGIN
                  IF to_regclass({old_table}) IS NOT NULL
                     AND to_regclass({new_table}) IS NULL THEN
                    ALTER TABLE {old_table_sql} RENAME TO identifier_evidence;
                  END IF;
                END
                $$;
                """
            ).format(
                old_table=sql.Literal(f'{schema}.identifier'),
                new_table=sql.Literal(f'{schema}.identifier_evidence'),
                old_table_sql=sql.SQL('{}.identifier').format(
                    sql.Identifier(schema)
                ),
            )
        )
        _drop_legacy_content_tables_if_needed(cur, schema)
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.identifier_evidence (
                  identifier_id uuid PRIMARY KEY,
                  identifier_type_id bigint NOT NULL
                    REFERENCES {}.vocab_identifier_type(identifier_type_id),
                  value text NOT NULL
                )
                """
            ).format(sql.Identifier(schema), sql.Identifier(schema))
        )
        _ensure_primary_key(cur, schema, 'identifier_evidence', ('identifier_id',))
        log_step('ensure identifier evidence shape')
        _ensure_identifier_evidence_key(cur, schema)
        log_step('ensure identifier evidence normalized-value column')
        _ensure_identifier_evidence_value_normalized(cur, schema)
        log_step('ensure identifier type value_pattern/normalization_rule')
        _ensure_identifier_type_value_pattern(cur, schema)
        log_step('create normalized dimension tables')
        _ensure_evidence_dimension_tables(cur, schema)
        log_step('create entity_evidence table')
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.entity_evidence (
                  source_id bigint NOT NULL
                    REFERENCES {}.data_source(source_id),
                  entity_evidence_id uuid NOT NULL,
                  dataset_id bigint NOT NULL
                    REFERENCES {}.dataset(dataset_id),
                  row_id bigint NOT NULL,
                  parent_entity_evidence_id uuid,
                  entity_role_id smallint NOT NULL
                    REFERENCES {}.vocab_entity_role(entity_role_id),
                  entity_type_id bigint
                    REFERENCES {}.vocab_entity_type(entity_type_id),
                  taxonomy_id bigint,
                  PRIMARY KEY (source_id, entity_evidence_id),
                  FOREIGN KEY (source_id, parent_entity_evidence_id)
                    REFERENCES {}.entity_evidence(source_id, entity_evidence_id)
                ) PARTITION BY LIST (source_id)
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
        _ensure_primary_key(
            cur, schema, 'entity_evidence', ('source_id', 'entity_evidence_id')
        )
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.entity_evidence_default
                PARTITION OF {}.entity_evidence DEFAULT
                """
            ).format(sql.Identifier(schema), sql.Identifier(schema))
        )
        log_step('create entity_evidence_identifier table')
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.entity_evidence_identifier (
                  source_id bigint NOT NULL
                    REFERENCES {}.data_source(source_id),
                  entity_evidence_id uuid NOT NULL,
                  identifier_id uuid NOT NULL
                    REFERENCES {}.identifier_evidence(identifier_id),
                  PRIMARY KEY (source_id, entity_evidence_id, identifier_id),
                  FOREIGN KEY (source_id, entity_evidence_id)
                    REFERENCES {}.entity_evidence(source_id, entity_evidence_id)
                    ON DELETE CASCADE
                ) PARTITION BY LIST (source_id)
                """
            ).format(
                sql.Identifier(schema),
                sql.Identifier(schema),
                sql.Identifier(schema),
                sql.Identifier(schema),
            )
        )
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.entity_evidence_identifier_default
                PARTITION OF {}.entity_evidence_identifier DEFAULT
                """
            ).format(sql.Identifier(schema), sql.Identifier(schema))
        )
        log_step('create relation_evidence table')
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.relation_evidence (
                  source_id bigint NOT NULL
                    REFERENCES {}.data_source(source_id),
                  relation_evidence_id uuid NOT NULL,
                  dataset_id bigint NOT NULL
                    REFERENCES {}.dataset(dataset_id),
                  row_id bigint NOT NULL,
                  subject_entity_evidence_id uuid,
                  subject_entity_id uuid,
                  predicate_id bigint NOT NULL
                    REFERENCES {}.vocab_relation_predicate(relation_predicate_id),
                  object_entity_evidence_id uuid,
                  object_entity_id uuid,
                  relation_category_id bigint NOT NULL
                    REFERENCES {}.vocab_relation_category(relation_category_id),
                  PRIMARY KEY (source_id, relation_evidence_id),
                  CHECK (
                    (subject_entity_evidence_id IS NOT NULL)::int
                    + (subject_entity_id IS NOT NULL)::int
                    = 1
                  ),
                  CHECK (
                    (object_entity_evidence_id IS NOT NULL)::int
                    + (object_entity_id IS NOT NULL)::int
                    = 1
                  ),
                  FOREIGN KEY (source_id, subject_entity_evidence_id)
                    REFERENCES {}.entity_evidence(source_id, entity_evidence_id),
                  FOREIGN KEY (source_id, object_entity_evidence_id)
                    REFERENCES {}.entity_evidence(source_id, entity_evidence_id)
                ) PARTITION BY LIST (source_id)
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
            )
        )
        _ensure_primary_key(
            cur, schema, 'relation_evidence', ('source_id', 'relation_evidence_id')
        )
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.relation_evidence_default
                PARTITION OF {}.relation_evidence DEFAULT
                """
            ).format(sql.Identifier(schema), sql.Identifier(schema))
        )
        log_step('create annotation table')
        cur.execute(
            sql.SQL(
                """
                CREATE TABLE IF NOT EXISTS {}.annotation (
                  annotation_key uuid PRIMARY KEY,
                  term text NOT NULL,
                  value text,
                  unit text
                )
                """
            ).format(sql.Identifier(schema))
        )
        _ensure_primary_key(cur, schema, 'annotation', ('annotation_key',))
        log_step('ensure annotation value schema')
        _ensure_annotation_value_schema(cur, schema)
        log_step('create evidence annotation tables')
        _ensure_evidence_annotation_tables(cur, schema)
        _ensure_resolution_schema(
            cur, schema, progress=progress, indexes=indexes
        )

    log_step('commit')
    conn.commit()
    if started is not None:
        import time

        log_step(f'ensure done elapsed={time.perf_counter() - started:.1f}s')


def ensure_content_primary_keys(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    progress: bool = False,
) -> None:
    """Restore content-table primary keys after constraint-free bulk loads."""

    def log_step(message: str) -> None:
        if progress:
            print(f'[schema] {message}', flush=True)

    log_step('ensure content primary keys')
    with conn.cursor() as cur:
        for table, columns in CONTENT_PRIMARY_KEYS:
            cur.execute(
                """
                SELECT to_regclass(%s) IS NOT NULL
                """,
                [f'{schema}.{table}'],
            )
            if not cur.fetchone()[0]:
                continue
            cur.execute(
                """
                SELECT 1
                FROM pg_constraint
                WHERE conrelid = %s::regclass
                  AND contype = 'p'
                """,
                [f'{schema}.{table}'],
            )
            if cur.fetchone() is not None:
                continue
            log_step(f'add primary key {table}')
            column_sql = sql.SQL(', ').join(
                sql.Identifier(column) for column in columns
            )
            cur.execute(
                sql.SQL(
                    'ALTER TABLE {}.{} ADD CONSTRAINT {} PRIMARY KEY ({})'
                ).format(
                    sql.Identifier(schema),
                    sql.Identifier(table),
                    sql.Identifier(f'{table}_pkey'),
                    column_sql,
                )
            )
    conn.commit()


def ensure_deferred_indexes(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    progress: bool = False,
) -> None:
    """Create high-volume evidence indexes used after bulk ingest."""

    def log_step(message: str) -> None:
        if progress:
            print(f'[schema] {message}', flush=True)

    log_step('ensure deferred indexes')
    with conn.cursor() as cur:
        _ensure_resolution_indexes(cur, schema)
    conn.commit()


def ensure_source_partitions(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    source: str,
) -> None:
    """Create source-specific partitions for partitioned evidence tables."""

    suffix = _source_partition_suffix(source)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {}.data_source (name)
                VALUES (%s)
                ON CONFLICT (name) DO UPDATE SET name = EXCLUDED.name
                RETURNING source_id
                """
            ).format(sql.Identifier(schema)),
            [source],
        )
        source_id = int(cur.fetchone()[0])
        for table in SOURCE_PARTITIONED_TABLES:
            cur.execute(
                sql.SQL(
                    """
                    CREATE TABLE IF NOT EXISTS {}.{}
                    PARTITION OF {}.{} FOR VALUES IN ({})
                    """
                ).format(
                    sql.Identifier(schema),
                    sql.Identifier(f'{table}_{suffix}'),
                    sql.Identifier(schema),
                    sql.Identifier(table),
                    sql.Literal(source_id),
                )
            )
    conn.commit()


def _source_partition_suffix(source: str) -> str:
    import re
    import hashlib

    slug = re.sub(r'[^a-z0-9]+', '_', source.lower()).strip('_')
    slug = slug[:40] or 'source'
    digest = hashlib.sha1(source.encode('utf-8')).hexdigest()[:8]
    return f'{slug}_{digest}'


def drop_deferred_content_indexes(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    progress: bool = False,
) -> list[str]:
    """Drop deferred evidence indexes to speed up a subsequent bulk reload."""

    names = [
        'entity_evidence_identifier_identifier_idx',
        'entity_evidence_source_dataset_row_idx',
        'entity_evidence_type_taxonomy_idx',
        'entity_resolution_status_idx',
        'entity_resolution_reason_idx',
        'entity_evidence_resolution_entity_idx',
        'entity_canonical_lookup_idx',
        'entity_status_idx',
        'relation_subject_idx',
        'relation_object_idx',
        'relation_subject_entity_idx',
        'relation_object_entity_idx',
        'relation_source_dataset_idx',
        'relation_evidence_predicate_category_idx',
        'relation_evidence_source_dataset_row_idx',
        'relation_canonical_subject_idx',
        'relation_canonical_object_idx',
        'entity_evidence_annotation_annotation_key_idx',
        'entity_annotation_relation_relation_id_idx',
        'relation_evidence_annotation_relation_evidence_idx',
        'relation_evidence_annotation_annotation_key_idx',
        'entity_annotation_annotation_key_idx',
        'relation_annotation_annotation_key_idx',
        'resources_build_status_idx',
        'identifier_evidence_value_lower_idx',
        'entity_identifier_entity_idx',
        'entity_identifier_identifier_idx',
        'entity_identifier_lookup_entity_idx',
        'entity_identifier_lookup_identifier_id_idx',
        'entity_identifier_lookup_value_idx',
        'entity_identifier_lookup_value_lower_idx',
        'entity_identifier_lookup_type_value_lower_idx',
    ]
    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        for name in names:
            if progress:
                print(f'[schema] drop index {name}', flush=True)
            cur.execute(
                sql.SQL('DROP INDEX IF EXISTS {}.{}').format(
                    schema_id,
                    sql.Identifier(name),
                )
            )
    conn.commit()
    return names


def _rename_legacy_vocab_tables(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    renames = (
        ('identifier_type', 'vocab_identifier_type'),
        ('entity_type', 'vocab_entity_type'),
        ('entity_role', 'vocab_entity_role'),
        ('relation_predicate', 'vocab_relation_predicate'),
        ('relation_category', 'vocab_relation_category'),
        ('annotation_scope', 'vocab_annotation_scope'),
        ('resolution_status', 'vocab_resolution_status'),
        ('resolution_reason', 'vocab_resolution_reason'),
    )
    for old_name, new_name in renames:
        cur.execute(
            """
            SELECT
              to_regclass(%s) IS NOT NULL AS old_exists,
              to_regclass(%s) IS NOT NULL AS new_exists
            """,
            [f'{schema}.{old_name}', f'{schema}.{new_name}'],
        )
        old_exists, new_exists = cur.fetchone()
        if old_exists and not new_exists:
            cur.execute(
                sql.SQL('ALTER TABLE {}.{} RENAME TO {}').format(
                    sql.Identifier(schema),
                    sql.Identifier(old_name),
                    sql.Identifier(new_name),
                )
            )


def _ensure_primary_key(
    cur: psycopg2.extensions.cursor,
    schema: str,
    table: str,
    columns: tuple[str, ...],
) -> None:
    """Retrofit a primary key ``CREATE TABLE IF NOT EXISTS`` cannot add.

    ``IF NOT EXISTS`` is a no-op on an already-existing table, so a table
    that predates its declared ``PRIMARY KEY`` (an older build, or schema
    drift) keeps missing it forever -- and every later ``ADD CONSTRAINT
    ... FOREIGN KEY REFERENCES <table>(<columns>)`` then fails with "no
    unique constraint matching given keys". Found live on 2026-09-09
    preparing to add spec 011 WP2's new FKs, on a database where several
    core tables had no primary key despite their key columns already being
    unique.
    """

    column_list = sql.SQL(', ').join(sql.Identifier(c) for c in columns)
    cur.execute(
        sql.SQL(
            """
            DO $$
            BEGIN
              IF NOT EXISTS (
                SELECT 1 FROM pg_constraint
                WHERE conrelid = {}::regclass AND contype = 'p'
              ) THEN
                ALTER TABLE {} ADD PRIMARY KEY ({});
              END IF;
            END $$
            """
        ).format(
            sql.Literal(f'{schema}.{table}'),
            sql.Identifier(schema, table),
            column_list,
        )
    )


def _drop_legacy_content_tables_if_needed(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """Drop empty generated-id content tables before recreating UUID keys."""

    expected_uuid_columns = {
        'identifier_evidence': 'identifier_id',
        'entity_evidence': 'entity_evidence_id',
        'relation_evidence': 'relation_evidence_id',
    }
    cur.execute(
        """
        SELECT table_name, column_name, data_type
        FROM information_schema.columns
        WHERE table_schema = %s
          AND (table_name, column_name) IN (
            ('identifier_evidence', 'identifier_id'),
            ('entity_evidence', 'entity_evidence_id'),
            ('relation_evidence', 'relation_evidence_id')
          )
        """,
        [schema],
    )
    column_types = {
        (table_name, column_name): data_type
        for table_name, column_name, data_type in cur.fetchall()
    }
    legacy_tables = [
        table
        for table, column in expected_uuid_columns.items()
        if column_types.get((table, column)) not in (None, 'uuid')
    ]
    cur.execute(
        """
        SELECT table_name, column_name
        FROM information_schema.columns
        WHERE table_schema = %s
          AND table_name IN (
            'entity_evidence',
            'relation_evidence',
            'entity_evidence_identifier',
            'entity_evidence_annotation',
            'relation_evidence_annotation',
            'entity_evidence_resolution',
            'relation',
            'relation_evidence_relation'
          )
        """,
        [schema],
    )
    columns_by_table: dict[str, set[str]] = {}
    for table_name, column_name in cur.fetchall():
        columns_by_table.setdefault(table_name, set()).add(column_name)
    incompatible_columns = {
        'entity_evidence': {
            'occurrence_id',
            'snapshot_id',
            'entity_role',
            'entity_type',
        },
        'relation_evidence': {
            'relation_occurrence_id',
            'snapshot_id',
            'predicate',
            'relation_category',
        },
        'entity_evidence_identifier': set(),
        'entity_evidence_annotation': {'scope'},
        'relation_evidence_annotation': {'scope'},
        'entity_evidence_resolution': set(),
        'relation': {'predicate', 'relation_category'},
        'relation_evidence_relation': set(),
    }
    source_partitioned_tables = {
        'entity_evidence',
        'relation_evidence',
        'entity_evidence_identifier',
        'entity_evidence_annotation',
        'relation_evidence_annotation',
        'entity_evidence_resolution',
        'relation_evidence_relation',
    }
    for table, columns in incompatible_columns.items():
        existing_columns = columns_by_table.get(table)
        if not existing_columns:
            continue
        if columns & existing_columns:
            legacy_tables.append(table)
        elif (
            table in source_partitioned_tables
            and 'source_id' not in existing_columns
        ):
            legacy_tables.append(table)
    legacy_tables = sorted(set(legacy_tables))
    if not legacy_tables:
        return

    existing_content_tables = _existing_content_tables(cur, schema)
    non_empty_tables: list[str] = []
    for table in existing_content_tables:
        cur.execute(
            sql.SQL('SELECT EXISTS (SELECT 1 FROM {}.{} LIMIT 1)').format(
                sql.Identifier(schema),
                sql.Identifier(table),
            )
        )
        if bool(cur.fetchone()[0]):
            non_empty_tables.append(table)

    if non_empty_tables:
        raise RuntimeError(
            'Cannot migrate omnipath_build content tables to deterministic UUID keys '
            'while content rows exist. Run `make reset-content` first. '
            f'Non-empty tables: {", ".join(non_empty_tables)}'
        )

    for table in existing_content_tables:
        cur.execute(
            sql.SQL('DROP TABLE IF EXISTS {}.{} CASCADE').format(
                sql.Identifier(schema),
                sql.Identifier(table),
            )
        )


def _existing_content_tables(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> list[str]:
    cur.execute(
        """
        SELECT table_name
        FROM information_schema.tables
        WHERE table_schema = %s
          AND table_type = 'BASE TABLE'
          AND table_name = ANY(%s)
        """,
        [schema, list(CONTENT_TABLES)],
    )
    existing = {row[0] for row in cur.fetchall()}
    return [table for table in CONTENT_TABLES if table in existing]


def reset_content_tables(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> list[str]:
    """Truncate refresh-loaded content while keeping resolver tables intact."""

    with conn.cursor() as cur:
        if _table_exists(cur, schema, 'annotation'):
            _drop_obsolete_annotation_indexes(cur, schema)
        tables_to_truncate = _existing_content_tables(cur, schema)
        if tables_to_truncate:
            cur.execute(
                sql.SQL('TRUNCATE {} RESTART IDENTITY').format(
                    sql.SQL(', ').join(
                        sql.SQL('{}.{}').format(
                            sql.Identifier(schema),
                            sql.Identifier(table),
                        )
                        for table in tables_to_truncate
                    )
                )
            )
        _drop_existing_source_partitions(cur, schema)
    conn.commit()
    return tables_to_truncate


def _drop_existing_source_partitions(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> list[str]:
    """Drop source partitions so source IDs can be reused after reset."""

    dropped: list[str] = []
    for parent_table in SOURCE_PARTITION_DROP_ORDER:
        cur.execute(
            """
            SELECT child.relname
            FROM pg_inherits i
            JOIN pg_class parent ON parent.oid = i.inhparent
            JOIN pg_namespace parent_ns ON parent_ns.oid = parent.relnamespace
            JOIN pg_class child ON child.oid = i.inhrelid
            JOIN pg_namespace child_ns ON child_ns.oid = child.relnamespace
            WHERE parent_ns.nspname = %s
              AND parent.relname = %s
              AND child_ns.nspname = %s
              AND child.relname <> %s
            ORDER BY child.relname
            """,
            [schema, parent_table, schema, f'{parent_table}_default'],
        )
        child_tables = [row[0] for row in cur.fetchall()]
        for child_table in child_tables:
            cur.execute(
                sql.SQL('ALTER TABLE {}.{} DETACH PARTITION {}.{}').format(
                    sql.Identifier(schema),
                    sql.Identifier(parent_table),
                    sql.Identifier(schema),
                    sql.Identifier(child_table),
                )
            )
            cur.execute(
                sql.SQL('DROP TABLE {}.{}').format(
                    sql.Identifier(schema),
                    sql.Identifier(child_table),
                )
            )
            dropped.append(child_table)
    return dropped


def _table_exists(
    cur: psycopg2.extensions.cursor,
    schema: str,
    table: str,
) -> bool:
    cur.execute('SELECT to_regclass(%s) IS NOT NULL', [f'{schema}.{table}'])
    return bool(cur.fetchone()[0])


def _is_partitioned_table(
    cur: psycopg2.extensions.cursor,
    schema: str,
    table: str,
) -> bool:
    cur.execute(
        """
        SELECT c.relkind = 'p'
        FROM pg_class c
        JOIN pg_namespace n ON n.oid = c.relnamespace
        WHERE n.nspname = %s
          AND c.relname = %s
        """,
        [schema, table],
    )
    row = cur.fetchone()
    return bool(row and row[0])


def _ensure_gene_anchored_schema(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """Gene-anchored entity model.

    Additive + idempotent. The base ``entity`` gains the stored label +
    resolution-mechanism columns; ``entity_evidence_resolution`` gains the cheap
    ``molecular_type_id`` (set even when only a symbol is known). The gene form
    lives in two tiers: a precomputed ``gene_protein_representative`` (1:1, the
    no-state-join representative UniProt) and the heavy opt-in ``state`` /
    ``state_component`` / ``evidence_state``. The base graph keeps referencing
    gene-level entities.
    """
    schema_id = sql.Identifier(schema)

    # vocab: molecular type
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.vocab_molecular_type (
              molecular_type_id smallint PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            INSERT INTO {0}.vocab_molecular_type (molecular_type_id, name)
            VALUES (1,'gene'),(2,'protein'),(3,'mrna'),(4,'mirna'),
                   (5,'lncrna'),(6,'ncrna')
            ON CONFLICT (molecular_type_id) DO NOTHING
            """
        ).format(schema_id)
    )

    # entity: the label is precomputed here and stored, never derived per
    # request; the resolution mechanism records what resolved the entity and
    # to what.
    for ddl in (
        'ALTER TABLE {0}.entity ADD COLUMN IF NOT EXISTS label text',
        'ALTER TABLE {0}.entity ADD COLUMN IF NOT EXISTS label_rule text',
        'ALTER TABLE {0}.entity ADD COLUMN IF NOT EXISTS '
        'resolution_mechanism text',
        'ALTER TABLE {0}.entity ADD COLUMN IF NOT EXISTS '
        'resolution_detail jsonb',
    ):
        cur.execute(sql.SQL(ddl).format(schema_id))

    # vocab_identity_level (spec 011 data-model.md section 3): how specific a
    # resolved chemical identity is. Reuses the same three-tier InChIKey-block
    # scheme chemical_resolution_level.py already established for the derive-
    # phase group tables (connectivity / stereo_isotope_tautomer / full), so
    # "how specific is this entity's own identity" and "what level would this
    # entity's structure group at" read the same way. A non-chemical entity's
    # identity_level_id stays NULL -- the concept is chemical-specific.
    from omnipath_build.chemical_resolution_level import LEVELS as _IDENTITY_LEVELS

    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.vocab_identity_level (
              identity_level_id smallint PRIMARY KEY,
              name text NOT NULL UNIQUE,
              specificity_rank smallint NOT NULL,
              description text NOT NULL
            )
            """
        ).format(schema_id)
    )
    cur.executemany(
        sql.SQL(
            """
            INSERT INTO {0}.vocab_identity_level
              (identity_level_id, name, specificity_rank, description)
            VALUES (%s, %s, %s, %s)
            ON CONFLICT (identity_level_id) DO UPDATE SET
              name = EXCLUDED.name,
              specificity_rank = EXCLUDED.specificity_rank,
              description = EXCLUDED.description
            """
        ).format(schema_id),
        [
            (level.level_id, level.name, level.specificity_rank, level.description)
            for level in _IDENTITY_LEVELS
        ],
    )

    # entity, WP2 additions (spec 011 T054): identity_level_id names how
    # specific the chosen representative's identity is; identity_variants
    # records what skeleton collapse (T057) folded into it -- the discarded
    # structures and which layer they differed in, so a consumer can see a
    # charge/stereo/tautomer variant was merged rather than dropped.
    for ddl in (
        'ALTER TABLE {0}.entity ADD COLUMN IF NOT EXISTS '
        'identity_level_id smallint '
        'REFERENCES {0}.vocab_identity_level(identity_level_id)',
        'ALTER TABLE {0}.entity ADD COLUMN IF NOT EXISTS '
        'identity_variants jsonb',
    ):
        cur.execute(sql.SQL(ddl).format(schema_id))

    # cheap-tier molecular type on the resolved evidence slot
    cur.execute(
        sql.SQL(
            'ALTER TABLE {0}.entity_evidence_resolution '
            'ADD COLUMN IF NOT EXISTS molecular_type_id smallint'
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            'CREATE INDEX IF NOT EXISTS '
            'entity_evidence_resolution_molecular_type_idx '
            'ON {0}.entity_evidence_resolution (molecular_type_id)'
        ).format(schema_id)
    )

    # precomputed representative UniProt (cheap tier, no state join)
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.gene_protein_representative (
              entity_id uuid PRIMARY KEY
                REFERENCES {0}.entity(entity_id) ON DELETE CASCADE,
              representative_uniprot text,
              is_reviewed boolean,
              uniprot_all text[]
            )
            """
        ).format(schema_id)
    )

    # heavy opt-in: full molecular state
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.state (
              state_id uuid PRIMARY KEY,
              gene_entity_id uuid NOT NULL
                REFERENCES {0}.entity(entity_id) ON DELETE CASCADE,
              molecular_type_id smallint
                REFERENCES {0}.vocab_molecular_type(molecular_type_id)
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.state_component (
              state_id uuid NOT NULL
                REFERENCES {0}.state(state_id) ON DELETE CASCADE,
              component_type text NOT NULL,
              value text NOT NULL,
              PRIMARY KEY (state_id, component_type, value)
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.evidence_state (
              source_id bigint NOT NULL
                REFERENCES {0}.data_source(source_id),
              entity_evidence_id uuid NOT NULL,
              state_id uuid NOT NULL
                REFERENCES {0}.state(state_id) ON DELETE CASCADE,
              PRIMARY KEY (source_id, entity_evidence_id, state_id),
              FOREIGN KEY (source_id, entity_evidence_id)
                REFERENCES {0}.entity_evidence(source_id, entity_evidence_id)
                ON DELETE CASCADE
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.evidence_state_default
            PARTITION OF {0}.evidence_state DEFAULT
            """
        ).format(schema_id)
    )

    # gene-centric output: the primary deliverable layout — one row per gene
    # entity with the human-readable `gene` (symbol, NCBI Gene id fallback),
    # the pure `ncbi_gene_id`, and the representative UniProt on demand from
    # gene_protein_representative. Reads only gene-level columns — NO
    # state/evidence_state join. A view (not a materialised table) so it stays
    # continuous with `entity` at zero cost.
    cur.execute(
        sql.SQL(
            """
            CREATE OR REPLACE VIEW {0}.gene_output AS
            SELECT
              e.entity_id,
              e.taxonomy_id,
              -- A gene product known only by a UniProt accession ("unknown gene")
              -- is keyed by that accession, not an NCBI Gene id — so ncbi_gene_id
              -- is populated only for Entrez-canonical genes; UniProt-keyed genes
              -- surface their accession via `uniprot` below and leave ncbi_gene_id
              -- NULL. The accession literals mirror the hardcoded 'Gene:MI:0250'
              -- entity type below (CV terms: Entrez MI:0477, UniProt MI:1097).
              CASE
                WHEN e.canonical_identifier_type_id = (
                  SELECT identifier_type_id
                  FROM {0}.vocab_identifier_type
                  WHERE name = 'Entrez:MI:0477'
                ) THEN e.canonical_identifier
              END AS ncbi_gene_id,
              coalesce(e.label, e.canonical_identifier) AS gene,
              coalesce(
                gpr.representative_uniprot,
                CASE
                  WHEN e.canonical_identifier_type_id = (
                    SELECT identifier_type_id
                    FROM {0}.vocab_identifier_type
                    WHERE name = 'Uniprot:MI:1097'
                  ) THEN e.canonical_identifier
                END
              ) AS uniprot,
              gpr.is_reviewed AS uniprot_is_reviewed,
              gpr.uniprot_all
            FROM {0}.entity e
            LEFT JOIN {0}.gene_protein_representative gpr
              ON gpr.entity_id = e.entity_id
            WHERE e.entity_type_id = (
              SELECT entity_type_id
              FROM {0}.vocab_entity_type
              WHERE name = 'Gene:MI:0250'
            )
            """
        ).format(schema_id)
    )


def _ensure_resolution_schema(
    cur: psycopg2.extensions.cursor,
    schema: str,
    *,
    progress: bool = False,
    indexes: bool = True,
) -> None:
    def log_step(message: str) -> None:
        if progress:
            print(f'[schema] {message}', flush=True)

    schema_id = sql.Identifier(schema)
    log_step('drop obsolete canonicalization tables')
    _drop_obsolete_canonicalization_tables(cur, schema)
    log_step('create identifier type table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_identifier_type (
              identifier_type_id bigint PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    log_step('create entity type table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.data_source (
              source_id bigserial PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.dataset (
              dataset_id bigserial PRIMARY KEY,
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              name text NOT NULL,
              UNIQUE (source_id, name)
            )
            """
        ).format(schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_entity_type (
              entity_type_id bigserial PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    log_step('create resolution status table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_resolution_status (
              resolution_status_id smallint PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    _ensure_static_resolution_statuses(cur, schema)
    log_step('create resolution reason table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_resolution_reason (
              resolution_reason_id smallint PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    _ensure_static_resolution_reasons(cur, schema)
    log_step('create entity table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity (
              entity_id uuid PRIMARY KEY,
              entity_type_id bigint NOT NULL
                REFERENCES {}.vocab_entity_type(entity_type_id),
              taxonomy_id bigint,
              canonical_identifier_type_id bigint
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              canonical_identifier text NOT NULL,
              resolution_status_id smallint NOT NULL
                REFERENCES {}.vocab_resolution_status(resolution_status_id),
              created_at timestamptz NOT NULL DEFAULT now()
            )
            """
        ).format(schema_id, schema_id, schema_id, schema_id)
    )
    _ensure_primary_key(cur, schema, 'entity', ('entity_id',))
    log_step('ensure identifier types')
    _ensure_static_identifier_types(cur, schema)
    log_step('ensure entity indexes')
    _ensure_entity_canonical_key(cur, schema)
    log_step('ensure classification vocab')
    _ensure_classification_vocab(cur, schema)
    log_step('create entity identifier table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_identifier (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id)
                ON DELETE CASCADE,
              identifier_id uuid NOT NULL
                REFERENCES {}.identifier_evidence(identifier_id),
              PRIMARY KEY (source_id, entity_id, identifier_id)
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id, schema_id, schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_identifier_default
            PARTITION OF {}.entity_identifier DEFAULT
            """
        ).format(schema_id, schema_id)
    )
    log_step('ensure annotation value schema')
    _ensure_annotation_value_schema(cur, schema)
    log_step('create evidence annotation tables')
    _ensure_evidence_annotation_tables(cur, schema)
    log_step('ensure relation evidence endpoints')
    _ensure_relation_evidence_entity_endpoints(cur, schema)
    log_step('drop obsolete entity_resolution_candidate table')
    cur.execute(
        sql.SQL(
            'DROP TABLE IF EXISTS {}.entity_resolution_candidate CASCADE'
        ).format(schema_id)
    )
    log_step('create entity_evidence_resolution table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_evidence_resolution (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              entity_evidence_id uuid NOT NULL,
              status_id smallint NOT NULL
                REFERENCES {}.vocab_resolution_status(resolution_status_id),
              entity_id uuid
                REFERENCES {}.entity(entity_id),
              reason_id smallint
                REFERENCES {}.vocab_resolution_reason(resolution_reason_id),
              resolved_at timestamptz NOT NULL DEFAULT now(),
              CHECK (
                (
                  status_id IN (1, 2, 3)
                  AND entity_id IS NOT NULL
                )
                OR
                (
                  status_id = 4
                  AND entity_id IS NULL
                )
              ),
              PRIMARY KEY (source_id, entity_evidence_id),
              FOREIGN KEY (source_id, entity_evidence_id)
                REFERENCES {}.entity_evidence(source_id, entity_evidence_id)
                ON DELETE CASCADE
            ) PARTITION BY LIST (source_id)
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
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_evidence_resolution_default
            PARTITION OF {}.entity_evidence_resolution DEFAULT
            """
        ).format(schema_id, schema_id)
    )
    log_step('create resolution_conflict table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.resolution_conflict (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              entity_evidence_id uuid NOT NULL,
              identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              value_normalized text NOT NULL,
              candidate_structure text NOT NULL,
              candidate_source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              candidate_role text NOT NULL CHECK (
                candidate_role IN ('authoritative', 'cross_reference')
              ),
              skeleton text NOT NULL,
              recorded_at timestamptz NOT NULL DEFAULT now(),
              PRIMARY KEY
                (source_id, entity_evidence_id, candidate_structure),
              FOREIGN KEY (source_id, entity_evidence_id)
                REFERENCES {}.entity_evidence(source_id, entity_evidence_id)
                ON DELETE CASCADE
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id, schema_id, schema_id, schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.resolution_conflict_default
            PARTITION OF {}.resolution_conflict DEFAULT
            """
        ).format(schema_id, schema_id)
    )
    log_step('ensure entity resolution reason')
    _ensure_entity_resolution_reason(cur, schema)
    log_step('ensure entity resolution check')
    _ensure_entity_resolution_entity_check(cur, schema)
    log_step('create resolver_mapping_policy table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.resolver_mapping_policy (
              entity_family text NOT NULL,
              resolver_source text,
              key_type text NOT NULL,
              action text NOT NULL CHECK (
                action IN ('accept', 'candidate_only', 'ignore')
              ),
              requires_taxonomy boolean NOT NULL DEFAULT false
            )
            """
        ).format(schema_id)
    )
    log_step('create identifier authority table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.identifier_authority (
              identifier_type_id bigint PRIMARY KEY
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              is_structure_authority boolean NOT NULL DEFAULT false
            )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    log_step('create identifier role table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.identifier_role (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              role text NOT NULL CHECK (
                role IN ('authoritative', 'cross_reference')
              ),
              mention_count bigint NOT NULL,
              PRIMARY KEY (source_id, identifier_type_id)
            )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    log_step('create chemical resolution coverage table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.chemical_resolution_coverage (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              role text NOT NULL CHECK (
                role IN ('authoritative', 'cross_reference')
              ),
              mentions bigint NOT NULL,
              entities bigint NOT NULL,
              reached_structure bigint NOT NULL,
              reached_name bigint NOT NULL,
              unresolved bigint NOT NULL,
              conflicted bigint NOT NULL,
              PRIMARY KEY (source_id, identifier_type_id, role)
            )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    log_step('seed the lipid name identifier type')
    # Purely build-computed -- unlike every other canonical identifier type
    # (InChIKey, Name, ChEBI, ...), no resource ever asserts a raw "Lipid
    # Name:OM:0209" identifier, so it never appears in entity_identifier_raw
    # for vocab_identifier_type's usual auto-discovery (schema.py's own
    # `missing`/`missing_new` INSERT below) to pick up. Seeded explicitly
    # here instead, before any data loads, so the auto-discovery step's own
    # MAX(identifier_type_id)+1 never collides with it.
    #
    # The MAX(...) must be computed in its own FROM-less scalar subquery, not
    # a `... FROM vocab_identifier_type WHERE NOT EXISTS (...)` on the outer
    # query -- an aggregate without GROUP BY always returns exactly one row
    # even when its WHERE excludes every row (MAX over zero rows is NULL, not
    # "no row"), so that shape still attempts an insert (with the bogus id
    # COALESCE(NULL, 0) + 1 = 1) and crashes on this sandbox's accumulated
    # vocab_identifier_type (which, across many builds, already has this row
    # at a real id from before it was added here) instead of being skipped as
    # intended. Found live running this cycle's final rebuild.
    cur.execute(
        sql.SQL(
            """
            INSERT INTO {}.vocab_identifier_type (identifier_type_id, name)
            SELECT
              (SELECT COALESCE(MAX(identifier_type_id), 0) + 1 FROM {}.vocab_identifier_type),
              'Lipid Name:OM:0209'
            WHERE NOT EXISTS (
              SELECT 1 FROM {}.vocab_identifier_type WHERE name = 'Lipid Name:OM:0209'
            )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    log_step('create lipid identity tables')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.lipid_name_node (
              lipid_name text NOT NULL,
              lipid_level text NOT NULL,
              chains_possible smallint NOT NULL DEFAULT 0,
              chains_listed smallint NOT NULL DEFAULT 0,
              lipid_category text,
              lipid_class text,
              total_carbon smallint,
              total_db smallint,
              sum_formula text,
              parser_version text,
              PRIMARY KEY (lipid_name, lipid_level, chains_listed, chains_possible)
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.lipid_name_edge (
              child_name text NOT NULL,
              child_level text NOT NULL,
              child_chains_listed smallint NOT NULL,
              child_chains_possible smallint NOT NULL,
              parent_name text NOT NULL,
              parent_level text NOT NULL,
              parent_chains_listed smallint NOT NULL,
              parent_chains_possible smallint NOT NULL,
              relation text NOT NULL CHECK (relation = 'is_a'),
              derivation text NOT NULL,
              PRIMARY KEY (
                child_name, child_level, child_chains_listed,
                child_chains_possible, parent_name, parent_level,
                parent_chains_listed, parent_chains_possible
              ),
              FOREIGN KEY (
                child_name, child_level, child_chains_listed,
                child_chains_possible
              ) REFERENCES {}.lipid_name_node (
                lipid_name, lipid_level, chains_listed, chains_possible
              ),
              FOREIGN KEY (
                parent_name, parent_level, parent_chains_listed,
                parent_chains_possible
              ) REFERENCES {}.lipid_name_node (
                lipid_name, lipid_level, chains_listed, chains_possible
              ),
              CHECK (
                child_name <> parent_name OR child_level <> parent_level
                OR child_chains_listed <> parent_chains_listed
                OR child_chains_possible <> parent_chains_possible
              )
            )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    log_step('create build capability table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.build_capability (
              capability text PRIMARY KEY,
              available boolean NOT NULL,
              provider text NOT NULL,
              reason text
            )
            """
        ).format(schema_id)
    )
    log_step('ensure gene-anchored schema')
    _ensure_gene_anchored_schema(cur, schema)
    log_step('create identifier type table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_identifier_type (
              identifier_type_id bigint PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    log_step('create resolver protein lookup table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.resolver_protein_identifier_lookup (
              key_identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              key_value text NOT NULL,
              taxonomy_id text,
              canonical_identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              canonical_identifier text NOT NULL
            )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    log_step('create ambiguous resolver protein lookup table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.resolver_protein_identifier_lookup_ambiguous (
              key_identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              key_value text NOT NULL,
              taxonomy_id text,
              canonical_identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              canonical_identifier text NOT NULL
            )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    log_step('create resolver chemical lookup table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.resolver_chemical_identifier_lookup (
              key_identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              key_value text NOT NULL,
              canonical_identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              canonical_identifier text NOT NULL
            )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    log_step('create ambiguous resolver chemical lookup table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.resolver_chemical_identifier_lookup_ambiguous (
              key_identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              key_value text NOT NULL,
              canonical_identifier_type_id bigint NOT NULL
                REFERENCES {}.vocab_identifier_type(identifier_type_id),
              canonical_identifier text NOT NULL
            )
            """
        ).format(schema_id, schema_id, schema_id)
    )
    log_step('create resources table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.resources (
              resource_id text PRIMARY KEY,
              resource_name text,
              resource_short text,
              resource_full text,
              synonyms text[] NOT NULL DEFAULT ARRAY[]::text[],
              description text,
              homepage_url text,
              license text,
              license_label text,
              pubmed_id text,
              resource_kind text,
              input_module text,
              input_module_commit text,
              input_module_dirty boolean NOT NULL DEFAULT false,
              categories jsonb,
              annotation_ontologies jsonb,
              entity_count bigint NOT NULL DEFAULT 0,
              interaction_count bigint NOT NULL DEFAULT 0,
              association_count bigint NOT NULL DEFAULT 0,
              identifier_count bigint NOT NULL DEFAULT 0,
              ontology_term_count bigint NOT NULL DEFAULT 0,
              total_size_bytes bigint NOT NULL DEFAULT 0,
              last_downloaded_at timestamptz,
              last_built_at timestamptz,
              build_status text
            )
            """
        ).format(schema_id)
    )
    _ensure_resources_table_schema(cur, schema)
    log_step('create ontology terms table')
    _ensure_ontology_terms_table(cur, schema)
    log_step('create entity ontology relation table')
    _ensure_entity_ontology_relation_table(cur, schema)
    log_step('create entity name table')
    _ensure_entity_name_table(cur, schema)
    log_step('create resolver policy index')
    cur.execute(
        sql.SQL(
            """
            CREATE UNIQUE INDEX IF NOT EXISTS resolver_mapping_policy_unique_idx
            ON {}.resolver_mapping_policy (
              entity_family,
              key_type,
              COALESCE(resolver_source, '')
            )
            """
        ).format(schema_id)
    )
    log_step('create relation table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.relation (
              relation_id uuid PRIMARY KEY,
              subject_entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id),
              predicate_id bigint NOT NULL
                REFERENCES {}.vocab_relation_predicate(relation_predicate_id),
              object_entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id),
              relation_category_id bigint
                REFERENCES {}.vocab_relation_category(relation_category_id),
              created_at timestamptz NOT NULL DEFAULT now()
            )
            """
        ).format(schema_id, schema_id, schema_id, schema_id, schema_id)
    )
    _ensure_primary_key(cur, schema, 'relation', ('relation_id',))
    log_step('create relation unique index')
    cur.execute(
        sql.SQL(
            """
            CREATE UNIQUE INDEX IF NOT EXISTS relation_unique_idx
            ON {}.relation (
              subject_entity_id,
              predicate_id,
              object_entity_id
            )
            """
        ).format(schema_id)
    )
    log_step('create relation_evidence_relation table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.relation_evidence_relation (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              relation_id uuid NOT NULL
                REFERENCES {}.relation(relation_id)
                ON DELETE CASCADE,
              relation_evidence_id uuid NOT NULL,
              PRIMARY KEY (source_id, relation_evidence_id),
              FOREIGN KEY (source_id, relation_evidence_id)
                REFERENCES {}.relation_evidence(source_id, relation_evidence_id)
                ON DELETE CASCADE
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id, schema_id, schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.relation_evidence_relation_default
            PARTITION OF {}.relation_evidence_relation DEFAULT
            """
        ).format(schema_id, schema_id)
    )
    log_step('create entity_annotation_relation table')
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_annotation_relation (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              entity_evidence_id uuid NOT NULL,
              annotation_key uuid NOT NULL,
              relation_id uuid NOT NULL
                REFERENCES {}.relation(relation_id)
                ON DELETE CASCADE,
              PRIMARY KEY (
                source_id,
                entity_evidence_id,
                annotation_key,
                relation_id
              ),
              FOREIGN KEY (
                source_id,
                entity_evidence_id,
                annotation_key
              )
                REFERENCES {}.entity_evidence_annotation(
                  source_id,
                  entity_evidence_id,
                  annotation_key
                )
                ON DELETE CASCADE
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id, schema_id, schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_annotation_relation_default
            PARTITION OF {}.entity_annotation_relation DEFAULT
            """
        ).format(schema_id, schema_id)
    )
    log_step('create canonical annotation tables')
    _ensure_canonical_annotation_tables(cur, schema)
    log_step('create resource license terms')
    _ensure_data_source_license(cur, schema)
    log_step('create interaction header and participant tables')
    _ensure_interaction_schema(cur, schema)
    log_step('drop relation_evidence_resolution table')
    cur.execute(
        sql.SQL('DROP TABLE IF EXISTS {}.relation_evidence_resolution').format(
            schema_id
        )
    )
    log_step('drop relation_annotation_evidence table')
    cur.execute(
        sql.SQL('DROP TABLE IF EXISTS {}.relation_annotation_evidence').format(
            schema_id
        )
    )
    log_step('ensure chemical resolution-level schema')
    from omnipath_build.chemical_resolution_level import (
        ensure_chemical_resolution_schema,
    )

    ensure_chemical_resolution_schema(cur, schema)
    if indexes:
        log_step('ensure resolution indexes')
        _ensure_resolution_indexes(cur, schema)
    else:
        log_step('defer resolution indexes')


def _ensure_identifier_evidence_value_normalized(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """``identifier_evidence.value_normalized``: the same canonical form
    ``normalize_identifier()`` (omnipath-utils, spec 011 T024-T026) computes.
    It sits beside the untouched raw ``value`` (T028/T038). WP1's resolver
    lookup joins on this column, not on ``value``. A source's raw
    prefix/case/padding must never matter to the join.
    """

    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.identifier_evidence
            ADD COLUMN IF NOT EXISTS value_normalized text
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS identifier_evidence_value_normalized_idx
            ON {}.identifier_evidence (identifier_type_id, value_normalized)
            """
        ).format(schema_id)
    )


#: One value_pattern per chemical namespace -- mirrors _VALUE_PATTERNS in
#: omnipath-utils' omnipath_utils/mapping/_id_types.py. Keep the two in
#: sync. The same duplication already exists in resolver_chemical.sql,
#: because the two repos share no runtime dependency.
#:
#: The build stores a bare accession as ``entity.canonical_identifier``,
#: with no namespace prefix -- the namespace is already carried by
#: ``canonical_identifier_type_id``. This differs from omnipath-utils'
#: CURIE-style resolver key (``CHEBI:123``, T024-T026). That prefix is a
#: join-key convention for the translation database, not the build's own
#: identity-storage convention, so this table keeps the two apart.
_CHEMICAL_VALUE_PATTERNS: dict[str, str] = {
    'pubchem': r'[1-9]\d*|0',
    'chembl': r'CHEMBL\d+',
    'chebi': r'\d+',
    'hmdb': r'HMDB\d{7}',
    'lipidmaps': r'LM[A-Z0-9]+',
    'swisslipids': r'SLM:\d+',
    'kegg': r'C\d{5}',
    'inchikey': r'[A-Z]{14}-[A-Z]{10}-[A-Z]',
}


def _ensure_identifier_type_value_pattern(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """``vocab_identifier_type.value_pattern``/``normalization_rule`` (T039).

    SC-007 checks every canonical chemical identity against its own
    namespace's syntax with ``value_pattern`` (see
    contracts/coverage-acceptance.sql). ``normalization_rule`` names which
    of omnipath-utils' normalizers applies. That rule itself lives in one
    place (omnipath-utils), not duplicated as code here -- this column is
    a pointer for a human reading the vocabulary table, nothing more.
    """

    from omnipath_build.duckdb_load import (
        RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE,
        STANDARD_INCHI_KEY_TYPE,
    )

    slug_to_type_name = dict(RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE)
    slug_to_type_name['inchikey'] = STANDARD_INCHI_KEY_TYPE

    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.vocab_identifier_type
            ADD COLUMN IF NOT EXISTS value_pattern text
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.vocab_identifier_type
            ADD COLUMN IF NOT EXISTS normalization_rule text
            """
        ).format(schema_id)
    )
    for slug, type_name in slug_to_type_name.items():
        pattern = _CHEMICAL_VALUE_PATTERNS.get(slug)
        if pattern is None:
            continue
        cur.execute(
            sql.SQL(
                """
                UPDATE {}.vocab_identifier_type
                SET value_pattern = %s, normalization_rule = %s
                WHERE name = %s
                """
            ).format(schema_id),
            [pattern, slug, type_name],
        )


def _ensure_identifier_evidence_key(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.identifier_evidence
            ADD COLUMN IF NOT EXISTS identifier_type_id bigint
            """
        ).format(schema_id)
    )
    cur.execute(
        """
        SELECT EXISTS (
          SELECT 1
          FROM information_schema.columns
          WHERE table_schema = %s
            AND table_name = 'identifier_evidence'
            AND column_name = 'type'
        )
        """,
        [schema],
    )
    if bool(cur.fetchone()[0]):
        cur.execute(
            sql.SQL(
                """
                WITH missing AS (
                  SELECT DISTINCT i.type AS name
                  FROM {}.identifier_evidence i
                  LEFT JOIN {}.vocab_identifier_type it
                    ON it.name = i.type
                  WHERE i.type IS NOT NULL
                    AND it.identifier_type_id IS NULL
                ),
                base AS (
                  SELECT COALESCE(MAX(identifier_type_id), 0) AS max_id
                  FROM {}.vocab_identifier_type
                )
                INSERT INTO {}.vocab_identifier_type (identifier_type_id, name)
                SELECT
                  base.max_id + row_number() OVER (ORDER BY missing.name),
                  missing.name
                FROM missing
                CROSS JOIN base
                ON CONFLICT (name) DO NOTHING
                """
            ).format(schema_id, schema_id, schema_id, schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                UPDATE {}.identifier_evidence i
                SET identifier_type_id = it.identifier_type_id
                FROM {}.vocab_identifier_type it
                WHERE i.identifier_type_id IS NULL
                  AND it.name = i.type
                """
            ).format(schema_id, schema_id)
        )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.identifier_evidence
            ALTER COLUMN identifier_type_id SET NOT NULL
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            DO $$
            BEGIN
              IF NOT EXISTS (
                SELECT 1
                FROM pg_constraint
                WHERE conname = 'identifier_evidence_type_id_fkey'
                  AND conrelid = {table_literal}::regclass
              ) THEN
                ALTER TABLE {table_sql}
                ADD CONSTRAINT identifier_evidence_type_id_fkey
                FOREIGN KEY (identifier_type_id)
                REFERENCES {type_table}(identifier_type_id);
              END IF;
            END
            $$;
            """
        ).format(
            table_literal=sql.Literal(f'{schema}.identifier_evidence'),
            table_sql=sql.SQL('{}.identifier_evidence').format(schema_id),
            type_table=sql.SQL('{}.vocab_identifier_type').format(schema_id),
        )
    )
    cur.execute(
        """
        SELECT conname
        FROM pg_constraint c
        JOIN pg_class rel ON rel.oid = c.conrelid
        JOIN pg_namespace ns ON ns.oid = rel.relnamespace
        WHERE ns.nspname = %s
          AND rel.relname = 'identifier_evidence'
          AND c.contype = 'u'
          AND pg_get_constraintdef(c.oid) IN (
            'UNIQUE (type, value)',
            'UNIQUE (type, value_hash)'
          )
        """,
        [schema],
    )
    for (constraint_name,) in cur.fetchall():
        cur.execute(
            sql.SQL(
                'ALTER TABLE {}.identifier_evidence DROP CONSTRAINT {}'
            ).format(
                schema_id,
                sql.Identifier(constraint_name),
            )
        )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.identifier_type_value_idx').format(
            schema_id
        )
    )
    cur.execute(
        sql.SQL(
            'DROP INDEX IF EXISTS {}.identifier_type_value_hash_idx'
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            'DROP INDEX IF EXISTS {}.identifier_evidence_type_value_hash_idx'
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE UNIQUE INDEX IF NOT EXISTS identifier_evidence_type_value_idx
            ON {}.identifier_evidence (identifier_type_id, value)
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.identifier_evidence
            DROP COLUMN IF EXISTS value_hash,
            DROP COLUMN IF EXISTS type
            """
        ).format(schema_id)
    )


def _ensure_evidence_dimension_tables(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.data_source (
              source_id bigserial PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.dataset (
              dataset_id bigserial PRIMARY KEY,
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              name text NOT NULL,
              UNIQUE (source_id, name)
            )
            """
        ).format(schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_entity_type (
              entity_type_id bigserial PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_entity_role (
              entity_role_id smallserial PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_relation_predicate (
              relation_predicate_id bigserial PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_relation_category (
              relation_category_id bigserial PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_annotation_scope (
              annotation_scope_id smallserial PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    _ensure_static_entity_roles(cur, schema)
    _ensure_static_annotation_scopes(cur, schema)


def _ensure_resources_table_schema(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.resources
            ADD COLUMN IF NOT EXISTS input_module text,
            ADD COLUMN IF NOT EXISTS input_module_commit text,
            ADD COLUMN IF NOT EXISTS input_module_dirty boolean
              NOT NULL DEFAULT false
            """
        ).format(schema_id)
    )


def _ensure_static_entity_roles(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    rows = ((1, 'parent'), (2, 'member'))
    cur.executemany(
        sql.SQL(
            """
            INSERT INTO {}.vocab_entity_role (entity_role_id, name)
            VALUES (%s, %s)
            ON CONFLICT (entity_role_id) DO UPDATE
            SET name = EXCLUDED.name
            """
        )
        .format(sql.Identifier(schema))
        .as_string(cur.connection),
        rows,
    )


def _ensure_static_annotation_scopes(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    rows = ((1, 'relation'), (2, 'subject'), (3, 'object'))
    cur.executemany(
        sql.SQL(
            """
            INSERT INTO {}.vocab_annotation_scope (annotation_scope_id, name)
            VALUES (%s, %s)
            ON CONFLICT (annotation_scope_id) DO UPDATE
            SET name = EXCLUDED.name
            """
        )
        .format(sql.Identifier(schema))
        .as_string(cur.connection),
        rows,
    )


def _annotation_legacy_columns(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> set[str]:
    cur.execute(
        """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = %s
          AND table_name = 'annotation'
        """,
        [schema],
    )
    return {row[0] for row in cur.fetchall()}


def _ensure_annotation_value_schema(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    legacy_columns = _annotation_legacy_columns(cur, schema)
    if legacy_columns and (
        'annotation_key' not in legacy_columns
        or legacy_columns
        & {'scope', 'entity_evidence_id', 'relation_evidence_id', 'entity_id'}
    ):
        for table in (
            'relation_evidence_annotation',
            'entity_evidence_annotation',
        ):
            cur.execute(
                sql.SQL('DROP TABLE IF EXISTS {}.{} CASCADE').format(
                    schema_id,
                    sql.Identifier(table),
                )
            )
        cur.execute(
            sql.SQL('DROP TABLE {}.annotation CASCADE').format(schema_id)
        )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.annotation (
              annotation_key uuid PRIMARY KEY,
              term text NOT NULL,
              value text,
              unit text
            )
            """
        ).format(schema_id)
    )
    _drop_obsolete_annotation_indexes(cur, schema)
    cur.execute(
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS annotation_term_idx
            ON {}.annotation (term)
            """
        ).format(schema_id)
    )


def _ensure_evidence_annotation_tables(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        """
        SELECT 1
        FROM information_schema.columns
        WHERE table_schema = %s
          AND table_name = 'relation_evidence_annotation'
          AND column_name = 'relation_id'
        """,
        [schema],
    )
    if cur.fetchone() is not None:
        cur.execute(
            sql.SQL(
                'DROP TABLE IF EXISTS {}.relation_evidence_annotation CASCADE'
            ).format(schema_id)
        )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_evidence_annotation (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              entity_evidence_id uuid NOT NULL,
              annotation_key uuid NOT NULL
                REFERENCES {}.annotation(annotation_key)
                ON DELETE CASCADE,
              PRIMARY KEY (source_id, entity_evidence_id, annotation_key),
              FOREIGN KEY (source_id, entity_evidence_id)
                REFERENCES {}.entity_evidence(source_id, entity_evidence_id)
                ON DELETE CASCADE
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id, schema_id, schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_evidence_annotation_default
            PARTITION OF {}.entity_evidence_annotation DEFAULT
            """
        ).format(schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.relation_evidence_annotation (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              relation_evidence_id uuid NOT NULL,
              annotation_key uuid NOT NULL
                REFERENCES {}.annotation(annotation_key)
                ON DELETE CASCADE,
              annotation_scope_id smallint NOT NULL
                REFERENCES {}.vocab_annotation_scope(annotation_scope_id),
              PRIMARY KEY (
                source_id,
                relation_evidence_id,
                annotation_key,
                annotation_scope_id
              ),
              FOREIGN KEY (source_id, relation_evidence_id)
                REFERENCES {}.relation_evidence(source_id, relation_evidence_id)
                ON DELETE CASCADE
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id, schema_id, schema_id, schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.relation_evidence_annotation_default
            PARTITION OF {}.relation_evidence_annotation DEFAULT
            """
        ).format(schema_id, schema_id)
    )


def _ensure_canonical_annotation_tables(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    cur.execute(
        sql.SQL('DROP TABLE IF EXISTS {}.entity_annotation CASCADE').format(
            sql.Identifier(schema)
        )
    )
    cur.execute(
        sql.SQL('DROP TABLE IF EXISTS {}.relation_annotation CASCADE').format(
            sql.Identifier(schema)
        )
    )


def _ensure_ontology_terms_table(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        """
        SELECT 1
        FROM information_schema.columns
        WHERE table_schema = %s
          AND table_name = 'ontology_terms'
          AND column_name = 'source_id'
        """,
        [schema],
    )
    has_source_id = cur.fetchone() is not None
    if _table_exists(cur, schema, 'ontology_terms') and (
        not has_source_id
        or not _is_partitioned_table(cur, schema, 'ontology_terms')
    ):
        cur.execute(
            sql.SQL('DROP TABLE {}.ontology_terms CASCADE').format(schema_id)
        )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.ontology_terms (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              term_entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id)
                ON DELETE CASCADE,
              term_id text NOT NULL,
              ontology_prefix text,
              label text NOT NULL,
              definition text,
              ontology_id text,
              synonyms text[] NOT NULL DEFAULT '{{}}',
              synonyms_text text NOT NULL DEFAULT '',
              sources text[] NOT NULL DEFAULT '{{}}',
              PRIMARY KEY (source_id, term_entity_id)
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id, schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.ontology_terms_default
            PARTITION OF {}.ontology_terms DEFAULT
            """
        ).format(schema_id, schema_id)
    )


def _ensure_entity_ontology_relation_table(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    if _table_exists(cur, schema, 'entity_ontology_relation') and (
        not _is_partitioned_table(cur, schema, 'entity_ontology_relation')
    ):
        cur.execute(
            sql.SQL('DROP TABLE {}.entity_ontology_relation CASCADE').format(
                schema_id
            )
        )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_ontology_relation (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              subject_entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id)
                ON DELETE CASCADE,
              predicate_id bigint NOT NULL
                REFERENCES {}.vocab_relation_predicate(
                  relation_predicate_id
                ),
              object_entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id)
                ON DELETE CASCADE,
              ontology_id text NOT NULL,
              PRIMARY KEY (
                source_id,
                subject_entity_id,
                predicate_id,
                object_entity_id,
                ontology_id
              )
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id, schema_id, schema_id, schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_ontology_relation_default
            PARTITION OF {}.entity_ontology_relation DEFAULT
            """
        ).format(schema_id, schema_id)
    )


def _ensure_entity_name_table(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """spec 011 data-model.md section 4: one row per name attached to an
    entity, replacing the implicit name/synonym identifier-type split with
    an explicit preference flag. Populated during ``derive``
    (``entity_name.py``), not per-source ingestion -- entity ids only
    stabilize after resolution, and the same candidate query
    ``chemical_labels.py`` already runs supplies the raw material.
    """

    schema_id = sql.Identifier(schema)
    if _table_exists(cur, schema, 'entity_name') and (
        not _is_partitioned_table(cur, schema, 'entity_name')
    ):
        cur.execute(
            sql.SQL('DROP TABLE {}.entity_name CASCADE').format(schema_id)
        )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_name (
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id)
                ON DELETE CASCADE,
              name text NOT NULL,
              name_normalized text NOT NULL,
              is_preferred boolean NOT NULL,
              name_kind text NOT NULL CHECK (
                name_kind IN (
                  'trivial', 'systematic', 'abbreviation', 'brand'
                )
              ),
              PRIMARY KEY (source_id, entity_id, name_normalized)
            ) PARTITION BY LIST (source_id)
            """
        ).format(schema_id, schema_id, schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.entity_name_default
            PARTITION OF {}.entity_name DEFAULT
            """
        ).format(schema_id, schema_id)
    )


def _ensure_static_identifier_types(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """Seed the small set of "core" identifier types this codebase references
    by a stable id (``identifier_type_id(CHEBI_TYPE)`` etc.), from
    ``omnipath_build.resolver.identifier_types.IDENTIFIER_TYPE_NAMES``.

    ``vocab_identifier_type`` is NOT a content table (absent from
    CONTENT_TABLES; reset-content never truncates it) -- it accumulates
    across every build this database has ever run, and the *dynamic*
    discovery step (duckdb_load.py, ``INSERT ... base.max_id + row_number()``)
    assigns its own sequential ids, alphabetically, to every OTHER identifier
    type name actually cited in loaded data. On a database with build
    history, a name that is now in this static list may already hold a
    *different* id from an earlier run/session that discovered it
    dynamically before it was added here (found live: ``Lipid Name:OM:0209``,
    ``Drugbank:MI:2002`` and ``Reactome Stable Id:OM:0130`` all landed in this
    list this cycle, spec 011 T092/T116, but this sandbox's accumulated table
    already had them at other ids from earlier partial runs -- an
    unconditional upsert-by-id would either crash on the id/name pair
    (`ON CONFLICT (identifier_type_id) DO UPDATE` hitting the separate
    `name` UNIQUE constraint) or, worse, silently rename an unrelated,
    already-FK-referenced row that happens to already hold the target id.

    So: purely additive. Only inserts a (id, name) pair when NEITHER the id
    NOR the name exists yet -- never updates or renames an existing row,
    whichever mechanism (this seed, or dynamic discovery) put it there.
    """
    from omnipath_build.resolver.identifier_types import identifier_type_rows

    cur.executemany(
        sql.SQL(
            """
            INSERT INTO {}.vocab_identifier_type (identifier_type_id, name)
            SELECT %s, %s
            WHERE NOT EXISTS (
              SELECT 1 FROM {}.vocab_identifier_type existing
              WHERE existing.name = %s OR existing.identifier_type_id = %s
            )
            """
        )
        .format(sql.Identifier(schema), sql.Identifier(schema))
        .as_string(cur.connection),
        [
            (row['identifier_type_id'], row['name'], row['name'], row['identifier_type_id'])
            for row in identifier_type_rows()
        ],
    )


def _ensure_static_resolution_statuses(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    rows = (
        (1, 'resolved'),
        (2, 'unresolved'),
        (3, 'ambiguous'),
        (4, 'unsupported'),
    )
    cur.executemany(
        sql.SQL(
            """
            INSERT INTO {}.vocab_resolution_status (resolution_status_id, name)
            VALUES (%s, %s)
            ON CONFLICT (resolution_status_id) DO UPDATE
            SET name = EXCLUDED.name
            """
        )
        .format(sql.Identifier(schema))
        .as_string(cur.connection),
        rows,
    )


def _ensure_static_resolution_reasons(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    rows = (
        (1, 'missing_entity_type'),
        (2, 'different_taxon'),
        (3, 'no_accepted_resolver_candidate'),
        (4, 'multiple_entity_candidates'),
        (5, 'legacy_unresolved'),
        # The record could not be resolved because a translation table it needed
        # was absent or empty in the connected identifier-resolution database —
        # distinct from a record for which the resolver had data but no match.
        (6, 'missing_translation_table'),
    )
    cur.executemany(
        sql.SQL(
            """
            INSERT INTO {}.vocab_resolution_reason (resolution_reason_id, name)
            VALUES (%s, %s)
            ON CONFLICT (resolution_reason_id) DO UPDATE
            SET name = EXCLUDED.name
            """
        )
        .format(sql.Identifier(schema))
        .as_string(cur.connection),
        rows,
    )


def _ensure_entity_resolution_reason(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.entity_evidence_resolution
            ADD COLUMN IF NOT EXISTS reason_id smallint
            """
        ).format(schema_id)
    )
    cur.execute(
        """
        SELECT EXISTS (
          SELECT 1
          FROM information_schema.columns
          WHERE table_schema = %s
            AND table_name = 'entity_evidence_resolution'
            AND column_name = 'reason'
        )
        """,
        [schema],
    )
    if bool(cur.fetchone()[0]):
        cur.execute(
            sql.SQL(
                """
                UPDATE {}.entity_evidence_resolution er
                SET reason_id = rr.resolution_reason_id
                FROM {}.vocab_resolution_reason rr
                WHERE er.reason_id IS NULL
                  AND rr.name = er.reason
                """
            ).format(schema_id, schema_id)
        )
    cur.execute(
        sql.SQL(
            """
            DO $$
            BEGIN
              IF NOT EXISTS (
                SELECT 1
                FROM pg_constraint
                WHERE conname = 'entity_evidence_resolution_reason_id_fkey'
                  AND conrelid = {table_literal}::regclass
              ) THEN
                ALTER TABLE {table_sql}
                ADD CONSTRAINT entity_evidence_resolution_reason_id_fkey
                FOREIGN KEY (reason_id)
                REFERENCES {reason_table}(resolution_reason_id);
              END IF;
            END
            $$;
            """
        ).format(
            table_literal=sql.Literal(f'{schema}.entity_evidence_resolution'),
            table_sql=sql.SQL('{}.entity_evidence_resolution').format(
                schema_id
            ),
            reason_table=sql.SQL('{}.vocab_resolution_reason').format(
                schema_id
            ),
        )
    )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.entity_evidence_resolution
            DROP COLUMN IF EXISTS candidate_count,
            DROP COLUMN IF EXISTS reason
            """
        ).format(schema_id)
    )


def _ensure_classification_vocab(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """Controlled vocabularies for build-time entity/predicate classification.

    Rows are seeded by the classifiers in ``omnipath_build.classify`` from their
    curated rule files; here we only ensure the tables exist.
    """
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_chemical_class (
              chemical_class_id smallint PRIMARY KEY,
              name text NOT NULL UNIQUE,
              precedence smallint NOT NULL
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_metabolic_domain (
              metabolic_domain_id smallint PRIMARY KEY,
              name text NOT NULL UNIQUE,
              precedence smallint NOT NULL
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_interaction_class (
              interaction_class_id smallint PRIMARY KEY,
              name text NOT NULL UNIQUE,
              label text,
              controlled_vocabulary_mapping text
            )
            """
        ).format(schema_id)
    )
    # `name` is the snake_case slug used for storage and filtering, `label` the
    # capitalised display form, `controlled_vocabulary_mapping` the nullable
    # ontology CURIE. Added here too, for databases created before the columns
    # existed.
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.vocab_interaction_class
            ADD COLUMN IF NOT EXISTS label text,
            ADD COLUMN IF NOT EXISTS controlled_vocabulary_mapping text
            """
        ).format(schema_id)
    )
    # The predicate vocabulary is created elsewhere (populated dynamically from
    # loaded data); attach the classification FK column here.
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.vocab_relation_predicate
            ADD COLUMN IF NOT EXISTS interaction_class_id smallint
            """
        ).format(schema_id)
    )


def _ensure_data_source_license(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """Resource license terms in ``data_source_license``, as ordinal levels.

    The build recorded no license anywhere: ``data_source`` held
    ``(source_id, name)`` and nothing else. This table adds the terms, one row
    per resource, in the ordinal form ``pypath/internals/license.py`` already
    models.

    **Levels, not a name.** ``License.enables(other)`` is ``self >= other``, so
    a license question is a range predicate over three small integers, and
    "everything usable commercially" is ``WHERE purpose_level >= 15`` over a
    44-row table. A name answers no question on its own — two resources under
    different names can permit the same use, and one name under two versions
    can not — so ``license_name`` is stored for attribution and is never what a
    filter compares. The vocabularies, ascending:

    * ``purpose_level``: ignore 0, academic 5, nonprofit 10, commercial 15,
      free 20, composite 25
    * ``sharing_level``: ignore 0, noshare 5, noderiv 10, alike 15, share 20,
      free 25
    * ``attrib_level``: ignore 0, attrib 5, free 10

    **Unknown terms are exclusions, never defaults.** A resource whose
    license could not be mapped carries ``is_known = false`` and must not
    appear in a license-filtered result. The flag defaults to ``false`` for the
    same reason: absence of terms is not permission, and an unmapped license
    admitted as permissive is the one failure mode of this table that cannot be
    detected downstream.

    Every other column is nullable, because a resource with no license record
    has no name, no url and no levels to store. Nullable levels are deliberate
    rather than lax: a range predicate silently drops a NULL, which is the
    right answer for the wrong reason, so the exclusion is carried by
    ``is_known`` and never by the numbers.

    The table is declared here; the load step fills it from
    ``pypath/resources/data/resources.json``.
    """
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.data_source_license (
              source_id bigint PRIMARY KEY
                REFERENCES {}.data_source(source_id),
              license_name text,
              license_full_name text,
              license_url text,
              purpose_level smallint,
              sharing_level smallint,
              attrib_level smallint,
              is_known boolean NOT NULL DEFAULT false
            )
            """
        ).format(schema_id, schema_id)
    )
    # A license filter resolves to a resource set before it touches an
    # interaction table, and that resolution is a range scan over the levels.
    # The table holds 44 rows, so the index is for the shape of the predicate
    # rather than for the scan; the partial clause states the rule the query
    # must follow, which is that an unknown license is never a candidate.
    cur.execute(
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS data_source_license_levels_idx
            ON {}.data_source_license (
              purpose_level,
              sharing_level,
              attrib_level
            )
            WHERE is_known
            """
        ).format(schema_id)
    )


def _ensure_interaction_schema(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """Interaction header, participants, role vocabulary and the fact tables.

    ``interaction``, ``interaction_party``, ``interaction_fact_resource`` and
    ``vocab_relation_role``. The header carries an endpoint-independent
    identity for one interaction or reaction event; ``interaction_party``
    generalises the two endpoints of ``relation_evidence`` and the
    reactant/product grain of the DuckDB ``reaction_member_signature`` to N
    participants, each in a role, on a side, at an ordinal, optionally with a
    stoichiometry, a compartment and its own organism.

    The binary projection is **one table**.
    ``interaction_fact_resource`` is the **record**: one row per ordered
    endpoint pair, class and contributing **resource**, plus the assertion
    signature that resource states. No scope is precomputed, the all-resources
    scope included — the collapse is the shape a query produces by folding the
    record for its own scope at request time, and it is declared nowhere here
    because it is not a table, and no database carries the materialisation that
    used to hold it: both names existed only between 2026-08-18 and 2026-08-24,
    and the one build made in that window no longer holds either.

    The tables are declared here; the derive step fills them.
    """
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.vocab_relation_role (
              relation_role_id bigserial PRIMARY KEY,
              name text NOT NULL UNIQUE
            )
            """
        ).format(schema_id)
    )
    _ensure_static_relation_roles(cur, schema)
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.interaction (
              interaction_id uuid PRIMARY KEY,
              interaction_class_id smallint
                REFERENCES {}.vocab_interaction_class(interaction_class_id),
              arity smallint NOT NULL,
              sources text[],
              dataset_tags text[],
              attributes jsonb,
              created_at timestamptz NOT NULL DEFAULT now()
            )
            """
        ).format(schema_id, schema_id)
    )
    # `interaction_id` is a deterministic content hash over the sorted
    # participant multiset and the interaction class, so it is stable across
    # builds and independent of which endpoint a caller asks from.
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.interaction_party (
              interaction_id uuid NOT NULL
                REFERENCES {}.interaction(interaction_id)
                ON DELETE CASCADE,
              entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id),
              role_id bigint NOT NULL
                REFERENCES {}.vocab_relation_role(relation_role_id),
              side smallint,
              ordinal smallint,
              stoichiometry numeric,
              organism bigint,
              compartment text,
              role_flag smallint
            )
            """
        ).format(schema_id, schema_id, schema_id, schema_id)
    )
    # `organism`, `compartment` and `role_flag` are hot per-participant
    # columns: they answer set-predicates and the per-entity compartment
    # `interaction_fact_resource` is the interaction **record**: one row per
    # ordered `(subject, object, class)` and contributing `source_id`, plus the
    # assertion signature that resource states. Keeping the resource on the row
    # is what makes every summary decomposable — a query restricted to a subset
    # of resources recollapses these rows for its own scope instead of reading
    # numbers computed over all of them.
    #
    # A resource that asserts two contradicting signs for the same endpoints
    # under different predicates keeps two rows, because the signature is part
    # of the key.
    #
    # `interaction_fact_resource_id` is a **surrogate primary key, and it is
    # not optional**. The rest of the key includes the nullable signature
    # columns, and Postgres foreign keys default to `MATCH SIMPLE`, under which
    # **a foreign key with any NULL column is not checked at all**. Drug-target
    # rows assert neither sign nor direction, so a composite foreign key would
    # pass silently on nearly every row of the largest detail table while
    # looking enforced. The surrogate is a deterministic content hash over the
    # full key — the same `md5(to_json(...))::uuid` scheme the header id uses
    # (`derived_tables.interaction_content_uuid_sql`, DuckDB macro
    # `content_uuid`) — so it is stable across builds and a rebuilt record does
    # not orphan its detail rows.
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {}.interaction_fact_resource (
              interaction_fact_resource_id uuid PRIMARY KEY,
              subject_entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id),
              object_entity_id uuid NOT NULL
                REFERENCES {}.entity(entity_id),
              interaction_class_id smallint NOT NULL
                REFERENCES {}.vocab_interaction_class(interaction_class_id),
              source_id bigint NOT NULL
                REFERENCES {}.data_source(source_id),
              is_directed boolean,
              is_stimulation boolean,
              is_inhibition boolean,
              subject_organism bigint,
              object_organism bigint,
              affinity double precision,
              pchembl double precision,
              score double precision,
              curation_flags text[],
              reference_pubmed_ids text[],
              reference_dois text[],
              attributes jsonb,
              interaction_id uuid
                REFERENCES {}.interaction(interaction_id)
                ON DELETE CASCADE
            )
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
    # Sign and direction are three-valued here for the same reason they are on
    # the collapse, and the statement is narrower: NULL means **this** resource
    # is silent, which is a different claim from an asserted `false`. They
    # therefore carry no NOT NULL and no DEFAULT. The numeric and reference
    # columns hold what this one resource asserts, not a cross-resource best or
    # union — that is what a collapse computes.

    # The full key of the record, which is the endpoint pair, the class, the
    # resource **and** the signature that resource asserts. `NULLS NOT
    # DISTINCT` is what makes it a key at all: the three signature columns are
    # nullable, and under the default a resource silent on both sign and
    # direction would be free to repeat its row without the index noticing.
    # The same rule `entity_canonical_key_idx` follows.
    cur.execute(
        sql.SQL(
            """
            CREATE UNIQUE INDEX IF NOT EXISTS
              interaction_fact_resource_key_idx
            ON {}.interaction_fact_resource (
              subject_entity_id,
              object_entity_id,
              interaction_class_id,
              source_id,
              is_directed,
              is_stimulation,
              is_inhibition
            ) NULLS NOT DISTINCT
            """
        ).format(schema_id)
    )
    specs = [
        (
            'interaction_class_idx',
            'interaction',
            ('interaction_class_id',),
        ),
        (
            'interaction_party_interaction_idx',
            'interaction_party',
            ('interaction_id',),
        ),
        (
            'interaction_party_entity_idx',
            'interaction_party',
            ('entity_id',),
        ),
        (
            'interaction_party_role_idx',
            'interaction_party',
            ('role_id',),
        ),
        # A collapse groups the record by the endpoint/class triple, so it
        # gets a B-tree on exactly that prefix and runs over sorted input
        # rather than through a hash over 14.7 million rows. The unique key
        # starts with the same three columns but continues into the resource
        # and the signature, so it does not serve a grouped scan as cheaply.
        (
            'interaction_fact_resource_collapse_idx',
            'interaction_fact_resource',
            ('subject_entity_id', 'object_entity_id', 'interaction_class_id'),
        ),
        # Every scoped query starts by restricting the resource set, so the
        # resource is the other selective column on the record.
        (
            'interaction_fact_resource_source_idx',
            'interaction_fact_resource',
            ('source_id',),
        ),
        # A query that groups one interaction whatever its arity keys on the
        # header instead of on the endpoint pair, and no other index on this
        # table leads with `interaction_id` — the collapse index starts at the
        # endpoints and the unique key at the endpoints too, so both leave such
        # a grouping to a full scan. This is also the index that follows a
        # single header to the rows that assert it, which is how a reaction's
        # star is read back.
        (
            'interaction_fact_resource_interaction_idx',
            'interaction_fact_resource',
            ('interaction_id',),
        ),
    ]
    for index_name, table, columns in specs:
        cur.execute(
            sql.SQL('CREATE INDEX IF NOT EXISTS {} ON {}.{} ({})').format(
                sql.Identifier(index_name),
                schema_id,
                sql.Identifier(table),
                sql.SQL(', ').join(
                    sql.Identifier(column) for column in columns
                ),
            )
        )
    # A scoped query restricts the resource and then groups the endpoint/class
    # triple, and the planner costs those two steps as if the resource and the
    # class varied independently. They do not: a resource contributes to one or
    # two classes, so multiplying a resource's row count by a class's *global*
    # rarity understates a single-resource scope by orders of magnitude — a
    # ligand-receptor drop of 44,455 rows planned at 213. The cost governor
    # prices every request from that same number, so an understated scope is an
    # underpriced one, and the governor admits the work it exists to refuse.
    # `dependencies` is what corrects it; `ndistinct` sharpens the group count
    # the DISTINCT over the triple asks for. Statistics objects survive the
    # truncate-and-refill the load performs, and the ANALYZE that follows it
    # fills these along with the column statistics, so they cost no extra scan.
    # The name is schema-qualified because a statistics object lives in a
    # schema of its own: left bare it would land wherever the search path
    # points, and the second schema to be ensured would find the name taken
    # and quietly get no statistics at all.
    cur.execute(
        sql.SQL(
            'CREATE STATISTICS IF NOT EXISTS {}.{} (dependencies, ndistinct) '
            'ON source_id, interaction_class_id, subject_entity_id, '
            'object_entity_id FROM {}.interaction_fact_resource'
        ).format(
            schema_id,
            sql.Identifier('interaction_fact_resource_scope_stx'),
            schema_id,
        )
    )
    # The aggregated provenance arrays are filtered with array containment;
    # the `attributes` GIN on the record is benchmark-gated and is not created
    # until its size and build cost are measured. The partial indexes on the
    # three-valued sign and direction columns wait on the same measurement.
    for index_name, table, column in (
        ('interaction_sources_gin_idx', 'interaction', 'sources'),
        ('interaction_dataset_tags_gin_idx', 'interaction', 'dataset_tags'),
    ):
        cur.execute(
            sql.SQL(
                'CREATE INDEX IF NOT EXISTS {} ON {}.{} USING gin ({})'
            ).format(
                sql.Identifier(index_name),
                schema_id,
                sql.Identifier(table),
                sql.Identifier(column),
            )
        )


def _ensure_static_relation_roles(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    """Seed the participant roles, by name.

    The id is a serial the database assigns; every consumer joins
    ``vocab_relation_role`` on ``name`` and never hardcodes an id. Inserting
    only the missing names keeps ids stable across reruns.
    """
    names = [
        'subject',
        'object',
        'reactant',
        'product',
        'enzyme',
        'cofactor',
        'regulator',
        'member',
    ]
    cur.execute(
        sql.SQL(
            """
            INSERT INTO {}.vocab_relation_role (name)
            SELECT candidate.name
            FROM unnest(%s::text[]) WITH ORDINALITY
              AS candidate(name, position)
            WHERE NOT EXISTS (
              SELECT 1
              FROM {}.vocab_relation_role existing
              WHERE existing.name = candidate.name
            )
            ORDER BY candidate.position
            """
        ).format(sql.Identifier(schema), sql.Identifier(schema)),
        [names],
    )


def _ensure_entity_canonical_key(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.entity
            ADD COLUMN IF NOT EXISTS entity_type_id bigint,
            ADD COLUMN IF NOT EXISTS resolution_status_id smallint,
            ADD COLUMN IF NOT EXISTS canonical_identifier_type_id bigint,
            ADD COLUMN IF NOT EXISTS canonical_identifier text,
            ADD COLUMN IF NOT EXISTS chemical_class_id smallint,
            ADD COLUMN IF NOT EXISTS metabolic_domain_id smallint
            """
        ).format(schema_id)
    )
    cur.execute(
        """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = %s
          AND table_name = 'entity'
          AND column_name IN ('entity_type', 'resolution_status')
        """,
        [schema],
    )
    old_entity_columns = {row[0] for row in cur.fetchall()}
    if 'entity_type' in old_entity_columns:
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {}.vocab_entity_type (name)
                SELECT DISTINCT entity_type
                FROM {}.entity
                WHERE entity_type IS NOT NULL
                ON CONFLICT (name) DO NOTHING
                """
            ).format(schema_id, schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                UPDATE {}.entity e
                SET entity_type_id = et.entity_type_id
                FROM {}.vocab_entity_type et
                WHERE et.name = e.entity_type
                  AND e.entity_type_id IS NULL
                """
            ).format(schema_id, schema_id)
        )
    if 'resolution_status' in old_entity_columns:
        cur.execute(
            sql.SQL(
                """
                UPDATE {}.entity e
                SET resolution_status_id = rs.resolution_status_id
                FROM {}.vocab_resolution_status rs
                WHERE rs.name = e.resolution_status
                  AND e.resolution_status_id IS NULL
                """
            ).format(schema_id, schema_id)
        )
    cur.execute(
        """
        SELECT column_name
        FROM information_schema.columns
        WHERE table_schema = %s
          AND table_name = 'entity'
          AND column_name IN ('id', 'id_type')
        """,
        [schema],
    )
    old_columns = {row[0] for row in cur.fetchall()}
    if {'id', 'id_type'} <= old_columns:
        cur.execute(
            sql.SQL(
                """
                UPDATE {}.entity e
                SET
                  canonical_identifier_type_id = it.identifier_type_id,
                  canonical_identifier = e.id
                FROM {}.vocab_identifier_type it
                WHERE e.canonical_identifier IS NULL
                  AND e.id IS NOT NULL
                  AND it.name = e.id_type
                """
            ).format(schema_id, schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                UPDATE {}.entity e
                SET canonical_identifier = e.id
                WHERE e.canonical_identifier IS NULL
                  AND e.id IS NOT NULL
                """
            ).format(schema_id)
        )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.entity
            ALTER COLUMN entity_type_id SET NOT NULL,
            ALTER COLUMN resolution_status_id SET NOT NULL,
            ALTER COLUMN canonical_identifier SET NOT NULL,
            DROP COLUMN IF EXISTS identifiers
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            DO $$
            BEGIN
              IF NOT EXISTS (
                SELECT 1
                FROM pg_constraint
                WHERE conname = 'entity_canonical_identifier_type_id_fkey'
                  AND conrelid = {entity_table_literal}::regclass
              ) THEN
                ALTER TABLE {entity_table_sql}
                ADD CONSTRAINT entity_canonical_identifier_type_id_fkey
                FOREIGN KEY (canonical_identifier_type_id)
                REFERENCES {identifier_type_table}(identifier_type_id);
              END IF;
            END
            $$;
            """
        ).format(
            entity_table_literal=sql.Literal(f'{schema}.entity'),
            entity_table_sql=sql.SQL('{}.entity').format(schema_id),
            identifier_type_table=sql.SQL('{}.vocab_identifier_type').format(
                schema_id
            ),
        )
    )
    for column, target_table, target_column, constraint_name in (
        (
            'entity_type_id',
            'vocab_entity_type',
            'entity_type_id',
            'entity_entity_type_id_fkey',
        ),
        (
            'resolution_status_id',
            'vocab_resolution_status',
            'resolution_status_id',
            'entity_resolution_status_id_fkey',
        ),
    ):
        cur.execute(
            sql.SQL(
                """
                DO $$
                BEGIN
                  IF NOT EXISTS (
                    SELECT 1
                    FROM pg_constraint
                    WHERE conname = {constraint_name}
                      AND conrelid = {entity_table_literal}::regclass
                  ) THEN
                    ALTER TABLE {entity_table_sql}
                    ADD CONSTRAINT {constraint_identifier}
                    FOREIGN KEY ({column_identifier})
                    REFERENCES {target_table_sql}({target_column_identifier});
                  END IF;
                END
                $$;
                """
            ).format(
                constraint_name=sql.Literal(constraint_name),
                entity_table_literal=sql.Literal(f'{schema}.entity'),
                entity_table_sql=sql.SQL('{}.entity').format(schema_id),
                constraint_identifier=sql.Identifier(constraint_name),
                column_identifier=sql.Identifier(column),
                target_table_sql=sql.SQL('{}.{}').format(
                    schema_id,
                    sql.Identifier(target_table),
                ),
                target_column_identifier=sql.Identifier(target_column),
            )
        )
    cur.execute(
        """
        SELECT conname
        FROM pg_constraint c
        JOIN pg_class rel ON rel.oid = c.conrelid
        JOIN pg_namespace ns ON ns.oid = rel.relnamespace
        WHERE ns.nspname = %s
          AND rel.relname = 'entity'
          AND c.contype = 'u'
          AND pg_get_constraintdef(c.oid) = 'UNIQUE (entity_type, id_type, id)'
        """,
        [schema],
    )
    for (constraint_name,) in cur.fetchall():
        cur.execute(
            sql.SQL('ALTER TABLE {}.entity DROP CONSTRAINT {}').format(
                schema_id,
                sql.Identifier(constraint_name),
            )
        )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.entity_id_idx').format(schema_id)
    )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.entity_type_id_hash_idx').format(
            schema_id
        )
    )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.entity_id_hash_lookup_idx').format(
            schema_id
        )
    )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.entity_type_taxonomy_idx').format(
            schema_id
        )
    )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.entity_status_idx').format(schema_id)
    )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.entity_canonical_key_idx').format(
            schema_id
        )
    )
    cur.execute(
        sql.SQL(
            """
            CREATE UNIQUE INDEX IF NOT EXISTS entity_canonical_key_idx
            ON {}.entity (
              entity_type_id,
              taxonomy_id,
              canonical_identifier_type_id,
              canonical_identifier
            )
            NULLS NOT DISTINCT
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.entity
            DROP COLUMN IF EXISTS id_hash,
            DROP COLUMN IF EXISTS canonical_identifier_hash,
            DROP COLUMN IF EXISTS id_type,
            DROP COLUMN IF EXISTS id,
            DROP COLUMN IF EXISTS entity_type,
            DROP COLUMN IF EXISTS resolution_status
            """
        ).format(schema_id)
    )


def _drop_obsolete_annotation_indexes(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.annotation_term_value_idx').format(
            schema_id
        )
    )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.annotation_dedupe_idx').format(
            schema_id
        )
    )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.annotation_value_idx').format(
            schema_id
        )
    )
    cur.execute(
        sql.SQL(
            'ALTER TABLE {}.annotation DROP COLUMN IF EXISTS value_hash'
        ).format(schema_id)
    )


def _ensure_relation_evidence_entity_endpoints(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    for column in ('subject_entity_id', 'object_entity_id'):
        constraint_name = f'relation_evidence_{column}_fkey'
        cur.execute(
            sql.SQL(
                """
                DO $$
                BEGIN
                  IF NOT EXISTS (
                    SELECT 1
                    FROM pg_constraint
                    WHERE conname = {constraint_name}
                      AND conrelid = {relation_table_literal}::regclass
                  ) THEN
                    ALTER TABLE {relation_table_sql}
                    ADD CONSTRAINT {constraint_identifier}
                    FOREIGN KEY ({column_identifier})
                    REFERENCES {entity_table}(entity_id)
                    ON DELETE CASCADE;
                  END IF;
                END
                $$;
                """
            ).format(
                constraint_name=sql.Literal(constraint_name),
                relation_table_literal=sql.Literal(
                    f'{schema}.relation_evidence'
                ),
                relation_table_sql=sql.SQL('{}.relation_evidence').format(
                    schema_id
                ),
                constraint_identifier=sql.Identifier(constraint_name),
                column_identifier=sql.Identifier(column),
                entity_table=sql.SQL('{}.entity').format(schema_id),
            )
        )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.relation_evidence
            DROP CONSTRAINT IF EXISTS relation_evidence_subject_endpoint_check,
            DROP CONSTRAINT IF EXISTS relation_evidence_object_endpoint_check,
            ADD CONSTRAINT relation_evidence_subject_endpoint_check
              CHECK (
                (subject_entity_evidence_id IS NOT NULL)::int
                + (subject_entity_id IS NOT NULL)::int
                = 1
              ),
            ADD CONSTRAINT relation_evidence_object_endpoint_check
              CHECK (
                (object_entity_evidence_id IS NOT NULL)::int
                + (object_entity_id IS NOT NULL)::int
                = 1
              )
            """
        ).format(schema_id)
    )


def _ensure_entity_resolution_entity_check(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.entity_evidence_resolution
            ADD COLUMN IF NOT EXISTS status_id smallint
            """
        ).format(schema_id)
    )
    cur.execute(
        """
        SELECT EXISTS (
          SELECT 1
          FROM information_schema.columns
          WHERE table_schema = %s
            AND table_name = 'entity_evidence_resolution'
            AND column_name = 'status'
        )
        """,
        [schema],
    )
    if bool(cur.fetchone()[0]):
        cur.execute(
            sql.SQL(
                """
                UPDATE {}.entity_evidence_resolution er
                SET status_id = rs.resolution_status_id
                FROM {}.vocab_resolution_status rs
                WHERE er.status_id IS NULL
                  AND rs.name = er.status
                """
            ).format(schema_id, schema_id)
        )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.entity_evidence_resolution
            ALTER COLUMN status_id SET NOT NULL
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            DO $$
            BEGIN
              IF NOT EXISTS (
                SELECT 1
                FROM pg_constraint
                WHERE conname = 'entity_evidence_resolution_status_id_fkey'
                  AND conrelid = {table_literal}::regclass
              ) THEN
                ALTER TABLE {table_sql}
                ADD CONSTRAINT entity_evidence_resolution_status_id_fkey
                FOREIGN KEY (status_id)
                REFERENCES {status_table}(resolution_status_id);
              END IF;
            END
            $$;
            """
        ).format(
            table_literal=sql.Literal(f'{schema}.entity_evidence_resolution'),
            table_sql=sql.SQL('{}.entity_evidence_resolution').format(
                schema_id
            ),
            status_table=sql.SQL('{}.vocab_resolution_status').format(
                schema_id
            ),
        )
    )
    cur.execute(
        """
        SELECT conname
        FROM pg_constraint c
        JOIN pg_class rel ON rel.oid = c.conrelid
        JOIN pg_namespace ns ON ns.oid = rel.relnamespace
        WHERE ns.nspname = %s
          AND rel.relname = 'entity_evidence_resolution'
          AND c.contype = 'c'
          AND pg_get_constraintdef(c.oid) LIKE 'CHECK %%entity_id%%'
          AND (
            pg_get_constraintdef(c.oid) LIKE '%%status%%'
            OR pg_get_constraintdef(c.oid) LIKE '%%status_id%%'
          )
        """,
        [schema],
    )
    for (constraint_name,) in cur.fetchall():
        cur.execute(
            sql.SQL(
                'ALTER TABLE {}.entity_evidence_resolution DROP CONSTRAINT {}'
            ).format(schema_id, sql.Identifier(constraint_name))
        )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.entity_evidence_resolution
            ADD CONSTRAINT entity_evidence_resolution_entity_check
            CHECK (
              (
                status_id IN (1, 2, 3)
                AND entity_id IS NOT NULL
              )
              OR
              (
                status_id = 4
                AND entity_id IS NULL
              )
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            ALTER TABLE {}.entity_evidence_resolution
            DROP COLUMN IF EXISTS status
            """
        ).format(schema_id)
    )


def _drop_obsolete_canonicalization_tables(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    cur.execute(
        """
        SELECT EXISTS (
          SELECT 1
          FROM information_schema.tables
          WHERE table_schema = %s
            AND table_name = 'canonical_entity'
        )
        """,
        [schema],
    )
    has_old_canonical_tables = bool(cur.fetchone()[0])
    if not has_old_canonical_tables:
        return

    schema_id = sql.Identifier(schema)
    for table in (
        'canonical_relation_evidence',
        'canonical_relation',
        'entity_evidence_resolution',
        'canonical_entity',
    ):
        cur.execute(
            sql.SQL('DROP TABLE IF EXISTS {}.{} CASCADE').format(
                schema_id,
                sql.Identifier(table),
            )
        )


def _ensure_resolution_indexes(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            'DROP INDEX IF EXISTS {}.entity_evidence_identifier_entity_idx'
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL('DROP INDEX IF EXISTS {}.entity_canonical_lookup_idx').format(
            schema_id
        )
    )
    specs = [
        (
            'entity_evidence_identifier_identifier_idx',
            'entity_evidence_identifier',
            ('identifier_id', 'entity_evidence_id'),
        ),
        (
            'entity_identifier_entity_idx',
            'entity_identifier',
            ('entity_id', 'identifier_id'),
        ),
        (
            'entity_identifier_identifier_idx',
            'entity_identifier',
            ('identifier_id', 'entity_id'),
        ),
        (
            'entity_resolution_status_idx',
            'entity_evidence_resolution',
            ('status_id',),
        ),
        (
            'entity_resolution_reason_idx',
            'entity_evidence_resolution',
            ('reason_id',),
        ),
        (
            'entity_evidence_resolution_entity_idx',
            'entity_evidence_resolution',
            ('entity_id',),
        ),
        (
            'entity_canonical_lookup_idx',
            'entity',
            (
                'canonical_identifier_type_id',
                'canonical_identifier',
            ),
        ),
        ('entity_status_idx', 'entity', ('resolution_status_id',)),
        (
            'relation_subject_idx',
            'relation_evidence',
            ('subject_entity_evidence_id',),
        ),
        (
            'relation_object_idx',
            'relation_evidence',
            ('object_entity_evidence_id',),
        ),
        (
            'relation_subject_entity_idx',
            'relation_evidence',
            ('subject_entity_id',),
        ),
        (
            'relation_object_entity_idx',
            'relation_evidence',
            ('object_entity_id',),
        ),
        (
            'relation_source_dataset_idx',
            'relation_evidence',
            ('source_id', 'dataset_id'),
        ),
        (
            'relation_evidence_predicate_category_idx',
            'relation_evidence',
            ('predicate_id', 'relation_category_id'),
        ),
        (
            'relation_canonical_subject_idx',
            'relation',
            ('subject_entity_id',),
        ),
        (
            'relation_canonical_object_idx',
            'relation',
            ('object_entity_id',),
        ),
        (
            'entity_evidence_annotation_annotation_key_idx',
            'entity_evidence_annotation',
            ('annotation_key',),
        ),
        (
            'entity_annotation_relation_relation_id_idx',
            'entity_annotation_relation',
            ('relation_id',),
        ),
        (
            'relation_evidence_annotation_relation_evidence_idx',
            'relation_evidence_annotation',
            ('relation_evidence_id',),
        ),
        (
            'relation_evidence_annotation_annotation_key_idx',
            'relation_evidence_annotation',
            ('annotation_key',),
        ),
        (
            'entity_ontology_relation_object_idx',
            'entity_ontology_relation',
            ('object_entity_id',),
        ),
        (
            'entity_ontology_relation_subject_idx',
            'entity_ontology_relation',
            ('subject_entity_id',),
        ),
        (
            'resolver_protein_lookup_key_tax_idx',
            'resolver_protein_identifier_lookup',
            ('key_identifier_type_id', 'key_value', 'taxonomy_id'),
        ),
        (
            'resolver_protein_lookup_key_idx',
            'resolver_protein_identifier_lookup',
            ('key_identifier_type_id', 'key_value'),
        ),
        (
            'resolver_protein_ambiguous_key_tax_idx',
            'resolver_protein_identifier_lookup_ambiguous',
            ('key_identifier_type_id', 'key_value', 'taxonomy_id'),
        ),
        (
            'resolver_chemical_lookup_key_idx',
            'resolver_chemical_identifier_lookup',
            ('key_identifier_type_id', 'key_value'),
        ),
        (
            'resolver_chemical_ambiguous_key_idx',
            'resolver_chemical_identifier_lookup_ambiguous',
            ('key_identifier_type_id', 'key_value'),
        ),
        (
            'resolver_chemical_lookup_canonical_idx',
            'resolver_chemical_identifier_lookup',
            ('canonical_identifier_type_id', 'canonical_identifier'),
        ),
        ('resources_build_status_idx', 'resources', ('build_status',)),
    ]
    for name, table, columns in specs:
        cur.execute(
            sql.SQL('CREATE INDEX IF NOT EXISTS {} ON {}.{} ({})').format(
                sql.Identifier(name),
                schema_id,
                sql.Identifier(table),
                sql.SQL(', ').join(
                    sql.Identifier(column) for column in columns
                ),
            )
        )
    cur.execute(
        sql.SQL(
            """
            CREATE INDEX IF NOT EXISTS identifier_evidence_value_lower_idx
            ON {}.identifier_evidence (lower(value), identifier_type_id)
            INCLUDE (identifier_id)
            """
        ).format(schema_id)
    )
