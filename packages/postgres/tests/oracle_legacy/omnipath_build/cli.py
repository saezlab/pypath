"""Administrative CLI for the DuckDB/PostgreSQL build pipeline.

Source loading is handled by ``omnipath_build.duckdb_direct_pipeline``. This
module keeps the supporting database, resolver-materialization, source-deletion,
and derived-table maintenance commands.
"""

from __future__ import annotations

import os
import sys
import argparse
import logging
import time

import psycopg2

from omnipath_build import configure_build_tmpdir

from omnipath_build.db.derived_tables import InteractionDeriveStats
from omnipath_build.db.resources import INTERACTION_RECORD_STEP
from omnipath_build.db import (
    ensure_schema,
    reset_content_tables,
    sync_resources_table,
    emit_build_manifest,
    delete_source_content,
    rebuild_bitmap_tables,
    rebuild_derived_tables,
    rebuild_interaction_tables,
    rebuild_resource_overlap_summary,
    sweep_staging_tables,
    ensure_deferred_indexes,
    create_secondary_indexes,
    ensure_content_primary_keys,
    drop_deferred_content_indexes,
)
from omnipath_build.chemical_resolution_level import (
    rebuild_chemical_ambiguous_name_candidates,
    rebuild_chemical_resolution_levels,
)
from omnipath_build.classify import (
    classify_chemical_class,
    classify_metabolic_domain,
    classify_interaction_class,
)
from omnipath_build.labels import (
    populate_chemical_labels,
    populate_entity_labels,
    populate_entity_name,
)
from omnipath_build.cosmos import (
    CosmosBuildStats,
    build_cosmos_projection,
    ensure_cosmos_edge_table,
)
from omnipath_build.metsigdb import build_metsigdb
from omnipath_build.network_views import NETWORKS, apply_all as apply_network_views
from omnipath_build.resources import discover_resources, populate_identifier_authority
from omnipath_build.resolver.mapping_tables import (
    SOURCE_NAMES as RESOLVER_SOURCE_NAMES,
    run_sources as build_resolver_sources,
)


# The derive orchestration logs through the package logger, which
# `omnipath_build._session` configured on import. Bare `print` call sites
# elsewhere in the build keep what they have; this cycle does not sweep them.
#
# The name is spelled out rather than taken from `__name__`: the build runs
# this module as `python -m omnipath_build.cli`, where `__name__` is
# `__main__` — a logger outside the package tree, left at the root level
# (WARNING by default), which would silence the derive progress output in
# exactly the way the build is actually invoked.
_logger = logging.getLogger('omnipath_build.cli')


# Build-phase Postgres session tuning. The build connection runs a few heavy
# one-off statements (large sorts/aggregations for derived tables and network
# views, plus bulk index builds), so it is given generous memory + parallelism
# for its lifetime only — these are SESSION GUCs, never global config, so the
# modest global settings of a Postgres instance shared with the web API/app are
# untouched.
#
# Deployment tuning: each value is overridable via the matching environment
# variable below; set one to an empty string to leave that GUC at the server
# default. Defaults assume the build's Postgres has ~10 GB of headroom (e.g.
# the lab's docker.service hard cap is 200 GB shared across all containers, so
# a single build at ~8-10 GB peak is comfortable). Lower them on smaller hosts.
_BUILD_SESSION_TUNING: tuple[tuple[str, str, str], ...] = (
    ('work_mem', 'OMNIPATH_BUILD_WORK_MEM', '512MB'),
    ('maintenance_work_mem', 'OMNIPATH_BUILD_MAINTENANCE_WORK_MEM', '2GB'),
    (
        'max_parallel_workers_per_gather',
        'OMNIPATH_BUILD_MAX_PARALLEL_WORKERS_PER_GATHER',
        '6',
    ),
    (
        'max_parallel_maintenance_workers',
        'OMNIPATH_BUILD_MAX_PARALLEL_MAINTENANCE_WORKERS',
        '4',
    ),
)


def _apply_build_session_tuning(conn) -> None:
    """Apply session-level memory/parallelism GUCs to the build connection."""
    applied = []
    with conn.cursor() as cur:
        for guc, env_var, default in _BUILD_SESSION_TUNING:
            value = os.environ.get(env_var, default)
            if not value:
                continue
            # GUC names are a fixed allow-list above; the value is parameterised.
            cur.execute(f'SET {guc} = %s', (value,))
            applied.append(f'{guc}={value}')
    conn.commit()
    if applied:
        print(
            '[omnipath_build] build-session tuning: ' + ', '.join(applied),
            flush=True,
        )


def main(argv: list[str] | None = None) -> int:
    """Run the omnipath_build administrative command line interface."""

    parser = argparse.ArgumentParser(
        prog='omnipath_build',
        description='DuckDB/PostgreSQL pipeline administration.',
    )
    parser.add_argument(
        '--database-url',
        default=os.environ.get('DATABASE_URL'),
        help='PostgreSQL connection URL. Defaults to DATABASE_URL.',
    )
    parser.add_argument('--schema', default='public')
    parser.add_argument(
        '--tmpdir',
        default=None,
        help=(
            'Base directory for all builder temporary files (DuckDB spill, '
            'staging). Defaults to the OMNIPATH_BUILD_TMPDIR environment '
            'variable, else the system temp dir.'
        ),
    )

    subparsers = parser.add_subparsers(dest='command', required=True)
    init_db = subparsers.add_parser('init-db')
    init_db.add_argument('--drop-existing', action='store_true')
    init_db.add_argument(
        '--indexes',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='Create deferrable indexes during schema setup.',
    )

    reset_content = subparsers.add_parser('reset-content')
    reset_content.add_argument(
        '--drop-indexes',
        action='store_true',
        help='Drop secondary content indexes after truncating for faster reload.',
    )

    drop_source = subparsers.add_parser('drop-source')
    drop_source.add_argument('--source', required=True)
    drop_source.add_argument(
        '--row-delete',
        action='store_true',
        help='Force row deletion instead of dropping source partitions.',
    )

    build_resolver = subparsers.add_parser('build-resolver')
    build_resolver.add_argument(
        'sources',
        nargs='*',
        choices=RESOLVER_SOURCE_NAMES,
        help='Resolver sources to materialize. Defaults to all sources.',
    )
    build_resolver.add_argument('--output-dir', default='data')
    build_resolver.add_argument(
        '--taxonomy-id',
        dest='taxonomy_ids',
        action='append',
        default=None,
        help='Optional UniProt taxonomy filter. Can be repeated.',
    )
    build_resolver.add_argument(
        '--max-records',
        type=int,
        default=None,
        help='Optional cap for chemical source rows in smoke tests.',
    )
    build_resolver.add_argument(
        '--pubchem-url',
        default=None,
        help=(
            'Optional single PubChem SDF .gz URL/path. '
            'Defaults to all current PubChem full-SDF shards.'
        ),
    )
    build_resolver.add_argument(
        '--pubchem-shards',
        type=int,
        default=None,
        help='Optional number of discovered PubChem SDF shards to stream.',
    )
    build_resolver.add_argument(
        '--jobs',
        type=int,
        default=1,
        help='Number of parallel jobs for resolver sources that support it.',
    )
    build_resolver.add_argument(
        '--skip-existing',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='Skip resolver sources already present in the output directory.',
    )

    derive = subparsers.add_parser('derive')
    derive.add_argument(
        '--indexes',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='Create secondary indexes.',
    )
    derive.add_argument(
        '--tables',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='Create and populate derived search/count tables.',
    )
    derive.add_argument(
        '--interactions',
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            'Project the canonical graph into interaction/interaction_party/'
            'interaction_fact_resource. Runs after the interaction-class '
            'classification, which fills the predicate map it reads. Nothing '
            'is folded: the collapse of a resource scope is what a query '
            'produces.'
        ),
    )
    derive.add_argument(
        '--bitmaps',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='Create and populate bitmap tables.',
    )
    derive.add_argument(
        '--network-views',
        action=argparse.BooleanOptionalAction,
        default=True,
        help='Apply + register the specialized network views (MetalinksDB, LIANA).',
    )
    derive.add_argument(
        '--metsigdb',
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            'Publish the MetSigDB membership substrate (KEGG, Reactome, '
            'WikiPathways, MACdb, ClassyFire). Runs after the build manifest, '
            'whose stamp every published row carries.'
        ),
    )
    derive.add_argument(
        '--cosmos',
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            'Project the metabolic reactions into the binary COSMOS network. '
            'Runs after the build manifest, and reads the reaction headers '
            'the interactions phase wrote.'
        ),
    )
    derive.add_argument('--inputs-package', default='pypath.inputs_v2')
    derive.add_argument('--database', default='omnipath')
    derive.add_argument(
        '--max-records',
        default=os.environ.get('MAX_RECORDS'),
        help=(
            'Per-source record cap the load ran with; a non-zero value flags '
            'build_manifest.partial_build. Defaults to the MAX_RECORDS env var.'
        ),
    )

    network_views = subparsers.add_parser('network-views')
    network_views.add_argument('--schema', default='public')
    network_views.add_argument(
        '--refresh',
        action='store_true',
        help='Refresh existing matviews instead of a full apply (faster).',
    )

    metsigdb = subparsers.add_parser('metsigdb')
    metsigdb.add_argument('--schema', default='public')
    metsigdb.add_argument(
        '--max-records',
        default=None,
        help=(
            'Cap each resource for the dev loop. Unlike the derive, this does '
            'not default to the MAX_RECORDS env var: there the cap describes a '
            'load that already happened, here it decides what gets written and '
            'permanently flags build_manifest.partial_build. A destructive '
            'default has to be typed, not inherited from the shell.'
        ),
    )

    args = parser.parse_args(argv)
    configure_build_tmpdir(args.tmpdir)
    if args.command == 'build-resolver':
        summary = build_resolver_sources(
            sources=args.sources or RESOLVER_SOURCE_NAMES,
            output_dir=args.output_dir,
            taxonomy_ids=args.taxonomy_ids,
            max_records=args.max_records,
            pubchem_url=args.pubchem_url,
            pubchem_shards=args.pubchem_shards,
            jobs=args.jobs,
            skip_existing=args.skip_existing,
        )
        for key, value in summary.items():
            print(f'{key}: {value}', flush=True)
        return 0

    if not args.database_url:
        print('DATABASE_URL is required.', file=sys.stderr)
        return 2

    print(
        f'[omnipath_build] connecting database schema={args.schema}', flush=True
    )
    with psycopg2.connect(args.database_url) as conn:
        print('[omnipath_build] database connected', flush=True)
        _apply_build_session_tuning(conn)
        if args.command == 'network-views':
            from omnipath_build.network_views import refresh_all
            runner = refresh_all if args.refresh else apply_network_views
            stats = runner(
                conn,
                NETWORKS,
                registry_schema=args.schema,
                log=lambda message: print(message, flush=True),
            )
            print(
                f'[network-views] {"refreshed" if args.refresh else "applied"}: '
                f'{", ".join(stats.applied)}',
                flush=True,
            )
            return 0
        if args.command == 'metsigdb':
            cap = _record_cap(args.max_records)
            stats = build_metsigdb(
                conn, database_url=args.database_url,
                max_records=cap, schema=args.schema,
            )
            print(
                '[metsigdb] '
                f'build={stats.build_id} '
                f'rows={stats.rows} '
                f'partial={stats.partial} '
                f'seconds={stats.seconds} '
                + ' '.join(
                    f'{loaded.resource}={loaded.rows}'
                    for loaded in stats.resources
                ),
                flush=True,
            )
            return 0
        if args.command == 'init-db':
            ensure_schema(
                conn,
                schema=args.schema,
                drop_existing=args.drop_existing,
                indexes=args.indexes,
            )
            return 0
        if args.command == 'reset-content':
            tables = reset_content_tables(conn, schema=args.schema)
            print(
                f'[reset-content] schema={args.schema} tables={len(tables)}',
                flush=True,
            )
            if tables:
                print(
                    '[reset-content] truncated=' + ','.join(tables),
                    flush=True,
                )
            if args.drop_indexes:
                dropped = drop_deferred_content_indexes(
                    conn,
                    schema=args.schema,
                    progress=True,
                )
                print(
                    f'[reset-content] dropped_indexes={len(dropped)}',
                    flush=True,
                )
            return 0
        if args.command == 'drop-source':
            stats = delete_source_content(
                conn,
                schema=args.schema,
                source=args.source,
                drop_partitions=not args.row_delete,
            )
            print(
                '[drop-source] '
                f'source={stats.source} '
                f'source_id={stats.source_id} '
                f'strategy={stats.strategy} '
                f'partitions_dropped={stats.partitions_dropped} '
                f'affected_relations={stats.affected_relations} '
                f'affected_entities={stats.affected_entities} '
                f'affected_identifiers={stats.affected_identifiers} '
                f'deleted_relations={stats.deleted_relations} '
                f'deleted_entities={stats.deleted_entities} '
                f'deleted_identifiers={stats.deleted_identifiers} '
                f'deleted_annotations={stats.deleted_annotations} '
                f'refreshed_relation_counts={stats.refreshed_relation_counts}',
                flush=True,
            )
            return 0
        if args.command == 'derive':
            derive_started = time.perf_counter()
            # Supplementary steps are caught rather than raised. Their names
            # land here so the derive can say, in its exit code, that it
            # published less than it was asked to.
            failed_steps: list[str] = []
            # `MAX_RECORDS` reaches the derive as a cap the load already
            # applied: the derived tables are projections of what was loaded,
            # so a capped load caps them by construction. What the derive owes
            # it is visibility — the progress output says the run was capped,
            # and the manifest flags the build partial rather than letting a
            # capped build read as authoritative.
            _derive_log(
                'start',
                schema=args.schema,
                max_records=args.max_records or '',
                partial_build=_is_partial_build(args.max_records),
            )
            step_started = time.perf_counter()
            _derive_log('primary_keys_start')
            ensure_content_primary_keys(conn, schema=args.schema, progress=True)
            _derive_log(
                'primary_keys_done',
                seconds=f'{time.perf_counter() - step_started:.3f}',
            )
            step_started = time.perf_counter()
            _derive_log('schema_start')
            ensure_schema(conn, schema=args.schema, indexes=False)
            _derive_log(
                'schema_done',
                seconds=f'{time.perf_counter() - step_started:.3f}',
            )
            if args.indexes:
                step_started = time.perf_counter()
                _derive_log('deferred_indexes_start')
                ensure_deferred_indexes(
                    conn,
                    schema=args.schema,
                    progress=True,
                )
                create_secondary_indexes(conn, schema=args.schema)
                _derive_log(
                    'deferred_indexes_done',
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
            if args.tables:
                step_started = time.perf_counter()
                _derive_log('tables_start')
                table_stats = rebuild_derived_tables(
                    conn,
                    schema=args.schema,
                    progress=True,
                    # The interaction projection is a step of its own below,
                    # after `classify_interaction_class` fills the
                    # predicate->class map it reads. Left at the default it
                    # would run here as well — a second pass over 14M relations
                    # for every build, out of order and for nothing.
                    interactions=False,
                )
                _derive_log(
                    'tables_done',
                    entity_identifier_lookup=(
                        table_stats.entity_identifier_lookup
                    ),
                    entity_relation_counts=(
                        table_stats.entity_relation_counts
                    ),
                    entity_ontology_terms=table_stats.entity_ontology_terms,
                    ontology_terms=table_stats.ontology_terms,
                    entity_source_count=table_stats.entity_source_count,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('chemical_resolution_levels_start')
                resolution_level_stats = rebuild_chemical_resolution_levels(
                    conn,
                    schema=args.schema,
                    progress=True,
                )
                _derive_log(
                    'chemical_resolution_levels_done',
                    chemical_entities=(
                        resolution_level_stats.chemical_entities_with_inchikey
                    ),
                    groups=resolution_level_stats.groups,
                    members=resolution_level_stats.group_members,
                    relations=resolution_level_stats.relations,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('chemical_ambiguous_name_candidates_start')
                ambiguous_name_stats = (
                    rebuild_chemical_ambiguous_name_candidates(
                        conn,
                        schema=args.schema,
                        progress=True,
                    )
                )
                _derive_log(
                    'chemical_ambiguous_name_candidates_done',
                    ambiguous_names=ambiguous_name_stats.ambiguous_names,
                    candidate_links=ambiguous_name_stats.candidate_links,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('classify_chemical_class_start')
                chemical_class_stats = classify_chemical_class(
                    conn,
                    schema=args.schema,
                )
                _derive_log(
                    'classify_chemical_class_done',
                    classified=chemical_class_stats.classified,
                    by_default=chemical_class_stats.by_default,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('classify_metabolic_domain_start')
                metabolic_domain_stats = classify_metabolic_domain(
                    conn,
                    schema=args.schema,
                )
                _derive_log(
                    'classify_metabolic_domain_done',
                    classified=metabolic_domain_stats.classified,
                    by_default=metabolic_domain_stats.by_default,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('classify_interaction_class_start')
                interaction_class_stats = classify_interaction_class(
                    conn,
                    schema=args.schema,
                )
                _derive_log(
                    'classify_interaction_class_done',
                    mapped=interaction_class_stats.mapped,
                    by_default=interaction_class_stats.by_default,
                    default_predicates=','.join(
                        interaction_class_stats.default_predicates
                    ),
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                interaction_stats = (
                    _derive_interactions(conn, schema=args.schema)
                    if args.interactions
                    else None
                )
                step_started = time.perf_counter()
                _derive_log('entity_labels_start')
                label_stats = populate_entity_labels(conn, schema=args.schema)
                _derive_log(
                    'entity_labels_done',
                    gene_symbol=label_stats.gene_symbol,
                    identifier_fallback=label_stats.identifier_fallback,
                    without_label=label_stats.without_label,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('entity_name_start')
                entity_name_stats = populate_entity_name(
                    conn, schema=args.schema,
                )
                _derive_log(
                    'entity_name_done',
                    rows_written=entity_name_stats.rows_written,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('chemical_labels_start')
                chem_label_stats = populate_chemical_labels(
                    conn,
                    schema=args.schema,
                )
                _derive_log(
                    'chemical_labels_done',
                    chemical_name=chem_label_stats.chemical_name,
                    chemical_iupac_name=chem_label_stats.chemical_iupac_name,
                    chemical_preferred_name=(
                        chem_label_stats.chemical_preferred_name
                    ),
                    chemical_identifier=chem_label_stats.chemical_identifier,
                    without_real_label=(
                        chem_label_stats.chemical_without_real_label
                    ),
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('discover_resources_start')
                discovered, configs = discover_resources(
                    database_name=args.database,
                    inputs_package=args.inputs_package,
                    progress=True,
                )
                _derive_log(
                    'discover_resources_done',
                    sources=len(discovered),
                    functions=sum(
                        len(functions) for functions in discovered.values()
                    ),
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('identifier_authority_start')
                with conn.cursor() as cur:
                    authority_rows = populate_identifier_authority(
                        cur, args.schema, configs, progress=True,
                    )
                conn.commit()
                _derive_log(
                    'identifier_authority_done',
                    rows=authority_rows,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                bitmap_stats = None
                if args.bitmaps:
                    step_started = time.perf_counter()
                    _derive_log('bitmaps_start')
                    bitmap_stats = rebuild_bitmap_tables(
                        conn,
                        schema=args.schema,
                        progress=True,
                    )
                    _derive_log(
                        'bitmaps_done',
                        annotation_term_entities=(
                            bitmap_stats.annotation_term_entities
                        ),
                        annotation_term_direct_relations=(
                            bitmap_stats.annotation_term_direct_relations
                        ),
                        entity_relations=bitmap_stats.entity_relations,
                        entity_facets=bitmap_stats.entity_facets,
                        relation_facets=bitmap_stats.relation_facets,
                        seconds=f'{time.perf_counter() - step_started:.3f}',
                    )
                    step_started = time.perf_counter()
                    _derive_log('resource_overlap_start')
                    overlap_pairs = rebuild_resource_overlap_summary(
                        conn,
                        schema=args.schema,
                        progress=True,
                    )
                    _derive_log(
                        'resource_overlap_done',
                        pairs=overlap_pairs,
                        seconds=f'{time.perf_counter() - step_started:.3f}',
                    )
                step_started = time.perf_counter()
                _derive_log('resources_start')
                resource_stats = sync_resources_table(
                    conn,
                    discovered,
                    schema=args.schema,
                    prefer_bitmaps=args.bitmaps,
                )
                _derive_log(
                    'resources_done',
                    resources=resource_stats.resources,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                step_started = time.perf_counter()
                _derive_log('build_manifest_start')
                manifest_stats = emit_build_manifest(
                    conn,
                    schema=args.schema,
                    inputs_package=args.inputs_package,
                    partial_build=_is_partial_build(args.max_records),
                    derive_cost=_interaction_derive_cost(interaction_stats),
                    deferral_cost=_interaction_deferral_cost(interaction_stats),
                    utils_db_url=os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL'),
                )
                _derive_log(
                    'build_manifest_done',
                    build_id=manifest_stats.build_id,
                    partial_build=manifest_stats.partial_build,
                    resources=manifest_stats.resources,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
                if args.network_views:
                    step_started = time.perf_counter()
                    _derive_log('network_views_start')
                    try:
                        network_stats = apply_network_views(
                            conn,
                            NETWORKS,
                            registry_schema=args.schema,
                            log=lambda message: print(message, flush=True),
                        )
                        _derive_log(
                            'network_views_done',
                            networks=','.join(network_stats.applied),
                            seconds=f'{time.perf_counter() - step_started:.3f}',
                        )
                    except Exception as exc:
                        # Network views are a supplementary layer; a failure here
                        # (e.g. a missing source) must not abort the core build.
                        conn.rollback()
                        _supplementary_step_failed(
                            failed_steps,
                            'network_views',
                            exc,
                            time.perf_counter() - step_started,
                        )
                if args.metsigdb:
                    step_started = time.perf_counter()
                    _derive_log('metsigdb_start')
                    metsigdb_conn = None
                    try:
                        # The step needs a session of its own. Its extractions
                        # are temp tables sized by `temp_buffers`, which
                        # Postgres refuses to raise once a session has touched
                        # one — and by now this session has touched many. On
                        # the derive's connection the widening is refused, the
                        # step keeps the 8 MB default and dies on KEGG with
                        # "no empty local buffer available". So commit what the
                        # derive owes, then hand the step a fresh connection.
                        conn.commit()
                        metsigdb_conn = psycopg2.connect(args.database_url)
                        # No record cap here. The substrate is a projection of
                        # what the load wrote, so a capped load caps it by
                        # construction; re-applying the cap would truncate a
                        # full build whose derive was handed the load's figure.
                        metsigdb_stats = build_metsigdb(
                            metsigdb_conn, database_url=args.database_url,
                            schema=args.schema,
                        )
                        _derive_log(
                            'metsigdb_done',
                            build=metsigdb_stats.build_id,
                            rows=metsigdb_stats.rows,
                            resources=','.join(
                                f'{loaded.resource}={loaded.rows}'
                                for loaded in metsigdb_stats.resources
                            ),
                            seconds=f'{time.perf_counter() - step_started:.3f}',
                        )
                    except Exception as exc:
                        # A published dataset over the core content, like the
                        # network views: a failure is loud but must not abort
                        # the build that produced the content it reads.
                        conn.rollback()
                        _supplementary_step_failed(
                            failed_steps,
                            'metsigdb',
                            exc,
                            time.perf_counter() - step_started,
                        )
                    finally:
                        if metsigdb_conn is not None:
                            metsigdb_conn.close()
                if args.cosmos:
                    step_started = time.perf_counter()
                    _derive_log('cosmos_start')
                    try:
                        ensure_cosmos_edge_table(conn, schema=args.schema)
                        # The identifier mappings live in the utils database,
                        # which the manifest step above already reads through
                        # the same variable. A build that has none still
                        # projects: the labels keep the identifiers the build
                        # canonicalised them to and say so.
                        cosmos_stats = build_cosmos_projection(
                            conn,
                            schema=args.schema,
                            utils_db_url=os.environ.get(
                                'OMNIPATH_BUILD_UTILS_PG_URL'
                            ),
                        )
                        _derive_log(
                            'cosmos_done',
                            edges=cosmos_stats.edges,
                            reactions=cosmos_stats.reactions,
                            orphan_reactions=cosmos_stats.orphan_reactions,
                            connectors=cosmos_stats.connectors,
                            skipped_parties=cosmos_stats.skipped_parties,
                            labels=_cosmos_label_counts(cosmos_stats),
                            seconds=f'{time.perf_counter() - step_started:.3f}',
                        )
                    except Exception as exc:
                        # A projection over the reaction headers the core phase
                        # wrote, like the network views and MetSigDB: loud on
                        # failure, but it must not abort the build that
                        # produced the content it reads.
                        conn.rollback()
                        _supplementary_step_failed(
                            failed_steps,
                            'cosmos',
                            exc,
                            time.perf_counter() - step_started,
                        )
                print(
                    '[derive] '
                    f'entity_identifier_lookup='
                    f'{table_stats.entity_identifier_lookup} '
                    f'entity_relation_counts='
                    f'{table_stats.entity_relation_counts} '
                    f'entity_ontology_terms='
                    f'{table_stats.entity_ontology_terms} '
                    f'ontology_terms={table_stats.ontology_terms} '
                    f'entity_source_count={table_stats.entity_source_count} '
                    f'resources={resource_stats.resources}',
                    flush=True,
                )
            elif args.bitmaps:
                step_started = time.perf_counter()
                _derive_log('bitmaps_start')
                bitmap_stats = rebuild_bitmap_tables(
                    conn,
                    schema=args.schema,
                    progress=True,
                )
                _derive_log(
                    'bitmaps_done',
                    annotation_term_entities=(
                        bitmap_stats.annotation_term_entities
                    ),
                    annotation_term_direct_relations=(
                        bitmap_stats.annotation_term_direct_relations
                    ),
                    entity_relations=bitmap_stats.entity_relations,
                    entity_facets=bitmap_stats.entity_facets,
                    relation_facets=bitmap_stats.relation_facets,
                    seconds=f'{time.perf_counter() - step_started:.3f}',
                )
            else:
                bitmap_stats = None
            if bitmap_stats is not None:
                print(
                    '[derive] '
                    f'annotation_term_entities='
                    f'{bitmap_stats.annotation_term_entities} '
                    f'annotation_term_direct_relations='
                    f'{bitmap_stats.annotation_term_direct_relations} '
                    f'entity_relations={bitmap_stats.entity_relations} '
                    f'entity_facets={bitmap_stats.entity_facets} '
                    f'relation_facets={bitmap_stats.relation_facets}',
                    flush=True,
                )
            sweep_started = time.perf_counter()
            _derive_log('sweep_staging_start')
            swept = sweep_staging_tables(conn, schema=args.schema, progress=True)
            _derive_log(
                'sweep_staging_done',
                dropped=swept,
                seconds=f'{time.perf_counter() - sweep_started:.3f}',
            )
            _derive_log(
                'done',
                seconds=f'{time.perf_counter() - derive_started:.3f}',
                failed_steps=','.join(failed_steps),
            )
            return _derive_exit_code(failed_steps)

    return 0


def _derive_interactions(
    conn,
    *,
    schema: str = 'public',
    progress: bool = True,
) -> InteractionDeriveStats:
    """The interaction projection as a registered derive step.

    Runs after ``classify_interaction_class``, which fills the
    ``vocab_relation_predicate.interaction_class_id`` map the projection reads,
    and once per build: ``rebuild_derived_tables`` is called with
    ``interactions=False`` because the step is registered here.

    A failure is logged at error level in the structured derive shape and then
    **re-raised**: unlike the network views this is a phase of the build, not a
    supplementary layer, so a broken projection must abort the derive with a
    non-zero exit rather than pass as a successful build. A build that reports
    success must not have skipped a phase in silence.
    """
    step_started = time.perf_counter()
    _derive_log('interactions_start', schema=schema)
    try:
        stats = rebuild_interaction_tables(
            conn,
            schema=schema,
            progress=progress,
        )
    except Exception as exc:
        conn.rollback()
        _derive_log(
            'interactions_failed',
            _level=logging.ERROR,
            error=repr(exc),
            seconds=f'{time.perf_counter() - step_started:.3f}',
        )
        raise
    _derive_log(
        'interactions_done',
        interactions=stats.interactions,
        parties=stats.parties,
        records=stats.records,
        rows_by_class=','.join(
            f'{name}:{count}' for name, count in sorted(stats.rows_by_class.items())
        ),
        seconds=f'{stats.seconds:.3f}',
    )
    return stats


def _interaction_derive_cost(
    stats: InteractionDeriveStats | None,
) -> dict[str, dict[str, object]] | None:
    """What the interaction projection cost this run, for the build manifest.

    Read off the ``InteractionDeriveStats`` the derive step already returned
    rather than re-derived. **Seconds and rows are reported per step**, each
    taken from ``step_seconds``, rather than the projection's whole wall clock
    being attributed to one table — these are the numbers the projection's
    share of the build-cost budget is argued against, and a step total cannot
    say which step overran.

    The projection writes ``interaction_fact_resource`` and stops, so the
    ``interaction_fact_combined`` line is gone with the table it priced. The
    ``scope_cost`` its neighbour used to hand the manifest went with it — that
    neighbour reported one scope, the all-resources collapse, and the fold is a
    query now and materialises nothing, so this build has no scope to report
    rather than a scope reported as free.

    Returns ``None`` when the projection did not run, so ``--no-interactions``
    records no cost instead of claiming zeros.
    """
    if stats is None:
        return None
    step_seconds = stats.step_seconds or {}
    return {
        INTERACTION_RECORD_STEP: {
            'seconds': step_seconds.get(INTERACTION_RECORD_STEP),
            'rows': stats.records,
        },
        'interaction_party': {'rows': stats.parties},
        'interaction_header': {
            'seconds': step_seconds.get('interaction_header'),
            'rows': stats.interactions,
        },
        # The two staging steps carry no rows of their own and are reported
        # anyway, so the step seconds add up to the projection's wall clock
        # rather than leaving an unattributed remainder. A ceiling argument
        # made against a partial accounting is an argument about the wrong
        # number.
        'interaction_class_evidence': {
            'seconds': step_seconds.get('interaction_class_evidence'),
        },
        'interaction_evidence_fold': {
            'seconds': step_seconds.get('interaction_evidence_fold'),
        },
    }


def _interaction_deferral_cost(
    stats: InteractionDeriveStats | None,
) -> dict[str, object] | None:
    """What deferring the constraints over the load bought, for the manifest.

    Read off the record the derive step already returned, and handed to
    ``emit_build_manifest`` unchanged: the normalising is that function's job,
    including the rule that an unmeasured field stays ``null`` rather than
    becoming a zero, so that a build which ran without the deferral and a
    deferral that saved nothing do not read alike.

    Returns ``None`` when the projection did not run at all, so
    ``--no-interactions`` records no deferral rather than claiming one that
    saved nothing.
    """
    if stats is None:
        return None
    return dict(stats.deferral) or None


def _cosmos_label_counts(stats: CosmosBuildStats) -> str:
    """What the label translation reached, as one field of the derive log.

    Both sides in one field rather than six of their own. What a reader wants
    from the line is the shape of the translation — whether a mapping database
    answered at all, and how much of each side reached the namespace COSMOS
    asks for against how much kept the build's own identifier — and six
    near-identical keys bury that rather than state it. The MetSigDB step
    reports its per-resource row counts the same way.
    """
    return (
        f'mappings={stats.identifier_translation},'
        f'gene={stats.gene_labels_translated}/{stats.gene_labels_fallback},'
        f'chemical={stats.chemical_labels_translated}'
        f'/{stats.chemical_labels_fallback}'
    )


def _supplementary_step_failed(
    failed: list[str],
    step: str,
    exc: BaseException,
    seconds: float,
) -> None:
    """Record a supplementary step's failure so the build cannot hide it.

    These steps are caught on purpose: a published dataset over the core
    content must not abort the build that produced the content it reads. What
    was missing is that the operator could not tell. The error went to the
    package logger, and `make` exports a quiet logging config, so a derive that
    published nothing printed its usual summary and exited 0. That is the
    silently skipped phase Principle V forbids, arrived at through the back
    door: the phase was not skipped, it failed and said so to nobody.

    So the failure goes to `print` as well, which is the channel the steps'
    own progress already uses and no logging config suppresses, and the name is
    recorded for the exit code.
    """
    failed.append(step)
    _derive_log(
        f'{step}_failed',
        _level=logging.ERROR,
        error=repr(exc),
        seconds=f'{seconds:.3f}',
    )
    print(f'[derive] STEP FAILED step={step} error={exc!r}', flush=True)


def _derive_exit_code(failed: list[str]) -> int:
    """0 when every step ran, 1 when a supplementary one did not.

    The build still completed, and everything that did run is committed. The
    non-zero code says the result is not the whole build, so a caller that
    treats exit 0 as "the database is ready to serve" is not misled.
    """
    if not failed:
        return 0
    print(
        f'[derive] incomplete: {len(failed)} step(s) failed: '
        + ', '.join(failed),
        flush=True,
    )
    return 1


def _record_cap(max_records: object) -> int | None:
    """``MAX_RECORDS`` as a cap the MetSigDB step can apply, or None.

    Same input as :func:`_is_partial_build` and the same tolerance for an unset
    or unparseable value: an absent cap is the authoritative case, so anything
    that is not a positive integer means "do not cap" rather than an error.
    """
    if max_records in (None, ''):
        return None
    try:
        cap = int(max_records)
    except (TypeError, ValueError):
        return None
    return cap if cap > 0 else None


def _is_partial_build(max_records: object) -> bool:
    """A non-zero MAX_RECORDS cap means the load was truncated (not authoritative)."""
    if max_records in (None, ''):
        return False
    try:
        return int(max_records) > 0
    except (TypeError, ValueError):
        return False


def _derive_log(event: str, _level: int = logging.INFO, **fields: object) -> None:
    """One structured derive-orchestration line.

    The ``event=… key=value`` shape is a contract, not a preference: the
    sign-conflict summary and the derive-cost figures the build manifest
    records are read out of this ``--progress`` output. Only the sink changed,
    from ``print`` to the package logger, which `omnipath_build._session`
    configures.
    """
    details = ' '.join(f'{key}={value}' for key, value in fields.items())
    _logger.log(
        _level,
        '[derive] event=%s%s',
        event,
        f' {details}' if details else '',
    )


if __name__ == '__main__':
    raise SystemExit(main())
