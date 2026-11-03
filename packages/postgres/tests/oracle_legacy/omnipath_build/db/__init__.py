"""PostgreSQL schema, indexes, and derived structures for omnipath_build."""

from omnipath_build.db.schema import (
    ensure_schema,
    ensure_content_primary_keys,
    reset_content_tables,
    ensure_deferred_indexes,
    ensure_source_partitions,
    drop_deferred_content_indexes,
)
from omnipath_build.db.bitmaps import BitmapStats, rebuild_bitmap_tables
from omnipath_build.db.indexes import create_secondary_indexes
from omnipath_build.db.refresh import (
    SourceContentDropStats,
    delete_source_content,
    source_has_content,
)
from omnipath_build.db.licenses import (
    LicenseTableStats,
    sync_data_source_licenses,
)
from omnipath_build.db.resources import (
    BuildManifestStats,
    ResourceTableStats,
    emit_build_manifest,
    sync_resources_table,
)
from omnipath_build.db.derived_tables import (
    DerivedTableStats,
    InteractionDeriveStats,
    interaction_content_uuid_sql,
    rebuild_derived_tables,
    rebuild_interaction_tables,
    rebuild_resource_overlap_summary,
    sweep_staging_tables,
)

__all__ = [
    'BitmapStats',
    'BuildManifestStats',
    'DerivedTableStats',
    'InteractionDeriveStats',
    'LicenseTableStats',
    'interaction_content_uuid_sql',
    'rebuild_interaction_tables',
    'ResourceTableStats',
    'SourceContentDropStats',
    'create_secondary_indexes',
    'delete_source_content',
    'drop_deferred_content_indexes',
    'ensure_deferred_indexes',
    'ensure_content_primary_keys',
    'emit_build_manifest',
    'ensure_schema',
    'ensure_source_partitions',
    'rebuild_bitmap_tables',
    'rebuild_derived_tables',
    'rebuild_resource_overlap_summary',
    'reset_content_tables',
    'sweep_staging_tables',
    'source_has_content',
    'sync_data_source_licenses',
    'sync_resources_table',
]
