"""MetSigDB unified membership substrate (cycle 010).

One build-time materialized table, ``metsigdb_membership``, holding canonical
metabolite memberships in resource-native sets from KEGG, Reactome,
WikiPathways, MACdb and ClassyFire.

The step reads the build database, not the upstream resources. The core build
already loaded all five through ``pypath.inputs_v2`` and already resolved their
entities, so the memberships sit in ``relation_evidence``,
``entity_evidence_resolution`` and ``entity_ontology_relation``. Reading those
needs no network access, and it keeps the metabolite side identical to every
other consumer of the canonical entity layer.
"""

from __future__ import annotations

from omnipath_build.metsigdb.build import (
    BuildStats,
    ResourceLoadStats,
    build_id,
    build_metsigdb,
    ensure_membership_table,
    load_resource,
)
from omnipath_build.metsigdb.mapping import RESOURCES, ResourceRule, rule_for

__all__ = [
    'RESOURCES',
    'BuildStats',
    'ResourceLoadStats',
    'ResourceRule',
    'build_id',
    'build_metsigdb',
    'ensure_membership_table',
    'load_resource',
    'rule_for',
]
