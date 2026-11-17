"""MetSigDB: canonical metabolite memberships in resource-native sets.

One table, ``metsigdb_membership``, built from the published release as Parquet
(``omnipath-subsets build``): memberships in sets from KEGG, Reactome, WikiPathways, MACdb
and ClassyFire. The step reads the release's per-resource Parquet. It parses no upstream
resource and downloads nothing: the resources already hold every membership, resolved
against the release's entities.
"""

from __future__ import annotations

from omnipath_subsets.metsigdb.mapping import RESOURCES, ResourceRule, rule_for
from omnipath_subsets.metsigdb.parquet_build import write_metsigdb

__all__ = ["RESOURCES", "ResourceRule", "rule_for", "write_metsigdb"]
