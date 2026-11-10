"""Facets and filter options, from the step files."""

from __future__ import annotations

import json
import time
from typing import Any

from omnipath_api.queries.filters import normalize_filters


class FacetsQueries:
    """Relation and entity filter options a selection leaves; relation counts."""

    def get_facets(
        self,
        filters: dict[str, Any] | None = None,
        resources: list[str] | None = None,
    ) -> dict[str, Any]:
        """The matching relations counted by category, predicate and sign: the query
        relation_counts."""
        from omnipath_api.queries.services.export import relation_params

        started = time.perf_counter()
        result = json.loads(self.run_query("relation_counts", relation_params(filters), resources))
        result["elapsed_ms"] = round((time.perf_counter() - started) * 1000, 2)
        return result

    def get_scoped_entity_facets(self, payload=None, resources=None):
        """The entity filter options a search leaves: the query entity_facets."""
        payload = dict(payload or {})
        filters = normalize_filters(payload)
        selection = resources or payload.get("resources") or None
        tokens = list(filters["entity_ids"]) + list(filters["entity_pks"])
        params = {
            "query": str(payload.get("query") or ""),
            "sources": list(filters["sources"]) or None,
            "entity_types": list(filters["entity_types"]) or None,
            "taxonomy_ids": [str(t) for t in filters["taxonomy_ids"]] or None,
            "reference_entity_keys": list(filters["reference_entity_keys"]) or None,
            "entity_ids": tokens or None,
            "facet_limit": int(payload.get("facetLimit") or 50),
            "selection": list(selection) if selection else None,
        }
        text = self.run_query("entity_facets", params)
        return json.loads(text) if text else []

    def get_scoped_relation_facets(self, payload=None, resources=None):
        """The relation filter options a selection leaves: the query relation_facets."""
        from omnipath_api.queries.services.export import relation_params

        payload = dict(payload or {})
        limit = payload.pop("taxonomyLimit", None)
        query = str(payload.pop("taxonomyQuery", "") or "")
        params = {
            **relation_params(_flattened(payload)),
            "taxonomy_limit": None if limit is None else int(limit),
            "taxonomy_query": query,
        }
        return json.loads(self.run_query("relation_facets", params, resources))

    def get_entity_filter_options(self, resources: list[str] | None = None) -> dict[str, Any]:
        return json.loads(self.run_query("entity_filter_options", {}, resources))

    def get_relation_filter_options(self, resources: list[str] | None = None) -> dict[str, Any]:
        return json.loads(self.run_query("relation_filter_options", {}, resources))


def _flattened(payload: dict[str, Any]) -> dict[str, Any]:
    """The payload with its nested ``filters`` merged into the top level."""
    out = dict(payload)
    for key, value in (out.pop("filters", None) or {}).items():
        outer = out.get(key)
        if isinstance(value, list) and isinstance(outer, list):
            out[key] = outer + value
        elif outer is None:
            out[key] = value
    return out
