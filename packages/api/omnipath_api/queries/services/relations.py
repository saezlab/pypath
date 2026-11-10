"""Relation queries: pages and single relations, from the step files."""

from __future__ import annotations

import json
import time
from typing import Any


class RelationsQueries:
    """Relation pages and relations by key."""

    def search_relations(
        self,
        filters: dict[str, Any] | None = None,
        resources: list[str] | None = None,
        limit: int = 50,
        offset: int = 0,
        include_details: bool = False,
    ) -> dict[str, Any]:
        """A page of matching relations as stored (with details, their annotations and
        evidence): the query relation_page's records."""
        from omnipath_api.queries.services.export import relation_params

        started = time.perf_counter()
        params = {
            **relation_params(filters), "limit": int(limit), "offset": int(offset),
            "details": bool(include_details), "records": True,
        }
        result = json.loads(self.run_query("relation_page", params, resources))
        return {
            "rows": result["records"],
            "total": int(result["total"]),
            "limit": limit,
            "offset": offset,
            "elapsed_ms": round((time.perf_counter() - started) * 1000, 2),
        }

    def search_relations_api(self, filters=None, resources=None, limit=20, offset=0):
        """A page of matching relations with both endpoints: the query relation_page,
        which the clients run too."""
        from omnipath_api.queries.services.export import relation_params

        started = time.perf_counter()
        params = {**relation_params(filters), "limit": int(limit), "offset": int(offset)}
        result = json.loads(self.run_query("relation_page", params, resources))
        result["elapsed_ms"] = round((time.perf_counter() - started) * 1000, 2)
        return result

    def get_relation_by_pk(
        self, relation_pk: str, resources: list[str] | None = None
    ) -> dict[str, Any] | None:
        text = self.run_query("relation", {"relation_key": str(relation_pk)}, resources)
        return json.loads(text) if text else None
