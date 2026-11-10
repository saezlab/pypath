"""Stats query service, including domain SQL and result shaping."""

from __future__ import annotations

import json
import logging
from typing import Any


logger = logging.getLogger(__name__)


class StatsQueries:
    """Stats queries over the release database and the shaping contract."""

    def get_stats_sources(self) -> list[dict[str, Any]]:
        return [
            {
                "slug": item["resource_id"],
                "short": item["resource_name"],
                "full": item["resource_name"],
                "entityCount": item.get("entity_count") or 0,
                "interactionCount": item.get("interaction_count") or 0,
                "associationCount": item.get("association_count") or 0,
                "identifierCount": item.get("identifier_count") or 0,
                "ontologyTermCount": item.get("ontology_term_count") or 0,
            }
            for item in self.list_resource_catalog()
            if item["queryable"]
        ]

    def get_stats_entity_types(self) -> list[dict[str, Any]]:
        return json.loads(self.run_query("stats_entity_types"))

    def get_stats_interaction_types(self) -> list[dict[str, Any]]:
        return json.loads(self.run_query("stats_interaction_types"))

    def get_stats_build_manifest(self) -> dict[str, Any]:
        catalog = self.list_resource_catalog()
        release = self._release_scope.get()
        return {
            "buildId": release["version"] if release else "latest",
            "builtAt": release.get("created_at", "") if release else "",
            "partialBuild": any(
                item.get("sample_build") or item.get("dataset_subset") for item in catalog
            ),
            "engine": "duckdb-release",
            "resources": [item["resource_id"] for item in catalog],
            "resourceVersions": {item["resource_id"]: item["version"] for item in catalog},
        }
