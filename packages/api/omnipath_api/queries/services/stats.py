"""Stats query service, including domain SQL and result shaping."""

from __future__ import annotations

import logging
from typing import Any


logger = logging.getLogger(__name__)


def entity_types_query(read_expr: str) -> str:
    return (
        f"SELECT entity_type, count(*) AS n "
        f"FROM {read_expr} "
        f"WHERE entity_type IS NOT NULL "
        f"GROUP BY 1 "
        f"ORDER BY n DESC"
    )


def interaction_types_query(read_expr: str) -> str:
    return (
        f"SELECT predicate, interaction_class, count(*) AS n "
        f"FROM {read_expr} "
        f"GROUP BY 1, 2 "
        f"ORDER BY n DESC"
    )


class StatsQueries:
    """Stats queries over ParquetQueries' storage and shaping contract."""

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
        if not self._selected_resource_infos():
            return []
        rows = self._db.execute(entity_types_query(self._table("entity"))).fetchall()
        return [{"entityType": row[0], "count": int(row[1])} for row in rows if row[0]]

    def get_stats_interaction_types(self) -> list[dict[str, Any]]:
        if not self._selected_resource_infos():
            return []
        rows = self._db.execute(interaction_types_query(self._table("relation"))).fetchall()
        return [
            {"interactionType": row[0], "interactionClass": row[1], "count": int(row[2])}
            for row in rows
        ]

    def get_stats_build_manifest(self) -> dict[str, Any]:
        catalog = self.list_resource_catalog()
        release = self._release_scope.get()
        return {
            "buildId": release["version"] if release else "latest",
            "builtAt": release.get("created_at", "") if release else "",
            "partialBuild": any(
                item.get("sample_build") or item.get("dataset_subset") for item in catalog
            ),
            "engine": "duckdb-parquet",
            "resources": [item["resource_id"] for item in catalog],
            "resourceVersions": {item["resource_id"]: item["version"] for item in catalog},
        }
