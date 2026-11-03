"""Selection query service, including domain SQL and result shaping."""

from __future__ import annotations

import logging
from typing import Any


logger = logging.getLogger(__name__)


class SelectionQueries:
    """Selection queries over ParquetQueries' storage and shaping contract."""

    def resolve_selection_scope(self, payload=None):
        return self._cached_facets(
            "selection-scope",
            dict(payload or {}),
            None,
            lambda p, _: self._resolve_selection_scope_uncached(p),
        )

    def _resolve_selection_scope_uncached(
        self, payload: dict[str, Any] | None = None
    ) -> dict[str, Any]:
        payload = payload or {}
        seeds = [str(v).strip() for v in (payload.get("entityPks") or []) if str(v).strip()]
        resolved = self.resolve_entity_keys(seeds) if seeds else []
        entity_pks = resolved or seeds
        include_associated = payload.get("includeAssociatedEntities") is not False
        if include_associated and entity_pks:
            # Neighbours: the other endpoint of each relation of the seeds.
            pairs = self._fetch_dicts(
                f"""SELECT DISTINCT resource, relation_id FROM {self._table("relation_endpoint")}
                WHERE key IN (SELECT unnest(?::VARCHAR[])) AND key_kind = 'entity'""",
                [entity_pks],
            )
            seen = set(entity_pks)
            for row in self._lookup(
                "relation",
                "relation_id",
                [(r["resource"], r["relation_id"]) for r in pairs],
                "subject_entity_key, object_entity_key",
            ):
                for key in (row["subject_entity_key"], row["object_entity_key"]):
                    if key and key not in seen:
                        seen.add(key)
                        entity_pks.append(key)
        return {
            "entityPks": entity_pks,
            "seedEntityPks": seeds,
            "termEntityPks": [],
            "ontologyTermIds": [str(v) for v in (payload.get("annotationTermIds") or [])],
            "criteriaCount": len(seeds),
            "expandedEntityCount": max(0, len(entity_pks) - len(seeds)),
        }
