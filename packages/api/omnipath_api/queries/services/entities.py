"""Entity queries: each runs a step file (omnipath_client/sql/queries), the same the
clients run."""

from __future__ import annotations

import json
import time
from typing import Any

from omnipath_api.queries.filters import EntitySearchCursor, normalize_filters


class EntitiesQueries:
    """Entity search, details, evidence, groups and molecular context."""

    def resolve_entity_keys(
        self,
        tokens: list[str] | None,
        resources: list[str] | None = None,
    ) -> list[str]:
        values = [str(t) for t in (tokens or [])]
        if not any(v.strip() for v in values):
            return []
        return json.loads(self.run_query("resolve_keys", {"tokens": values}, resources))

    def search_entity_groups(self, *, strategy="chemical_connectivity", **kwargs):
        from omnipath_api.queries.services.groups import search_groups

        if strategy not in {"auto", "chemical_connectivity", "gene_reference"}:
            raise ValueError(f"Unknown grouping strategy: {strategy}")
        return search_groups(self, strategy=strategy, **kwargs)

    def search_connectivity_groups(self, **kwargs):
        return self.search_entity_groups(strategy="chemical_connectivity", **kwargs)

    def get_molecular_context(
        self, entity_id, resources=None, *, view="reference", isoform_identifier=None, limit=20, offset=0
    ):
        """An entity's molecular context: the query molecular_context."""
        params = {
            "entity_id": str(entity_id), "view": view, "isoform_identifier": isoform_identifier,
            "limit": int(limit), "offset": int(offset),
        }
        text = self.run_query("molecular_context", params, resources)
        return json.loads(text) if text else None

    def get_entity_examples(self):
        from omnipath_api.queries.examples import get_examples

        return get_examples(self)

    def search_entities_api(
        self,
        query: str = "",
        filters: dict[str, Any] | None = None,
        resources: list[str] | None = None,
        limit: int = 20,
        cursor: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """A page of matching entities: the query ``entity_search``, which the clients
        run too."""
        started = time.perf_counter()
        filters = normalize_filters(filters)
        phase = offset = key = None
        if cursor:
            cursor = EntitySearchCursor.model_validate(cursor).model_dump()
            phase, offset, key = cursor["phase"], cursor["offset"], cursor["queryKey"]
        tokens = list(filters["entity_ids"]) + list(filters["entity_pks"])
        params = {
            "query": query or "",
            "sources": list(filters["sources"]) or None,
            "entity_types": list(filters["entity_types"]) or None,
            "taxonomy_ids": [str(t) for t in filters["taxonomy_ids"]] or None,
            "reference_entity_keys": list(filters["reference_entity_keys"]) or None,
            "entity_ids": tokens or None,
            "limit": int(limit),
            "cursor_phase": phase,
            "cursor_offset": offset,
            "cursor_key": key,
        }
        result = json.loads(self.run_query("entity_search", params, resources))
        result["elapsed_ms"] = round((time.perf_counter() - started) * 1000, 2)
        return result

    def get_entities_by_pks(
        self,
        pks: list[str],
        resources: list[str] | None = None,
        *,
        slim: bool = False,
    ) -> dict[str, Any]:
        values = [str(v).strip() for v in pks if str(v).strip()]
        if not values:
            return {"entities": []}
        text = self.run_query("entities_by_keys", {"tokens": values, "slim": bool(slim)}, resources)
        return json.loads(text)

    def get_entities_by_public_ids(
        self, public_ids: list[str], resources: list[str] | None = None
    ) -> dict[str, Any]:
        return self.get_entities_by_pks(public_ids, resources)

    def get_entity_core(self, public_id, resources=None):
        text = self.run_query("entity_core", {"entity_key": str(public_id)}, resources)
        return json.loads(text) if text else None

    def get_entity_evidence(self, public_id, resources=None, limit=20, offset=0):
        text = self.run_query(
            "entity_evidence",
            {"entity_key": str(public_id), "limit": int(limit), "offset": int(offset)},
            resources,
        )
        return json.loads(text) if text else None

    def get_entity_relationships(self, public_id, resources=None, limit=50, offset=0):
        keys = public_id if isinstance(public_id, list) else [public_id]
        text = self.run_query(
            "entity_relationships",
            {"entity_keys": [str(k) for k in keys], "limit": int(limit), "offset": int(offset)},
            resources,
        )
        return json.loads(text)

    def get_entity_details(
        self, public_id: str, resources: list[str] | None = None
    ) -> dict[str, Any] | None:
        # Details are addressed solely by the exact entity key. Identifier and
        # alias discovery belongs to search, never to single-entity hydration.
        text = self.run_query("entity_details", {"entity_key": str(public_id)}, resources)
        return json.loads(text) if text else None

    def resolve_identifiers(
        self, identifiers: list[str], resources: list[str] | None = None
    ) -> dict[str, Any]:
        values = [str(v).strip() for v in identifiers if str(v).strip()]
        if not values:
            return {"matches": [], "entities": []}
        found = self.get_entities_by_pks(values, resources).get("entities") or []
        if not found:
            for value in values[:20]:
                extra = (
                    self.search_entities_api(query=value, limit=5, resources=resources).get(
                        "entities"
                    )
                    or []
                )
                found.extend(extra)
        unique: dict[str, dict[str, Any]] = {}
        for entity in found:
            unique[entity["entityPk"]] = entity
        entities = list(unique.values())
        matches = []
        for identifier in values:
            needle = identifier.lower()
            entity_ids = []
            for entity in entities:
                public_id = f"{entity['canonicalIdentifierType']}|{entity['canonicalIdentifier']}"
                haystacks = [
                    entity["entityPk"].lower(),
                    entity["canonicalIdentifier"].lower(),
                    public_id.lower(),
                    str(entity.get("label") or "").lower(),
                ]
                if any(needle == h or needle in h for h in haystacks if h):
                    entity_ids.append(public_id)
            matches.append({"identifier": identifier, "entityIds": entity_ids})
        return {"matches": matches, "entities": entities}
