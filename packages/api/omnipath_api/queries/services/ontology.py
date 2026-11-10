"""Ontology terms, children and trees, from the step files."""

from __future__ import annotations

import json
from typing import Any


class OntologyQueries:
    """Ontology search, term metadata, children and ancestry trees."""

    def get_ontology_tree(
        self,
        term_ids: list[str],
        ontology_id: str | None = None,
        resources: list[str] | None = None,
    ) -> dict[str, Any]:
        """An ontology term's ancestry as a tree: the query ontology_tree."""
        if not term_ids:
            return {"root": None}
        params = {"term_id": str(term_ids[0]).strip(), "ontology_id": ontology_id}
        return json.loads(self.run_query("ontology_tree", params, resources))

    def get_ontology_children(
        self,
        term_id: str,
        ontology_id: str | None = None,
        resources: list[str] | None = None,
        limit: int = 100,
        offset: int = 0,
    ) -> dict[str, Any]:
        """A page of an ontology term's children: the query ontology_children."""
        params = {"term_id": str(term_id), "ontology_id": ontology_id, "limit": int(limit), "offset": int(offset)}
        return json.loads(self.run_query("ontology_children", params, resources))

    def search_ontology(self, payload=None, resources=None, counts=False):
        """Ontology terms with their annotation counts: the query ontology_search."""
        payload = dict(payload or {})
        scope = payload.get("selectionScope")
        params = {
            "query": str(payload.get("query") or payload.get("q") or ""),
            "limit": int(payload.get("limit") or 24),
            "offset": int(payload.get("offset") or 0),
            "ontology_ids": list(payload.get("ontologyIds") or []) or None,
            "prefixes": list(payload.get("prefixes") or []) or None,
            "counts": bool(counts),
            "scoped": bool(payload.get("scoped", False)),
            "entity_pks": [str(v) for v in (payload.get("entityPks") or [])] or None,
            "term_ids": [str(v) for v in (payload.get("termIds") or payload.get("annotationTermIds") or [])] or None,
            "has_selection": scope is not None,
        }
        if scope is not None:
            params.update(
                selection_entity_pks=[str(v) for v in (scope.get("entityPks") or [])] or None,
                selection_include_associated=scope.get("includeAssociatedEntities") is not False,
                selection_term_ids=[str(v) for v in (scope.get("annotationTermIds") or [])] or None,
            )
        return json.loads(self.run_query("ontology_search", params, resources))

    def get_terms_info(
        self, term_ids: list[str], resources: list[str] | None = None
    ) -> dict[str, Any]:
        """Ontology term metadata by id: the query ontology_terms."""
        if not term_ids:
            return {"terms": {}}
        return json.loads(self.run_query("ontology_terms", {"term_ids": [str(t) for t in term_ids]}, resources))
