"""Entity selections, from the step files."""

from __future__ import annotations

import json


class SelectionQueries:
    """An entity selection with its neighbours."""

    def resolve_selection_scope(self, payload=None):
        """An entity selection with its neighbours: the query selection_scope."""
        payload = dict(payload or {})
        params = {
            "entity_pks": [str(v) for v in (payload.get("entityPks") or [])] or None,
            "include_associated": payload.get("includeAssociatedEntities") is not False,
            "annotation_term_ids": [str(v) for v in (payload.get("annotationTermIds") or [])] or None,
        }
        return json.loads(self.run_query("selection_scope", params))
