"""Facets query service, including domain SQL and result shaping."""

from __future__ import annotations

import logging
import json
import time
from typing import Any


from omnipath_api.queries.store.query_cache import normalized


import pyarrow as pa
from omnipath_core.schema import RELATION_FACET_COLUMNS, RELATION_FACET_TABLE

from omnipath_api.queries.services.constants import RELATION_QUALIFIER_FILTERS
from omnipath_api.queries.filters import normalize_filters
from omnipath_api.queries.molecular import has_form_filters

logger = logging.getLogger(__name__)


class FacetsQueries:
    """Facets queries over ParquetQueries' storage and shaping contract."""

    def get_facets(
        self,
        filters: dict[str, Any] | None = None,
        resources: list[str] | None = None,
    ) -> dict[str, Any]:
        """Compute aggregated facet counts in a single vectorized pass."""
        t0 = time.perf_counter()
        selection, params = self._relation_selection(filters, resources)
        facet_sql = f"""
            SELECT category, predicate, sign, count(*) AS n
            FROM ({selection})
            GROUP BY category, predicate, sign
        """
        rows = self._db.execute(facet_sql, params).fetchall()

        categories: dict[str, int] = {}
        predicates: dict[str, int] = {}
        signs: dict[str, int] = {}

        for cat, pred, sign, count in rows:
            if cat:
                categories[cat] = categories.get(cat, 0) + count
            if pred:
                predicates[pred] = predicates.get(pred, 0) + count
            if sign is not None:
                s_key = str(sign)
                signs[s_key] = signs.get(s_key, 0) + count

        # Sort predicates by count descending
        sorted_predicates = dict(sorted(predicates.items(), key=lambda item: -item[1])[:25])
        sorted_categories = dict(sorted(categories.items(), key=lambda item: -item[1]))

        elapsed_ms = round((time.perf_counter() - t0) * 1000, 2)

        return {
            "categories": sorted_categories,
            "predicates": sorted_predicates,
            "signs": signs,
            "elapsed_ms": elapsed_ms,
        }

    def _cached_facets(self, kind, payload, resources, compute):
        # Paths identify pinned releases; fingerprint also invalidates moving latest.
        infos = self._selected_resource_infos(None)
        key = (
            kind,
            self._inventory_fingerprint,
            self._taxonomy_cache_version(),
            tuple(sorted(i["key"] for i in infos)),
            json.dumps(
                normalized({"payload": payload or {}, "resources": resources}), sort_keys=True
            ),
            json.dumps(self._release_scope.get(), sort_keys=True, default=str),
            json.dumps(
                [m.get("references") for m in self.releases.list()], sort_keys=True, default=str
            ),
        )
        return self._facet_cache.get(key, lambda: compute(payload, resources))

    def get_scoped_entity_facets(self, payload=None, resources=None):
        return self._cached_facets("entities", payload, resources, self._compute_entity_facets)

    def get_scoped_relation_facets(self, payload=None, resources=None):
        """The relation filter options a selection leaves: for each filter group, the values
        of the relations that the scope (an entity selection, molecular forms) and the
        other groups' filters keep, most relations first. A group's own filters do not
        narrow it, so its other values stay listed (choosing one adds it)."""
        payload = dict(payload or {})
        limit = payload.pop("taxonomyLimit", None)
        query = str(payload.pop("taxonomyQuery", "") or "").strip().casefold()
        items = self._cached_facets("relations", payload, resources, self._compute_relation_facets)
        if limit is None:
            return items  # Existing API callers can still request the complete list.
        limit = max(1, min(int(limit), 10000))
        selected = set(map(str, payload.get("taxonomyIds") or payload.get("ncbi_tax_id") or []))
        taxa = [i for i in items if i["facetName"] == "taxonomy_id"]
        present = {i["facetValue"] for i in taxa}
        missing = self._label_taxonomy_facets(
            [
                dict(facetName="taxonomy_id", facetValue=value, facetCategory=None)
                for value in sorted(selected - present)
            ]
        )
        items.extend(missing)
        taxa.extend(missing)
        matching = [
            i
            for i in taxa
            if not query
            or query in i["facetValue"].casefold()
            or query in str(i.get("facetLabel") or "").casefold()
        ]
        visible = {i["facetValue"] for i in matching[:limit]} | selected
        return [i for i in items if i["facetName"] != "taxonomy_id" or i["facetValue"] in visible]

    def _selected_combinations(self, scope_payload, resources):
        """The ``relation_facet`` rows (with ``resource``) of the relations an endpoint or
        molecular form scope selects, as columns of plain values."""
        scope = normalize_filters(_flattened(scope_payload))
        selection, params = self._relation_selection(scope, resources, qualifiers=())
        return (
            self._db.execute(
                f"""SELECT {", ".join(SELECTED_COMBINATIONS.names)}
                FROM (SELECT {", ".join(RELATION_FACET_COLUMNS)}, resource,
                      max(evidence_count) AS evidence_count, count(*) AS relation_count
                      FROM ({selection}) GROUP BY ALL)""",
                params,
            )
            .fetch_arrow_table()
            .to_pydict()
        )

    def _compute_relation_facets(
        self,
        payload: dict[str, Any] | None = None,
        resources: list[str] | None = None,
    ) -> list[dict[str, Any]]:
        payload = payload or {}
        filters = normalize_filters(_flattened(payload))
        scope_payload = _without(payload, RELATION_FILTER_KEYS)
        scope = normalize_filters(_flattened(scope_payload))
        if not self._selected_resource_infos(resources):
            return []
        if any(scope[key] for key in _ENDPOINT_SCOPE) or has_form_filters(scope):
            # The selection's combinations are kept: a filter change reads no relations.
            columns = self._cached_facets(
                "relation-combinations", scope_payload, resources, self._selected_combinations
            )
            selected = pa.Table.from_pydict(columns, schema=SELECTED_COMBINATIONS)
            self._db.register("selected_combinations", selected)
            try:
                rows = self._facet_rows("SELECT * FROM selected_combinations", [], filters)
            finally:
                self._db.unregister("selected_combinations")
        else:
            where, params = self._relation_where(scope, qualifiers=())
            table = self._table("relation_facet", resources)
            rows = self._facet_rows(f"SELECT * FROM {table} WHERE {where}", params, filters)
        out = [
            {"facetName": name, "facetValue": value, "facetCategory": category}
            for name, value, category, _ in sorted(
                rows, key=lambda row: (-row[3], row[0], row[1], row[2] or "")
            )
        ]
        return self._label_taxonomy_facets(out)

    def _facet_rows(self, combinations, params, filters):
        """(facet, value, category, relations) of each option over the ``combinations``
        SQL, each facet under the filters of the other groups."""
        parts, params = [], list(params)
        for name, (keys, value, category) in RELATION_FACETS.items():
            where, values = self._relation_where(
                {key: item for key, item in filters.items() if key not in keys}
            )
            parts.append(
                f"""SELECT '{name}' AS facet, value, category, sum(relation_count) AS n
                FROM (SELECT {value} AS value, {category} AS category, relation_count
                      FROM combinations WHERE {where})
                WHERE value IS NOT NULL AND value <> '' GROUP BY ALL"""
            )
            params += values
        return self._db.execute(
            f"WITH combinations AS MATERIALIZED ({combinations}) {' UNION ALL '.join(parts)}",
            params,
        ).fetchall()

    def _entity_scoped_facet_where(
        self,
        payload: dict[str, Any],
        *,
        apply_sources: bool = True,
    ) -> tuple[str, list[Any], list[str]]:
        payload = payload or {}
        query = str(payload.get("query") or "").strip()
        filters = normalize_filters(payload)
        scope = (
            filters["sources"] if apply_sources and filters["sources"] else payload.get("resources")
        )
        if not apply_sources:
            filters["sources"] = []
        where_clauses, params = self._entity_filter_clauses(filters, resources=scope)
        where, params = self._entity_match_where(query, where_clauses, params, scope)
        return where, params, scope

    def _compute_entity_facets(
        self,
        payload: dict[str, Any] | None = None,
        resources: list[str] | None = None,
    ) -> list[dict[str, Any]]:
        payload = dict(payload or {})
        if resources:
            payload["resources"] = resources
        where_sql, params, scope = self._entity_scoped_facet_where(payload, apply_sources=True)
        if not self._selected_resource_infos(scope):
            return []
        groups = self._db.execute(
            f"""SELECT CASE WHEN GROUPING(entity_type)=0 THEN 'entity_type' ELSE 'taxonomy_id' END,
                       coalesce(entity_type, taxon), count(*)
                FROM {self._table("entity", scope)} WHERE {where_sql}
                GROUP BY GROUPING SETS ((entity_type), (taxon))""",
            params,
        ).fetchall()
        # Source counts ignore the source filter, so other sources stay selectable.
        source_where, source_params, _ = self._entity_scoped_facet_where(
            payload, apply_sources=False
        )
        infos = self._selected_resource_infos(None)
        source_counts = dict(
            self._db.execute(
                f"SELECT resource, count(*) FROM {self._table('entity')} WHERE {source_where} GROUP BY resource",
                source_params,
            ).fetchall()
        )
        source_items = [
            dict(
                facetName="source",
                facetValue=info["resource"],
                scopedCount=source_counts.get(info["key"], 0),
            )
            for info in infos
        ]
        limit = int(payload.get("facetLimit") or 50)
        items = [
            dict(facetName=name, facetValue=str(value), scopedCount=count)
            for name, value, count in groups
            if value
        ]
        items.extend(source_items)
        items.sort(key=lambda item: (-item["scopedCount"], item["facetName"], item["facetValue"]))
        type_tax_limit = max(limit * 2, 20)
        typed = [item for item in items if item["facetName"] != "source"][:type_tax_limit]
        sources = sorted(source_items, key=lambda item: (-item["scopedCount"], item["facetValue"]))
        return self._label_taxonomy_facets(typed + sources)

    def get_entity_filter_options(self, resources: list[str] | None = None) -> dict[str, Any]:
        if not self._selected_resource_infos(resources):
            return {"entity_types": [], "sources": [], "taxonomy_ids": []}
        read_expr = self._table("entity", resources)
        types = [
            r[0]
            for r in self._db.execute(
                f"SELECT DISTINCT entity_type FROM {read_expr} WHERE entity_type IS NOT NULL ORDER BY 1"
            ).fetchall()
        ]
        taxons = [
            str(r[0])
            for r in self._db.execute(
                f"SELECT DISTINCT taxon FROM {read_expr} WHERE taxon IS NOT NULL ORDER BY 1"
            ).fetchall()
        ]
        sources = sorted({info["resource"] for info in self._selected_resource_infos(resources)})
        return {"entity_types": types, "sources": sources, "taxonomy_ids": taxons}

    def get_relation_filter_options(self, resources: list[str] | None = None) -> dict[str, Any]:
        facets = self.get_scoped_relation_facets({}, resources=resources)
        predicates_by_category: dict[str, list[str]] = {}
        sources: list[str] = []
        interaction_types: list[str] = []
        for item in facets:
            if item["facetName"] == "predicate":
                category = item.get("facetCategory") or "uncategorized"
                predicates_by_category.setdefault(category, []).append(item["facetValue"])
            elif item["facetName"] == "source":
                sources.append(item["facetValue"])
            elif item["facetName"] == "participant_type":
                interaction_types.append(item["facetValue"])
        return {
            "predicatesByCategory": {k: sorted(set(v)) for k, v in predicates_by_category.items()},
            "sources": sorted(set(sources)),
            "interactionTypes": sorted(set(interaction_types)),
        }


# The payload keys of the relation filter panels, with their ``filters.normalize_filters``
# aliases; the relation counts leave them out (the scope keys stay).
RELATION_FILTER_KEYS = (
    "sources",
    "resources",
    "taxonomy_ids",
    "taxonomyIds",
    "ncbi_tax_id",
    "predicates",
    "relation_categories",
    "categories",
    "object_aspect_qualifier",
    "object_direction_qualifier",
    "causal_mechanism_qualifier",
    "interaction_types",
    "interactionTypes",
    "entity_types",
    "entityTypes",
    "partner_types",
    "partnerTypes",
    "signs",
    "min_evidence_count",
)


# The filters that select relations by their endpoints rather than by their columns.
_ENDPOINT_SCOPE = ("entity_ids", "entity_pks", "scope_entity_ids", "reference_entity_keys")
# The combinations of selected relations, with the resource of each.
SELECTED_COMBINATIONS = RELATION_FACET_TABLE.append(pa.field("resource", pa.string()))

# Each relation facet: the filters of its group (normalized names), which do not narrow
# it, and the SQL of its value and category over a ``relation_facet`` row. A relation
# counts once for each type at its ends (what choosing the type keeps) and once for its
# pair of types in either direction ("chemical_entity|protein": what choosing both
# keeps, "protein|protein" for a type chosen twice).
_RELATION_TYPE = ("relation_categories", "predicates")
_PARTICIPANT = ("interaction_types", "partner_types")
RELATION_FACETS = {
    "category": (_RELATION_TYPE, "category", "NULL"),
    "predicate": (_RELATION_TYPE, "predicate", "category"),
    "participant_type": (
        _PARTICIPANT,
        "unnest(list_distinct([subject_type, object_type]))",
        "NULL",
    ),
    "participant_pair": (
        _PARTICIPANT,
        "CASE WHEN subject_type <> '' AND object_type <> '' "
        "THEN array_to_string(list_sort([subject_type, object_type]), '|') END",
        "NULL",
    ),
    "source": (("sources",), "unnest(list_distinct(sources))", "NULL"),
    "taxonomy_id": (("taxonomy_ids",), "taxon", "NULL"),
    "sign": (("signs", "sign"), "sign::VARCHAR", "NULL"),
    **{
        term: ((term,), f"unnest(list_distinct({term}))", "NULL")
        for term in RELATION_QUALIFIER_FILTERS
    },
}


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


def _without(payload: dict[str, Any], keys: tuple[str, ...]) -> dict[str, Any]:
    """The payload without a selection, at the top level and under ``filters``."""
    out = {key: value for key, value in payload.items() if key not in keys or not value}
    nested = out.get("filters")
    if isinstance(nested, dict) and any(nested.get(key) for key in keys):
        out["filters"] = {key: value for key, value in nested.items() if key not in keys}
    return out
