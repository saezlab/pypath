"""Relations query service, including domain SQL and result shaping."""

from __future__ import annotations

import logging
import time
from typing import Any


from omnipath_core.biolink import (
    predicate as biolink_predicate,
    entity_type as biolink_entity_type,
    annotation_value,
)

from omnipath_api.queries.molecular import form_match_sql, has_form_filters, matching_evidence
from omnipath_api.queries.filters import normalize_filters


from omnipath_api.queries.services.constants import RELATION_QUALIFIER_FILTERS

logger = logging.getLogger(__name__)

# Relation pages are ordered by these columns, then by resource and position.
RELATION_ORDER = "category, predicate, subject_label, object_label, resource, relation_id"


class RelationsQueries:
    """Relations queries over ParquetQueries' storage and shaping contract."""

    def _relation_where(self, filters, *, qualifiers=RELATION_QUALIFIER_FILTERS):
        """Conditions on the relation table's own columns, and their parameters."""
        where_clauses = ["TRUE"]
        params: list[Any] = []

        categories = self._as_list(filters.get("categories") or filters.get("relation_categories"))
        if categories:
            placeholders = ", ".join("?" for _ in categories)
            where_clauses.append(f"category IN ({placeholders})")
            params.extend(categories)

        predicates = self._as_list(filters.get("predicates"))
        if predicates:
            norm_predicates = [biolink_predicate(p) for p in predicates]
            placeholders = ", ".join("?" for _ in norm_predicates)
            where_clauses.append(f"predicate IN ({placeholders})")
            params.extend(norm_predicates)

        for term in qualifiers:
            values = self._as_list(filters.get(term))
            if values:
                where_clauses.append(f"list_has_any({term}, ?::VARCHAR[])")
                params.append([annotation_value(term, value) for value in values])

        if "is_directed" in filters and filters["is_directed"] is not None:
            where_clauses.append("is_directed = ?")
            params.append(bool(filters["is_directed"]))

        signs = self._as_list(filters.get("signs"))
        if "sign" in filters and filters["sign"] is not None and not signs:
            signs = [filters["sign"]]
        if signs:
            placeholders = ", ".join("?" for _ in signs)
            where_clauses.append(f"sign IN ({placeholders})")
            params.extend(int(s) for s in signs)

        min_ev = filters.get("min_evidence_count")
        if min_ev is not None and int(min_ev) > 1:
            where_clauses.append("evidence_count >= ?")
            params.append(int(min_ev))

        interaction_types = self._as_list(
            filters.get("interaction_types")
            or filters.get("interactionTypes")
            or filters.get("entity_types")
        )
        partner_types = self._as_list(filters.get("partner_types") or filters.get("partnerTypes"))
        if interaction_types:
            norm_types = [biolink_entity_type(t) for t in interaction_types]
            placeholders = ", ".join("?" for _ in norm_types)
            if partner_types:
                # One endpoint of a chosen type and the other of a chosen partner type, in
                # either direction: protein and partner protein is protein-protein.
                partners = [biolink_entity_type(t) for t in partner_types]
                partner_marks = ", ".join("?" for _ in partners)
                where_clauses.append(
                    f"((subject_type IN ({placeholders}) AND object_type IN ({partner_marks}))"
                    f" OR (object_type IN ({placeholders}) AND subject_type IN ({partner_marks})))"
                )
                params.extend([*norm_types, *partners, *norm_types, *partners])
            else:
                where_clauses.append(
                    f"(subject_type IN ({placeholders}) OR object_type IN ({placeholders}))"
                )
                params.extend(norm_types)
                params.extend(norm_types)

        taxons = self._as_list(
            filters.get("ncbi_tax_id") or filters.get("taxonomy_ids") or filters.get("taxonomyIds")
        )
        if taxons:
            placeholders = ", ".join("?" for _ in taxons)
            where_clauses.append(f"taxon IN ({placeholders})")
            params.extend(str(t) for t in taxons)

        sources = self._as_list(filters.get("sources"))
        if sources:
            where_clauses.append("list_has_any(sources, ?::VARCHAR[])")
            params.append([str(s) for s in sources])

        return " AND ".join(where_clauses), params

    def _endpoint_candidates(self, filters, resources=None):
        """{resource: relation ids} of relations matching the endpoint filters (entity
        keys, scope keys and reference keys, each with its endpoint mode), from the
        endpoint table; None without endpoint filters."""
        endpoint = self._table("relation_endpoint", resources)
        sets, params = [], []

        def add(keys, kind, mode):
            side, having = "TRUE", ""
            if mode in {"source", "source_only"}:
                side = "side = 'subject'"
            elif mode in {"target", "target_only"}:
                side = "side = 'object'"
            elif mode == "both":
                having = "HAVING bool_or(side = 'subject') AND bool_or(side = 'object')"
            sets.append(
                f"""SELECT resource, relation_id FROM {endpoint}
                WHERE key IN (SELECT unnest(?::VARCHAR[])) AND key_kind = ? AND {side}
                GROUP BY resource, relation_id {having}"""
            )
            params.extend([keys, kind])

        tokens = filters["entity_ids"] + filters["entity_pks"]
        if tokens:
            add(self.resolve_entity_keys(tokens, resources), "entity", filters["gene_mode"])
        if filters["scope_entity_ids"]:
            keys = self.resolve_entity_keys(filters["scope_entity_ids"], resources)
            add(keys, "entity", filters["scope_endpoint_mode"])
        if filters["reference_entity_keys"]:
            add(filters["reference_entity_keys"], "reference", "any")
        if not sets:
            return None
        rows = self._fetch_dicts(
            f"""SELECT resource, list(relation_id ORDER BY relation_id) AS ids
            FROM ({" INTERSECT ".join(sets)}) GROUP BY resource""",
            params,
        )
        return {row["resource"]: row["ids"] for row in rows}

    def _relation_selection(
        self, filters, resources=None, *, qualifiers=RELATION_QUALIFIER_FILTERS
    ):
        """SQL selecting the relation rows (with ``resource``) that match ``filters``.

        Endpoint filters select candidates from the endpoint table first, then read only
        those relations (and their evidence, for molecular forms) per resource.
        """
        filters = normalize_filters(filters)
        if not resources and filters["sources"]:
            # A relation lists only its own resource among its sources, so a source filter
            # needs the files of those resources only; the WHERE below keeps the filter. A
            # source that is no queryable resource of the release reads every file, as before,
            # and matches nothing.
            queryable = {info["resource"] for info in self._selected_resource_infos(None)}
            resources = [source for source in filters["sources"] if source in queryable] or None
        where, params = self._relation_where(filters, qualifiers=qualifiers)
        form = form_match_sql(filters) if has_form_filters(filters) else None
        candidates = self._endpoint_candidates(filters, resources)
        if candidates is None:
            sql = f"SELECT * FROM {self._table('relation', resources)} WHERE {where}"
            if form:
                sql += (
                    " AND (resource, relation_id) IN (SELECT resource, relation_id FROM "
                    f"{self._table('relation_evidence', resources)} WHERE {form[0]})"
                )
                params += form[1]
            return sql, params
        pairs = [(resource, i) for resource, ids in candidates.items() for i in ids]
        if form:
            # Relations with a matching evidence row, read by relation id per resource.
            sql, values = self._lookup_sql("relation_evidence", "relation_id", pairs, *form)
            pairs = [
                (r["resource"], r["relation_id"])
                for r in self._fetch_dicts(
                    f"SELECT DISTINCT resource, relation_id FROM ({sql})", values
                )
            ]
        return self._lookup_sql("relation", "relation_id", pairs, where, params)

    def search_relations(
        self,
        filters: dict[str, Any] | None = None,
        resources: list[str] | None = None,
        limit: int = 50,
        offset: int = 0,
        include_details: bool = False,
    ) -> dict[str, Any]:
        """A page of matching relations; with details, their annotations and evidence."""
        t0 = time.perf_counter()
        filters = normalize_filters(filters)
        selection, params = self._relation_selection(filters, resources)
        # The page is ordered on its sort columns only; then its rows are read in full.
        page = self._fetch_dicts(
            f"SELECT resource, relation_id FROM ({selection}) ORDER BY {RELATION_ORDER} LIMIT ? OFFSET ?",
            [*params, int(limit), int(offset)],
        )
        pairs = [(r["resource"], r["relation_id"]) for r in page]
        found = {
            (r["resource"], r["relation_id"]): r
            for r in self._lookup("relation", "relation_id", pairs)
        }
        rows = [found[pair] for pair in pairs]
        if rows and (offset or len(rows) == limit) or not rows and offset:
            total = self._fetch_dicts(f"SELECT count(*) AS n FROM ({selection})", params)[0]["n"]
        else:
            total = offset + len(rows)
        if include_details:
            self._children("relation_annotation", rows, "annotations")
            self._children("relation_evidence", rows, "evidence")
            if has_form_filters(filters):
                for row in rows:
                    row["evidence"] = matching_evidence(row["evidence"], filters)
                    row["evidence_count"] = len(row["evidence"])
        return {
            "rows": rows,
            "total": int(total),
            "limit": limit,
            "offset": offset,
            "elapsed_ms": round((time.perf_counter() - t0) * 1000, 2),
            "files_scanned": len(self._selected_resource_infos(resources)),
        }

    def search_relations_api(self, filters=None, resources=None, limit=20, offset=0):
        started = time.perf_counter()
        payload = dict(filters=filters or {}, limit=limit, offset=offset)
        result = self._cached_facets(
            "relation-page",
            payload,
            resources,
            lambda p, r: self._search_relations_api_uncached(resources=r, **p),
        )
        result["elapsed_ms"] = round((time.perf_counter() - started) * 1000, 2)
        return result

    def _search_relations_api_uncached(
        self,
        filters: dict[str, Any] | None = None,
        resources: list[str] | None = None,
        limit: int = 20,
        offset: int = 0,
    ) -> dict[str, Any]:
        t0 = time.perf_counter()
        result = self.search_relations(
            filters=filters,
            resources=resources,
            limit=limit,
            offset=offset,
            include_details=has_form_filters(normalize_filters(filters)),
        )
        rows = self._relation_items(result["rows"])
        next_offset = offset + len(rows)
        return {
            "relations": [row["relation"] for row in rows],
            "rows": rows,
            "total": result["total"],
            "nextCursor": str(next_offset) if next_offset < int(result["total"]) else None,
            "elapsed_ms": round((time.perf_counter() - t0) * 1000, 2),
        }

    def get_relation_by_pk(
        self, relation_pk: str, resources: list[str] | None = None
    ) -> dict[str, Any] | None:
        rows = self._fetch_dicts(
            f"SELECT * FROM {self._table('relation', resources)} WHERE relation_key = ? "
            "ORDER BY resource LIMIT 1",
            [relation_pk],
        )
        return self._relation_items(rows)[0] if rows else None
