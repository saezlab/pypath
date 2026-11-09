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
from omnipath_api.queries.store.connection import sql_literal

logger = logging.getLogger(__name__)

# Relation pages are ordered by these columns, then by resource and position.
RELATION_ORDER = "category, predicate, subject_label, object_label, resource, relation_id"
# Endpoint selections up to this many relations keep their page order between pages.
ORDER_CACHE_LIMIT = 100_000


class RelationsQueries:
    """Relations queries over the release database and the shaping contract."""

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
        resources, where, params, form, pairs = self._selection_parts(
            filters, resources, qualifiers
        )
        if pairs is None:
            sql = f"SELECT * FROM {self._table('relation', resources)} WHERE {where}"
            if form:
                sql += (
                    " AND (resource, relation_id) IN (SELECT resource, relation_id FROM "
                    f"{self._table('relation_evidence', resources)} WHERE {form[0]})"
                )
                params += form[1]
            return sql, params
        return self._lookup_sql("relation", "relation_id", pairs, where, params)

    def _selection_parts(self, filters, resources, qualifiers):
        """The resources, relation conditions (and parameters), molecular form condition
        and, under endpoint filters, the candidate (resource, relation id) pairs (else
        None) of normalized ``filters``."""
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
            return resources, where, params, form, None
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
        return resources, where, params, form, pairs

    def _entity_relation_sql(self, filters, resources=None, *, qualifiers=RELATION_QUALIFIER_FILTERS):
        """SQL and parameters of the relation rows (one per relation and resource, with
        the page-order, filter and facet columns, ``resource`` and ``_rid``) of a plain
        entity selection (entity tokens, or reference keys), from the database's endpoint
        table sorted by key: the selection is one range of it. None for other selections
        (both kinds, a scope, molecular forms) or without such a table."""
        filters = normalize_filters(filters)
        tokens = filters["entity_ids"] + filters["entity_pks"]
        references = filters["reference_entity_keys"]
        if (
            not (tokens or references)
            or (tokens and references)
            or filters["scope_entity_ids"]
            or has_form_filters(filters)
        ):
            return None
        if not resources and filters["sources"]:
            # As in _selection_parts: a source filter reads only those resources.
            queryable = {info["resource"] for info in self._selected_resource_infos(None)}
            resources = [source for source in filters["sources"] if source in queryable] or None
        # The rewritten endpoint table (database formats 7 and 8), or entity_relation
        # (formats 2 to 6, entity tokens only).
        columns = self._served().columns
        endpoint = columns.get("relation_endpoint", ())
        if "page_rank" in endpoint or "subject_label" in endpoint:
            table, key, kind = "relation_endpoint", "key", "reference" if references else "entity"
        elif "entity_relation" in columns and not references:
            table, key, kind = "entity_relation", "entity_key", None
        else:
            return None
        source, in_resources = self._derived(table, resources)
        keys = references if references else self._entity_keys(tokens, resources)
        if not keys:
            return f"SELECT * FROM {source} WHERE FALSE", []
        where, params = self._relation_where(filters, qualifiers=qualifiers)
        # A reference key selects relations at either end, as before.
        mode = "any" if references else filters["gene_mode"]
        side = "TRUE"
        if mode in {"source", "source_only"}:
            side = "side = 'subject'"
        elif mode in {"target", "target_only"}:
            side = "side = 'object'"
        elif mode == "both":
            side = (
                "_rid IN (SELECT _rid FROM er GROUP BY _rid "
                "HAVING bool_or(side = 'subject') AND bool_or(side = 'object'))"
            )
        literals = ", ".join(sql_literal(k) for k in dict.fromkeys(keys))
        kind_condition = f" AND key_kind = '{kind}'" if kind else ""
        # The keys alone first (the table is sorted by them); a relation with both ends
        # among them is one relation.
        return (
            f"""WITH er AS MATERIALIZED (SELECT * FROM {source} WHERE {key} IN ({literals}){kind_condition})
            SELECT DISTINCT ON (_rid) * FROM er WHERE {in_resources} AND {side} AND {where}""",
            params,
        )

    def _entity_keys(self, tokens, resources=None):
        """resolve_entity_keys, kept with the page orders: every page of a gene uses them."""
        return self._cached_facets(
            "entity-keys",
            {"tokens": list(tokens)},
            resources,
            lambda payload, resources: self.resolve_entity_keys(payload["tokens"], resources),
            cache=self._order_cache,
        )

    def _relation_order(self, filters, resources=None):
        """The (resource, relation id) pairs of an endpoint selection in page order, or
        None without endpoint filters or above ORDER_CACHE_LIMIT relations.

        Kept between requests: paging through a gene's relations finds and sorts them
        once, and every page has the total.
        """

        def compute(payload, resources):
            fast = self._entity_relation_sql(filters, resources)
            if fast is not None:
                sql, values = fast
                # The endpoint table's place of each relation in page order (format 8).
                order = (
                    "page_rank"
                    if "page_rank" in self._served().columns.get("relation_endpoint", ())
                    else RELATION_ORDER
                )
                rows = self._db.execute(
                    f"SELECT resource, relation_id FROM ({sql}) ORDER BY {order}", values
                ).fetchall()
                return None if len(rows) > ORDER_CACHE_LIMIT else [list(row) for row in rows]
            resources, where, params, _, pairs = self._selection_parts(
                filters, resources, RELATION_QUALIFIER_FILTERS
            )
            if pairs is None or len(pairs) > ORDER_CACHE_LIMIT:
                return None
            sql, values = self._lookup_sql("relation", "relation_id", pairs, where, params)
            return [
                [r["resource"], r["relation_id"]]
                for r in self._fetch_dicts(
                    f"SELECT resource, relation_id FROM ({sql}) ORDER BY {RELATION_ORDER}", values
                )
            ]

        return self._cached_facets(
            "relation-order", {"filters": filters}, resources, compute, cache=self._order_cache
        )

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
        order = self._relation_order(filters, resources)
        if order is not None:
            pairs = [tuple(pair) for pair in order[int(offset) : int(offset) + int(limit)]]
            total = len(order)
        else:
            selection, params = self._relation_selection(filters, resources)
            # The page is ordered on its sort columns only; then its rows are read in full.
            page = self._fetch_dicts(
                f"SELECT resource, relation_id FROM ({selection}) ORDER BY {RELATION_ORDER} LIMIT ? OFFSET ?",
                [*params, int(limit), int(offset)],
            )
            pairs = [(r["resource"], r["relation_id"]) for r in page]
            total = None
        found = {
            (r["resource"], r["relation_id"]): r
            for r in self._lookup("relation", "relation_id", pairs)
        }
        rows = [found[pair] for pair in pairs]
        if total is None:
            if rows and (offset or len(rows) == limit) or not rows and offset:
                total = self._fetch_dicts(f"SELECT count(*) AS n FROM ({selection})", params)
                total = total[0]["n"]
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
        found = self._relations_by_keys([relation_pk], resources)[:1]
        rows = self._lookup(
            "relation", "relation_id", [(r["resource"], r["relation_id"]) for r in found]
        )
        return self._relation_items(rows)[0] if rows else None
