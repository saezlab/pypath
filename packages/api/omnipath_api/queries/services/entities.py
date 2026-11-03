"""Entities query service: search, lookup by key and merging an entity's resource copies."""

from __future__ import annotations

import hashlib
import logging
import json
import re
import time
from collections import defaultdict
from typing import Any


from omnipath_core.biolink import entity_type as biolink_entity_type

from omnipath_api.queries.filters import normalize_filters, EntitySearchCursor
from omnipath_api.queries.shape import shape_entity_summary, shape_relation_summary


logger = logging.getLogger(__name__)

# Term hits rank exact label, exact identifier, label prefix, identifier prefix;
# label substring matches (the fallback phase) come after all of them.
_TERM_RANK = (
    "CASE WHEN term = ? THEN (CASE WHEN kind = 'label' THEN 1 ELSE 2 END) "
    "ELSE (CASE WHEN kind = 'label' THEN 3 ELSE 4 END) END"
)


def _prefix_range(text):
    """[low, high) bounds of the strings starting with ``text``."""
    return text, text[:-1] + chr(ord(text[-1]) + 1)


def _curie(query):
    """(lowercase namespace, identifier) of a ``ns:id`` or ``ns|id`` query, else None."""
    match = re.fullmatch(r"([A-Za-z][\w.-]*)[:|](.+)", query)
    return (match.group(1).lower(), match.group(2)) if match else None


# Group keys of the grouped explorer: a gene's reference, a structure, or one entity.
_GROUP_PREFIXES = ("gene:", "connectivity:", "entity:")


class EntitiesQueries:
    """Entities queries over ParquetQueries' storage and shaping contract."""

    def _term_hits(self, query, resources=None):
        """SQL and parameters of the prefix hits of a query: resource, entity_id,
        match_rank, term and the namespace an entity must have (NULL for any).

        Terms are lowercase labels and identifiers (``entity_term``). A ``ns:id`` query
        also matches exactly that identifier within that namespace, and an entity by
        its reference key (``entity_group``).
        """
        q = query.strip().lower()
        terms = self._table("entity_term", resources)
        parts = [
            f"SELECT resource, entity_id, {_TERM_RANK} AS match_rank, term, NULL::VARCHAR AS namespace "
            f"FROM {terms} WHERE term >= ? AND term < ?"
        ]
        params = [q, *_prefix_range(q)]
        if curie := _curie(query.strip()):
            namespace, identifier = curie
            parts.append(
                f"SELECT resource, entity_id, {_TERM_RANK}, term, ? FROM {terms} WHERE term = ?"
            )
            params += [identifier.lower(), namespace, identifier.lower()]
            parts.append(
                f"SELECT resource, entity_id, 2, ?, NULL FROM {self._table('entity_group', resources)} "
                "WHERE group_key = ? AND kind = 'reference'"
            )
            params += [q, f"{namespace}:{identifier}"]
        return " UNION ALL ".join(parts), params

    def _search_phases(self, query):
        q = (query or "").strip()
        if not q:
            return ("all",)
        if self._is_entity_key(q):
            return ("key",)
        return ("prefix", "contains")

    def _phase_where(self, query, phase, resources=None):
        """WHERE clause over the entity table selecting a query's matches in one phase."""
        q = (query or "").strip()
        if phase == "all":
            return "TRUE", []
        if phase == "key":
            return "entity_key = ?", [q]
        if phase == "contains":
            return "contains(lower(label), ?)", [q.lower()]
        hits, params = self._term_hits(q, resources)
        return (
            "((resource, entity_id) IN (SELECT resource, entity_id FROM hits WHERE namespace IS NULL) "
            "OR (resource, entity_id, lower(namespace)) IN (SELECT resource, entity_id, namespace FROM hits))"
        ).replace("hits", f"({hits})"), [*params, *params]

    def _entity_match_where(self, query, clauses, params, resources=None):
        """The first search phase with a match among filtered entities, as a WHERE clause.

        Grouped views and facets select the same entities as the ungrouped search.
        """
        where = " AND ".join(clauses) or "TRUE"
        phases = self._search_phases(query)
        for phase in phases:
            text, values = self._phase_where(query, phase, resources)
            if len(phases) == 1 or self._fetch_dicts(
                f"SELECT 1 FROM {self._table('entity', resources)} WHERE ({where}) AND {text} LIMIT 1",
                [*params, *values],
            ):
                return f"({where}) AND {text}", [*params, *values]
        return "FALSE", []

    def resolve_entity_keys(
        self,
        tokens: list[str] | None,
        resources: list[str] | None = None,
    ) -> list[str]:
        """Map labels, identifiers, ``ns:id`` / ``ns|id`` CURIEs, reference keys, entity
        keys or the explorer's group keys (``gene:entrez:7157``, ``connectivity:<block>``,
        ``entity:<key>``) to entity keys, by exact (case-insensitive) match."""
        values = [str(t).strip() for t in (tokens or []) if str(t).strip()]
        if not values:
            return []
        groups = [v for v in values if v.startswith(_GROUP_PREFIXES)]
        keys, other = self._split_entity_tokens([v for v in values if v not in groups])
        found = list(keys) + self._group_members(groups, resources)
        if other:
            terms = self._table("entity_term", resources)
            # Exact keys (files are sorted by key), lowercase terms, and CURIEs below.
            parts = [
                f"SELECT resource, entity_id, NULL::VARCHAR AS namespace FROM {self._table('entity', resources)} "
                f"WHERE entity_key IN ({','.join('?' for _ in other)})",
                f"SELECT resource, entity_id, NULL FROM {terms} "
                "WHERE term IN (SELECT unnest(?::VARCHAR[]))",
            ]
            params: list[Any] = [*other, [v.lower() for v in other]]
            curies = [c for v in other if (c := _curie(v))]
            if curies:
                parts.append(
                    f"SELECT resource, entity_id, namespace FROM {terms} JOIN "
                    "(SELECT unnest(?::VARCHAR[]) AS namespace, unnest(?::VARCHAR[]) AS term) USING (term)"
                )
                params += [[ns for ns, _ in curies], [i.lower() for _, i in curies]]
                parts.append(
                    f"SELECT resource, entity_id, NULL FROM {self._table('entity_group', resources)} "
                    "WHERE kind = 'reference' AND group_key IN (SELECT unnest(?::VARCHAR[]))"
                )
                params.append([f"{ns}:{i}" for ns, i in curies])
            hits = self._fetch_dicts(" UNION ALL ".join(parts), params)
            rows = self._lookup(
                "entity",
                "entity_id",
                [(h["resource"], h["entity_id"]) for h in hits],
                "entity_id, entity_key, namespace",
            )
            namespaces = defaultdict(set)
            for hit in hits:
                namespaces[(hit["resource"], hit["entity_id"])].add(hit["namespace"])
            found += sorted(
                row["entity_key"]
                for row in rows
                if namespaces[(row["resource"], row["entity_id"])]
                & {None, str(row["namespace"] or "").lower()}
            )
        return list(dict.fromkeys(found))

    def _group_members(self, groups: list[str], resources=None) -> list[str]:
        """Entity keys of group keys as the grouped explorer forms them."""
        found = [g.removeprefix("entity:") for g in groups if g.startswith("entity:")]
        wanted = [
            ("reference" if g.startswith("gene:") else "connectivity", g.split(":", 1)[1])
            for g in groups
            if not g.startswith("entity:")
        ]
        if wanted:
            hits = self._fetch_dicts(
                f"""SELECT g.resource, g.entity_id FROM {self._table("entity_group", resources)} g
                JOIN (SELECT unnest(?::VARCHAR[]) AS kind, unnest(?::VARCHAR[]) AS group_key) w
                USING (kind, group_key)""",
                [[kind for kind, _ in wanted], [key for _, key in wanted]],
            )
            rows = self._lookup(
                "entity", "entity_id", [(h["resource"], h["entity_id"]) for h in hits], "entity_key"
            )
            found += sorted({row["entity_key"] for row in rows})
        return found

    def _merge_by_key(self, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """One merged row per entity key, in the order keys first appear."""
        groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            groups[str(row["entity_key"])].append(row)
        return [self._merge_entity_row_group(members) for members in groups.values()]

    def _merge_entity_row_group(self, members: list[dict[str, Any]]) -> dict[str, Any]:
        """One entity from its resource copies: the richest copy's scalars, the union of
        their identifiers, annotations and evidence, and the keys and resources merged."""

        def count(row, field):
            nested = row.get(field)
            return (
                len(nested)
                if isinstance(nested, list)
                else int(row.get(f"{field[:-1]}_count") or 0)
            )

        def rank(row: dict[str, Any]) -> tuple[Any, ...]:
            label = str(row.get("label") or "")
            ident = str(row.get("identifier") or "")
            gene_like = int(bool(re.fullmatch(r"[A-Z][A-Z0-9-]{1,14}", label)))
            return (
                count(row, "annotations"),
                count(row, "identifiers"),
                gene_like,
                int(bool(label) and label != ident and not label.isdigit()),
                # Equal-ranked resource copies must not depend on Parquet scan order.
                label.casefold(),
                label,
                str(row.get("entity_key") or ""),
                str(row.get("resource") or ""),
                bool(row.get("has_hierarchy")),
                int(row.get("parent_count") or 0),
                int(row.get("child_count") or 0),
            )

        merged = dict(max(members, key=rank))
        merged["relation_count"] = sum(int(row.get("relation_count") or 0) for row in members)
        merged["_resources"] = list(
            dict.fromkeys(row["resource"] for row in members if row.get("resource"))
        )
        merged["_source_entity_keys"] = list(
            dict.fromkeys(str(row["entity_key"]) for row in members if row.get("entity_key"))
        )
        merged["gene_reference_keys"] = sorted(
            {ref for row in members for ref in row.get("gene_reference_keys") or []}
        )
        # Reference selection is independent of display-label ranking. A native
        # fallback in any source copy must survive catalogue gene enrichment.
        references = {
            row.get("reference_entity_key") for row in members if row.get("reference_entity_key")
        }
        native = f"{merged.get('namespace')}:{merged.get('identifier')}"
        merged["reference_entity_key"] = (
            native
            if native in references or len(references) > 1
            else next(iter(references))
            if references
            else None
        )
        if len(members) == 1:
            return merged
        for field, identity in (
            ("identifiers", lambda item: (item.get("ns"), item.get("id"))),
            ("annotations", lambda item: json.dumps(item, sort_keys=True, default=str)),
            ("evidence", lambda item: json.dumps(item, sort_keys=True, default=str)),
        ):
            if any(field in row for row in members):
                items = {}
                for row in members:
                    for item in row.get(field) or []:
                        items.setdefault(identity(item), item)
                merged[field] = list(items.values())
        return merged

    def _taxonomy_names(self) -> dict[str, str]:
        scoped = self._taxonomy_scope.get()
        return scoped if scoped is not None else self._resolve_taxonomy_names()

    def _taxonomy_cache_version(self):
        from omnipath_core.taxonomy_lookup import lookup_path

        return lookup_path(self.data_root) if self._release_scope.get() is None else None

    def _resolve_taxonomy_names(self) -> dict[str, str]:
        from omnipath_core.taxonomy import read_reference

        manifest = self._release_scope.get()
        if manifest is None:
            from omnipath_core.taxonomy_lookup import lookup_path, read_lookup

            if path := lookup_path(self.data_root):
                return read_lookup(path)
            # "latest" is a moving view; named releases never inherit newer labels.
            manifest = next(
                (m for m in self.releases.list() if m.get("references", {}).get("taxonomy")), {}
            )
        reference = manifest.get("references", {}).get("taxonomy")
        if not reference:
            return {}
        path = self.data_root / "references/taxonomy" / reference["version"] / "taxonomy.parquet"
        return read_reference(str(path))

    def _label_taxonomy_facets(self, items: list[dict]) -> list[dict]:
        names = self._taxonomy_names()
        return [
            {**item, "facetLabel": names.get(item["facetValue"])}
            if item["facetName"] == "taxonomy_id"
            else dict(item)
            for item in items
        ]

    def _to_entity_summary(
        self, row: dict[str, Any], resource_hint: str | None = None
    ) -> dict[str, Any]:
        result = shape_entity_summary(row, resource_hint=resource_hint)
        result["taxonomyName"] = self._taxonomy_names().get(result.get("taxonomyId"))
        return result

    def _to_entity_relation(self, row: dict[str, Any]) -> dict[str, Any]:
        if "annotations" not in row:
            self._with_qualifiers(row)
        return shape_relation_summary(row)

    def _endpoint_summary(self, row: dict[str, Any], side: str, endpoints) -> dict[str, Any]:
        """A relation endpoint's entity from its resource, or from the relation's columns."""
        found = endpoints.get((row.get("resource"), row.get(f"{side}_entity_key")))
        if found is not None:
            return self._to_entity_summary(found)
        return self._to_entity_summary(
            {
                "entity_key": row.get(f"{side}_entity_key"),
                "namespace": row.get(f"{side}_type") or "unknown",
                "identifier": row.get(f"{side}_label") or row.get(f"{side}_entity_key"),
                "label": row.get(f"{side}_label"),
                "entity_type": row.get(f"{side}_type"),
                "taxon": row.get("taxon"),
                "reference_entity_key": row.get(f"{side}_reference_entity_key"),
            }
        )

    def _relation_items(self, rows: list[dict[str, Any]], **extra) -> list[dict[str, Any]]:
        """Relation rows with their endpoint entities, as the API returns them.

        ``extra`` maps response fields to functions of the relation row.
        """
        endpoints = self._endpoint_rows(rows)
        return [
            {
                "relation": self._to_entity_relation(row),
                "subjectEntity": self._endpoint_summary(row, "subject", endpoints),
                "objectEntity": self._endpoint_summary(row, "object", endpoints),
                **{name: value(row) for name, value in extra.items()},
            }
            for row in rows
        ]

    def _entity_filter_clauses(
        self,
        filters: dict[str, Any] | None = None,
        *,
        resources: list[str] | None = None,
    ) -> tuple[list[str], list[Any]]:
        filters = normalize_filters(filters)
        clauses: list[str] = []
        params: list[Any] = []
        entity_types = self._as_list(filters.get("entity_types"))
        if entity_types:
            norm_types = [biolink_entity_type(t) for t in entity_types]
            placeholders = ", ".join("?" for _ in norm_types)
            clauses.append(f"entity_type IN ({placeholders})")
            params.extend(norm_types)
        taxons = self._as_list(filters.get("ncbi_tax_id") or filters.get("taxonomy_ids"))
        if taxons:
            placeholders = ", ".join("?" for _ in taxons)
            clauses.append(f"taxon IN ({placeholders})")
            params.extend(str(t) for t in taxons)
        if filters["reference_entity_keys"]:
            clauses.append("reference_entity_key IN (SELECT unnest(?::VARCHAR[]))")
            params.append(filters["reference_entity_keys"])
        tokens = filters["entity_ids"] + filters["entity_pks"]
        if tokens:
            keys = self.resolve_entity_keys(tokens, resources or filters["sources"] or None)
            clauses.append("entity_key IN (SELECT unnest(?::VARCHAR[]))")
            params.append(keys)
        return clauses, params

    def search_entity_groups(self, *, strategy="chemical_connectivity", **kwargs):
        from omnipath_api.queries.services.groups import search_groups

        if strategy not in {"auto", "chemical_connectivity", "gene_reference"}:
            raise ValueError(f"Unknown grouping strategy: {strategy}")
        return search_groups(self, strategy=strategy, **kwargs)

    def search_connectivity_groups(self, **kwargs):
        return self.search_entity_groups(strategy="chemical_connectivity", **kwargs)

    def get_molecular_context(self, entity_id, resources=None, **kwargs):
        from omnipath_api.queries.services.molecular_context import context

        return context(self, entity_id, resources, **kwargs)

    def get_entity_examples(self):
        from omnipath_api.queries.examples import get_examples

        return get_examples(self)

    def _search_page(self, query, phase, where, params, resources, count):
        """Entity rows of the first ``count`` matching entity keys, grouped by key in
        match order: rank, then shorter terms (or labels), then more relations, then key."""
        q = (query or "").strip()
        if phase == "prefix":
            return self._prefix_page(q, where, params, resources, count)
        text, values = self._phase_where(q, phase, resources)
        order = (
            "length(sort_label), relations DESC, sort_label, entity_key"
            if phase == "contains"
            else "entity_key"
        )
        return self._fetch_dicts(
            f"""WITH matched AS (
                SELECT * FROM {self._table("entity", resources)} WHERE ({where}) AND {text}
            ), page AS (
                SELECT entity_key, row_number() OVER (ORDER BY {order}) AS position
                FROM (SELECT entity_key, coalesce(min(label), '') AS sort_label,
                    sum(relation_count) AS relations FROM matched GROUP BY entity_key)
                ORDER BY position LIMIT ?
            ) SELECT matched.* FROM matched JOIN page USING (entity_key)
            ORDER BY page.position, matched.resource""",
            [*params, *values, count],
        )

    def _prefix_page(self, query, where, params, resources, count):
        """Term hits in match order, read per resource until ``count`` keys pass the
        filters. Hits are cheap (``entity_term`` is sorted by term); entity rows are read
        only for the hits needed, at most doubling the number read each round."""
        hits, values = self._term_hits(query, resources)
        ordered = f"""SELECT resource, entity_id, min(match_rank) AS match_rank,
              min(length(term)) AS term_length, min(term) AS term,
              list(DISTINCT namespace) FILTER (WHERE namespace IS NOT NULL) AS namespaces,
              bool_or(namespace IS NULL) AS any_namespace
            FROM ({hits}) GROUP BY resource, entity_id
            ORDER BY match_rank, term_length, term, resource, entity_id LIMIT ?"""
        limit = max(4 * count, 64)
        while True:
            found = self._fetch_dicts(ordered, [*values, limit])
            rows = self._lookup(
                "entity",
                "entity_id",
                [(h["resource"], h["entity_id"]) for h in found],
                extra=f"AND ({where})",
                params=params,
            )
            by_hit = {(r["resource"], r["entity_id"]): r for r in rows}
            best, relations, matched = {}, defaultdict(int), defaultdict(list)
            for hit in found:
                row = by_hit.get((hit["resource"], hit["entity_id"]))
                if row is None or not (
                    hit["any_namespace"]
                    or str(row["namespace"] or "").lower() in (hit["namespaces"] or [])
                ):
                    continue
                key = row["entity_key"]
                rank = (hit["match_rank"], hit["term_length"], hit["term"])
                best[key] = min(best.get(key, rank), rank)
                if row not in matched[key]:
                    matched[key].append(row)
                    relations[key] += int(row["relation_count"] or 0)
            # Equally good matches: the more connected entity first.
            keys = sorted(
                best, key=lambda k: (best[k][0], best[k][1], -relations[k], best[k][2], k)
            )
            # Stop once the hits read cover every tie of the last key needed.
            if len(found) < limit or (
                len(keys) >= count
                and (found[-1]["match_rank"], found[-1]["term_length"]) > best[keys[count - 1]][:2]
            ):
                page = keys[:count]
                return self._with_other_copies(page, matched, where, params, resources)
            limit *= 4

    def _with_other_copies(self, keys, matched, where, params, resources):
        """The rows of ``keys`` in key order: the matched rows and every other resource's
        row of the same entity that passes the filters.

        An entity's rows whose own terms rank beyond the hits read would otherwise be
        left out of some pages, so its relation count and label would depend on the page
        size. The other rows are located by key reading few columns (files are sorted
        by key), then read whole by id.
        """
        rows = {key: list(matched[key]) for key in keys}
        if not keys:
            return []
        known = {(r["resource"], r["entity_id"]) for found in rows.values() for r in found}
        condition, values = self._in("entity_key", keys)
        located = self._fetch_dicts(
            f"SELECT resource, entity_id FROM {self._table('entity', resources)} "
            f"WHERE {condition} AND ({where})",
            [*values, *params],
        )
        others = [pair for r in located if (pair := (r["resource"], r["entity_id"])) not in known]
        for row in self._lookup("entity", "entity_id", others):
            rows[row["entity_key"]].append(row)
        return [row for key in keys for row in sorted(rows[key], key=lambda r: r["resource"])]

    def search_entities_api(
        self,
        query: str = "",
        filters: dict[str, Any] | None = None,
        resources: list[str] | None = None,
        limit: int = 20,
        cursor: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        t0 = time.perf_counter()
        filters = normalize_filters(filters)
        scope = resources or filters.get("sources") or None
        infos = self._selected_resource_infos(scope)
        if not infos:
            return {"entities": [], "nextCursor": None, "total": 0, "elapsed_ms": 0.0}
        q = (query or "").strip()
        query_key = hashlib.sha256(
            json.dumps([q, filters, sorted(i["key"] for i in infos)], sort_keys=True).encode()
        ).hexdigest()
        offset = 0
        phases = self._search_phases(q)
        if cursor:
            cursor = EntitySearchCursor.model_validate(cursor).model_dump()
            if cursor["queryKey"] != query_key:
                raise ValueError("Entity cursor belongs to another query")
            offset, phases = cursor["offset"], (cursor["phase"],)
        clauses, params = self._entity_filter_clauses(filters, resources=resources)
        where = " AND ".join(clauses) or "TRUE"
        merged: list[dict[str, Any]] = []
        phase = phases[0]
        for phase in phases:
            rows = self._search_page(q, phase, where, params, scope, offset + int(limit) + 1)
            merged = self._merge_by_key(rows)
            if merged:
                break
        page = merged[offset : offset + int(limit)]
        more = len(merged) > offset + int(limit)
        next_cursor = (
            {"phase": phase, "offset": offset + int(limit), "queryKey": query_key} if more else None
        )
        entities = [self._to_entity_summary(row) for row in page]
        return {
            "entities": entities,
            "nextCursor": next_cursor,
            "total": offset + len(entities) + int(more),
            "elapsed_ms": round((time.perf_counter() - t0) * 1000, 2),
        }

    def _fetch_entities_by_keys(
        self,
        keys: list[str],
        resources: list[str] | None = None,
        *,
        slim: bool = False,
        evidence: bool = False,
    ) -> list[dict[str, Any]]:
        """Merged entities by key; ``slim`` leaves out identifiers and annotations."""
        keys = list(dict.fromkeys(str(k) for k in keys if k))
        if not keys:
            return []
        rows = self._entity_rows(
            "entity_key IN (SELECT unnest(?::VARCHAR[]))",
            [keys],
            resources,
            nested=not slim,
            evidence=evidence,
        )
        merged = {row["entity_key"]: row for row in self._merge_by_key(rows)}
        return [merged[key] for key in keys if key in merged]

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
        keys = self.resolve_entity_keys(values, resources)
        rows = self._fetch_entities_by_keys(keys, resources, slim=slim)
        return {"entities": [self._to_entity_summary(row) for row in rows]}

    def get_entities_by_public_ids(
        self, public_ids: list[str], resources: list[str] | None = None
    ) -> dict[str, Any]:
        return self.get_entities_by_pks(public_ids, resources)

    def get_entity_core(self, public_id, resources=None):
        from omnipath_api.queries.entity_details import core

        return core(self, public_id, resources)

    def get_entity_evidence(self, public_id, resources=None, limit=20, offset=0):
        from omnipath_api.queries.entity_details import evidence

        return evidence(self, public_id, resources, limit, offset)

    def get_entity_relationships(self, public_id, resources=None, limit=50, offset=0):
        from omnipath_api.queries.entity_details import relationships

        return relationships(self, public_id, resources, limit, offset)

    def get_entity_details(
        self, public_id: str, resources: list[str] | None = None
    ) -> dict[str, Any] | None:
        # Details are addressed solely by the exact entity key. Identifier and
        # alias discovery belongs to search, never to single-entity hydration.
        from omnipath_api.queries.entity_details import details

        return details(self, public_id, resources)

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
