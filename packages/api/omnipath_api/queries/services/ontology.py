"""Ontology query service, including domain SQL and result shaping."""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Any
from omnipath_api.queries.store.connection import sql_literal


logger = logging.getLogger(__name__)


class OntologyQueries:
    """Ontology queries over the release database and the shaping contract."""

    def _hierarchy_read(self, rel_read: str) -> str:
        """Normalize only hierarchical predicates into child-to-parent edges."""
        from omnipath_core.biolink import hierarchy_direction

        context = "statement_kind = 'ontology'"
        cases = []
        allowed = []
        for row in self._fetch_dicts(f"SELECT DISTINCT predicate FROM {rel_read} WHERE {context}"):
            predicate = row["predicate"]
            reverse = hierarchy_direction(predicate)
            if reverse is not None:
                literal = sql_literal(predicate)
                allowed.append(literal)
                if reverse:
                    cases.append(literal)
        reversed_sql = "predicate IN (" + (", ".join(cases) or "NULL") + ")"
        allowed_sql = "predicate IN (" + (", ".join(allowed) or "NULL") + ")"
        return f"""(SELECT DISTINCT
            CASE WHEN ({reversed_sql}) THEN object_entity_key ELSE subject_entity_key END AS subject_entity_key,
            CASE WHEN ({reversed_sql}) THEN object_label ELSE subject_label END AS subject_label,
            predicate,
            CASE WHEN ({reversed_sql}) THEN subject_entity_key ELSE object_entity_key END AS object_entity_key,
            CASE WHEN ({reversed_sql}) THEN subject_label ELSE object_label END AS object_label
            FROM {rel_read} WHERE {context} AND ({allowed_sql}))"""

    def _hierarchy_edges(self, resources=None) -> str:
        """Child-to-parent edges of the ontologies: the ``ontology_edge`` table of the
        database, else from the relations."""
        derived = self._derived("ontology_edge", resources)
        if derived is None:
            return self._hierarchy_read(self._table("relation", resources))
        source, in_resources = derived
        return f"""(SELECT DISTINCT child_key AS subject_entity_key, child_label AS subject_label,
            predicate, parent_key AS object_entity_key, parent_label AS object_label
            FROM {source} WHERE {in_resources})"""

    def _ontology_seed(self, ent_read, term_id, ontology_id=None, resources=None):
        namespace = ontology_id
        if "|" in term_id:
            qualified_namespace, term_id = term_id.split("|", 1)
            if namespace and namespace.casefold() != qualified_namespace.casefold():
                raise ValueError("Term namespace contradicts ontologyId")
            namespace = qualified_namespace
        if self._has_term_keys():
            rows = self._seed_rows(term_id, namespace, resources)
        else:
            where = "(entity_key = ? OR identifier = ? OR label = ?)"
            params = [term_id, term_id, term_id]
            if namespace:
                where += " AND lower(namespace) = lower(?)"
                params.append(namespace)
            rows = self._fetch_dicts(
                f"SELECT DISTINCT entity_key, label, identifier, namespace FROM {ent_read} WHERE {where} ORDER BY entity_key",
                params,
            )
        keys = {row["entity_key"] for row in rows}
        if len(keys) > 1:
            raise ValueError(
                "Ambiguous ontology term; use ontologyId and an exact identifier or entity key"
            )
        return rows[0] if rows else None

    def _seed_rows(self, term_id, namespace, resources):
        """The entities whose key, identifier or label is ``term_id`` (in ``namespace``),
        found among the lowercase terms (sorted) and read by id or by key, instead of
        a scan of every entity row."""
        hits = self._fetch_dicts(
            f"SELECT DISTINCT resource, entity_id FROM {self._table('entity_term', resources)} WHERE term = ?",
            [term_id.lower()],
        )
        fields = "entity_id, entity_key, label, identifier, namespace"
        rows = self._lookup(
            "entity", "entity_id", [(h["resource"], h["entity_id"]) for h in hits], fields, narrow=True
        )
        rows += self._fetch_dicts(
            f"SELECT {fields}, resource FROM {self._table('entity', resources)} WHERE entity_key = ?",
            [term_id],
        )
        found = {
            (r["entity_key"], r["label"], r["identifier"], r["namespace"])
            for r in rows
            if term_id in (r["entity_key"], r["identifier"], r["label"])
            and (not namespace or str(r["namespace"] or "").lower() == namespace.lower())
        }
        return [
            dict(entity_key=k, label=l, identifier=i, namespace=n)
            for k, l, i, n in sorted(found, key=lambda row: tuple("" if v is None else v for v in row))
        ]

    def _in_namespace(self, keys, namespace, resources):
        """The ``keys`` with an entity row in ``namespace``, read by key."""
        if not keys:
            return set()
        literals = ", ".join(sql_literal(k) for k in sorted(keys))
        return {
            r["entity_key"]
            for r in self._fetch_dicts(
                f"SELECT DISTINCT entity_key FROM {self._table('entity', resources)} "
                f"WHERE entity_key IN ({literals}) AND namespace = ?",
                [namespace],
            )
        }

    def _edge_range(self, column, keys, resources):
        """Distinct hierarchy edges with ``column`` (child_key or parent_key) among ``keys``:
        a few ranges of ``ontology_edge`` (sorted by child) for children."""
        source, in_resources = self._derived("ontology_edge", resources)
        literals = ", ".join(sql_literal(k) for k in sorted(keys))
        return self._fetch_dicts(
            f"""SELECT DISTINCT child_key, child_label, predicate, parent_key, parent_label
            FROM {source} WHERE {column} IN ({literals}) AND {in_resources}"""
        )

    def _tree_edges(self, seed_k, namespace, resources, *, ancestors=True):
        """(the seed's child rows, its ancestor edges) within one namespace, read level by
        level from ``ontology_edge``; None without that table. The same rows as the
        recursive query over the scoped hierarchy, which rescanned it at each level."""
        if self._derived("ontology_edge", resources) is None:
            return None
        # Both ends of an edge in the namespace (the seed's own rows may be in another).
        known = self._in_namespace({seed_k}, namespace, resources)
        children = self._edge_range("parent_key", [seed_k], resources) if known else []
        known |= self._in_namespace({e["child_key"] for e in children}, namespace, resources)
        child_rows = sorted(
            (
                dict(subject_entity_key=e["child_key"], subject_label=e["child_label"])
                for e in children
                if e["child_key"] in known
            ),
            key=lambda r: r["subject_entity_key"],
        )
        edges, seen, frontier = set(), {seed_k}, {seed_k}
        for _ in range(10 if ancestors else 0):
            frontier &= known
            if not frontier:
                break
            found = self._edge_range("child_key", frontier, resources)
            known |= self._in_namespace({e["parent_key"] for e in found} - known, namespace, resources)
            found = [e for e in found if e["parent_key"] in known]
            edges |= {(e["parent_key"], e["parent_label"], e["child_key"], e["child_label"]) for e in found}
            frontier = {e["parent_key"] for e in found} - seen
            seen |= frontier
        ancestor_edges = [
            dict(parent_key=p, parent_label=pl, child_key=c, child_label=cl) for p, pl, c, cl in edges
        ]
        return child_rows, ancestor_edges

    def _scoped_hierarchy(self, ent_read, rel_read, namespace):
        # Keep both endpoints in the resolved seed namespace; cross-ontology links
        # are associations, not children in this namespace's tree.
        literal = sql_literal(namespace)
        keys = f"SELECT entity_key FROM {ent_read} WHERE namespace = {literal}"
        return f"(SELECT * FROM {rel_read} WHERE subject_entity_key IN ({keys}) AND object_entity_key IN ({keys}))"

    def _tree_edges_sql(self, seed_k, rel_read):
        """_tree_edges by a recursive query over a hierarchy read."""
        # 1. Fetch direct children of seed
        child_rows = self._fetch_dicts(
            f"""
            SELECT subject_entity_key, subject_label
            FROM {rel_read}
            WHERE object_entity_key = ?
            ORDER BY subject_entity_key
            LIMIT 101
            """,
            [seed_k],
        )

        # 2. Fetch all ancestor edges (walking parents recursively)
        ancestor_edges = self._fetch_dicts(
            f"""
            WITH RECURSIVE ancestors AS (
                SELECT
                    subject_entity_key AS child_key,
                    subject_label AS child_label,
                    predicate,
                    object_entity_key AS parent_key,
                    object_label AS parent_label,
                    1 AS depth
                FROM {rel_read}
                WHERE subject_entity_key = ?


                UNION ALL

                SELECT
                    r.subject_entity_key AS child_key,
                    r.subject_label AS child_label,
                    r.predicate,
                    r.object_entity_key AS parent_key,
                    r.object_label AS parent_label,
                    a.depth + 1 AS depth
                FROM {rel_read} r
                JOIN ancestors a ON r.subject_entity_key = a.parent_key
                WHERE a.depth < 10
            )
            SELECT DISTINCT parent_key, parent_label, child_key, child_label
            FROM ancestors;
            """,
            [seed_k],
        )

        return child_rows, ancestor_edges

    def get_ontology_tree(
        self,
        term_ids: list[str],
        ontology_id: str | None = None,
        resources: list[str] | None = None,
    ) -> dict[str, Any]:
        """Build ontology tree hierarchy node with all ancestor paths down to the seed term."""
        if not term_ids:
            return {"root": None}
        seed_id = str(term_ids[0]).strip()
        if not self._selected_resource_infos(resources):
            return {"root": {"id": seed_id, "name": seed_id, "children": [], "distance": 0}}

        ent_read = self._table("entity", resources)
        rel_read = self._hierarchy_edges(resources)

        seed_ent = self._ontology_seed(ent_read, seed_id, ontology_id, resources)
        if seed_ent is None:
            return {"root": {"id": seed_id, "name": seed_id, "children": [], "distance": 0}}
        seed_k = seed_ent["entity_key"]
        seed_public_id = seed_ent.get("identifier") or seed_ent.get("label") or seed_id
        found = self._tree_edges(seed_k, seed_ent["namespace"], resources)
        if found is not None:
            child_rows, ancestor_edges = found
        else:
            child_rows, ancestor_edges = self._tree_edges_sql(
                seed_k, self._scoped_hierarchy(ent_read, rel_read, seed_ent["namespace"])
            )
        children_has_more = len(child_rows) > 100
        child_rows = child_rows[:100]

        if not ancestor_edges and not child_rows:
            return {"root": None}

        # 3. Resolve metadata (id, label) for all involved keys
        all_keys = list(
            set(
                [seed_k]
                + [e["parent_key"] for e in ancestor_edges]
                + [e["child_key"] for e in ancestor_edges]
                + [c["subject_entity_key"] for c in child_rows]
            )
        )
        ent_meta: dict[str, dict[str, str]] = {}
        if all_keys:
            placeholders = ", ".join("?" for _ in all_keys)
            for r in self._fetch_dicts(
                f"SELECT entity_key, label, identifier FROM {ent_read} WHERE entity_key IN ({placeholders})",
                all_keys,
            ):
                ident = str(r.get("identifier") or r.get("label") or r.get("entity_key"))
                lbl = str(r.get("label") or ident)
                ent_meta[str(r.get("entity_key"))] = {"id": ident, "name": lbl}

        # 4. Build adjacency graph (parent_key -> child_keys)
        parent_to_children: dict[str, set[str]] = {}
        has_parents: set[str] = set()
        for edge in ancestor_edges:
            p_key = str(edge["parent_key"])
            c_key = str(edge["child_key"])
            parent_to_children.setdefault(p_key, set()).add(c_key)
            has_parents.add(c_key)

        # Attach children under seed node
        for c in child_rows:
            c_key = str(c["subject_entity_key"])
            parent_to_children.setdefault(seed_k, set()).add(c_key)

        roots = [p_k for p_k in parent_to_children if p_k not in has_parents]
        if not roots:
            roots = [seed_k]

        def build_node(key: str, visited: set[str]) -> dict[str, Any]:
            meta = ent_meta.get(key, {"id": key, "name": key})
            node_id = meta["id"]
            node_name = meta["name"]
            children = []
            if key not in visited:
                next_visited = visited | {key}
                for child_k in sorted(
                    parent_to_children.get(key, set()),
                    key=lambda k: ent_meta.get(k, {}).get("name", ""),
                ):
                    children.append(build_node(child_k, next_visited))
            return {
                "id": node_id,
                "name": node_name,
                "children": children,
                "distance": 0 if node_id == seed_public_id or key == seed_k else 1,
                "hasMoreChildren": children_has_more if key == seed_k else False,
                "nextCursor": 100 if key == seed_k and children_has_more else None,
                "ontologyId": seed_ent["namespace"],
            }

        if len(roots) == 1:
            tree_root = build_node(roots[0], set())
        else:
            tree_root = {
                "id": "root",
                "name": "Available ancestry (display root)",
                "children": [
                    build_node(r, set())
                    for r in sorted(roots, key=lambda k: ent_meta.get(k, {}).get("name", ""))
                ],
                "distance": 1,
            }

        return {"root": tree_root}

    def get_ontology_children(
        self,
        term_id: str,
        ontology_id: str | None = None,
        resources: list[str] | None = None,
        limit: int = 100,
        offset: int = 0,
    ) -> dict[str, Any]:
        """Return an explicit page of unique child entities in one namespace."""
        if not 1 <= limit <= 1000 or offset < 0:
            raise ValueError("Invalid ontology child page")
        empty = {"termId": term_id, "children": [], "total": 0, "nextCursor": None}
        if not self._selected_resource_infos(resources):
            return empty
        ent_read = self._table("entity", resources)
        seed = self._ontology_seed(ent_read, term_id, ontology_id, resources)
        if seed is None:
            return empty
        found = self._tree_edges(seed["entity_key"], seed["namespace"], resources, ancestors=False)
        if found is not None:
            # The children's keys from the edges, then their entity rows by key.
            keys = sorted({row["subject_entity_key"] for row in found[0]})
            grouped = (
                self._fetch_dicts(
                    f"""SELECT entity_key, min(identifier) AS identifier, min(label) AS label,
                        min(namespace) AS namespace FROM {ent_read}
                    WHERE entity_key IN ({", ".join(sql_literal(k) for k in keys)})
                    GROUP BY entity_key ORDER BY entity_key"""
                )
                if keys
                else []
            )
            total, rows = len(grouped), grouped[offset : offset + limit]
        else:
            rel_read = self._scoped_hierarchy(ent_read, self._hierarchy_edges(resources), seed["namespace"])
            children_sql = f"""SELECT e.entity_key, min(e.identifier) AS identifier, min(e.label) AS label, min(e.namespace) AS namespace
                FROM {ent_read} e WHERE e.entity_key IN
                    (SELECT subject_entity_key FROM {rel_read} WHERE object_entity_key = ?)
                GROUP BY e.entity_key"""
            total = self._db.execute(
                f"SELECT count(*) FROM ({children_sql})", [seed["entity_key"]]
            ).fetchone()[0]
            rows = self._fetch_dicts(
                f"{children_sql} ORDER BY entity_key LIMIT ? OFFSET ?",
                [seed["entity_key"], limit, offset],
            )
        children = [
            dict(
                id=row["identifier"],
                termId=row["identifier"],
                name=row["label"] or row["identifier"],
                label=row["label"] or row["identifier"],
                entityPk=row["entity_key"],
                ontologyId=row["namespace"],
            )
            for row in rows
        ]
        return {
            "termId": term_id,
            "children": children,
            "total": total,
            "nextCursor": offset + limit if offset + limit < total else None,
        }

    def _compute_ontology_items(self, payload, resources=None):
        if not self._selected_resource_infos(resources):
            return []
        scope = payload.get("selectionScope")
        scoped = (
            scope is not None
            or payload.get("scoped", False)
            or bool(payload.get("entityPks") or payload.get("termIds"))
        )
        keys = list(map(str, payload.get("entityPks") or []))
        terms = list(map(str, payload.get("termIds") or payload.get("annotationTermIds") or []))
        if scope is not None:
            resolved = self.resolve_selection_scope(scope)
            keys = resolved["entityPks"]
            terms = resolved["ontologyTermIds"]
        if scoped and not keys and not terms:
            return []
        entities = self._table("entity", resources)
        relations = self._table("relation", resources)
        # Each endpoint can carry the ontology term; count unique neighbours and relations.
        predicate = "WHERE entity_key IN (SELECT unnest(?::VARCHAR[]))" if scoped else ""
        derived = self._derived("term_entity", resources)
        if derived is not None:
            # Precomputed from the term's side, for terms only.
            source, in_resources = derived
            # Entities and relations by their numbers, where the table has them: the
            # counts below count distinct integers instead of keys.
            numbered = "entity_n" in self._served().columns.get("term_entity", ())
            counted = "entity_n AS entity_key, relation_n AS relation_key" if numbered else "entity_key, relation_key"
            edges = f"""SELECT term_key, {counted} FROM {source}
                WHERE {in_resources}{" AND entity_key IN (SELECT unnest(?::VARCHAR[]))" if scoped else ""}"""
        else:
            edges = f"""SELECT * FROM (
                SELECT object_entity_key term_key, subject_entity_key entity_key, relation_key FROM {relations}
                UNION ALL
                SELECT subject_entity_key term_key, object_entity_key entity_key, relation_key FROM {relations}
            ) {predicate}"""
        params = [keys] if scoped else []
        where = ""
        if scoped:
            where = "WHERE a.term_key IS NOT NULL OR t.entity_key IN (SELECT unnest(?::VARCHAR[])) OR t.identifier IN (SELECT unnest(?::VARCHAR[]))"
            params.extend([keys + terms, terms])
        rows = self._fetch_dicts(
            f"""
            WITH terms AS (
                SELECT entity_key, min(identifier) AS identifier, min(label) AS label, min(namespace) AS namespace
                FROM {entities} WHERE entity_type = 'ontology_class' GROUP BY entity_key
            ), edges AS ({edges}), counts AS (
                SELECT e.term_key, count(DISTINCT e.entity_key) entities, count(DISTINCT e.relation_key) relations
                FROM edges e JOIN terms t ON t.entity_key = e.term_key GROUP BY e.term_key
            )
            SELECT t.*, coalesce(a.entities, 0) entity_count, coalesce(a.relations, 0) relation_count
            FROM terms t LEFT JOIN counts a ON a.term_key = t.entity_key {where}
            ORDER BY entity_count DESC, t.identifier, t.entity_key
        """,
            params,
        )
        return [
            dict(
                entityPk=r["entity_key"],
                termId=r["identifier"] or r["entity_key"],
                ontologyPrefix=(r["identifier"] or "").split(":")[0]
                if ":" in (r["identifier"] or "")
                else r["namespace"],
                ontologyId=r["namespace"] or "unknown",
                label=r["label"],
                definition=None,
                synonyms=[],
                sources=[],
                annotatedEntityCount=r["entity_count"],
                annotatedRelationCount=r["relation_count"],
                annotatedItemCount=r["entity_count"] + r["relation_count"],
            )
            for r in rows
        ]

    def search_ontology(self, payload=None, resources=None, counts=False):
        payload = dict(payload or {})
        base = {
            k: v
            for k, v in payload.items()
            if k not in ("query", "q", "limit", "offset", "ontologyIds", "prefixes")
        }
        items = self._cached_facets("ontology", base, resources, self._compute_ontology_items)
        query = str(payload.get("query") or payload.get("q") or "").strip().casefold()
        if query:
            items = [
                i
                for i in items
                if query in str(i["label"] or "").casefold() or query in i["termId"].casefold()
            ]
        if counts:
            totals = defaultdict(int)
            for item in items:
                totals[item["ontologyId"]] += 1
            return [
                dict(ontologyId=k, scopedCount=v)
                for k, v in sorted(totals.items(), key=lambda kv: (-kv[1], kv[0]))
            ]
        if payload.get("ontologyIds"):
            items = [i for i in items if i["ontologyId"] in payload["ontologyIds"]]
        if payload.get("prefixes"):
            items = [i for i in items if i["ontologyPrefix"] in payload["prefixes"]]
        offset = max(0, int(payload.get("offset") or 0))
        limit = max(1, min(100, int(payload.get("limit") or 24)))
        return items[offset : offset + limit]

    def get_terms_info(
        self, term_ids: list[str], resources: list[str] | None = None
    ) -> dict[str, Any]:
        """Fetch metadata for ontology terms."""
        terms_map: dict[str, Any] = {}
        if not term_ids:
            return {"terms": {}}
        if not self._selected_resource_infos(resources):
            return {"terms": {t: None for t in term_ids}}
        placeholders = ", ".join("?" for _ in term_ids)
        rows = self._fetch_dicts(
            f"""SELECT resource, entity_id, entity_key, label, identifier, namespace
            FROM {self._table("entity", resources)}
            WHERE identifier IN ({placeholders}) OR entity_key IN ({placeholders})""",
            term_ids + term_ids,
        )
        self._children("entity_annotation", rows, "annotations")
        for r in rows:
            ident = str(r.get("identifier") or "")
            label = str(r.get("label") or ident)
            definition = next(
                (
                    str(a.get("value"))
                    for a in (r.get("annotations") or [])
                    if a.get("value")
                    and (
                        str(a.get("term", ""))
                        in {"description", "biolink:description", "IAO:0000115"}
                        or "0801" in str(a.get("term", ""))
                        or "def" in str(a.get("term", "")).lower()
                        or "comment" in str(a.get("term", "")).lower()
                    )
                ),
                None,
            )
            item = {
                "id": ident or label,
                "termId": ident or label,
                "label": label,
                "name": label,
                "namespace": r.get("namespace"),
                "definition": definition,
            }
            if ident:
                terms_map[ident] = item
            if label:
                terms_map[label] = item

        for t in term_ids:
            if t not in terms_map:
                terms_map[t] = None
        return {"terms": terms_map}
