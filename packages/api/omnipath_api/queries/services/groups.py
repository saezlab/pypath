"""Grouped entity views: chemicals by structure, gene products by their gene.

A query-time view; canonical identities and relation endpoints stay intact. Chemicals
group by connectivity (the first InChIKey block of their reference), other entities by
a shared Entrez gene reference. An entity whose resource copies disagree, or that has
no group, is a singleton. ``auto`` applies both, the other strategies one of them.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
import json
import time

from omnipath_api.queries.filters import normalize_filters
from omnipath_api.queries.shape.display_name import display_name, preferred_name
from omnipath_core.display_names import name_rank


def _agreed(column):
    """The value all non-null copies share, else NULL (compared by hash)."""
    present = f"FILTER (WHERE {column} IS NOT NULL)"
    return (
        f"CASE WHEN min(hash({column})) {present} = max(hash({column})) {present} "
        f"THEN any_value({column}) END"
    )


def _group_members(queries, scope, groups):
    """(resource, entity_id) of the entities listed under (group_key, kind) pairs."""
    if not groups:
        return []
    rows = queries._fetch_dicts(
        f"""SELECT DISTINCT resource, entity_id FROM {queries._table("entity_group", scope)}
        WHERE (group_key, kind) IN (SELECT unnest(?::VARCHAR[]), unnest(?::VARCHAR[]))""",
        [[key for key, _ in groups], [kind for _, kind in groups]],
    )
    return [(r["resource"], r["entity_id"]) for r in rows]


def _candidate_keys(queries, group_key, scope):
    """Entity keys that may belong to one group: those listed under it."""
    prefix, _, value = group_key.partition(":")
    if prefix == "entity":
        return [value]
    kind = {"connectivity": "connectivity", "gene": "reference"}.get(prefix)
    if kind is None:
        return []
    pairs = _group_members(queries, scope, [(value, kind)])
    return sorted(
        {r["entity_key"] for r in queries._lookup("entity", "entity_id", pairs, "entity_key")}
    )


def _member_rows(queries, page, scope):
    """Every copy in scope of the page's member entities, found by group or key."""
    groups = [
        (group["connectivity"], "connectivity") if group["connectivity"]
        else (group["reference_entity_key"], "reference")
        for group in page
        if group["connectivity"] or group["reference_entity_key"]
    ]  # fmt: skip
    rows = queries._lookup("entity", "entity_id", _group_members(queries, scope, groups))
    singletons = [
        key for group in page if not (group["connectivity"] or group["reference_entity_key"])
        for key in group["member_keys"]
    ]  # fmt: skip
    if singletons:
        rows += queries._entity_rows(
            "entity_key IN (SELECT unnest(?::VARCHAR[]))", [singletons], scope, nested=False
        )
    members = {key for group in page for key in group["member_keys"]}
    return [row for row in rows if row["entity_key"] in members]


@lru_cache(maxsize=1)
def chemical_types():
    from omnipath_core.biolink import entity_type, schema

    return tuple(
        entity_type(name) for name in schema().class_descendants("chemical entity", reflexive=True)
    )


def search_groups(queries, *, strategy="chemical_connectivity", **kwargs):
    """Cached per release and inventory: the explorer's first page is the same for everyone."""
    from omnipath_api.queries.entity_details import cache_key, page_details

    start = time.perf_counter()
    detail_limit = kwargs.pop("detail_limit", 20)
    detail_offset = kwargs.pop("detail_offset", 0)
    detail_field = kwargs.pop("detail_field", None)
    resources = kwargs.get("resources") or (kwargs.get("filters") or {}).get("sources")
    key = cache_key(
        queries, "groups", json.dumps([strategy, kwargs], sort_keys=True, default=str), resources
    )
    result = queries._detail_cache.get(key, lambda: _search_groups(queries, strategy, **kwargs))
    groups = [dict(group) for group in result["groups"]]
    if kwargs.get("include_details"):
        for group in groups:
            group["entity"] = page_details(
                group["entity"], detail_limit, detail_offset, detail_field
            )
    return dict(result, groups=groups, elapsed_ms=round((time.perf_counter() - start) * 1000, 2))


def _search_groups(
    queries,
    strategy,
    *,
    query="",
    filters=None,
    resources=None,
    limit=20,
    cursor=None,
    group_key=None,
    member_cursor=None,
    member_limit=5,
    include_details=False,
    include_member_keys=False,
    **kwargs,
):
    filters = normalize_filters(filters)
    scope = resources or filters["sources"] or None
    if not queries._selected_resource_infos(scope):
        return {"groups": [], "nextCursor": None}
    clauses, params = queries._entity_filter_clauses(filters, resources=resources)
    if group_key:
        # One group: its candidate entities, with all their copies in scope.
        clauses.append("entity_key IN (SELECT unnest(?::VARCHAR[]))")
        params.append(_candidate_keys(queries, group_key, scope))
    where, params = queries._entity_match_where(query, clauses, params, scope)
    types = chemical_types()
    marks = ",".join("?" for _ in types)
    chemical = strategy in {"auto", "chemical_connectivity"}
    gene = strategy in {"auto", "gene_reference"}
    chemical_sql = f"bool_or(coalesce(entity_type IN ({marks}), FALSE))"
    # One row per entity key (its copies agree on a group, or it is a singleton), then
    # one per group. Release-wide pages aggregate millions of keys: hashes keep the
    # aggregation on fixed-width values, and only the page's singletons get a key.
    cte = f"""WITH keys AS MATERIALIZED (
        SELECT any_value(entity_key) AS entity_key,
            CASE WHEN {chemical} AND {chemical_sql} THEN {_agreed("group_connectivity")} END
                AS connectivity,
            CASE WHEN {gene} AND NOT {chemical_sql}
                AND bool_and(coalesce(starts_with(reference_entity_key, 'entrez:'), FALSE))
                THEN {_agreed("reference_entity_key")} END AS reference_entity_key,
            sum(coalesce(relation_count, 0))::BIGINT AS relations
        FROM (SELECT entity_key, entity_type, group_connectivity, reference_entity_key,
                relation_count
            FROM {queries._table("entity", scope)} WHERE {where})
        GROUP BY hash(entity_key)
    )"""
    params = [*types, *types, *params]
    extra, values = "TRUE", []
    if group_key:
        extra, values = "group_key = ?", [group_key]
    elif cursor:
        try:
            relations, members, key = json.loads(cursor)
            if (
                not isinstance(relations, int)
                or not isinstance(members, int)
                or relations < 0
                or members < 1
                or not isinstance(key, str)
            ):
                raise ValueError("Invalid group cursor")
        except (TypeError, ValueError) as exc:
            raise ValueError("Invalid group cursor") from exc
        extra, values = (
            "(relations < ? OR (relations = ? AND (member_count < ?"
            " OR (member_count = ? AND group_key > ?))))",
            [relations, relations, members, members, key],
        )
    size = 1 if group_key else int(limit) + 1
    page = queries._fetch_dicts(
        cte
        + f""", counts AS (
            SELECT 'connectivity:' || connectivity AS group_key, connectivity,
                NULL::VARCHAR AS reference_entity_key, member_count, relations
            FROM (SELECT any_value(connectivity) AS connectivity, count(*) AS member_count,
                    sum(relations) AS relations
                FROM keys WHERE connectivity IS NOT NULL GROUP BY hash(connectivity))
            UNION ALL SELECT 'gene:' || reference_entity_key, NULL, reference_entity_key,
                member_count, relations
            FROM (SELECT any_value(reference_entity_key) AS reference_entity_key,
                    count(*) AS member_count, sum(relations) AS relations
                FROM keys WHERE reference_entity_key IS NOT NULL GROUP BY hash(reference_entity_key))
            UNION ALL (SELECT * FROM (
                SELECT 'entity:' || entity_key AS group_key, NULL, NULL, 1 AS member_count,
                    relations FROM keys
                WHERE connectivity IS NULL AND reference_entity_key IS NULL)
                WHERE {extra} ORDER BY relations DESC, member_count DESC, group_key LIMIT {size})
        ), page AS (
            -- Most studied first, as the ungrouped search ranks equal matches: a search
            -- for EGFR leads with the human gene, not the group with the most members.
            -- Member count only breaks ties.
            SELECT * FROM counts WHERE {extra} ORDER BY relations DESC, member_count DESC, group_key LIMIT {size}
        ), members AS (
            SELECT group_key, entity_key FROM page JOIN keys USING (connectivity)
            UNION ALL SELECT group_key, entity_key FROM page JOIN keys USING (reference_entity_key)
            UNION ALL SELECT p.group_key, k.entity_key FROM page p JOIN keys k
                ON p.group_key = 'entity:' || k.entity_key
        ) SELECT page.*, list(m.entity_key ORDER BY m.entity_key) AS member_keys
        FROM page JOIN members m USING (group_key)
        GROUP BY ALL ORDER BY relations DESC, member_count DESC, group_key""",
        [*params, *values, *values],
    )
    last = page[int(limit) - 1] if not group_key and len(page) > int(limit) else None
    next_cursor = (
        json.dumps([int(last["relations"]), int(last["member_count"]), last["group_key"]])
        if last
        else None
    )
    page = page[: int(limit)]
    rows = _member_rows(queries, page, scope)
    if include_details:
        queries._with_entity_children(rows)
    by_key = defaultdict(list)
    for row in rows:
        by_key[row["entity_key"]].append(row)
    groups = []
    for group in page:
        keys = group.pop("member_keys")
        members = [queries._merge_entity_row_group(by_key[key]) for key in keys if key in by_key]
        summaries = [queries._to_entity_summary(row) for row in members]
        visible = [s for s in summaries if s["entityPk"] > (member_cursor or "")]
        if group["connectivity"]:
            entity = _chemical_card(queries, group, members, summaries, include_details)
        elif group["reference_entity_key"]:
            entity = _gene_card(group, members, summaries)
        else:
            entity = summaries[0]
        is_group = bool(group["connectivity"] or group["reference_entity_key"])
        next_member_cursor = (
            visible[member_limit - 1]["entityPk"] if len(visible) > member_limit else None
        )
        if is_group:
            entity.update(
                # the group's relations, not its representative member's
                relationCount=sum(int(m.get("relation_count") or 0) for m in members),
                entityPk=group["group_key"],
                groupMemberKeys=[
                    s["entityPk"]
                    for s in (summaries if include_member_keys else visible[:member_limit])
                ],
                groupMemberCount=len(summaries),
                groupDetailsLoaded=include_details,
                groupQuery=query,
                groupFilters=filters,
                groupResources=scope or [],
            )
        entity["groupMemberCursor"] = next_member_cursor
        group.update(
            entity=entity,
            members=visible[:member_limit],
            is_group=is_group,
            group_label=group["connectivity"] or group["reference_entity_key"],
            nextMemberCursor=next_member_cursor,
        )
        groups.append(group)
    return {"groups": groups, "nextCursor": next_cursor}


def _union(summaries, field, identity):
    return list(
        {
            identity(item): item for summary in summaries for item in summary.get(field) or []
        }.values()
    )


def _chemical_card(queries, group, members, summaries, include_details):
    """A structure group: every member's identifiers, labels and sources; its name the
    preferred one of the member labels."""
    entity = queries._to_entity_summary(
        queries._merge_entity_row_group(members if include_details else members[:1])
    )
    identifiers = _union(summaries, "identifiers", lambda i: (i["identifierType"], i["identifier"]))
    entity.update(
        identifiers=identifiers,
        identifiersTotal=len(identifiers),
        sources=sorted({s for summary in summaries for s in summary["sources"]}),
        displayName=_card_name(members, summaries) or display_name(entity),
        canonicalIdentifier=group["connectivity"],
        canonicalIdentifierType="connectivity",
        primaryNamespace="connectivity",
        primaryIdentifier=group["connectivity"],
        groupStrategy="chemical_connectivity",
        # A stereoisomer's depiction cannot represent the whole group.
        ontologyHierarchy=None,
    )
    if include_details:
        entity["entityAttributes"] = (
            _union(
                summaries, "entityAttributes", lambda a: json.dumps(a, sort_keys=True, default=str)
            )
            or None
        )
    if len({m.get("taxon") for m in members}) > 1:
        entity.update(taxonomyId=None, taxonomyName=None, taxon=None)
    return entity


def _card_name(members, summaries):
    """A group's name: the label of its most connected well-named member, so a structure
    group of 52 hexoses is 'D-Glucose', not the shortest name among them ('Allose')."""
    named = [
        (-int(member.get("relation_count") or 0), name_rank(summary["label"]), summary["label"])
        for member, summary in zip(members, summaries)
        if summary.get("label") and name_rank(summary["label"])[0] < 8
    ]
    return min(named)[2] if named else preferred_name(s["label"] for s in summaries)


def _gene_card(group, members, summaries):
    """A gene group: annotations, identifiers and sources of every member; a gene's
    function text sits on the protein record resolved to it, not on the first member.
    It is named by its gene, not by a product's entry name (Q53GA5_HUMAN)."""
    genes = [i for i, s in enumerate(summaries) if s["entityType"] == "gene"]
    entity = dict(summaries[genes[0]] if genes else summaries[0])
    entity["displayName"] = (
        summaries[genes[0]]["displayName"] if genes else _card_name(members, summaries)
    ) or entity.get("displayName")
    identifiers = _union(summaries, "identifiers", lambda i: (i["identifierType"], i["identifier"]))
    attributes = _union(
        summaries, "entityAttributes", lambda a: json.dumps(a, sort_keys=True, default=str)
    )
    entity.update(
        entityAttributes=attributes or None,
        identifiers=identifiers,
        identifiersTotal=len(identifiers),
        sources=sorted({s for summary in summaries for s in summary["sources"]}),
        referenceEntityKey=group["reference_entity_key"],
        entityType="gene",
        memberEntityTypes=sorted({s["entityType"] for s in summaries}),
        canonicalIdentifier=group["reference_entity_key"].split(":", 1)[-1],
        canonicalIdentifierType="entrez",
        groupStrategy="gene_reference",
    )
    return entity
