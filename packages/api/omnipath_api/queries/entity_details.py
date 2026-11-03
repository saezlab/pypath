"""Release-aware entity hydration and paged, resource-local relationship reads."""

import json

# Structural relationships shown in an entity's dialog; other relations are paged
# from relation search.
STRUCTURAL = ("has_input", "has_output", "enabled_by", "has_member", "has_part")


def cache_key(queries, kind, public_id, resources, *page):
    return json.dumps(
        [
            kind,
            public_id,
            page,
            [info["key"] for info in queries._selected_resource_infos(resources)],
            queries._inventory_fingerprint,
            queries._release_scope.get(),
            queries._taxonomy_cache_version(),
            [r.get("references") for r in queries.releases.list()],
        ],
        sort_keys=True,
        default=str,
    )


def _copies(queries, public_id, resources):
    """The entity's rows in each selected resource, by key (files are sorted by it)."""
    return queries._entity_rows("entity_key = ?", [public_id], resources, nested=False)


def core(queries, public_id, resources=None):
    def load():
        copies = queries._with_entity_children(_copies(queries, public_id, resources))
        if not copies:
            return None
        entity = queries._to_entity_summary(queries._merge_entity_row_group(copies))
        # Evidence is paged from its own endpoint; the entity carries only its count.
        entity["molecularEvidenceTotal"] = sum(int(c["evidence_count"] or 0) for c in copies)
        return {"entity": entity}

    return queries._detail_cache.get(cache_key(queries, "entity", public_id, resources), load)


def evidence(queries, public_id, resources=None, limit=20, offset=0):
    """A page of an entity's evidence items, by resource and position."""

    def load():
        copies = _copies(queries, public_id, resources)
        if not copies:
            return None
        total = sum(int(c["evidence_count"] or 0) for c in copies)
        items, start = [], 0
        for copy in copies:
            count = int(copy["evidence_count"] or 0)
            first, last = max(offset - start, 0), min(offset + limit - start, count)
            if first < last:
                rows = queries._lookup(
                    "entity_evidence",
                    "entity_id",
                    [(copy["resource"], copy["entity_id"])],
                    extra="ORDER BY ordinal LIMIT ? OFFSET ?",
                    params=[last - first, first],
                )
                for row in rows:
                    for name in ("entity_id", "ordinal", "resource"):
                        row.pop(name)
                items += rows
            start += count
        end = offset + len(items)
        return dict(
            evidence=items,
            evidenceTotal=total,
            nextCursor=str(end) if items and end < total else None,
        )

    return queries._detail_cache.get(
        cache_key(queries, "evidence", public_id, resources, limit, offset), load
    )


def _endpoint_page(queries, keys, resources, where="TRUE", params=(), limit=50, offset=0):
    """A page of distinct (resource, relation_id, predicate) of relations with an
    endpoint among ``keys``, ordered by predicate; each row carries the total."""
    return queries._fetch_dicts(
        f"""SELECT resource, relation_id, predicate, count(*) OVER () AS total FROM (
            SELECT DISTINCT resource, relation_id, predicate
            FROM {queries._table("relation_endpoint", resources)}
            WHERE key IN (SELECT unnest(?::VARCHAR[])) AND key_kind = 'entity' AND {where})
        ORDER BY predicate, resource, relation_id LIMIT ? OFFSET ?""",
        [keys, *params, limit, offset],
    )


def _relations(queries, pairs):
    """Relation rows for (resource, relation_id) pairs, in the pairs' order."""
    found = {
        (r["resource"], r["relation_id"]): r
        for r in queries._lookup("relation", "relation_id", pairs)
    }
    return [found[pair] for pair in pairs if pair in found]


def relationships(queries, public_id, resources=None, limit=50, offset=0):
    def load():
        keys = list(dict.fromkeys(public_id if isinstance(public_id, list) else [public_id]))
        empty = dict(relationships=[], relationshipsTotal=0, nextCursor=None)
        if not keys:
            return empty
        marks = ",".join("?" for _ in STRUCTURAL)
        rows = _endpoint_page(
            queries, keys, resources, f"predicate IN ({marks})", STRUCTURAL, limit, offset
        )
        if not rows:
            return empty
        records = _relations(queries, [(r["resource"], r["relation_id"]) for r in rows])
        queries._children("relation_annotation", records, "annotations")
        result = queries._relation_items(records, annotations=lambda r: r["annotations"])
        total = int(rows[0]["total"])
        return dict(
            relationships=result,
            relationshipsTotal=total,
            nextCursor=str(offset + len(rows)) if offset + len(rows) < total else None,
        )

    return queries._relationship_cache.get(
        cache_key(queries, "relationships", public_id, resources, limit, offset), load
    )


def details(queries, public_id, resources=None):
    """The entity, its structural relationships, relation count and associations."""
    found = core(queries, public_id, resources)
    if not found:
        return None
    entity = found["entity"]
    keys = list(dict.fromkeys([*(entity.get("sourceEntityPks") or []), entity["entityPk"]]))
    counts = queries._fetch_dicts(
        f"""SELECT count(*) AS n, count(*) FILTER (WHERE category = 'association') AS a FROM (
            SELECT DISTINCT resource, relation_id, category
            FROM {queries._table("relation_endpoint", resources)}
            WHERE key IN (SELECT unnest(?::VARCHAR[])) AND key_kind = 'entity')""",
        [keys],
    )[0]
    associations = []
    if counts["a"]:
        rows = _endpoint_page(queries, keys, resources, "category = 'association'", limit=100)
        associations = [
            {"relationPk": r["relation_key"], "predicate": r["predicate"]}
            for r in sorted(
                _relations(queries, [(r["resource"], r["relation_id"]) for r in rows]),
                key=lambda r: (r["relation_key"], r["resource"]),
            )
        ]
    page = relationships(queries, public_id, resources)
    return {
        "entity": entity,
        "relationships": page["relationships"],
        "relationshipsTotal": page["relationshipsTotal"],
        "summary": {"interactionCount": int(counts["n"] or 0)},
        "annotations": associations,
    }


# Descriptive annotations, in the order the explorer shows them; they lead the
# attribute pages, so the first page carries the entity's description.
DESCRIPTIVE_TERMS = (
    "up:function_annotation",
    "function",
    "description",
    "iao:0000115",
    "up:disease_annotation",
    "disease_involvement",
    "up:subcellular_location_annotation",
    "subcellular_location",
    "up:pathway_annotation",
    "up:activity_regulation_annotation",
    "up:ptm_annotation",
    "up:mutagenesis_annotation",
)
_DESCRIPTIVE_RANK = {term: rank for rank, term in enumerate(DESCRIPTIVE_TERMS)}


def _descriptions_first(attributes):
    """Descriptions lead; the rest is grouped by term so a page continues a group."""

    def rank(item):
        term = str(item.get("term") or "").lower().removeprefix("biolink:")
        position = _DESCRIPTIVE_RANK.get(term, len(DESCRIPTIVE_TERMS))
        return position, term if position == len(DESCRIPTIVE_TERMS) else ""

    return sorted(attributes, key=rank)  # stable: otherwise in source order


PAGED_FIELDS = (
    "identifiers",
    "entityAttributes",
    "sources",
    "resources",
    "sourceEntityPks",
    "molecularEvidence",
)


def page_details(entity, limit=20, offset=0, field=None):
    """Bound each nested collection in an entity response.

    With ``field``, only that collection is paged from ``offset``; the other
    collections are left out, so a tab can load more of its own list.
    """
    if entity.get("entityAttributes"):
        entity = dict(entity, entityAttributes=_descriptions_first(entity["entityAttributes"]))
    result = dict(entity)
    longest = 0
    for name in PAGED_FIELDS:
        if name == "molecularEvidence" and "molecularEvidenceTotal" in entity:
            continue  # counted only; the items are paged from /entities/{id}/evidence
        if field and name != field:
            result.pop(name, None)
            continue
        rows = entity.get(name) or []
        longest = max(longest, len(rows))
        result[name] = rows[offset : offset + limit]
        result[name + "Total"] = len(rows)
        result[name + "NextCursor"] = str(offset + limit) if offset + limit < len(rows) else None
    if not field:
        result["detailNextCursor"] = str(offset + limit) if offset + limit < longest else None
    return result
