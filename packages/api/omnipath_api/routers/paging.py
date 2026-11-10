"""Paging the nested collections of an entity response (the web app's detail tabs)."""

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
