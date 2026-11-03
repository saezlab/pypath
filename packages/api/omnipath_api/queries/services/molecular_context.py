"""Catalogue products and evidence-supported forms have distinct navigation."""

import json
from typing import Literal

from omnipath_api.queries.molecular import matching_evidence

_ISOFORM = "molecular_form.isoform_identifier.ns || ':' || molecular_form.isoform_identifier.id = ?"


def _group_members(queries, resources, group_key, kind):
    """(resource, entity_id) of the entities in one group of ``entity_group``."""
    rows = queries._fetch_dicts(
        f"SELECT DISTINCT resource, entity_id FROM {queries._table('entity_group', resources)} "
        "WHERE group_key = ? AND kind = ? ORDER BY resource, entity_id",
        [group_key, kind],
    )
    return [(r["resource"], r["entity_id"]) for r in rows]


def _product_rows(queries, resources, reference):
    """Catalogue product rows of a gene: the entities linked to it, except Entrez records."""
    rows = queries._lookup(
        "entity",
        "entity_id",
        _group_members(queries, resources, reference, "gene"),
        extra="AND namespace <> 'entrez'",
    )
    return queries._with_entity_children(rows)


def _occurrences(queries, sql, params, keys, limit, offset):
    """A page of entity evidence rows as {entityPk, occurrence}, by resource and position."""
    rows = queries._fetch_dicts(
        f"SELECT * FROM ({sql}) ORDER BY resource, entity_id, ordinal LIMIT ? OFFSET ?",
        [*params, limit + 1, offset],
    )
    if keys is None:
        keys = {
            (r["resource"], r["entity_id"]): r["entity_key"]
            for r in queries._lookup(
                "entity", "entity_id", [(r["resource"], r["entity_id"]) for r in rows],
                "entity_id, entity_key",
            )
        }  # fmt: skip
    result = []
    for row in rows:
        entity = keys[(row.pop("resource"), row.pop("entity_id"))]
        row.pop("ordinal")
        result.append({"entityPk": entity, "occurrence": row})
    return result


def _entity_evidence(queries, resources, pairs, isoform_identifier, limit, offset):
    """Evidence items of the given entity rows, optionally of one isoform."""
    keys = {
        (r["resource"], r["entity_id"]): r["entity_key"]
        for r in queries._lookup("entity", "entity_id", pairs, "entity_id, entity_key")
    }
    if not keys:
        return []
    condition, values = "TRUE", []
    if isoform_identifier:
        condition, values = _ISOFORM, [isoform_identifier]
    sql, params = queries._lookup_sql("entity_evidence", "entity_id", list(keys), condition, values)
    return _occurrences(queries, sql, params, keys, limit, offset)


def _product_evidence(
    queries, resources, product_kind, entity_id, isoform_identifier, limit, offset
):
    # An asserted native product can be named by observations on a gene row
    # without having any catalogue gene links: every row's evidence is read.
    conditions = [f"molecular_form.{product_kind}_entity_key = ?"]
    values = [entity_id]
    if isoform_identifier:
        conditions.append(_ISOFORM)
        values.append(isoform_identifier)
    sql = f"""SELECT * FROM {queries._table("entity_evidence", resources)}
        WHERE {" AND ".join(conditions)}"""
    return _occurrences(queries, sql, values, None, limit, offset)


def _product_kind(summary):
    entity_type = (summary.get("entityType") or "").lower()
    # Reusable writer records have these literal types. A namespace describes
    # an identifier, not whether this navigation selected a product record.
    return entity_type if entity_type in {"protein", "transcript"} else None


def _referenced_products(queries, relations, standalone, known_keys, resources):
    """Hydrate at most 100 additional exact product records in one scalar lookup."""
    keys = {}

    def remember(form):
        for field in ("protein_entity_key", "transcript_entity_key"):
            key = (form or {}).get(field)
            if key and key not in known_keys:
                keys.setdefault(key, None)
                if len(keys) == 100:
                    return True
        return False

    # Give standalone observations their labels even when a relation has many
    # different paired forms. Both endpoint products of each returned pair count.
    full = False
    for record in standalone:
        if remember(record["occurrence"].get("molecular_form")):
            full = True
            break
    if not full:
        for row in relations:
            for evidence in row["evidence"]:
                for side in ("subject", "object"):
                    if remember(evidence.get(side + "_molecular_form")):
                        full = True
                        break
                if full:
                    break
            if full:
                break
    if not keys:
        return []
    rows = queries._fetch_entities_by_keys(list(keys), resources, slim=True)
    return [queries._to_entity_summary(row) for row in rows]


def context(
    queries,
    entity_id,
    resources=None,
    isoform_identifier=None,
    limit=20,
    offset=0,
    view: Literal["reference", "product"] = "reference",
):
    if view not in {"reference", "product"}:
        raise ValueError("Context view must be reference or product.")
    if isoform_identifier and view != "product":
        raise ValueError("Isoform selection requires product view.")
    gene_group = entity_id.startswith("gene:")
    entity = None if gene_group else queries.get_entity_core(entity_id, resources)
    summary = entity["entity"] if entity else {}
    reference = entity_id[5:] if gene_group else summary.get("referenceEntityKey")
    if not gene_group and not entity:
        return None
    product_kind = _product_kind(summary)
    is_product = view == "product"
    if is_product and product_kind is None:
        raise ValueError("Product view requires a protein or transcript record.")
    gene_reference = bool(reference and reference.startswith("entrez:"))
    filters = (
        {"reference_entity_keys": [reference]} if gene_reference else {"entity_pks": [entity_id]}
    )
    if is_product:
        field = product_kind + "_entity_keys"
        filters = {field: [entity_id]}
    if isoform_identifier:
        filters["isoform_identifiers"] = [isoform_identifier]
    page = queries.search_relations(filters, resources, limit, offset, include_details=True)
    relations = queries._relation_items(page["rows"], evidence=lambda row: row["evidence"])
    forms = {}
    for row in page["rows"]:
        for evidence in row.get("evidence") or []:
            for side in ("subject", "object"):
                if not is_product:
                    field, selected = (
                        ("_reference_entity_key", reference)
                        if gene_reference
                        else ("_entity_key", entity_id)
                    )
                    if row.get(side + field) != selected:
                        continue
                form = evidence.get(side + "_molecular_form")
                if (
                    form
                    and is_product
                    and not matching_evidence(
                        [{"subject_molecular_form": form}],
                        {**filters, "molecular_endpoint_mode": "source"},
                    )
                ):
                    continue
                if form:
                    forms.setdefault(json.dumps(form, sort_keys=True), form)
    products = []
    if reference:
        products = [
            queries._to_entity_summary(r)
            for r in queries._merge_by_key(_product_rows(queries, resources, reference))
        ]
    # Standalone evidence can predate stored references; native typed-key
    # selection still applies, and product view keeps its exact form filter.
    if is_product:
        standalone = _product_evidence(
            queries, resources, product_kind, entity_id, isoform_identifier, limit, offset
        )
    else:
        pairs = (
            _group_members(queries, resources, reference, "reference")
            if gene_reference
            else [
                (r["resource"], r["entity_id"])
                for r in queries._entity_rows("entity_key = ?", [entity_id], resources, nested=False)
            ]
        )
        standalone = _entity_evidence(queries, resources, pairs, isoform_identifier, limit, offset)
    standalone_more = len(standalone) > limit
    standalone = standalone[:limit]
    for record in standalone:
        form = record["occurrence"].get("molecular_form")
        if form:
            forms.setdefault(json.dumps(form, sort_keys=True), form)
    referenced_products = _referenced_products(
        queries,
        relations,
        standalone,
        {r[side]["entityPk"] for r in relations for side in ("subjectEntity", "objectEntity")}
        | {product["entityPk"] for product in products},
        resources,
    )
    return {
        "referenceEntityKey": reference,
        "standaloneEvidence": standalone,
        "catalogueProducts": products,
        "referencedProducts": referenced_products,
        "observedForms": list(forms.values()),
        "relations": relations,
        "relationsTotal": page["total"],
        "nextCursor": str(offset + limit)
        if offset + len(relations) < page["total"] or standalone_more
        else None,
        "filters": filters,
    }
