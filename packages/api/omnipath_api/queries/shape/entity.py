"""Shape raw DuckDB entity rows into API and Svelte domain models."""

from __future__ import annotations

from typing import Any
from omnipath_core.naming import normalize_namespace
from omnipath_api.queries.shape.display_name import display_name


def shape_entity_identifier(
    entity_pk: str,
    raw_type: str,
    raw_value: str,
    source: str = "",
    is_canonical: bool = False,
) -> dict[str, Any]:
    """Shape an identifier dictionary with slug type, id, and optional curie."""
    slug = normalize_namespace(raw_type)
    val = str(raw_value or "").strip()
    curie = None
    return {
        "entityPk": entity_pk,
        "identifier": val,
        "identifierType": slug,
        "type": slug,
        "id": val,
        "curie": curie,
        "source": source,
        "isCanonical": is_canonical,
    }


def shape_entity_summary(row: dict[str, Any], resource_hint: str | None = None) -> dict[str, Any]:
    """Convert DuckDB entity row into standard API entity summary."""
    entity_pk = str(row.get("entity_key") or "")
    raw_ns = str(row.get("namespace") or "unknown").strip()
    identifier = str(row.get("identifier") or entity_pk).strip()
    label = str(row.get("label") or identifier).strip()
    slug_ns = normalize_namespace(raw_ns)

    identifiers: list[dict[str, Any]] = []
    sources: list[str] = []
    seen: set[tuple[str, str]] = set()

    def add_ident(id_type: str, val: str, src: str = "", is_canon: bool = False) -> None:
        item = shape_entity_identifier(entity_pk, id_type, val, src, is_canon)
        key = (item["type"], item["id"])
        if not key[0] or not key[1] or key in seen:
            return
        seen.add(key)
        identifiers.append(item)
        if src and src not in sources:
            sources.append(src)

    add_ident(slug_ns, identifier, src=resource_hint or "", is_canon=True)
    if label and label != identifier:
        add_ident("genesymbol" if slug_ns in ("uniprot", "entrez", "ensg") else "name", label)

    nested_ids = row.get("identifiers") or []
    if isinstance(nested_ids, list):
        for item in nested_ids:
            if not isinstance(item, dict):
                continue
            ns = item.get("ns") or item.get("type") or ""
            val = item.get("id") or item.get("identifier") or ""
            src = item.get("source") or ""
            canon = bool(item.get("is_canonical", False))
            add_ident(ns, val, src=src, is_canon=canon)

    if resource_hint and resource_hint not in sources:
        sources.append(resource_hint)

    annotations: list[dict[str, Any]] = []
    nested_anns = row.get("annotations") or []
    if isinstance(nested_anns, list):
        for ann in nested_anns:
            if not isinstance(ann, dict):
                continue
            term = str(ann.get("term") or "").strip()
            val = str(ann.get("value") or "").strip()
            src = str(ann.get("source") or "").strip()
            dataset = str(ann.get("dataset") or "").strip()
            if src and src not in sources:
                sources.append(src)
            annotations.append(
                {
                    "term": term,
                    "value": val,
                    **({"quantity": ann["quantity"]} if ann.get("quantity") is not None else {}),
                    "source": src,
                    "dataset": dataset,
                }
            )

    has_hierarchy = bool(row.get("has_hierarchy", False))
    parent_count = int(row.get("parent_count") or 0)
    child_count = int(row.get("child_count") or 0)

    ontology_hierarchy = None
    if has_hierarchy or (parent_count + child_count) > 0:
        prefix = slug_ns
        if slug_ns in {"OM:0204", "OM:0012", "cvterm", "cv_term"} or not prefix:
            if "-" in identifier:
                prefix = identifier.split("-")[0]
            elif ":" in identifier:
                prefix = identifier.split(":")[0]
            else:
                prefix = "KW" if identifier.startswith("KW") else slug_ns

        ontology_hierarchy = {
            "termId": identifier,
            "ontologyPrefix": prefix,
            "label": label,
            "definition": next(
                (
                    a["value"]
                    for a in annotations
                    if a["term"] in {"description", "biolink:description", "IAO:0000115"}
                    and a["value"]
                ),
                None,
            ),
            "ontologyId": None,
            "childCount": child_count,
            "parentCount": parent_count,
        }

    # Keys of the per-resource rows folded into this summary by the queries'
    # cross-resource merge; relation hydration and detail counts fan out over them.
    source_pks = [str(k) for k in (row.get("_source_entity_keys") or []) if str(k).strip()]
    if entity_pk and entity_pk not in source_pks:
        source_pks.append(entity_pk)

    result = {
        "entityPk": entity_pk,
        "referenceEntityKey": row.get("reference_entity_key"),
        "geneReferenceKeys": row.get("gene_reference_keys") or [],
        "molecularEvidence": row.get("evidence") or [],
        "resolutionStatus": "resolved",
        "sourceEntityPks": source_pks,
        "entityFacetHints": {
            "chemicalClasses": [],
            "metabolicDomains": [],
            "structuralSpecificities": [],
        },
        "entityType": str(row.get("entity_type") or "unknown").strip().lower(),
        "primaryNamespace": slug_ns,
        "canonicalIdentifier": identifier,
        "canonicalIdentifierType": slug_ns,
        "primaryIdentifier": identifier,
        "label": label,
        "taxon": row.get("taxon"),
        "taxonomyId": str(row.get("taxon") or "") or None,
        "hasHierarchy": has_hierarchy,
        "parentCount": parent_count,
        "childCount": child_count,
        "relationCount": int(row.get("relation_count") or 0),
        "sources": sources,
        "sourceCount": len(sources),
        "resources": sources,
        "entityAttributes": annotations if annotations else None,
        "ontologyHierarchy": ontology_hierarchy,
        "identifiersTotal": len(identifiers),
        "identifiers": identifiers,
        "annotations": annotations,
    }
    result["displayName"] = row.get("_display_name") or display_name(result)
    return result
