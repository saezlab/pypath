"""API display names: the resolver's label, with identifier fallbacks."""

import re
from omnipath_core.display_names import (
    NAME_NAMESPACES,
    CHEMICAL_TYPES,
    preferred_name,
)


def display_name(entity):
    kind = re.sub(r"[\s_]", "", str(entity.get("entityType") or "")).lower()
    label = str(entity.get("label") or "").strip()
    identifiers = entity.get("identifiers") or []
    names = [label] if label else []
    for item in identifiers:
        ns = str(item.get("identifierType") or "").lower()
        if ns in NAME_NAMESPACES or ns.endswith(":name") or " entry name" in ns:
            names.append(item["identifier"])
    if kind in CHEMICAL_TYPES:
        # The resolver's label, unless it is only a number or an InChI string.
        if label and not label.isdigit() and not label.lower().startswith("inchi="):
            return label
        fallbacks = [
            item["identifier"]
            for ns in ("chebi", "chembl", "hmdb", "pubchem")
            for item in identifiers
            if item["identifierType"] == ns
        ]
        return (
            next((value for value in fallbacks if not value.isdigit()), None)
            or next(iter(fallbacks), None)
            or label
            or entity.get("canonicalIdentifier")
            or entity["entityPk"]
        )
    if kind in {"protein", "gene"}:
        symbols = [
            item["identifier"]
            for item in identifiers
            if item["identifierType"]
            in {"genesymbol", "gene_symbol", "gene_name_primary", "symbol", "om:0200"}
        ]
        uniprot = [
            item["identifier"]
            for item in identifiers
            if item["identifierType"] in {"uniprot", "uniprotkb"}
        ]
        return (
            label
            or next(iter(symbols), None)
            or next(iter(uniprot), None)
            or next(iter(names), None)
            or entity.get("canonicalIdentifier")
            or entity["entityPk"]
        )
    return (
        label or next(iter(names), None) or entity.get("canonicalIdentifier") or entity["entityPk"]
    )


__all__ = [
    "display_name",
    "preferred_name",
    "NAME_NAMESPACES",
    "CHEMICAL_TYPES",
]
