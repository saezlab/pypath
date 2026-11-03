"""Authoritative label assignment for canonical entities."""

from __future__ import annotations

import re
from typing import Mapping, Sequence

from .policy import EntityPolicy

_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z0-9]$")

_UNIPROT_KEYWORD_LABELS: dict[str, str] = {
    "KW-9990": "Technical term",
    "KW-9991": "PTM",
    "KW-9992": "Molecular function",
    "KW-9993": "Ligand",
    "KW-9994": "Domain",
    "KW-9995": "Disease",
    "KW-9996": "Developmental stage",
    "KW-9997": "Coding sequence diversity",
    "KW-9998": "Cellular component",
    "KW-9999": "Biological process",
}


def assign_label(
    canonical_namespace: str,
    canonical_identifier: str,
    ids_by_ns: Mapping[str, Sequence[str]],
    policy: EntityPolicy,
) -> str:
    """Assign display label strictly from policy priorities, filters, or fallback."""
    for ns in policy.label_priority:
        if ns in ("glossary", "uniprot_keyword"):
            if canonical_identifier in _UNIPROT_KEYWORD_LABELS:
                return _UNIPROT_KEYWORD_LABELS[canonical_identifier]
            continue

        candidates = ids_by_ns.get(ns, [])
        if not candidates:
            continue

        valid: list[str] = []
        for cand in candidates:
            s = str(cand).strip()
            if not s:
                continue
            if all(f(s) for f in policy.label_filters):
                valid.append(s)

        if not valid:
            continue
        valid = sorted(set(valid))

        if ns == "genesymbol":
            return sorted(valid, key=lambda s: (len(s) > 15, len(s), s))[0]
        if ns == "name":
            if policy.entity_class == "chemical":
                return sorted(valid, key=lambda s: (len(s) > 80, len(s), s))[0]
            return sorted(valid, key=lambda s: (len(s) > 120, len(s), s))[0]
        if ns == "chebi":
            c = valid[0]
            return c if str(c).upper().startswith("CHEBI") else f"CHEBI:{c}"
        if ns == "pubchem":
            return f"CID:{valid[0]}"
        return valid[0]

    return canonical_identifier
