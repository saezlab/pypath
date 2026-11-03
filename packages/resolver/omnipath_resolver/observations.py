"""Stable namespace encoding and normalized observation evidence policy."""

from __future__ import annotations

from omnipath_core.biolink import entity_type as normalize_entity_type
import os
from functools import lru_cache

NS = (
    "inchikey",
    "pubchem",
    "chebi",
    "chembl",
    "hmdb",
    "lipidmaps",
    "swisslipids",
    "bigg",
    "metanetx",
    "uniprot",
    "uniprot-sec",
    "uniprot_entry",
    "entrez",
    "ensg",
    "enst",
    "ensp",
    "hgnc",
    "refseq",
    "refseq_protein",
    "genesymbol",
    "genesymbol-syn",
    "inchi",
    "kegg",
    "cas",
    "drugbank",
    "goslin",
    "refmet",
    "genbank",
    "ramp",
    "ramp_gene",
    "kegg_gene",
    # reactions (appended: existing codes must not move)
    "rhea",
    "kegg_reaction",
    "reactome",
    "metacyc_reaction",
    "ecocyc_reaction",
    "macie",
    "bigg_reaction",
    "vmh_reaction",
    "seed_reaction",
    "sabiork_reaction",
    "metanetx_reaction",
)
CODES = {n: i + 1 for i, n in enumerate(NS)}
# Key target byte per reference library.
TARGETS = {"chemical": 1, "gene_protein": 2, "reaction": 3}

_LIPID_CACHE = None
_LIPID_PENDING = []


def key(kind, route, ns, scope, value):
    return (
        bytes([1, kind, route])
        + CODES[ns].to_bytes(2, "big")
        + (b"\x01" + int(scope).to_bytes(4, "big") if scope else b"\x00")
        + value.encode()
    )


def flush_lipid_cache():
    if _LIPID_PENDING:
        _LIPID_CACHE.store(_LIPID_PENDING)
        _LIPID_PENDING.clear()


@lru_cache(maxsize=16384)
def lipid_name_result(name):
    """Reuse the same Goslin normalizer and persistent cache as the reference."""
    global _LIPID_CACHE
    from .goslin_cache import NormalizationCache
    from .goslin import normalize

    if _LIPID_CACHE is None:
        _LIPID_CACHE = NormalizationCache(
            os.environ.get("OMNIPATH_GOSLIN_CACHE", "data/reference/.goslin-cache"), normalize
        )
    cached = _LIPID_CACHE.lookup([name])
    if cached:
        return cached[0]
    result = normalize(name)
    _LIPID_PENDING.append(result)
    if len(_LIPID_PENDING) >= 512:
        flush_lipid_cache()
    return result


def lipid_fallback_votes(normalized, observed, library):
    """Use exact lipid descriptors only when no stronger identifier can vote."""
    # A native RefMet ID can describe an underspecified lipid; it does not
    # replace a structural descriptor (and older reference versions lack it).
    if library != "chemical" or any(v.ns != "refmet" for v in normalized):
        return normalized
    from .canonical.match import Vote

    names = set(observed.get("name", []) + observed.get("synonym", []))
    parsed = [lipid_name_result(n) for n in sorted(names) if n.strip()]
    parsed = [p for p in parsed if p["status"] == "parsed"]
    if not parsed:
        return normalized
    specificity = max(p["specificity"] for p in parsed)
    return normalized + [
        Vote("goslin", value, None, "id", False)
        for value in sorted({p["goslin"] for p in parsed if p["specificity"] == specificity})
    ]


def observation_bundle(input_id, obs, normalized, observed, library, source="runtime"):
    """Encode normalized evidence for both serving builds and offline replay."""
    from .canonical.policy import INCHIKEY_RE, RNA_ENTITY_TYPES

    target = TARGETS[library]
    normalized = lipid_fallback_votes(normalized, observed, library)
    symbol_scoped = target == 2 and any(v.kind == "symbol" for v in normalized)
    bundle = set()
    for v in normalized:
        if v.ns not in CODES:
            raise ValueError("Unhandled eligible namespace " + v.ns)
        scope = (v.taxon or "") if symbol_scoped else ""
        route = 2 if target == 2 and v.ns in {"entrez", "ramp_gene", "kegg_gene"} else 1
        # Only a stated full InChIKey decides alone; one derived from SMILES is an ordinary vote.
        anchor = (
            v.id
            if target == 1
            and v.ns == "inchikey"
            and not getattr(v, "derived", False)
            and INCHIKEY_RE.fullmatch(v.id)
            else ""
        )
        bundle.add((v.ns, v.id, scope, route, anchor, v.primary))
    query = dict(
        input_id=input_id,
        source=source,
        entity_key=input_id,
        entity_type=obs.entity_type,
        namespace=str(obs.namespace),
        identifier=obs.identifier,
        target=target,
    )
    rows = [
        dict(
            input_id=input_id,
            ns=ns,
            identifier=value,
            scope=scope,
            anchor=anchor,
            target=target,
            route=route,
            ordinal=j,
            primary=primary,
            lookup_key=key(target, route, ns, scope, value),
        )
        for j, (ns, value, scope, route, anchor, primary) in enumerate(sorted(bundle))
    ]
    if target == 2 and normalize_entity_type(obs.entity_type) in RNA_ENTITY_TYPES:
        for row in rows:
            row["gene_only"] = True
    return query, rows
