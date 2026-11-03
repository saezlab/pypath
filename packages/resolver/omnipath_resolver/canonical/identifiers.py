"""Identifier normalization shared by the library builder (SQL) and the matcher (Python).

Both sides must agree on the written form of an identifier, otherwise a vote
never finds its library row.  Rules are deliberately few:

* namespace: glossary slug (``UniProt`` → ``uniprot``, ``OM:0202`` → ``name`` …)
* chebi:     ``CHEBI:<digits>`` (bare numbers get the prefix, case normalised)
* hgnc:      ``HGNC:<digits>``
* inchikey:  upper case, ``InChIKey=`` prefix removed
* ensg:      trailing ``.<version>`` removed
* uniprot:   isoform suffix ``-<n>`` is *kept* on the vote, and the base
             accession is emitted as an additional vote
"""

from __future__ import annotations

import re

from omnipath_core.naming import normalize_namespace

_ISOFORM_RE = re.compile(r"^([A-Z][0-9][A-Z0-9]{3}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)-\d+$")
_ENSG_VERSION_RE = re.compile(r"\.\d+$")
_REACTOME_VERSION_RE = re.compile(r"\.[0-9]+$")


def normalize_ns(ns: str | None) -> str:
    slug = str(normalize_namespace(str(ns or "")))
    # Input modules name BiGG metabolites explicitly; chemical hubs use `bigg`.
    # Keep reaction namespaces distinct and preserve compartment suffixes.
    return {"bigg_metabolite": "bigg", "uniprot_trembl": "uniprot"}.get(slug, slug)


def normalize_id(ns: str, value: str | None) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    if ns == "chebi":
        digits = text.upper().removeprefix("CHEBI:").strip()
        return f"CHEBI:{digits}" if digits.isdigit() else text
    if ns == "hgnc":
        digits = text.upper().removeprefix("HGNC:").strip()
        return f"HGNC:{digits}" if digits.isdigit() else text
    if ns == "hmdb":
        match = re.fullmatch(r"(?:HMDB:?)?([0-9]+)", text.upper())
        if match:
            return "HMDB" + (match.group(1).lstrip("0") or "0").zfill(7)
    if ns == "inchikey":
        return text.removeprefix("InChIKey=").removeprefix("INCHIKEY=").upper()
    if ns in {"ensg", "enst", "ensp"}:
        return _ENSG_VERSION_RE.sub("", text)
    if ns == "chembl":
        return text.upper()
    if ns == "rhea":
        digits = re.sub(r"^RHEA:", "", text, flags=re.I).strip()
        return digits if digits.isdigit() else text
    if ns == "reactome":
        return _REACTOME_VERSION_RE.sub("", re.sub(r"^reactome:", "", text, flags=re.I).strip())
    return text


# KEGG names the genes of these organisms by their NCBI Gene ID: hsa:7157 is Entrez 7157.
# No reference source carries KEGG gene IDs, so this is how a KEGG-only gene resolves.
KEGG_NCBI_GENE_ORGANISMS = frozenset({"hsa", "mmu", "rno"})


def normalize_identifier(ns: str | None, value: str | None) -> list[tuple[str, str]]:
    """All ``(ns, id)`` votes an observed identifier contributes (usually one)."""
    slug = normalize_ns(ns)
    if slug == "ensembl":
        match = re.fullmatch(r"ENS[A-Z]*([GTP])\d+(?:\.\d+)?", str(value or "").strip())
        if match:
            slug = {"G": "ensg", "T": "enst", "P": "ensp"}[match.group(1)]
    ident = normalize_id(slug, value)
    if not slug or not ident:
        return []
    out = [(slug, ident)]
    if slug == "uniprot":
        m = _ISOFORM_RE.match(ident)
        if m:
            out.append((slug, m.group(1)))
    if slug == "kegg_gene":
        organism, _, number = ident.partition(":")
        if organism.lower() in KEGG_NCBI_GENE_ORGANISMS and number.isdigit():
            out.append(("entrez", number))
    return out


def normalize_id_sql(ns_expr: str, id_expr: str) -> str:
    """DuckDB expression producing the same written form as :func:`normalize_id`."""
    return f"""
    CASE {ns_expr}
      WHEN 'chebi' THEN
        CASE WHEN regexp_matches(upper(trim({id_expr})), '^(CHEBI:)?[0-9]+$')
             THEN 'CHEBI:' || regexp_replace(upper(trim({id_expr})), '^CHEBI:', '')
             ELSE trim({id_expr}) END
      WHEN 'hgnc' THEN
        CASE WHEN regexp_matches(upper(trim({id_expr})), '^(HGNC:)?[0-9]+$')
             THEN 'HGNC:' || regexp_replace(upper(trim({id_expr})), '^HGNC:', '')
             ELSE trim({id_expr}) END
      WHEN 'hmdb' THEN
        CASE WHEN regexp_full_match(upper(trim({id_expr})), '(HMDB:?)?[0-9]+')
             THEN 'HMDB' || lpad(coalesce(nullif(ltrim(regexp_extract(trim({id_expr}), '[0-9]+'), '0'), ''), '0'),
                  greatest(7, length(coalesce(nullif(ltrim(regexp_extract(trim({id_expr}), '[0-9]+'), '0'), ''), '0')))::INTEGER, '0')
             ELSE trim({id_expr}) END
      WHEN 'inchikey' THEN upper(regexp_replace(trim({id_expr}), '^(?i)inchikey=', ''))
      WHEN 'ensg' THEN regexp_replace(trim({id_expr}), '\\.[0-9]+$', '')
      WHEN 'enst' THEN regexp_replace(trim({id_expr}), '\\.[0-9]+$', '')
      WHEN 'ensp' THEN regexp_replace(trim({id_expr}), '\\.[0-9]+$', '')
      WHEN 'chembl' THEN upper(trim({id_expr}))
      WHEN 'rhea' THEN
        CASE WHEN regexp_full_match(trim({id_expr}), '(?i)(rhea:)?[0-9]+')
             THEN regexp_replace(trim({id_expr}), '(?i)^rhea:', '')
             ELSE trim({id_expr}) END
      WHEN 'reactome' THEN regexp_replace(regexp_replace(trim({id_expr}), '(?i)^reactome:', ''), '\\.[0-9]+$', '')
      ELSE trim({id_expr})
    END"""
