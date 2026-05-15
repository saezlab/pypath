"""KEGG REST bulk-download parser for metabolic reactions, enzymes and compounds."""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any


_DATA_CACHE: dict[str, list[dict[str, Any]]] = {}

_MOUSE_TAXON = '10090'

_ARROW_RE = re.compile(r'<=>|=>|<=')
_STOICH_RE = re.compile(r'^(?:>\s*)?(?:\d+|[a-zA-Z]|\([^)]+\))\s+')


def _raw(
    opener,
    data_type: str,
    organism: str = 'mmu',
    force_refresh: bool = False,
    max_records: int | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield normalised KEGG records for the requested dataset.

    Parameters
    ----------
    opener:
        A pypath curl opener whose ``result`` attribute is a dict of
        file-like objects keyed by filename, matching the multi-file
        download pattern used by dlmachine.
    data_type:
        One of ``'reactions'``.  Additional types (e.g. ``'compounds'``)
        may be added here as the module evolves.
    organism:
        KEGG organism code, e.g. ``'mmu'`` (mouse) or ``'hsa'`` (human).
    force_refresh:
        Bypass the in-process cache and re-parse from the opener.
    max_records:
        If set, yield at most this many records (useful during development).
    """
    records = _load_records(opener, organism=organism, force_refresh=force_refresh)
    rows = records.get(data_type, [])

    if max_records is not None:
        rows = rows[:max_records]

    yield from rows


def _load_records(
    opener,
    organism: str = 'mmu',
    force_refresh: bool = False,
) -> dict[str, list[dict[str, Any]]]:

    cache_key = organism
    if _DATA_CACHE.get(cache_key) and not force_refresh:
        return _DATA_CACHE[cache_key]

    handles = _get_handles(opener)
    records = _build_records(handles, organism=organism)

    _DATA_CACHE[cache_key] = records
    return records


def _get_handles(opener) -> dict[str, Any]:
    """Extract named file handles from the opener's result."""
    if not opener or not opener.result:
        return {}
    if isinstance(opener.result, dict):
        return opener.result
    # Single-file opener — caller passed the wrong opener type; return empty
    return {}


def _iter_tsv(handle) -> Generator[tuple[str, str], None, None]:
    """Yield (col1, col2) pairs from a two-column KEGG TSV response."""
    if not handle:
        return
    if hasattr(handle, 'seek'):
        handle.seek(0)
    for raw_line in handle:
        if isinstance(raw_line, bytes):
            raw_line = raw_line.decode('utf-8', 'ignore')
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split('\t', 1)
        if len(parts) == 2:
            yield parts[0].strip(), parts[1].strip()


def _strip_prefix(value: str, prefix: str) -> str:
    """Remove a leading 'prefix:' string if present."""
    if value.startswith(prefix + ':'):
        return value[len(prefix) + 1:]
    return value


# ---------------------------------------------------------------------------
# Translation table builders
# ---------------------------------------------------------------------------

def _build_uniprot_map(handle) -> dict[str, list[str]]:
    """
    Build gene-to-UniProt mapping from /conv/uniprot/{organism}.

    Returns dict mapping bare gene ID (e.g. 'mmu:12385') to list of
    UniProt accessions.
    """
    result: dict[str, list[str]] = {}
    for gene_id, uniprot_raw in _iter_tsv(handle):
        uniprot = _strip_prefix(uniprot_raw, 'up')
        result.setdefault(gene_id, []).append(uniprot)
    return result


def _build_ec_map(handle) -> dict[str, list[str]]:
    """
    Build gene-to-EC mapping from /link/enzyme/{organism}.

    Returns dict mapping bare gene ID to list of EC number strings.
    """
    result: dict[str, list[str]] = {}
    for gene_id, ec_raw in _iter_tsv(handle):
        result.setdefault(gene_id, []).append(ec_raw)
    return result


def _build_rxn_map(handle) -> dict[str, list[str]]:
    """
    Build EC-to-reaction mapping from /link/reaction/enzyme.

    Returns dict mapping EC string to list of bare reaction IDs.
    """
    result: dict[str, list[str]] = {}
    for ec_id, rxn_raw in _iter_tsv(handle):
        rxn_id = _strip_prefix(rxn_raw, 'rn')
        result.setdefault(ec_id, []).append(rxn_id)
    return result


def _build_compound_map(
    list_handle,
    chebi_handle,
) -> dict[str, dict[str, Any]]:
    """
    Build compound name-to-identifiers mapping.

    Combines /list/compound (names) and /conv/chebi/compound (ChEBI IDs).

    Returns dict mapping each compound name to::

        {
            'kegg_id': 'cpd:C00001',
            'chebi':   'chebi:15377',
        }
    """
    # cpd:CXXXXX -> ChEBI ID
    chebi_lookup: dict[str, str] = {}
    for cpd_raw, chebi_raw in _iter_tsv(chebi_handle):
        cpd_id = _strip_prefix(cpd_raw, 'cpd')
        chebi_id = chebi_raw  # keep full "chebi:XXXXX" form
        chebi_lookup[cpd_id] = chebi_id

    # cpd:CXXXXX -> first (canonical) name
    result: dict[str, dict[str, Any]] = {}
    for cpd_raw, names_raw in _iter_tsv(list_handle):
        cpd_id = _strip_prefix(cpd_raw, 'cpd')
        for name in names_raw.split(';'):
            name = name.strip()
            if not name:
                continue
            if name not in result:
                result[name] = {
                    'kegg_id': f'cpd:{cpd_id}',
                    'chebi':   chebi_lookup.get(cpd_id),
                }
    return result


def _build_rxn_to_uniprot(
    uniprot_map: dict[str, list[str]],
    ec_map: dict[str, list[str]],
    rxn_map: dict[str, list[str]],
) -> dict[str, list[str]]:
    """
    Chain gene→EC→reaction to produce reaction→UniProt mapping.

    Returns dict mapping bare reaction ID to sorted list of UniProt accessions.
    """
    # EC → UniProt (via gene)
    ec_to_uniprot: dict[str, set[str]] = {}
    for gene_id, ec_list in ec_map.items():
        uniprots = uniprot_map.get(gene_id, [])
        for ec in ec_list:
            ec_to_uniprot.setdefault(ec, set()).update(uniprots)

    # Reaction → UniProt
    rxn_to_uniprot: dict[str, list[str]] = {}
    for ec, rxn_list in rxn_map.items():
        uniprots = ec_to_uniprot.get(ec, set())
        for rxn_id in rxn_list:
            rxn_to_uniprot.setdefault(rxn_id, set()).update(uniprots)  # type: ignore[arg-type]

    return {rxn_id: sorted(ups) for rxn_id, ups in rxn_to_uniprot.items()}


# ---------------------------------------------------------------------------
# Equation parser
# ---------------------------------------------------------------------------

def _strip_stoich(name: str) -> str:
    """Remove leading stoichiometric coefficients (2, n, (n+1), > 2, etc.)."""
    return _STOICH_RE.sub('', name).strip()


def _parse_equation(description: str) -> tuple[list[str], list[str]]:
    """
    Extract substrate and product name lists from a KEGG reaction description.

    KEGG descriptions follow the pattern::

        "Enzyme name; possibly more names; A + B <=> C + D"

    The equation is always after the last semicolon.  A greedy ``rsplit``
    on ``;`` correctly handles multi-semicolon descriptions without the
    non-greedy regex bug present in earlier versions of this pipeline.

    Returns
    -------
    tuple[list[str], list[str]]
        (substrates, products) — each a list of cleaned compound name strings.
    """
    equation = description.rsplit(';', 1)[-1].strip()

    match = _ARROW_RE.search(equation)
    if not match:
        return [], []

    sub_str  = equation[:match.start()].strip()
    prod_str = equation[match.end():].strip()

    substrates = [_strip_stoich(s.strip()) for s in sub_str.split(' + ') if s.strip()]
    products   = [_strip_stoich(p.strip()) for p in prod_str.split(' + ') if p.strip()]

    return substrates, products


def _resolve_compound(
    name: str,
    compound_map: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    """Look up a compound name and return its identifiers dict."""
    entry = compound_map.get(name)
    if entry:
        return {'name': name, 'kegg_id': entry['kegg_id'], 'chebi': entry['chebi']}
    return {'name': name, 'kegg_id': None, 'chebi': None}


# ---------------------------------------------------------------------------
# Record builder
# ---------------------------------------------------------------------------

def _build_records(
    handles: dict[str, Any],
    organism: str = 'mmu',
) -> dict[str, list[dict[str, Any]]]:
    """
    Parse all KEGG bulk-download handles and return a dict of record lists.

    Expected keys in ``handles``
    ----------------------------
    ``'conv_uniprot'``
        Response from ``/conv/uniprot/{organism}``
    ``'link_enzyme'``
        Response from ``/link/enzyme/{organism}``
    ``'link_reaction'``
        Response from ``/link/reaction/enzyme``
    ``'list_compound'``
        Response from ``/list/compound``
    ``'conv_chebi'``
        Response from ``/conv/chebi/compound``
    ``'list_reaction'``
        Response from ``/list/reaction``
    """
    uniprot_map  = _build_uniprot_map(handles.get('conv_uniprot'))
    ec_map       = _build_ec_map(handles.get('link_enzyme'))
    rxn_map      = _build_rxn_map(handles.get('link_reaction'))
    compound_map = _build_compound_map(
        handles.get('list_compound'),
        handles.get('conv_chebi'),
    )
    rxn_to_up    = _build_rxn_to_uniprot(uniprot_map, ec_map, rxn_map)

    reactions: list[dict[str, Any]] = []

    for rxn_raw, description in _iter_tsv(handles.get('list_reaction')):
        rxn_id = _strip_prefix(rxn_raw, 'rn')

        uniprot_ids = rxn_to_up.get(rxn_id)
        if not uniprot_ids:
            continue

        substrates, products = _parse_equation(description)

        reactions.append({
            'Reaction':  rxn_id,
            'UniProt':   ';'.join(uniprot_ids),
            'Substrate': [_resolve_compound(s, compound_map) for s in substrates],
            'Product':   [_resolve_compound(p, compound_map) for p in products],
        })

    return {'reactions': reactions}
