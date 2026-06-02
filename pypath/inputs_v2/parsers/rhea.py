"""Parser for the Rhea reaction database GUI API TSV export.

Yields one dict per master reaction row from the Rhea web API.

Role assignment uses the equation structure: the ``ChEBI identifier`` column
lists ChEBI IDs in left-to-right equation order (reactants first, products
after the separator), so counting `` + ``-separated terms on each side of the
equation separator gives an exact reactant/product boundary in the ChEBI list.
"""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv


_LIST_DELIMITER = '||'

# Equation separators tried in decreasing specificity order.
# Master reactions use ' = '; directional variants use '=>' or '<=>'
# but those are not returned by the API.  Handled defensively.
_EQ_SEPARATORS = (' <=> ', ' => ', ' <= ', ' = ')

_GO_ID_RE = re.compile(r'GO:\d+')
_UNIPROT_RE = re.compile(
    r'\b[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:-\d+)?\b'
    r'|\b[A-NR-Z][0-9][A-Z0-9]{3}[0-9][A-Z0-9]{4}[0-9](?:-\d+)?\b'
)
_RHEA_ID_RE = re.compile(r'(?:RHEA:)?(\d+)')

_TRANSPORT_COMPARTMENT_RE = re.compile(r'\((in|out)\)$')
_COEFFICIENT_RE = re.compile(r'^\d+ ')
_EXCLUDED_TERMS = frozenset({'H(+)', 'H2O'})


def _count_equation_sides(equation: str) -> tuple[int, int]:
    """
    Count reactants and products by parsing the equation string.

    Args:
        equation: Reaction equation string, e.g.
            ``ATP + H2O = ADP + phosphate``.

    Returns:
        A ``(n_reactants, n_products)`` tuple.  Compounds within each side
        are separated by `` + `` and sides are split on the equation
        separator.  Returns ``(0, 0)`` if no recognised separator is found.
    """

    for sep in _EQ_SEPARATORS:

        if sep in equation:

            left, right = equation.split(sep, 1)

            return len(left.split(' + ')), len(right.split(' + '))

    return 0, 0


def _split_semicolons(raw: str) -> list[str]:

    return [x.strip() for x in raw.split(';') if x.strip()]


def _uniprot_ids(raw: str) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for match in _UNIPROT_RE.finditer(raw):
        accession = match.group(0)
        if accession in seen:
            continue
        seen.add(accession)
        out.append(accession)
    return out


def _normalize_rhea_id(raw: str) -> str:
    match = _RHEA_ID_RE.search(raw or '')
    return f'RHEA:{match.group(1)}' if match else (raw or '').strip()


def _uniprot_by_master_id(opener: Any | None) -> dict[str, list[str]]:
    if opener is None:
        return {}

    mapping: dict[str, list[str]] = {}
    seen: set[tuple[str, str]] = set()

    for row in iter_tsv(opener):
        rhea_id = _normalize_rhea_id(row.get('MASTER_ID') or '')
        uniprot_id = (row.get('ID') or '').strip()
        if not rhea_id or not uniprot_id:
            continue
        key = (rhea_id, uniprot_id)
        if key in seen:
            continue
        seen.add(key)
        mapping.setdefault(rhea_id, []).append(uniprot_id)

    return mapping


def _go_ids(raw: str) -> list[str]:
    """
    Extract GO accessions from a free-text GO annotation string.

    Args:
        raw: Raw Gene Ontology column value, e.g.
            ``GO:0006520 cellular amino acid metabolic process``.

    Returns:
        List of GO accession strings (e.g. ``['GO:0006520']``).
    """

    return _GO_ID_RE.findall(raw)


def _strip_term(term: str) -> str:
    """Strip leading numeric coefficient and trailing (in)/(out) compartment annotation."""
    term = _COEFFICIENT_RE.sub('', term.strip())
    return _TRANSPORT_COMPARTMENT_RE.sub('', term)


def _is_transport(equation: str) -> bool:
    """
    Return True if the equation represents a membrane transport reaction.

    Two conditions must both hold:

    1. At least one compound term ends with ``(in)`` or ``(out)``.
    2. At least one non-proton, non-water compound name appears unchanged
       on both sides of the equation (i.e. the molecule is translocated,
       not chemically transformed).

    Condition 2 correctly excludes oxidative-phosphorylation reactions such
    as cytochrome c oxidase, where protons appear on both sides but the
    electron-carrier compounds change oxidation state.  It also excludes
    pure proton pumps (e.g. V-type H⁺ ATPase) where only ``H(+)`` is shared.

    Args:
        equation: Reaction equation string from the Rhea API.

    Returns:
        True if the equation encodes a membrane transport event.
    """

    if '(in)' not in equation and '(out)' not in equation:
        return False

    for sep in _EQ_SEPARATORS:
        if sep in equation:
            left, right = equation.split(sep, 1)
            left_terms = [t.strip() for t in left.split(' + ')]
            right_terms = [t.strip() for t in right.split(' + ')]
            break
    else:
        return False

    if not any(_TRANSPORT_COMPARTMENT_RE.search(t) for t in left_terms + right_terms):
        return False

    left_names = {_strip_term(t) for t in left_terms} - _EXCLUDED_TERMS
    right_names = {_strip_term(t) for t in right_terms} - _EXCLUDED_TERMS
    return bool(left_names & right_names)


def _compartments_from_equation(equation: str, n_participants: int) -> list[str]:
    """
    Extract the membrane-side label for each participant in equation order.

    Args:
        equation: Reaction equation string with ``(in)``/``(out)`` annotations.
        n_participants: Expected number of participants (``len(chebi_ids)``).
            Used to pad or truncate the result to the correct length.

    Returns:
        List of compartment strings in equation order (reactants first, then
        products).  Entries are empty strings where no annotation is found.
    """

    for sep in _EQ_SEPARATORS:
        if sep in equation:
            left, right = equation.split(sep, 1)
            terms = (
                [t.strip() for t in left.split(' + ')] +
                [t.strip() for t in right.split(' + ')]
            )
            break
    else:
        return [''] * n_participants

    compartments = []

    for term in terms:
        m = _TRANSPORT_COMPARTMENT_RE.search(term)
        compartments.append(m.group(1) if m else '')

    if len(compartments) < n_participants:
        compartments += [''] * (n_participants - len(compartments))

    return compartments[:n_participants]


def _participant_fields(
    chebi_ids: list[str],
    chebi_names: list[str],
    n_reactants: int,
    compartments: list[str] | None = None,
) -> dict[str, str]:
    """
    Build ``||``-delimited participant record fields.

    Args:
        chebi_ids: Ordered list of ChEBI identifiers (reactants first,
            then products), as returned by the API.
        chebi_names: Display names in the same order as ``chebi_ids``.
            May be shorter; missing names are filled with empty strings.
        n_reactants: Number of leading entries in ``chebi_ids`` that are
            reactants; the remainder are treated as products.
        compartments: Compartment label per participant in the same order
            as ``chebi_ids``.  ``None`` or a short list is padded with
            empty strings.

    Returns:
        Dict with keys ``participant_chebi``, ``participant_display_name``,
        ``participant_role``, and ``participant_compartment``, each a
        ``||``-joined string aligned by position.
    """

    n_products = max(len(chebi_ids) - n_reactants, 0)
    roles = ['reactant'] * n_reactants + ['product'] * n_products

    names = chebi_names + [''] * max(len(chebi_ids) - len(chebi_names), 0)

    if compartments is None:
        compartments = [''] * len(chebi_ids)
    elif len(compartments) < len(chebi_ids):
        compartments = compartments + [''] * (len(chebi_ids) - len(compartments))

    return {
        'participant_chebi': _LIST_DELIMITER.join(chebi_ids),
        'participant_display_name': _LIST_DELIMITER.join(names[:len(chebi_ids)]),
        'participant_role': _LIST_DELIMITER.join(roles),
        'participant_compartment': _LIST_DELIMITER.join(compartments[:len(chebi_ids)]),
    }


def _raw(
    opener: Any,
    data_type: str = 'reactions',
    uniprot_opener: Any | None = None,
    **_kwargs: Any,
) -> Generator[dict, None, None]:
    """
    Parse the Rhea GUI API TSV and yield one dict per master reaction.

    Args:
        opener: File-like object for the downloaded TSV, as supplied by the
            ``Dataset`` machinery.
        data_type: Controls which reactions are yielded.  One of
            ``'reactions'`` (all), ``'metabolic_reactions'`` (non-transport
            only), or ``'transport_reactions'`` (transport only).
        **_kwargs: Ignored; present for framework compatibility.

    Yields:
        One dict per master reaction with keys: ``rhea_id``, ``equation``,
        ``ec``, ``pubmed``, ``go``, ``ecocyc``, ``metacyc``, ``kegg``,
        ``reactome``, ``participant_chebi``, ``participant_display_name``,
        ``participant_role``, and ``participant_compartment``.

    Note:
        Expected TSV columns (set by the configured API URL):
        ``Reaction identifier``, ``Equation``, ``ChEBI name``,
        ``ChEBI identifier``, ``EC number``, ``Enzymes``,
        ``Gene Ontology``, ``PubMed``, ``Cross-reference (EcoCyc)``,
        ``Cross-reference (MetaCyc)``, ``Cross-reference (KEGG)``,
        ``Cross-reference (Reactome)``, ``Cross-reference (M-CSA)``.
    """

    uniprot_by_rhea_id = _uniprot_by_master_id(uniprot_opener)

    for row in iter_tsv(opener):

        rhea_id = (row.get('Reaction identifier') or '').strip()
        equation = (row.get('Equation') or '').strip()

        if not rhea_id or not equation:
            continue

        is_transport = _is_transport(equation)

        if data_type == 'metabolic_reactions' and is_transport:
            continue
        if data_type == 'transport_reactions' and not is_transport:
            continue

        chebi_ids = _split_semicolons(row.get('ChEBI identifier') or '')
        chebi_names = _split_semicolons(row.get('ChEBI name') or '')

        ec_raw = (row.get('EC number') or '').strip()
        go_raw = (row.get('Gene Ontology') or '').strip()
        pubmed = (row.get('PubMed') or '').strip()
        uniprot = ';'.join(
            uniprot_by_rhea_id.get(rhea_id)
            or _uniprot_ids(row.get('UniProt') or '')
        )

        ecocyc = (row.get('Cross-reference (EcoCyc)') or '').strip()
        metacyc = (row.get('Cross-reference (MetaCyc)') or '').strip()
        kegg = (row.get('Cross-reference (KEGG)') or '').strip()
        reactome = (row.get('Cross-reference (Reactome)') or '').strip()

        n_reactants, _ = _count_equation_sides(equation)

        compartments = (
            _compartments_from_equation(equation, len(chebi_ids))
            if is_transport else None
        )

        yield {
            'rhea_id': rhea_id,
            'equation': equation,
            'ec': ec_raw.removeprefix('EC:'),
            'pubmed': pubmed,
            'uniprot': uniprot,
            'go': ';'.join(_go_ids(go_raw)),
            'ecocyc': ecocyc,
            'metacyc': metacyc,
            'kegg': kegg,
            'reactome': reactome,
            **_participant_fields(chebi_ids, chebi_names, n_reactants, compartments),
        }
