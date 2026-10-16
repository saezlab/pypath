"""Parser for the Rhea reaction database GUI API TSV export.

Yields one dict per master reaction row from the Rhea web API.

Participants follow the equation. The ``ChEBI identifier`` column lists each
compound once, in order of its first appearance in the equation, so a compound
on both sides (the cargo of every transport) appears in it once. Each equation
term is matched to its compound by that order, which restores the second
occurrence with its own role, membrane side and coefficient.
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
# UniProt's accession format: O, P and Q accessions have their own shape.
_UNIPROT_RE = re.compile(
    r'\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+)?\b'
)
_RHEA_ID_RE = re.compile(r'(?:RHEA:)?(\d+)')

_TRANSPORT_COMPARTMENT_RE = re.compile(r'\((in|out)\)$')
_COEFFICIENT_RE = re.compile(r'^\d+ ')
# A variable coefficient, e.g. ``n ATP`` or ``(n+1) H(+)``: no number to keep.
_VARIABLE_COEFFICIENT_RE = re.compile(r'^\(?[0-9+-]*n[0-9n+-]*\)? ')
_EXCLUDED_TERMS = frozenset({'H(+)', 'H2O'})


def _equation_terms(equation: str) -> list[tuple[str, str]]:
    """
    The terms of an equation in order, each with its side.

    Args:
        equation: Reaction equation string, e.g.
            ``sulfate(out) + ATP + H2O = sulfate(in) + ADP + phosphate + H(+)``.

    Returns:
        ``(role, term)`` pairs, ``role`` being ``'reactant'`` left of the
        separator and ``'product'`` right of it. Empty when no recognised
        separator is found.
    """

    for sep in _EQ_SEPARATORS:

        if sep in equation:

            left, right = equation.split(sep, 1)

            return [
                (role, term.strip())
                for role, side in (('reactant', left), ('product', right))
                for term in side.split(' + ')
            ]

    return []


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


def _participant_fields(
    equation: str,
    chebi_ids: list[str],
    chebi_names: list[str],
) -> dict[str, str]:
    """
    Build ``||``-delimited participant record fields, one entry per equation
    term.

    Args:
        equation: Reaction equation string.
        chebi_ids: ChEBI identifiers, each compound once, in order of first
            appearance in the equation.
        chebi_names: Display names in the same order as ``chebi_ids``.

    Returns:
        Dict with keys ``participant_chebi``, ``participant_display_name``,
        ``participant_role``, ``participant_compartment`` (``in``/``out`` for
        a transport, else empty) and ``participant_stoichiometry`` (the
        coefficient, ``1`` when the equation writes none, empty when it is
        variable), each a ``||``-joined string aligned by position.
    """

    names = chebi_names + [''] * max(len(chebi_ids) - len(chebi_names), 0)
    terms = _equation_terms(equation)
    compounds = list(dict.fromkeys(_strip_term(term) for _, term in terms))
    participants = []

    if terms and len(compounds) == len(chebi_ids):
        index = {compound: i for i, compound in enumerate(compounds)}
        for role, term in terms:
            i = index[_strip_term(term)]
            side = _TRANSPORT_COMPARTMENT_RE.search(term)
            coefficient = _COEFFICIENT_RE.match(term.strip())
            participants.append((
                chebi_ids[i],
                names[i],
                role,
                side.group(1) if side else '',
                coefficient.group(0).strip() if coefficient
                else '' if _VARIABLE_COEFFICIENT_RE.match(term.strip()) else '1',
            ))
    else:
        # The compound list does not match the equation (e.g. a term the
        # column leaves out): keep the compounds once, with roles by count.
        n_reactants = sum(role == 'reactant' for role, _ in terms)
        for i, chebi in enumerate(chebi_ids):
            role = 'reactant' if i < n_reactants else 'product'
            participants.append((chebi, names[i], role, '', ''))

    fields = (
        'participant_chebi',
        'participant_display_name',
        'participant_role',
        'participant_compartment',
        'participant_stoichiometry',
    )
    return {
        field: _LIST_DELIMITER.join(values)
        for field, values in zip(fields, zip(*participants) if participants else [()] * 5)
    }


def _raw(
    opener: Any,
    uniprot_opener: Any | None = None,
    **_kwargs: Any,
) -> Generator[dict, None, None]:
    """
    Parse the Rhea GUI API TSV and yield one dict per master reaction.

    Args:
        opener: File-like object for the downloaded TSV, as supplied by the
            ``Dataset`` machinery.
        **_kwargs: Ignored; present for framework compatibility.

    Yields:
        One dict per master reaction with keys: ``rhea_id``, ``equation``,
        ``ec``, ``pubmed``, ``go``, ``ecocyc``, ``metacyc``, ``kegg``,
        ``reactome``, ``is_transport``, ``participant_chebi``,
        ``participant_display_name``, ``participant_role``,
        ``participant_compartment`` and ``participant_stoichiometry``.

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
            'is_transport': is_transport,
            **_participant_fields(equation, chebi_ids, chebi_names),
        }
