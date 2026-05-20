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


def _participant_fields(
    chebi_ids: list[str],
    chebi_names: list[str],
    n_reactants: int,
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

    Returns:
        Dict with keys ``participant_chebi``, ``participant_display_name``,
        and ``participant_role``, each a ``||``-joined string aligned by
        position.
    """

    n_products = max(len(chebi_ids) - n_reactants, 0)
    roles = ['reactant'] * n_reactants + ['product'] * n_products

    names = chebi_names + [''] * max(len(chebi_ids) - len(chebi_names), 0)

    return {
        'participant_chebi': _LIST_DELIMITER.join(chebi_ids),
        'participant_display_name': _LIST_DELIMITER.join(names[:len(chebi_ids)]),
        'participant_role': _LIST_DELIMITER.join(roles),
    }


def _raw(opener: Any, **_kwargs: Any) -> Generator[dict, None, None]:
    """
    Parse the Rhea GUI API TSV and yield one dict per master reaction.

    Args:
        opener: File-like object for the downloaded TSV, as supplied by the
            ``Dataset`` machinery.
        **_kwargs: Ignored; present for framework compatibility.

    Yields:
        One dict per master reaction with keys: ``rhea_id``, ``equation``,
        ``ec``, ``pubmed``, ``go``, ``ecocyc``, ``metacyc``, ``kegg``,
        ``reactome``, ``participant_chebi``, ``participant_display_name``,
        and ``participant_role``.

    Note:
        Expected TSV columns (set by the configured API URL):
        ``Reaction identifier``, ``Equation``, ``ChEBI name``,
        ``ChEBI identifier``, ``EC number``, ``Enzymes``,
        ``Gene Ontology``, ``PubMed``, ``Cross-reference (EcoCyc)``,
        ``Cross-reference (MetaCyc)``, ``Cross-reference (KEGG)``,
        ``Cross-reference (Reactome)``, ``Cross-reference (M-CSA)``.
    """

    for row in iter_tsv(opener):

        rhea_id = (row.get('Reaction identifier') or '').strip()
        equation = (row.get('Equation') or '').strip()

        if not rhea_id or not equation:
            continue

        chebi_ids = _split_semicolons(row.get('ChEBI identifier') or '')
        chebi_names = _split_semicolons(row.get('ChEBI name') or '')

        ec_raw = (row.get('EC number') or '').strip()
        go_raw = (row.get('Gene Ontology') or '').strip()
        pubmed = (row.get('PubMed') or '').strip()

        ecocyc = (row.get('Cross-reference (EcoCyc)') or '').strip()
        metacyc = (row.get('Cross-reference (MetaCyc)') or '').strip()
        kegg = (row.get('Cross-reference (KEGG)') or '').strip()
        reactome = (row.get('Cross-reference (Reactome)') or '').strip()

        # Count reactants to assign participant roles
        n_reactants, _ = _count_equation_sides(equation)

        yield {
            'rhea_id': rhea_id,
            'equation': equation,
            'ec': ec_raw.removeprefix('EC:'),
            'pubmed': pubmed,
            'go': ';'.join(_go_ids(go_raw)),
            'ecocyc': ecocyc,
            'metacyc': metacyc,
            'kegg': kegg,
            'reactome': reactome,
            **_participant_fields(chebi_ids, chebi_names, n_reactants),
        }
