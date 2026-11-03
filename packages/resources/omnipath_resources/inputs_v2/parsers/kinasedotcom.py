"""Read the human kinome table of kinase.com.

``Kincat_Hsap.08.02.xls`` is the table behind Manning et al. 2002, kept up to
date on kinase.com: one row per human protein kinase gene or kinase pseudogene,
with the kinase's group, family and subfamily (``AGC`` > ``Akt``;
``AGC`` > ``DMPK`` > ``GEK``). A kinase without a subfamily sits in its family. A kinase with two kinase domains has
the classification of the second domain in the ``KD2_*`` columns.

The reader returns the classification as a tree, one dictionary per group:

    level      ``group``, ``family`` or ``subfamily``
    name       the name as the table writes it
    children   the dictionaries of the next level inside this one
    members    the kinases classified at this level, each with
               ``sk`` (the kinase.com accession, ``SK018``), ``name`` (the
               kinase.com name), ``entrez`` and ``symbol`` (Entrez Gene ID and
               symbol; None for the two genes the table gives none), and
               ``kinase_domain`` (``1`` or ``2``: which domain's
               classification places the kinase here)

Rows that the ``Pseudogene?`` column marks ``Y`` (pseudogenes, 102 rows) or
``R`` (2 rows) are not kinase genes and are left out.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

NAME = 'Name'
ACCESSION = 'SK'
PSEUDOGENE = 'Pseudogene?'
ENTREZ = 'Entrez_GeneID'
SYMBOL = 'Entrez_Symbol'
# Column names of the classification of the first and second kinase domain.
CLASSIFICATION = {
    '1': ('Group', 'Family', 'Subfamily'),
    '2': ('KD2_Group', 'KD2_Family', 'KD2_Subfamily'),
}
COLUMNS = (
    NAME,
    ACCESSION,
    PSEUDOGENE,
    ENTREZ,
    SYMBOL,
    *CLASSIFICATION['1'],
    *CLASSIFICATION['2'],
)
LEVELS = ('group', 'family', 'subfamily')
# ``Pseudogene?``: N is a kinase gene; Y a pseudogene; R is not explained in
# the table (its two rows are named ``-rs``).
GENES = 'N'
NOT_GENES = frozenset({'Y', 'R'})


def _entrez(value: Any, row: dict[str, Any]) -> str | None:
    """The Entrez Gene ID, which the sheet stores as a number (207.0)."""
    if value is None:
        return None
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    if not str(value).isdigit():
        raise ValueError(f'Malformed kinase.com Entrez Gene ID in {row!r}')
    return str(value)


def _path(row: dict[str, Any], domain: str) -> tuple[str, ...]:
    """Group, family and subfamily of a row, as far as the table gives them."""
    names = tuple(row[column] for column in CLASSIFICATION[domain])
    path = tuple(name for name in names if name is not None)
    # A level is not given below one that is missing.
    if path != names[: len(path)] or (domain == '1' and not path):
        raise ValueError(f'Malformed kinase.com classification in {row!r}')
    return path


def iter_kinase_tree(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one tree per kinase group, in the order of the table.

    A pseudogene row is left out; an unknown ``Pseudogene?`` value, a
    classification with a gap and a repeated kinase.com accession are errors.
    """
    roots: dict[str, dict[str, Any]] = {}
    seen = set()
    for row in iter_excel(opener, sheet='Sheet1', required=COLUMNS):
        flag = row[PSEUDOGENE]
        if flag in NOT_GENES:
            continue
        if flag != GENES:
            raise ValueError(f'Unknown kinase.com Pseudogene? value in {row!r}')
        if row[ACCESSION] in seen:
            raise ValueError(f'Repeated kinase.com accession in {row!r}')
        seen.add(row[ACCESSION])
        for domain in CLASSIFICATION:
            path = _path(row, domain) if domain == '1' or row[CLASSIFICATION[domain][0]] else ()
            if not path:
                continue
            level_nodes = roots
            node = None
            for level, name in zip(LEVELS, path):
                node = level_nodes.setdefault(
                    name,
                    {'level': level, 'name': name, 'children': {}, 'members': []},
                )
                level_nodes = node['children']
            node['members'].append(
                {
                    'sk': row[ACCESSION],
                    'name': row[NAME],
                    'entrez': _entrez(row[ENTREZ], row),
                    'symbol': row[SYMBOL],
                    'kinase_domain': domain,
                }
            )
    if not roots:
        raise ValueError('kinase.com table has no kinase genes')

    def listed(node: dict[str, Any]) -> dict[str, Any]:
        return {**node, 'children': [listed(c) for c in node['children'].values()]}

    for node in roots.values():
        yield listed(node)
