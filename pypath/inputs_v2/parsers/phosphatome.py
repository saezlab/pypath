"""Read the human phosphatases of Chen et al. 2017 (Phosphatome.Net).

The supplement of "Genomics and evolution of protein phosphatases" (Science
Signaling 10, eaag1796) is a zip archive with one workbook, ``Tables S1 to
S23``. Two of its sheets are read:

- ``S1-HumanPhosphatases``: one row per human phosphatase gene (264), with its
  fold, family and, for 189 of them, subfamily, and the columns of the
  "Substrate type" group: ``Protein substrates?``, ``Non-protein substrates?``
  and ``Catalytic activity?``. ``Included in Protein Phosphatome?`` tells
  whether the gene is one of the protein phosphatases.
- ``S2-Fold``: one row per fold with its SCOP identifier and name.

The readers below give the classification as a tree, one dictionary per fold
(the same shape as the one of the kinase.com reader), the substrate-type
columns of each gene, and the values that those columns take.
"""

from __future__ import annotations

from collections.abc import Generator
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel, read_bytes

SYMBOL = 'Entrez Gene Symbol'
ENTREZ = 'Entrez Gene ID'
FOLD = 'Fold'
FAMILY = 'Family'
SUBFAMILY = 'Subfamily'
PROTEIN_PHOSPHATOME = 'Included in Protein Phosphatome?'
PROTEIN_SUBSTRATES = 'Protein substrates?'
NON_PROTEIN_SUBSTRATES = 'Non-protein substrates?'
CATALYTIC_ACTIVITY = 'Catalytic activity?'
GENE_COLUMNS = (
    SYMBOL,
    ENTREZ,
    FOLD,
    FAMILY,
    SUBFAMILY,
    PROTEIN_PHOSPHATOME,
    PROTEIN_SUBSTRATES,
    NON_PROTEIN_SUBSTRATES,
    CATALYTIC_ACTIVITY,
)
FOLD_COLUMNS = (FOLD, 'SCOP ID and name')
# The columns of yes/no flags, and the columns whose values are terms.
FLAGS = (PROTEIN_PHOSPHATOME, PROTEIN_SUBSTRATES, NON_PROTEIN_SUBSTRATES)
CLASS_COLUMNS = (*FLAGS, CATALYTIC_ACTIVITY)
LEVELS = ('fold', 'family', 'subfamily')


def _sheets(opener) -> tuple[list[dict[str, Any]], dict[str, str | None]]:
    """The gene rows of S1 and the SCOP names of the folds of S2.

    S1's first row groups the columns ("Classification", "Substrate type"),
    so its header is the second row; S2's header is its first row. A sheet
    without a column that is read, or without rows, is an error.
    """
    data = read_bytes(opener)
    if not data:
        raise ValueError('The Phosphatome supplement is empty')
    genes = list(
        iter_excel(
            SimpleNamespace(result=data),
            sheet='S1-HumanPhosphatases',
            header_row=2,
            required=GENE_COLUMNS,
        )
    )
    folds = {
        row[FOLD]: None if row['SCOP ID and name'] == 'N/A' else row['SCOP ID and name']
        for row in iter_excel(
            SimpleNamespace(result=data),
            sheet='S2-Fold',
            header_row=1,
            required=FOLD_COLUMNS,
        )
    }
    return genes, folds


def _entrez(row: dict[str, Any]) -> str:
    """The Entrez Gene ID; the sheet stores it as a number."""
    value = row[ENTREZ]
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    if not str(value).isdigit() or not row[SYMBOL]:
        raise ValueError(f'Malformed Phosphatome gene row {row!r}')
    return str(value)


def iter_phosphatase_tree(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one tree per fold, in the order of the sheet.

    A tree has ``level`` (``fold``, ``family`` or ``subfamily``), ``name``,
    ``description`` (the SCOP name, for a fold), ``children`` (the trees of the
    next level) and ``members``: the genes classified at this level, each with
    ``entrez`` and ``symbol``. A gene without a subfamily is a member of its
    family. A fold of the gene sheet that S2 does not list is an error.
    """
    genes, folds = _sheets(opener)
    roots: dict[str, dict[str, Any]] = {}
    for row in genes:
        names = (row[FOLD], row[FAMILY], row[SUBFAMILY])
        path = tuple(name for name in names if name is not None)
        if len(path) < 2 or path != names[: len(path)]:
            raise ValueError(f'Malformed Phosphatome classification in {row!r}')
        if row[FOLD] not in folds:
            raise ValueError(f'Phosphatome fold {row[FOLD]!r} is not in sheet S2')
        level_nodes = roots
        node = None
        for level, name in zip(LEVELS, path):
            node = level_nodes.setdefault(
                name,
                {
                    'level': level,
                    'name': name,
                    'description': folds[name] if level == 'fold' else None,
                    'children': {},
                    'members': [],
                },
            )
            level_nodes = node['children']
        node['members'].append({'entrez': _entrez(row), 'symbol': row[SYMBOL]})

    def listed(node: dict[str, Any]) -> dict[str, Any]:
        return {**node, 'children': [listed(c) for c in node['children'].values()]}

    for node in roots.values():
        yield listed(node)


def _flag(row: dict[str, Any], column: str) -> str:
    value = row[column]
    if value not in {'Yes', 'No'}:
        raise ValueError(f'Unknown Phosphatome value {value!r} in {column!r}')
    return value


def iter_phosphatase_classes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per gene with the columns of its classes.

    ``entrez`` and ``symbol``, and the values as the sheet writes them:
    ``protein_phosphatome``, ``protein_substrates`` and
    ``non_protein_substrates`` (``Yes`` or ``No``), and ``catalytic_activity``.
    """
    genes, _ = _sheets(opener)
    for row in genes:
        yield {
            'entrez': _entrez(row),
            'symbol': row[SYMBOL],
            'protein_phosphatome': _flag(row, PROTEIN_PHOSPHATOME),
            'protein_substrates': _flag(row, PROTEIN_SUBSTRATES),
            'non_protein_substrates': _flag(row, NON_PROTEIN_SUBSTRATES),
            'catalytic_activity': row[CATALYTIC_ACTIVITY],
        }


def iter_class_values(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield each value that a class column takes, once, in the order seen.

    A dictionary has ``column`` (the sheet's name of the column) and ``value``.
    """
    genes, _ = _sheets(opener)
    seen = set()
    for row in genes:
        for column in CLASS_COLUMNS:
            if (column, row[column]) not in seen:
                seen.add((column, row[column]))
                yield {'column': column, 'value': row[column]}
