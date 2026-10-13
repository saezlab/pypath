"""
Parsers of the matrisome master lists.

The Matrisome Portal (Naba lab, University of Illinois Chicago) offers the list of
the genes of the matrisome, the ensemble of the extracellular matrix (ECM) proteins
and the proteins associated with it, for each species as a spreadsheet. The sheet
has a title row, then a header row, and a row for each gene:

    Matrisome Division  Matrisome Category  Gene Symbol  Gene Name  Synonyms
    HGNC_IDs  HGNC_IDs Links  UniProt_IDs  Refseq_IDs  Notes

The human list names its genes by HGNC id, the mouse list by MGI id (``MGI_IDs``,
``MGI_IDs Links``). Synonyms are separated by ``|``, UniProt and RefSeq accessions
by ``:``. A division has several categories: the core matrisome has collagens, ECM
glycoproteins and proteoglycans; the matrisome-associated proteins are ECM
regulators, ECM-affiliated proteins and secreted factors. ``Notes`` holds a free
comment on some genes (``Semaphorin``, ``Basement Membrane``).
"""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

# The species of the lists: the column of the gene id, the pattern the id follows
# and the NCBI taxon.
SPECIES = {
    'human': ('HGNC_IDs', re.compile(r'\d+'), '9606'),
    'mouse': ('MGI_IDs', re.compile(r'MGI:\d+'), '10090'),
}
COLUMNS = (
    'Matrisome Division',
    'Matrisome Category',
    'Gene Symbol',
    'Gene Name',
    'Synonyms',
    'UniProt_IDs',
    'Refseq_IDs',
    'Notes',
)
# The title row is the first row of the sheet.
HEADER_ROW = 2


def _text(value: Any) -> str | None:
    # A whole number is a numeric cell of the workbook (the HGNC ids of the
    # human list), not a decimal.
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    text = str(value).strip() if value is not None else ''
    return text or None


def _list(value: Any, separator: str) -> list[str]:
    """
    The distinct items of a cell, in order.

    ``-`` stands for none, and so does the text ``None``, which the sheets write in
    the synonyms of the genes that have none (72 human and 78 mouse genes), as they
    do in the notes.
    """
    text = _text(value)
    if text is None or text in ('-', 'None'):
        return []
    return list(dict.fromkeys(item.strip() for item in text.split(separator) if item.strip()))


def _note(value: Any) -> str | None:
    """The comment of a gene; the sheet writes the text ``None`` for no comment."""
    text = _text(value)
    return None if text in (None, 'None') else text


def iter_genes(
    opener, species: str = 'human', **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """
    One dictionary per gene of a master list.

    - ``division``, ``category``: the two levels of the matrisome classification.
    - ``symbol``, ``name``: the gene symbol and name; ``synonyms``: a list (empty
      where the sheet writes ``None``).
    - ``gene_id``: the HGNC id (digits) of the human list or the MGI id
      (``MGI:`` and digits) of the mouse list.
    - ``taxon``: the NCBI taxon id of the list.
    - ``uniprot``, ``refseq``: lists of accessions.
    - ``note``: the free comment, or None. The sheet writes the text ``None`` in
      the rows without a comment, and some comments with a trailing space
      (``FACIT``, ``FACIT ``): they are the same comment.

    A gene without a division, category, symbol or well-formed gene id, a symbol
    that is listed twice, or an empty sheet is an error.
    """
    id_column, id_pattern, taxon = SPECIES[species]
    seen = set()
    for row in iter_excel(
        opener, header_row=HEADER_ROW, required=(*COLUMNS, id_column)
    ):
        symbol = _text(row['Gene Symbol'])
        division = _text(row['Matrisome Division'])
        category = _text(row['Matrisome Category'])
        gene_id = _text(row[id_column])
        if not (symbol and division and category and gene_id) or not id_pattern.fullmatch(
            gene_id
        ):
            raise ValueError(f'Matrisome {species} row not understood: {row!r}')
        if symbol in seen:
            raise ValueError(f'Matrisome {species} list has {symbol} twice')
        seen.add(symbol)
        yield {
            'division': division,
            'category': category,
            'symbol': symbol,
            'name': _text(row['Gene Name']),
            'synonyms': _list(row['Synonyms'], '|'),
            'gene_id': gene_id,
            'taxon': taxon,
            'uniprot': _list(row['UniProt_IDs'], ':'),
            'refseq': _list(row['Refseq_IDs'], ':'),
            'note': _note(row['Notes']),
        }
    if not seen:
        raise ValueError(f'Matrisome {species} list has no rows')


def iter_terms(
    opener, mouse_download, force_refresh: bool = False, **kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """
    The divisions and categories of the human and mouse lists.

    A division is yielded as ``{'kind': 'division', 'term': ..., 'parent': None}``
    and a category as ``{'kind': 'category', 'term': ..., 'parent': <division>}``,
    divisions first, in the order of the sheet. A category that the lists put in
    two divisions is an error.
    """
    parents = {}
    divisions = {}
    for species, source in (('human', opener), ('mouse', None)):
        if source is None:
            source = mouse_download.open(force_refresh=force_refresh, **kwargs)
        for gene in iter_genes(source, species=species):
            divisions.setdefault(gene['division'], None)
            if parents.setdefault(gene['category'], gene['division']) != gene['division']:
                raise ValueError(
                    f'Matrisome category {gene["category"]!r} is in two divisions'
                )
    for division in divisions:
        yield {'kind': 'division', 'term': division, 'parent': None}
    for category, division in parents.items():
        yield {'kind': 'category', 'term': category, 'parent': division}
