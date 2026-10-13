"""Read Supplementary Table S1 of Zhong et al. 2015, the human cell adhesion molecules.

Zhong et al. updated the set of human genes that encode cell adhesion
molecules (CAMs) by searching databases, reading the literature and annotating,
and classified them by function. The table is a tab-separated file (OmniPath's
copy of the supplement) with one row per gene:

    symbol  entrezID  CAM type  Inter Pro  GO  NCBI  Sum  name
    CDH1    999       i         1          1   1     3    cadherin 1

``CAM type`` is a code of the functional subclass. ``Inter Pro``, ``GO`` and
``NCBI`` are 1 where the gene was found by the search of that database and
``Sum`` is their sum. Four more columns, without a name, are empty. The name
has trailing spaces.
"""

from __future__ import annotations

from collections.abc import Generator
from functools import partial
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = ('symbol', 'entrezID', 'CAM type', 'Inter Pro', 'GO', 'NCBI', 'Sum', 'name')
SEARCHES = ('Inter Pro', 'GO', 'NCBI')

# Written in place of an Entrez Gene ID for three genes (CD24, KIR2DS2, KIR3DL3).
# No other Entrez Gene ID of the table is above 728577, and the same number stands for
# three different genes, so it identifies none of them: these genes keep the symbol alone.
PLACEHOLDER_ENTREZ = '100000000'

_table = partial(iter_tsv, quotechar=None, required=COLUMNS)


def _flag(row: dict[str, Any], column: str) -> bool:
    if row[column] not in ('0', '1'):
        raise ValueError(f'{row["symbol"]}: {column} is {row[column]!r}')
    return row[column] == '1'


def iter_genes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One record per gene: ``symbol``, ``entrez`` (None for the placeholder),
    ``cam_type`` (the code), ``name``, ``interpro``, ``go``, ``ncbi`` (booleans:
    found by the search) and ``sources`` (their number, as the file's ``Sum``)."""
    seen: set[str] = set()
    for row in _table(opener):
        symbol = row['symbol']
        if symbol in seen:
            raise ValueError(f'Gene {symbol!r} occurs twice')
        seen.add(symbol)
        entrez = row['entrezID']
        if not entrez.isdigit():
            raise ValueError(f'{symbol}: Entrez Gene ID {entrez!r}')
        flags = {column: _flag(row, column) for column in SEARCHES}
        if int(row['Sum']) != sum(flags.values()):
            raise ValueError(f'{symbol}: Sum {row["Sum"]!r} is not the sum of the searches')
        yield {
            'symbol': symbol,
            'entrez': None if entrez == PLACEHOLDER_ENTREZ else entrez,
            'cam_type': row['CAM type'],
            'name': row['name'].strip() or None,
            'interpro': flags['Inter Pro'],
            'go': flags['GO'],
            'ncbi': flags['NCBI'],
            'sources': int(row['Sum']),
        }
    if not seen:
        raise ValueError('The CAM table has no genes')


def iter_cam_types(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The codes of the CAM types of the file, each once, in order of appearance."""
    seen: set[str] = set()
    for gene in iter_genes(opener):
        if gene['cam_type'] not in seen:
            seen.add(gene['cam_type'])
            yield {'cam_type': gene['cam_type']}
