"""
Parsers of the PROGENy model matrix.

``model_human_full.rda`` is an R data file with one data frame,
``model_human_full``, of the full human model: a row for each pair of a gene
and one of the 14 pathways, with the columns

    gene       HGNC symbol
    pathway    name of the PROGENy pathway
    weight     "z-scores for a given gene" (the package documentation)
    p.value    "significance of gene in pathway" (the package documentation)

The model has 22,479 genes; seven pathways have all of them, VEGF 22,472, TGFb
22,470, PI3K 21,733, Trail 20,514 and Androgen, Estrogen and WNT 9,867 each, so
the data frame has 274,143 rows.
"""

from __future__ import annotations

from collections.abc import Generator
import functools
from typing import Any

from pypath.inputs_v2.parsers.base import iter_rdata

MODEL = 'model_human_full'
COLUMNS = ('gene', 'pathway', 'weight', 'p.value')

_read = functools.partial(iter_rdata, name=MODEL)


def iter_model(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The rows of the model: ``gene``, ``pathway``, ``weight``, ``p_value``.

    A table without one of the columns, with a missing value or without rows
    raises; the file has no missing values.
    """
    count = 0
    for row in _read(opener):
        missing = [column for column in COLUMNS if row.get(column) is None]
        if missing:
            raise ValueError(f'PROGENy model row lacks {missing}: {row}')
        count += 1
        yield {
            'gene': row['gene'],
            'pathway': row['pathway'],
            'weight': row['weight'],
            'p_value': row['p.value'],
        }
    if not count:
        raise ValueError('PROGENy model has no rows')


def iter_pathways(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The pathways of the model, each once, in order of first appearance."""
    seen: dict[str, None] = {}
    for row in iter_model(opener):
        seen.setdefault(row['pathway'])
    for pathway in seen:
        yield {'pathway': pathway}
