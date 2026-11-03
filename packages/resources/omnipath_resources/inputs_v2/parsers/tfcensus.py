"""Read the census of human transcription factors (Vaquerizas et al. 2009).

Supplementary Table 2 of "A census of human transcription factors: function,
expression and evolution" (Nature Reviews Genetics 10, 252) is a tab-separated
text file. A preamble of title, legend and download notes comes first, then the
line ``--Table start--``, the header line and one line per Ensembl gene. The
legend says: "List of TF-encoding loci classified as 'a', 'b', 'c', 'x', or
selected from other databases (classified as 'other'). [...] For TFs included
in the repertoire, the list also contains accompanying information including:
Ensembl gene IDs (release 51), HGNC identifiers, IPI IDs, associated
DNA-binding Interpro domains and families, and tissue specificity if any."

Each row is returned as a dictionary:

    census_class   ``a``, ``b``, ``c``, ``x`` or ``other``
    ensembl        Ensembl gene ID (``ENSG`` and eleven digits)
    ipi            IPI protein IDs
    interpro_dbd   InterPro IDs of the DNA-binding domains
    interpro_family  InterPro IDs of the DNA-binding families
    symbol         HGNC symbol; None where the table has none (153 rows)
    tissues        tissue-specificity labels as the table writes them
                   (``general``, ``whole.blood.JJV``, ``fetal.brain`` ...)
"""

from __future__ import annotations

from collections.abc import Generator
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv, read_text

MARKER = '--Table start--'
COLUMNS = (
    'Class',
    'Ensembl ID',
    'IPI ID',
    'Interpro DBD',
    'Interpro DNA-binding family',
    'HGNC symbol',
    'Tissue-specificity',
)


def _list(value: str | None) -> list[str]:
    """A cell of values separated by semicolons; an empty cell is no values."""
    return [item.strip() for item in (value or '').split(';') if item.strip()]


def iter_census(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per gene of the table.

    A file without the table marker or with another header is an error; so is
    a gene row without class or Ensembl gene ID.
    """
    text = read_text(opener)
    start = text.find(MARKER)
    if start < 0:
        raise ValueError(f'TF census file lacks the {MARKER!r} line')
    table = text[start + len(MARKER):].lstrip('\r\n')
    for row in iter_tsv(SimpleNamespace(result=table), required=COLUMNS):
        if not row['Class'] or not row['Ensembl ID']:
            raise ValueError(f'Malformed TF census row {row!r}')
        yield {
            'census_class': row['Class'],
            'ensembl': row['Ensembl ID'],
            'ipi': _list(row['IPI ID']),
            'interpro_dbd': _list(row['Interpro DBD']),
            'interpro_family': _list(row['Interpro DNA-binding family']),
            'symbol': row['HGNC symbol'] or None,
            'tissues': _list(row['Tissue-specificity']),
        }


def iter_census_list(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary for the whole table: its ``classes``, sorted.

    The table is read, so that a changed or empty file is an error.
    """
    yield {'classes': sorted({row['census_class'] for row in iter_census(opener)})}


def iter_census_interpro(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield a dictionary per gene and InterPro entry of its two InterPro columns.

    The dictionary has the fields of ``iter_census`` (except the two lists of
    InterPro entries), ``interpro`` (the entry) and ``interpro_column``
    (``DNA-binding domain`` or ``DNA-binding family``: the column of the table
    that lists it).
    """
    columns = {
        'interpro_dbd': 'DNA-binding domain',
        'interpro_family': 'DNA-binding family',
    }
    for row in iter_census(opener):
        gene = {k: v for k, v in row.items() if k not in columns}
        for field, column in columns.items():
            for entry in row[field]:
                yield {**gene, 'interpro': entry, 'interpro_column': column}
