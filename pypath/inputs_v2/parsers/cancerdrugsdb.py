"""
Reader for the CancerDrugs_DB table.

The tab separated table is read by the shared TSV reader;
:func:`iter_drug_targets` is the CancerDrugs_DB interpretation of it: it decodes
the identifier cells, which hold HTML links, and the lists, and gives one record
per drug and target.
"""

from __future__ import annotations

from collections.abc import Generator
import html
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['iter_drug_targets']


# Columns of the table, in source order. The header row ends with a tab, so it
# has one more, nameless column, which is empty on every line.
COLUMNS = (
    'Product',
    'EMA',
    'FDA',
    'EN',
    'Other',
    'WHO',
    'Year',
    'Generic',
    'DrugBank ID',
    'ATC',
    'ChEMBL',
    'Indications',
    'Targets',
    'Last Update',
)
# An identifier cell: "<a href="URL">ID</a>", several separated by "; ". The
# link text of a few ChEMBL cells has a leading space.
_LINK_REGEX = re.compile(r'^<a href="([^"]*)">([^<]*)</a>$')
_DRUGBANK_REGEX = re.compile(r'^DB\d{5}$')
_CHEMBL_REGEX = re.compile(r'^CHEMBL\d+$')
# The text the DrugBank cell holds, in a link, for a drug DrugBank does not list.
_NO_DRUGBANK = 'Not found in DrugBank'


def _identifiers(cell: str, pattern: re.Pattern[str], line: int, column: str) -> list[str]:
    """The identifiers of a cell of links, in cell order.

    The placeholder for a drug that DrugBank does not list is dropped. Any
    other text that is not an identifier of the column raises ``ValueError``,
    as does a link whose address does not end in its text.
    """
    result: list[str] = []
    for part in cell.split(';'):
        part = part.strip()
        if not part:
            continue
        link = _LINK_REGEX.match(part)
        if not link:
            raise ValueError(f'CancerDrugsDB line {line}: {column} is not a link: {part!r}')
        text = link[2].strip()
        if text == _NO_DRUGBANK and column == 'DrugBank ID':
            continue
        if not pattern.match(text) or not link[1].strip().endswith(text):
            raise ValueError(f'CancerDrugsDB line {line}: unexpected {column}: {part!r}')
        if text not in result:
            result.append(text)
    return result


def _list(cell: str) -> list[str]:
    """Items of a cell separated by semicolons, in cell order, without repeats."""
    result: list[str] = []
    for item in cell.split(';'):
        item = item.strip()
        if item and item not in result:
            result.append(item)
    return result


def iter_drug_targets(opener, **_kwargs: object) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per drug and target listed for it.

    A drug without targets yields nothing. Output fields; the cells the source
    gives as text are kept as written, trimmed:

    - ``product``: the drug's name.
    - ``ema``, ``fda``, ``en``, ``other``, ``who``, ``year``, ``generic``,
      ``indications``, ``last_update``: the columns of those names (approval by
      the EMA, the FDA, a European national agency, an agency named in
      ``other``, WHO essential medicine, the first approval year, generic
      flag, indications with HTML character references decoded, and the date
      of the last update).
    - ``drugbank_ids``, ``chembl_ids``: the identifiers of the DrugBank and
      ChEMBL cells; a list, empty where the cell has none or says the drug is
      not in DrugBank.
    - ``atc_codes``: the codes of the ATC cell.
    - ``target``: one gene symbol of the targets cell, as written.
    - ``targets``: all gene symbols of that cell.

    A header without the expected columns and a table without records raise
    ``ValueError``.
    """
    found = False
    for number, row in enumerate(iter_tsv(opener, quotechar=None), start=2):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(f'CancerDrugsDB table lacks columns {missing}')
        if any(row.get(column) is None for column in COLUMNS):
            raise ValueError(f'CancerDrugsDB line {number} has fewer fields than the header')
        cells = {column: row[column].strip() for column in COLUMNS}
        found = True
        targets = _list(cells['Targets'])
        record = {
            'product': cells['Product'],
            'ema': cells['EMA'],
            'fda': cells['FDA'],
            'en': cells['EN'],
            'other': cells['Other'],
            'who': cells['WHO'],
            'year': cells['Year'],
            'generic': cells['Generic'],
            'indications': html.unescape(cells['Indications']),
            'last_update': cells['Last Update'],
            'drugbank_ids': _identifiers(cells['DrugBank ID'], _DRUGBANK_REGEX, number, 'DrugBank ID'),
            'chembl_ids': _identifiers(cells['ChEMBL'], _CHEMBL_REGEX, number, 'ChEMBL'),
            'atc_codes': _list(cells['ATC']),
            'targets': targets,
        }
        for target in targets:
            yield {**record, 'target': target}
    if not found:
        raise ValueError('CancerDrugsDB table holds no drugs')
