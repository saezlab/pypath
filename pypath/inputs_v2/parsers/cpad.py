"""
Reader for the CPAD download (``CPAD-data.txt``).

The file is tab separated, in the GBK encoding (Greek letters of the names
are two-byte characters), with a header row that ends in three empty column
names. Each line is one entry: a publication that shows how a regulator
(a gene, a microRNA, a lncRNA, a drug...) acts on a cancer pathway and what
follows for a cancer. The shared TSV reader reads the file; :func:`iter_entries`
is the CPAD interpretation: it checks the header and the values of the two
state columns, decodes the lists and gives one record per pair of things the
entry relates.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['iter_entries', 'iter_terms']

COLUMNS = (
    'Regulator',
    'Regulator_Type',
    'Pathway',
    'Pathway_Category',
    'KEGG_ID',
    'Regulation_Type',
    'Cancer',
    'ID',
    'Outcome_Description',
    'Description',
    'Method',
    'Species',
    'Tissue',
    'Cell line',
    'PMID',
)
# The pairs an entry relates. A pair needs both of its parts.
PAIRS = (
    ('regulator_pathway', 'regulator', 'pathway'),
    ('regulator_cancer', 'regulator', 'cancer'),
    ('pathway_cancer', 'pathway', 'cancer'),
)

# The file gives the expression pattern of the regulator (``Regulator_Type``,
# "Expression pattern" on the website) and the state of the pathway
# (``Regulation_Type``, "Pathway state"). An entry with several regulators may
# give the value of each, named: ``menin(Upregulation);Dnmt1(Downregulation)``.
_EXPRESSION = r'Upregulation|Downregulation'
_STATE = r'Activating|Inhibiting'


def _value_pattern(values: str) -> re.Pattern[str]:
    named = rf'[^;/()]+\((?:{values})\)'
    return re.compile(rf'(?:{values})(?:/(?:{values}))?|{named}(?:[;/]{named})*')


EXPRESSION = _value_pattern(_EXPRESSION)
STATE = _value_pattern(_STATE)
# The one cell of the snapshot that is a defect: the expression pattern ``0``
# (entry of the regulator 8-Bromo-7-methoxychrysin). It is read as no value.
DEFECTIVE_EXPRESSION = '0'


def _cell(value: str | None) -> str | None:
    value = None if value is None else value.strip()
    return None if value in (None, '', 'NULL') else value


def _list(value: str | None) -> list[str]:
    cell = _cell(value)
    return [] if cell is None else [item.strip() for item in cell.split(';') if item.strip()]


def _checked(value: str | None, pattern: re.Pattern[str], column: str, line: int) -> str | None:
    cell = _cell(value)
    if cell is not None and not pattern.fullmatch(cell):
        raise ValueError(f'CPAD line {line}: unknown {column}: {cell!r}')
    return cell


def _rows(opener):
    # The header repeats its trailing empty names. The descriptions that hold
    # commas or quotes are quoted.
    return iter_tsv(opener, allow_repeated=True, required=COLUMNS)


def iter_entries(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one record per pair of the entry's regulator, pathway and cancer.

    An entry without a regulator (``NULL`` or empty) relates the pathway and
    the cancer only. Output fields, the columns of the file with ``NULL`` as
    None:

    - ``pair``: ``regulator_pathway``, ``regulator_cancer`` or ``pathway_cancer``.
    - ``regulator``, ``pathway``, ``cancer``: the names as written. A regulator
      of several molecules (``AXIN and APC``) is one name, as on the website.
    - ``regulator_expression``: ``Regulator_Type``, the expression pattern of
      the regulator (``Upregulation``, ``Downregulation``).
    - ``pathway_state``: ``Regulation_Type``, the state of the pathway
      (``Activating``, ``Inhibiting``).
    - ``pathway_categories``, ``kegg_pathways``: the lists of the
      ``Pathway_Category`` column (names of KEGG pathways, ``Others``) and of
      the ``KEGG_ID`` column (``path:04151``).
    - ``disease_id``: the ``ID`` column (``H00031``, which the website labels
      ICD-10; they are accessions of KEGG DISEASE).
    - ``outcome``, ``description``: ``Outcome_Description`` and ``Description``.
    - ``methods``, ``cell_lines``: the lists of the ``Method`` and ``Cell line``
      columns.
    - ``species``, ``tissue``: as written. ``pmid``: the PubMed identifier.

    The ``0`` of the expression pattern column (one entry) is read as no
    value. A value of the expression pattern or the state that is not one of
    the source's, a PubMed identifier of another form and a table without
    entries raise ``ValueError``.
    """
    found = False
    for line, row in enumerate(_rows(opener), start=2):
        found = True
        expression = _cell(row['Regulator_Type'])
        if expression == DEFECTIVE_EXPRESSION:
            expression = None
        expression = _checked(expression, EXPRESSION, 'expression pattern', line)
        state = _checked(row['Regulation_Type'], STATE, 'pathway state', line)
        pmid = _cell(row['PMID'])
        if pmid is not None and not is_pmid(pmid):
            raise ValueError(f'CPAD line {line}: unexpected PubMed identifier {pmid!r}')
        record = {
            'regulator': _cell(row['Regulator']),
            'pathway': _cell(row['Pathway']),
            'cancer': _cell(row['Cancer']),
            'regulator_expression': expression,
            'pathway_state': state,
            'pathway_categories': _list(row['Pathway_Category']),
            'kegg_pathways': _list(row['KEGG_ID']),
            'disease_id': _cell(row['ID']),
            'outcome': _cell(row['Outcome_Description']),
            'description': _cell(row['Description']),
            'methods': _list(row['Method']),
            'cell_lines': _list(row['Cell line']),
            'species': _cell(row['Species']),
            'tissue': _cell(row['Tissue']),
            'pmid': pmid,
        }
        for pair, first, second in PAIRS:
            if record[first] is not None and record[second] is not None:
                yield {'pair': pair, **record}
    if not found:
        raise ValueError('CPAD table has no entries')


def iter_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield each pathway, then each cancer of the table once, in order of first appearance.

    Output fields: ``kind`` (``pathway`` or ``cancer``) and ``name``, as written.
    """
    names: dict[str, dict[str, None]] = {'pathway': {}, 'cancer': {}}
    for row in _rows(opener):
        for kind, column in (('pathway', 'Pathway'), ('cancer', 'Cancer')):
            name = _cell(row[column])
            if name is not None:
                names[kind].setdefault(name)
    for kind, found in names.items():
        for name in found:
            yield {'kind': kind, 'name': name}
