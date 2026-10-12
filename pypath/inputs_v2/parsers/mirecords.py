"""
Read the miRecords table of experimentally validated miRNA - target interactions.

Source file (``miRecords_v4.xls``, Excel 97-2003, one sheet ``miRecords_20130427``
with 26 columns; the first row is the header, the release is 27 April 2013):

    Pubmed_id | Target gene_species_scientific | Target gene_name | Target gene_Refseq_acc |
    Target site_number | miRNA_species | miRNA_mature_ID | miRNA_regulation | ...

Every row is one literature account of one miRNA - target pair at the level of
the target gene, region or site (the columns of the experiments, which the
database describes as "target gene level", "target region level" and "target
site level" evidence, are filled in accordingly; 1,825 rows have ``NA`` in them).

The header has two columns named ``Reporter link element``, so the columns are
read by position and the header is checked against ``COLUMNS``.

Defects of the source that the parser repairs, each only where the repair is
exact and unambiguous:

* spaces around values, and runs of spaces inside a species name;
* ``Xenopus. tropicalis`` (a stray period) for ``Xenopus tropicalis``;
* text that was UTF-8 but was read as Mac Roman before it was saved (``Œ∫B-Ras2``
  for ``κB-Ras2``); the text is converted back only if it is then valid UTF-8;
* a trailing ``|`` after a RefSeq accession;
* ``unknown`` (and its misspelling ``unkown``) where the source has no gene name or
  accession: it is a missing value, not a name.

What the parser does not repair, but reports: a gene name that Excel turned into
a date (``39692``, ``39693``: ``SEPT1``/``SEPT2`` typed as ``1-Sep``/``2-Sep``) is
given as ``target_gene_excel_date`` and not as a gene name.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

# Header of the sheet, in source order (the two ``Reporter link element``
# columns are the 10th and the 20th).
COLUMNS = (
    'Pubmed_id',
    'Target gene_species_scientific',
    'Target gene_name',
    'Target gene_Refseq_acc',
    'Target site_number',
    'miRNA_species',
    'miRNA_mature_ID',
    'miRNA_regulation',
    'Reporter_target gene/region',
    'Reporter link element',
    'Test_method_inter',
    'Target gene mRNA_level',
    'Original description',
    'Mutation_target region',
    'Post mutation_method',
    'Original description_mutation_region',
    'Target site_position',
    'miRNA_regulation_site',
    'Reporter_target site',
    'Reporter link element',
    'Test_method_inter_site',
    'Original description_inter_site',
    'Mutation_target site',
    'Post mutation_method_site',
    'Original description_mutation_site',
    'Additional note',
)

# Names for the columns, in the order of ``COLUMNS``; those starting with an
# underscore are not part of the parsed rows.
_NAMES = (
    'pmid',
    'target_species',
    'target_gene',
    'target_refseq',
    'site_number',
    'mirna_species',
    'mirna',
    'mirna_regulation',
    'reporter',
    '_reporter_link_element',
    'test_method',
    'mrna_level',
    '_original_description',
    '_mutation_target_region',
    '_post_mutation_method',
    '_original_description_mutation_region',
    '_target_site_position',
    '_mirna_regulation_site',
    '_reporter_target_site',
    '_reporter_link_element_site',
    '_test_method_site',
    '_original_description_site',
    '_mutation_target_site',
    '_post_mutation_method_site',
    '_original_description_mutation_site',
    '_additional_note',
)

# The cell values with which the source marks a missing RefSeq accession.
_MISSING_ACCESSION = frozenset({'unknown', 'unkown'})

# The cell value with which the source marks a missing gene name.
_MISSING_NAME = frozenset({'unknown'})

_WHITESPACE = re.compile(r'\s+')
_PMID = re.compile(r'^(\d+)(?:\.0+)?$')


def _text(value: Any) -> str | None:
    """A text cell with the spaces around it removed; None if it is empty."""
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def repair_encoding(text: str) -> str:
    """
    Undo text that was UTF-8 but was read as Mac Roman (``Œ∫`` for ``κ``).

    Text that is not made of Mac Roman characters or that does not then form
    valid UTF-8 is returned unchanged.
    """
    if text.isascii():
        return text
    try:
        return text.encode('mac_roman').decode('utf-8')
    except (UnicodeEncodeError, UnicodeDecodeError):
        return text


def _species(value: Any) -> str | None:
    text = _text(value)
    if text is None:
        return None
    text = _WHITESPACE.sub(' ', text)
    # A stray period after the genus.
    return text.replace('Xenopus. tropicalis', 'Xenopus tropicalis')


def _brace_list(value: Any) -> list[str]:
    """
    The labels of a cell written as ``{a}{b}`` (or ``a}{b``) as a list.

    ``NA`` (not applicable) and an empty cell give an empty list.
    """
    text = _text(value)
    if text is None or text == 'NA':
        return []
    return [
        label.strip()
        for label in text.strip('{}').split('}{')
        if label.strip()
    ]


def _record(values: dict[str, Any], number: int) -> dict[str, Any]:
    pmid = values['pmid']
    pmid = _PMID.match(str(pmid).strip()) if pmid is not None else None
    if not pmid:
        raise ValueError(f'Row {number}: Pubmed_id is not a number: {values["pmid"]!r}')

    mirna = _text(values['mirna'])
    if mirna is None:
        raise ValueError(f'Row {number}: no miRNA name')
    mirna = repair_encoding(mirna)
    # A name between square brackets is the name as the paper gives it; see the
    # module that maps the rows.
    bracketed = mirna.startswith('[') and mirna.endswith(']')
    if bracketed:
        mirna = mirna[1:-1].strip()

    gene = values['target_gene']
    excel_date = None
    if isinstance(gene, (int, float)):
        # A cell that Excel has turned into a date: its serial number.
        if gene != int(gene):
            raise ValueError(f'Row {number}: gene name is a fraction: {gene!r}')
        excel_date = int(gene)
        gene = None
    else:
        gene = _text(gene)
        if gene is not None:
            gene = repair_encoding(gene)
        if gene in _MISSING_NAME:
            gene = None

    refseq = _text(values['target_refseq'])
    if refseq is not None:
        refseq = refseq.rstrip('|').strip()
        if refseq in _MISSING_ACCESSION:
            refseq = None

    site_number = values['site_number']
    if isinstance(site_number, (int, float)):
        site_number = int(site_number)
    elif _text(site_number) in (None, 'NA'):
        site_number = None

    return {
        'pmid': pmid.group(1),
        'target_species': _species(values['target_species']),
        'target_gene': gene,
        'target_gene_excel_date': excel_date,
        'target_refseq': refseq,
        'site_number': site_number,
        'mirna_species': _species(values['mirna_species']),
        'mirna': mirna,
        'mirna_bracketed': bracketed,
        'mirna_regulation': _brace_list(values['mirna_regulation']),
        'reporter': None if _text(values['reporter']) == 'NA' else _text(values['reporter']),
        'test_method': _brace_list(values['test_method']),
        'mrna_level': _brace_list(values['mrna_level']),
    }


def iter_interactions(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the literature accounts of miRNA - target pairs of the sheet.

    Keys:

    ``pmid``
        PubMed ID as digits (the sheet has it as a float); ``0`` for two rows
        without a publication.
    ``target_species``, ``mirna_species``
        The scientific names of the two species, as written apart from the
        repairs of the module docstring. A few cells do not name a species
        (``[Virus]``, a RefSeq accession, ``Xenopus``, ``us musculus``) and are
        kept as they are.
    ``target_gene``
        The gene name as written (spaces removed, encoding repaired); None where
        the cell is empty or was turned into a date by Excel, in which case
        ``target_gene_excel_date`` holds the date's serial number.
    ``target_refseq``
        The accession of the ``Target gene_Refseq_acc`` column; None where the
        source writes ``unknown``. It is not always a RefSeq accession (GenBank
        accessions and ``NM213137`` occur).
    ``site_number``
        The number of target sites of the account (None for ``NA``).
    ``mirna``
        The miRNA name as written (spaces removed, encoding repaired), without
        the square brackets ``mirna_bracketed`` tells about.
    ``mirna_regulation``, ``test_method``, ``mrna_level``, ``reporter``
        The experimental columns of the target gene or region level: how the
        miRNA level was changed, the methods that measured the target's
        protein level, whether the target's mRNA level changed, and the
        reporter. None of them gives the direction of the effect on the target.

    The columns of the target site level evidence are not returned. The header
    must be that of the source; an empty table is an error.
    """
    count = 0
    # header_row=0: the header is the first row of the data, which is checked here.
    rows = iter_excel(
        opener,
        sheet=0,
        header_row=0,
        fieldnames=[f'c{index}' for index in range(len(COLUMNS))],
        **kwargs,
    )
    for number, row in enumerate(rows, start=1):
        values = [row.get(f'c{index}') for index in range(len(COLUMNS))]
        if number == 1:
            if tuple(values) != COLUMNS:
                raise ValueError(f'Unexpected miRecords columns: {values}')
            continue
        count += 1
        yield _record(
            {
                name: value
                for name, value in zip(_NAMES, values)
                if not name.startswith('_')
            },
            number,
        )
    if not count:
        raise ValueError('miRecords table has no rows')

