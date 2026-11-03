"""
Read the ORegAnno combined export and select the transcription factor - gene records.

Source file (tab separated, with a header; 18 columns):

    ORegAnno_ID  Species  Outcome  Type  Gene_Symbol  Gene_ID  Gene_Source
    Regulatory_Element_Symbol  Regulatory_Element_ID  Regulatory_Element_Source
    dbSNP_ID  PMID  Dataset  Build  Strand  Chr  Start  End

One ORegAnno record (``ORegAnno_ID``) takes one line for each genome build it
is mapped to (hg18, hg19, hg38, ...), so a record is on three or more lines
that differ only in ``Build``, ``Strand``, ``Chr``, ``Start`` and ``End``. A
record can also have lines for several genes or several factors.

A line is a transcription factor - target gene record when its ``Type`` is
``TRANSCRIPTION FACTOR BINDING SITE``, its ``Outcome`` is ``POSITIVE OUTCOME``,
and the line names both the factor (``Regulatory_Element_Symbol``) and the gene
(``Gene_Symbol``). ``exclusion_reason`` says why any other line is not one.
The rest of the file holds binding sites without a gene (mostly ChIP-seq
peaks), regulatory regions, polymorphisms and haplotypes, miRNA binding sites
and sites whose outcome was negative, neutral, dual or not given.
"""

from __future__ import annotations

from collections.abc import Generator
import logging
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

_log = logging.getLogger(__name__)

COLUMNS = (
    'ORegAnno_ID',
    'Species',
    'Outcome',
    'Type',
    'Gene_Symbol',
    'Gene_ID',
    'Gene_Source',
    'Regulatory_Element_Symbol',
    'Regulatory_Element_ID',
    'Regulatory_Element_Source',
    'dbSNP_ID',
    'PMID',
    'Dataset',
    'Build',
    'Strand',
    'Chr',
    'Start',
    'End',
)
BINDING_SITE = 'TRANSCRIPTION FACTOR BINDING SITE'
POSITIVE = 'POSITIVE OUTCOME'
# The value of the regulatory element symbol of sites whose factor is not known.
UNSPECIFIED = 'Unspecified'
# Regulatory_Element_Source of a binding site of a microRNA, not a factor.
MIRNA_SOURCE = 'miRBase'
_NULLS = frozenset({'', 'N/A', 'NA'})


def _clean(value: Any) -> str | None:
    """A cell without surrounding blanks; None for the file's missing values."""
    value = (value or '').strip()
    return None if value in _NULLS else value


def exclusion_reason(row: dict[str, Any]) -> str | None:
    """
    Why a line of the file is not a transcription factor - gene record.

    None for a record line. The reasons are checked in this order, so each
    other line has exactly one.
    """
    if row['Type'].strip() != BINDING_SITE:
        return 'not a transcription factor binding site'
    if row['Outcome'].strip() != POSITIVE:
        return 'outcome is not positive'
    tf = _clean(row['Regulatory_Element_Symbol'])
    if tf is None:
        return 'no transcription factor'
    if tf == UNSPECIFIED:
        return 'transcription factor unspecified'
    if _clean(row['Regulatory_Element_Source']) == MIRNA_SOURCE:
        return 'regulatory element is a microRNA'
    if _clean(row['Gene_Symbol']) is None:
        return 'no target gene'
    return None


def _malformed(row: dict[str, Any]) -> bool:
    """A line with fewer or more fields than the header has columns."""
    return None in row or any(value is None for value in row.values())


def iter_lines(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    All lines of the file keyed by the header.

    The rescued copy of the file ends in the middle of a line. A line with a
    wrong number of fields is an error unless it is the last one, which is
    left out with a warning. A header without one of the columns and a file
    without lines are errors as well.
    """
    previous = None
    count = 0
    for row in iter_tsv(opener, quotechar=None, **kwargs):
        if not count:
            missing = [c for c in COLUMNS if c not in row]
            if missing:
                raise ValueError(f'ORegAnno file lacks the columns {missing}')
        count += 1
        if previous is not None:
            if _malformed(previous):
                raise ValueError(
                    f'ORegAnno line with another number of fields than the '
                    f'header: {previous!r}'
                )
            yield previous
        previous = row
    if not count:
        raise ValueError('ORegAnno file has no rows')
    if _malformed(previous):
        _log.warning(
            'ORegAnno file ends in an incomplete line, left out: %r', previous
        )
    else:
        yield previous


def iter_binding_sites(
    opener,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the transcription factor - target gene records of the file.

    The lines of one record that differ only in the genome build are yielded
    once. Keys: ``oreganno_id``, ``species`` (as the file names it), ``tf_symbol``,
    ``tf_id``, ``tf_source``, ``gene_symbol``, ``gene_id``, ``gene_source``,
    ``pmid`` and ``dataset``. A missing value is None. An identifier whose
    source column is missing is None as well: the file then repeats the
    symbol, or says ``Unspecified``, and it is not an accession. A header
    without one of the columns, a line with another number of fields than the
    header (except an incomplete last line) and a file without lines are
    errors.
    """
    seen = set()
    for row in iter_lines(opener, **kwargs):
        if exclusion_reason(row) is not None:
            continue
        record = {
            'oreganno_id': row['ORegAnno_ID'].strip(),
            'species': row['Species'].strip(),
            'tf_symbol': _clean(row['Regulatory_Element_Symbol']),
            'tf_source': _clean(row['Regulatory_Element_Source']),
            'gene_symbol': _clean(row['Gene_Symbol']),
            'gene_source': _clean(row['Gene_Source']),
            'pmid': _clean(row['PMID']),
            'dataset': _clean(row['Dataset']),
        }
        record['tf_id'] = (
            _clean(row['Regulatory_Element_ID']) if record['tf_source'] else None
        )
        record['gene_id'] = (
            _clean(row['Gene_ID']) if record['gene_source'] else None
        )
        key = tuple(record.values())
        if key in seen:
            continue
        seen.add(key)
        yield record
