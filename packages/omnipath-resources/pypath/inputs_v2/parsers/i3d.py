"""
Reader for the Interactome3D interaction table (``interactions.dat``).

Interactome3D annotates protein-protein interactions with structures: a row
describes one experimental structure or one model of the complex of two
proteins. The human representative set has one row, the best ranked structure
or model, for each interaction. The file is a tab separated table with a header
row and 22 columns (the meaning of each is in the help of the web site):

    PROT1 PROT2 RANK_MAJOR RANK_MINOR TYPE PDB_ID BIO_UNIT
    CHAIN1 MODEL1 SEQ_IDENT1 COVERAGE1 SEQ_BEGIN1 SEQ_END1 DOMAIN1
    CHAIN2 MODEL2 SEQ_IDENT2 COVERAGE2 SEQ_BEGIN2 SEQ_END2 DOMAIN2
    FILENAME

``TYPE`` is ``Structure`` (an experimental structure), ``Model`` (a homology
model on a template of the PDB) or ``Dom_dom_model`` (a model on a
domain-domain template of 3did). For a model, ``PDB_ID`` is the template.
``SEQ_IDENT`` and ``COVERAGE`` are percentages and ``SEQ_BEGIN`` and
``SEQ_END`` are residues of the UniProt sequence. The columns ``BIO_UNIT``,
``MODEL`` and ``DOMAIN`` hold ``-`` where they do not apply.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Generator
import logging
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['COLUMNS', 'iter_interactions']

_log = logging.getLogger(__name__)

COLUMNS = (
    'PROT1',
    'PROT2',
    'RANK_MAJOR',
    'RANK_MINOR',
    'TYPE',
    'PDB_ID',
    'BIO_UNIT',
    'CHAIN1',
    'MODEL1',
    'SEQ_IDENT1',
    'COVERAGE1',
    'SEQ_BEGIN1',
    'SEQ_END1',
    'DOMAIN1',
    'CHAIN2',
    'MODEL2',
    'SEQ_IDENT2',
    'COVERAGE2',
    'SEQ_BEGIN2',
    'SEQ_END2',
    'DOMAIN2',
    'FILENAME',
)
_ACCESSION_RE = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
)
_PERCENT_RE = re.compile(r'\d{1,3}(?:\.\d+)?')


def _none(value: str) -> str | None:
    """``-`` stands for no value."""
    return None if value == '-' else value


def _percent(value: str, column: str) -> str:
    if not _PERCENT_RE.fullmatch(value) or float(value) > 100:
        raise ValueError(f'Interactome3D: {column} is not a percentage: {value!r}')
    return value


def _position(value: str, column: str) -> int:
    if not value.isdigit() or int(value) < 1:
        raise ValueError(f'Interactome3D: {column} is not a residue number: {value!r}')
    return int(value)


def iter_interactions(
    opener,
    stats: Counter | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per structure or model of the table.

    The fields are those of :data:`COLUMNS`, in lower case and with the suffix
    ``_a`` and ``_b`` for the two proteins (``uniprot_a`` for ``PROT1``,
    ``chain_a`` for ``CHAIN1``): ``uniprot``, ``rank_major``, ``rank_minor``,
    ``type`` (as written), ``pdb_id``, ``bio_unit``, ``chain``, ``model``,
    ``seq_identity`` and ``coverage`` (percentages, as text),
    ``seq_begin`` and ``seq_end`` (numbers), ``domain`` and ``filename``.
    ``bio_unit``, ``model`` and ``domain`` are None where the file has ``-``.

    A file without the expected columns, a UniProt accession of another form,
    a percentage or residue number that is not one, a range that ends before it
    begins, or no rows at all is an error. ``stats``, if given, is filled with
    the counts of rows by ``TYPE``, which are also logged.
    """
    stats = Counter() if stats is None else stats
    rows = iter_tsv(opener, quotechar=None, strict=True, required=COLUMNS)
    for row in rows:
        stats[row['TYPE']] += 1
        record: dict[str, Any] = {
            'rank_major': row['RANK_MAJOR'],
            'rank_minor': row['RANK_MINOR'],
            'type': row['TYPE'],
            'pdb_id': row['PDB_ID'],
            'bio_unit': _none(row['BIO_UNIT']),
            'filename': row['FILENAME'],
        }
        for number, side in (('1', 'a'), ('2', 'b')):
            accession = row[f'PROT{number}']
            if not _ACCESSION_RE.fullmatch(accession):
                raise ValueError(f'Interactome3D: not a UniProt accession: {accession!r}')
            begin = _position(row[f'SEQ_BEGIN{number}'], f'SEQ_BEGIN{number}')
            end = _position(row[f'SEQ_END{number}'], f'SEQ_END{number}')
            if end < begin:
                raise ValueError(
                    f'Interactome3D: {accession} range {begin}-{end} ends before it begins'
                )
            record.update(
                {
                    f'uniprot_{side}': accession,
                    f'chain_{side}': row[f'CHAIN{number}'],
                    f'model_{side}': _none(row[f'MODEL{number}']),
                    f'seq_identity_{side}': _percent(
                        row[f'SEQ_IDENT{number}'], f'SEQ_IDENT{number}'
                    ),
                    f'coverage_{side}': _percent(
                        row[f'COVERAGE{number}'], f'COVERAGE{number}'
                    ),
                    f'seq_begin_{side}': begin,
                    f'seq_end_{side}': end,
                    f'domain_{side}': _none(row[f'DOMAIN{number}']),
                }
            )
        yield record
    _log.info('Interactome3D table, rows by type: %s', dict(sorted(stats.items())))
