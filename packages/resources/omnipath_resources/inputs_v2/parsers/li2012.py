"""Read the Li et al. 2012 supplement: tyrosine kinase circuits and families.

The supplement is one ``.xls`` workbook. File S1 lists the 583 curated
"tyrosine kinase circuits", each a tyrosine kinase family (the writer), a
substrate with a phosphotyrosine position and a family of SH2 or PTB domain
proteins that binds the phosphorylated site (the reader):

    substrate family  protein/pY position  TK writers  SH2/PTB Readers
    BCR               BCR/177              ABL         GRB2 SH2

The paper assembled the circuits from curated kinase-substrate and
substrate-reader interactions, so a site with two writers and two readers
gives four rows. The readers below read the two kinds of interaction back
out: one row per distinct writer and site, and one per distinct site and
reader. File S3 lists the members of the writer and reader families. The
predicted circuits (Files S2 and S4) are not read.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

# Column names of File S1 (circuits) and File S3 (families).
SITE = 'protein/pY position'
WRITER = 'TK writers'
READER = 'SH2/PTB Readers'
CIRCUIT_COLUMNS = (
    'substrate family',
    SITE,
    WRITER,
    READER,
    'signaling route',
    'tissue specificity',
    'evo-patterns',
)
FAMILY = 'TK/SH2/PTB family'
MEMBERS = 'gene members'
MEMBER_IDS = 'member ID (Gene ID)'
FAMILY_COLUMNS = (
    FAMILY,
    'family type',
    MEMBERS,
    MEMBER_IDS,
    'tissue specificity',
    'cellular localization',
)

# "BCR/177": the substrate's gene symbol and the phosphotyrosine position.
_SITE_REGEX = re.compile(r'^([A-Za-z0-9-]+)/([1-9]\d*)$')
# "GRB2 SH2": the reader family and its phosphotyrosine binding domain.
_READER_REGEX = re.compile(r'^(\S+) (SH2|PTB)$')


def _sheet(
    opener, sheet: str, columns: tuple[str, ...]
) -> Generator[dict[str, Any], None, None]:
    """Rows of a supplemental sheet, whose first row is its title.

    A changed header is an error, not an empty dataset.
    """
    for row in iter_excel(opener, sheet=sheet, header_row=2):
        missing = [column for column in columns if column not in row]
        if missing:
            raise ValueError(f'Li 2012 sheet {sheet!r} lacks columns {missing}')
        yield {column: row[column] for column in columns}


def iter_circuits(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the circuit rows of File S1.

    Below the circuits the sheet lists the references of the data set; those
    rows have neither a writer nor a reader and are skipped. A circuit with
    a malformed site or reader is an error.
    """
    for row in _sheet(opener, 'File S1', CIRCUIT_COLUMNS):
        if row[WRITER] is None and row[READER] is None:
            continue
        if (
            row[WRITER] is None
            or not _SITE_REGEX.match(str(row[SITE] or ''))
            or not _READER_REGEX.match(str(row[READER] or ''))
        ):
            raise ValueError(f'Malformed Li 2012 circuit: {row!r}')
        yield row


def _distinct(
    rows, key: tuple[str, str]
) -> Generator[dict[str, Any], None, None]:
    """One row per distinct value of the two key columns, in file order.

    Every other column holds the list of its distinct values over the
    circuits that share the key.
    """
    groups: dict[tuple[Any, Any], dict[str, list[Any]]] = {}
    for row in rows:
        group = groups.setdefault(
            (row[key[0]], row[key[1]]),
            {column: [] for column in row},
        )
        for column, value in row.items():
            if value not in group[column]:
                group[column].append(value)
    for group in groups.values():
        yield {
            column: values[0] if column in key else values
            for column, values in group.items()
        }


def _site_fields(row: dict[str, Any]) -> dict[str, Any]:
    substrate, position = _SITE_REGEX.match(row[SITE]).groups()
    return {'substrate': substrate, 'position': int(position)}


def iter_kinase_substrate(
    opener, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per distinct tyrosine kinase family and site.

    ``TK writers`` and ``protein/pY position`` are as written; ``substrate``
    and ``position`` are read from the latter. The other File S1 columns are
    lists of their values over the circuits of this kinase and site.
    """
    for row in _distinct(iter_circuits(opener), (WRITER, SITE)):
        yield {**row, **_site_fields(row)}


def iter_phosphosite_binding(
    opener, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per distinct site and reader.

    ``protein/pY position`` and ``SH2/PTB Readers`` are as written;
    ``substrate`` and ``position`` are read from the former,
    ``reader_family`` and ``reader_domain`` from the latter. The other File
    S1 columns are lists of their values over the circuits of this pair.
    """
    for row in _distinct(iter_circuits(opener), (SITE, READER)):
        family, domain = _READER_REGEX.match(row[READER]).groups()
        yield {
            **row,
            **_site_fields(row),
            'reader_family': family,
            'reader_domain': domain,
        }


def _list(cell: object) -> list[str]:
    """Items of a comma separated cell; Excel stores a single gene ID as a number."""
    if isinstance(cell, float) and cell.is_integer():
        cell = int(cell)
    return [item.strip() for item in str(cell or '').split(',') if item.strip()]


def iter_families(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the kinase and reader families of File S3.

    ``members`` pairs each gene symbol of ``gene members`` with its Entrez
    Gene ID from ``member ID (Gene ID)``; the two lists must align.
    """
    for row in _sheet(opener, 'File S3', FAMILY_COLUMNS):
        symbols, gene_ids = _list(row[MEMBERS]), _list(row[MEMBER_IDS])
        if len(symbols) != len(gene_ids) or not all(i.isdigit() for i in gene_ids):
            raise ValueError(f'Li 2012 family members do not align: {row!r}')
        row[MEMBER_IDS] = ','.join(gene_ids)
        yield {
            **row,
            'members': [
                {'symbol': symbol, 'entrez': gene_id}
                for symbol, gene_id in zip(symbols, gene_ids)
            ],
        }
