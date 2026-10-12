"""
Parser of Supplementary Data 2 of Ramilowski et al. 2015 (Nature Communications).

The workbook has the sheet `All.Pairs`, one row for each candidate ligand-receptor
pair, and the sheets `HPMR`, `IUPHAR` and `HPRD`, which carry the PubMed IDs those
databases give for a pair (one row for each `Pair.Name`). The parser yields the
rows of `All.Pairs`, keyed by its header (`Pair.Name`, `Ligand.ApprovedSymbol`,
`Ligand.Name`, `Receptor.ApprovedSymbol`, `Receptor.Name`, `DLRP`, `HPMR`, `IUPHAR`,
`HPRD`, `STRING.binding`, `STRING.experiment`, `HPMR.Ligand`, `HPMR.Receptor`,
`PMID.Manual`, `Pair.Source`, `Pair.Evidence`; empty cells are None), and adds the
cell of each database sheet as `HPMR.PMID`, `IUPHAR.PMID` and `HPRD.PMID`
(None for a pair the sheet does not list).

The six database columns (`DLRP`, `HPMR`, `IUPHAR`, `HPRD`, `STRING.binding`,
`STRING.experiment`) hold the column's own name when the database lists the pair
and are empty otherwise; any other value raises.

A PubMed cell is a comma separated list of numbers; a few carry a stray closing
bracket, which `split_pubmed_cell` removes.
"""

from __future__ import annotations

from collections.abc import Generator
import io
import re
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_excel, read_bytes

# Sheets holding the PubMed IDs a source database gives for a pair.
PMID_SHEETS = ('HPMR', 'IUPHAR', 'HPRD')
# Columns of `All.Pairs` flagging the databases a pair was found in: the cell
# repeats the column's name or is empty.
DATABASE_COLUMNS = (
    'DLRP',
    'HPMR',
    'IUPHAR',
    'HPRD',
    'STRING.binding',
    'STRING.experiment',
)
# A PubMed ID as the sheets write it; a few carry a stray closing bracket.
_pmid_re = re.compile('(\\d+)\\)?')


def iter_pairs(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Rows of `All.Pairs`, each with its pair's PubMed IDs from the database sheets.

    The joined cells are added unchanged as `HPMR.PMID`, `IUPHAR.PMID` and
    `HPRD.PMID`; every sheet has one row per `Pair.Name`.
    """
    content = read_bytes(opener)
    if not content:
        return

    def sheet(
        name: str, required: tuple[str, ...] = ()
    ) -> Generator[dict[str, Any], None, None]:
        return iter_excel(
            SimpleNamespace(result=io.BytesIO(content)),
            sheet=name,
            required=required or None,
        )

    pmids = {
        name: {row['Pair.Name']: row.get('PMID') for row in sheet(name)}
        for name in PMID_SHEETS
    }
    for row in sheet('All.Pairs', required=DATABASE_COLUMNS):
        for column in DATABASE_COLUMNS:
            if row.get(column) not in (None, column):
                raise ValueError(
                    f'Unknown Ramilowski 2015 database flag {row[column]!r} '
                    f'in {column} (pair {row.get("Pair.Name")!r})'
                )
        for name in PMID_SHEETS:
            row[f'{name}.PMID'] = pmids[name].get(row.get('Pair.Name'))
        yield row


def split_pubmed_cell(value: Any, column: str) -> list[str]:
    """The PubMed IDs of one cell in source order; a token that is not a number
    (with at most a stray closing bracket) raises."""
    out: list[str] = []
    for token in str(value or '').split(','):
        if not (token := token.strip()):
            continue
        match = _pmid_re.fullmatch(token)
        if not match:
            raise ValueError(
                f'Unknown Ramilowski 2015 PubMed ID {token!r} in {column}'
            )
        # A number that is not the form of a PubMed ID is left out.
        if is_pmid(match.group(1)):
            out.append(match.group(1))
    return out
