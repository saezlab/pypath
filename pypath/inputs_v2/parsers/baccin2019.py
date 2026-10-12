"""
Parser of Supplementary Table 3 of Baccin et al. 2019 (Nature Cell Biology).

The table is the sheet `SuppTable3` of the supplementary workbook: the
ligand-receptor database of the RNA-Magnet method (version 1.0.0 of the
RNAMagnet package). The sheet's third row is the header, the columns are

    Pair.Name, Ligand.Mouse, Receptor.Mouse, Source, ManualAnnotation,
    Ligand.CC, Ligand.GO, Reference (PMID / KEGG)

A partner cell names one gene, several genes joined by `&` that form a complex,
or several genes or complexes joined by `|` of which any one is the partner.
The parser yields one row for each combination of alternatives: the source
columns unchanged, plus `ligand` and `receptor`, the partner of this row (the
whole cell unless it lists alternatives). A table without the expected columns,
or without rows, raises.

The reference cell holds PubMed IDs, which the spreadsheet damaged, and the names
of KEGG pathway maps; `reference_tokens` and `pubmed_ids` decode it.
"""

from __future__ import annotations

from collections.abc import Generator
import itertools
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

SHEET = 'SuppTable3'
# Columns of `SuppTable3`, in source order. The header is the sheet's third row.
COLUMNS = (
    'Pair.Name',
    'Ligand.Mouse',
    'Receptor.Mouse',
    'Source',
    'ManualAnnotation',
    'Ligand.CC',
    'Ligand.GO',
    'Reference (PMID / KEGG)',
)
REFERENCE = 'Reference (PMID / KEGG)'
# Stricter than the shared PubMed identifier rule (1 to 8 digits), because this
# file has a known defect. Every intact PMID in the sheet has seven or eight
# digits. The spreadsheet fused many comma-separated lists into one number
# (`1578147310751170`) and cut a few PMIDs short; the cut ones are left out. The
# five numbers of five or six digits were looked up at NCBI (esummary, October
# 2026) and are 1975-1978 papers on other subjects (`939660` of Hspg2-Sdc1 is
# "The operative report", `62193` of Il1b-Il1r1 is "Assessment of paracetamol
# poisoning"), and three of them are the tail of a PMID that the same cell lists
# in full (`62193` of `9062193`, the structure of the IL-1 receptor complex;
# `414473` of `9414473`; `807467` of `8074647`). Of the fused ones only a token
# of exactly two PMIDs is repaired (see `split_fused_pmids`). Two PMIDs carry a
# stray closing parenthesis.
_pmid_re = re.compile(r'(\d{7,8})\)?')


def _alternatives(value: Any) -> list[str]:
    """The partners a cell names: one, or each gene or complex of a `|` list."""
    return [
        token for part in str(value or '').split('|') if (token := part.strip())
    ]


def iter_pairs(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Rows of `SuppTable3`, one for each pair of alternative partners.

    Every source column is kept under its own name. `ligand` and `receptor`
    are added: the partner of this record, which is the whole cell unless the
    cell lists alternatives separated by `|`. A row whose cells list
    alternatives is yielded once for each combination.
    """
    count = 0
    for row in iter_excel(opener, sheet=SHEET, header_row=3):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f'Baccin 2019 table lacks expected columns: {", ".join(missing)}'
            )
        count += 1
        for ligand, receptor in itertools.product(
            _alternatives(row['Ligand.Mouse']),
            _alternatives(row['Receptor.Mouse']),
        ):
            yield {**row, 'ligand': ligand, 'receptor': receptor}
    if not count:
        raise ValueError('Baccin 2019 table has no rows')


def reference_tokens(value: Any) -> list[str]:
    """The comma or semicolon separated entries of a reference cell; a number
    cell (an integer stored as a float) is read as its integer."""
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return [
        token
        for part in re.split('[,;]', str(value if value is not None else ''))
        if (token := part.strip())
    ]


def split_fused_pmids(token: str) -> list[str] | None:
    """The two PMIDs a token is made of, or None.

    The spreadsheet ran the two PMIDs of a list together into one number. The
    token is split only when it has one reading as two well-formed PMIDs of
    seven or eight digits, none starting with a zero: 14 digits are 7 + 7 and 16
    digits are 8 + 8. A token of 15 digits reads as 7 + 8 and as 8 + 7, so it is
    not split, and neither is one of any other length (13, or 21 and more, where
    the spreadsheet also turned the trailing digits into zeros).
    """
    if not token.isascii() or not token.isdigit() or len(token) not in (14, 16):
        return None
    half = len(token) // 2
    first, second = token[:half], token[half:]
    if first[0] == '0' or second[0] == '0':
        return None
    return [first, second]


def pubmed_ids(value: Any) -> list[str]:
    """PubMed IDs of a reference cell in source order, each once."""
    out: list[str] = []
    for token in reference_tokens(value):
        match = _pmid_re.fullmatch(token)
        pmids = [match.group(1)] if match else split_fused_pmids(token) or []
        for pmid in pmids:
            if pmid not in out:
                out.append(pmid)
    return out
