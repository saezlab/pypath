"""Read the DEPOD phosphatase-substrate worksheets.

DEPOD-2019 publishes one worksheet of human phosphatases with their protein
substrates (``PSprots``, in ``PPase_protSubtrates_201903.xls``) and one with
their non-protein substrates (``PSnonprots``, in
``PPase_nonprotSubtrates_201903.xls``). A protein substrate row lists the
dephosphorylated residues of one phosphatase-substrate pair
("Ser-339, Ser-340"), so a row becomes one parsed row per residue. Rows
without a site ("N/A") and rows whose site cell holds free text
("Ser-660 in Isoform Beta-II") stay single rows: the first has no site, the
second keeps the source text without a position.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_excel

# Column names of the two worksheets.
ENZYME = 'Phosphatase accession numbers'
SUBSTRATE = 'Substrate accession numbers'
SITES = 'Dephosphosites'
WINDOWS = '5 amino acid window around the dephosphosite (small letters)'
ASSAY = 'Assay/s to infer dephosphorylation'
REFERENCES = 'Literature source/s'
SCORE = 'Reliability score'
PROTEIN_COLUMNS = (
    ENZYME,
    'Phosphatase entry names',
    SUBSTRATE,
    'Substrate entry names',
    SITES,
    WINDOWS,
    ASSAY,
    REFERENCES,
    SCORE,
)
NON_PROTEIN_COLUMNS = (
    ENZYME,
    'Phosphatase entry name',
    'Substrate ChEBI ID',
    'Substrate KEGG ID',
    'Substrate ChEBI name',
    'Dephosphosite',
    ASSAY,
    REFERENCES,
    SCORE,
)

# The source writes a missing value both ways.
_MISSING = {'N/A', 'NA'}
# "Ser-88"; a few cells lack the hyphen ("Ser211").
_SITE_REGEX = re.compile(r'^(Ser|Thr|Tyr|His)-?([1-9]\d*)$')
_RESIDUES = {'Ser': 'S', 'Thr': 'T', 'Tyr': 'Y', 'His': 'H'}
# A sequence window: residues in capitals, the dephosphorylated one in lower case.
_WINDOW_REGEX = re.compile(r'^(?=[A-Za-z-]*[A-Z])(?=[A-Za-z-]*[a-z])[A-Za-z-]+$')


def _text(value: object) -> str | None:
    """A cell as text; Excel stores a single PubMed ID as a number."""
    if value is None:
        return None
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return str(value).strip() or None


def _checked_rows(
    opener, sheet: str, columns: tuple[str, ...]
) -> Generator[dict[str, Any], None, None]:
    """Rows of a DEPOD worksheet; a changed header is an error, not an empty dataset."""
    for row in iter_excel(opener, sheet=sheet):
        missing = [column for column in columns if column not in row]
        if missing:
            raise ValueError(f'DEPOD worksheet {sheet!r} lacks columns {missing}')
        yield row


def _windows(cell: str | None, sites: list[tuple[str, str, int]]) -> list[str | None]:
    """Sequence windows of the sites, where the cell lists one per site.

    A window is assigned only if the two lists have the same length and each
    window marks its site's residue in lower case; repeat motifs given once
    for many sites, and comments in this column, stay in the source cell.
    """
    windows = [w for w in re.split(r'[,;\s]+', cell or '') if w]
    aligned = len(windows) == len(sites) and all(
        [c for c in window if c.islower()] == [residue.lower()]
        for window, (_, residue, _) in zip(windows, sites)
    )
    return windows if aligned else [None] * len(sites)


def _window_comment(cell: str | None) -> str | None:
    """Text of the windows column that is not a list of windows.

    One row comments on its site there ("1148 is not a Tyr in HsEGFR"). A
    bare number, which a few cells hold, says nothing about the site.
    """
    if cell is None or cell in _MISSING or cell.isdigit():
        return None
    tokens = [token for token in re.split(r'[,;\s]+', cell) if token]
    return None if all(_WINDOW_REGEX.match(token) for token in tokens) else cell


def iter_protein_substrates(
    opener, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per dephosphorylation site of the protein substrate sheet.

    Each dictionary has the source row's columns as written, and:

    - ``site``: the site as written (``Ser-88``); for a cell with free text
      the whole cell, followed in brackets by a comment the source left in
      the windows column; None if the row names no site.
    - ``residue``, ``position``: one-letter residue and integer position of a
      site written as ``Ser-88``; None for free text and commented sites.
    - ``window``: the site's sequence window from the windows column, if that
      column lists one window per site.
    - ``pmids``: the PubMed IDs of the literature column.
    """
    for row in _checked_rows(opener, 'PSprots', PROTEIN_COLUMNS):
        row[REFERENCES] = _text(row[REFERENCES])
        row[WINDOWS] = _text(row[WINDOWS])
        row['pmids'] = pmids(row[REFERENCES])
        cell = _text(row[SITES])
        tokens = [token.strip() for token in re.split(r'[,;]', cell or '')]
        matches = [_SITE_REGEX.match(token) for token in tokens]
        comment = _window_comment(row[WINDOWS])
        if cell is None or cell in _MISSING:
            yield {**row, 'site': None, 'residue': None, 'position': None, 'window': None}
        elif all(matches) and comment is None:
            sites = [
                (token, _RESIDUES[match[1]], int(match[2]))
                for token, match in zip(tokens, matches)
            ]
            windows = _windows(row[WINDOWS], sites)
            for (token, residue, position), window in zip(sites, windows):
                yield {
                    **row,
                    'site': token,
                    'residue': residue,
                    'position': position,
                    'window': window,
                }
        else:
            # Isoform numbering, alternatives and comments: no position is read.
            text = cell if comment in (None, cell) else f'{cell} ({comment})'
            yield {**row, 'site': text, 'residue': None, 'position': None, 'window': None}


def iter_non_protein_substrates(
    opener, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield the rows of the non-protein substrate sheet, with ``pmids`` added."""
    for row in _checked_rows(opener, 'PSnonprots', NON_PROTEIN_COLUMNS):
        row[REFERENCES] = _text(row[REFERENCES])
        row['pmids'] = pmids(row[REFERENCES])
        yield row


def reference_tokens(cell: str | None) -> list[str]:
    return [token for token in re.split(r'[,;\s]+', cell or '') if token]


def pmids(cell: str | None) -> list[str]:
    """PubMed IDs of a literature cell; the label "UniProt" is not one."""
    return [token for token in reference_tokens(cell) if is_pmid(token)]


def reference_labels(cell: str | None) -> list[str]:
    """Tokens of a literature cell that are not PubMed IDs (the label "UniProt")."""
    return [token for token in reference_tokens(cell) if not token.isdigit()]
