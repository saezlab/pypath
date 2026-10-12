"""
Parser of Supplementary Table S12 of Kirouac et al. 2010 (Molecular Systems Biology).

The table is the sheet `S12` of the workbook `msb201071-s3.xls`. Its first row is
the title, its second the header (`LIGAND`, `RECEPTOR(S)`), and every row after
it describes one ligand: a free-text ligand cell, usually `description (SYMBOL /
alias)`, and a receptor cell in the authors' shorthand, where `BMPR 1A/1B/2`
stands for BMPR1A, BMPR1B and BMPR2 and `FZD 1-10` for FZD1 to FZD10.

The parser expands the cells into one row for each ligand-receptor pair. A row
has the fields

    LIGAND, RECEPTOR(S)  the two cells the pair comes from, unchanged
    ligand, receptor     one symbol each, spelled as the table spells it
    receptor_synonym     the second designation a receptor cell gives in brackets
                         (`IL8RB (CXCR2)` has the symbol IL8RB and the synonym
                         CXCR2), else None

A ligand with no receptor (the cell is empty or `???`, `????` or `xxx`) gives one
row whose `ligand`, `receptor` and `receptor_synonym` are None. A cell in a notation the grammar does
not know raises.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

# Receptor cells of ligands the table gives no receptor for.
NO_RECEPTOR = frozenset({'???', '????', 'xxx'})
# A symbol in capitals, as the table writes the chemokines.
_capitals_re = re.compile('[A-Z0-9]+')
# A symbol among the bracketed designations of a ligand: TNFSF10, VEGF-D, WNT7b.
_symbol_re = re.compile('[A-Z0-9]+(?:-[A-Z0-9]+)*[a-z]?')
# `Attractin (ATRN)`, `IL8RB (CXCR2)`: a receptor with a second designation.
_aliased_re = re.compile('([A-Za-z0-9]+) \\(([A-Z][A-Z0-9]+)\\)')
# `NOTCH 1-4`, `FGFR1-4`: a numbered series.
_series_re = re.compile('([A-Z]+) ?(\\d+)-(\\d+)')
# `BMPR 1A/1B/2`, `TNFRSF1 A/B`: a stem and its endings.
_endings_re = re.compile('([A-Z0-9]+) ([A-Za-z0-9]+(?:/[A-Za-z0-9]+)+)')
# `IL1R1/2`: the same without the space, with numbers after a stem ending in a letter.
_numbers_re = re.compile('([A-Z0-9]*[A-Z])(\\d+(?:/\\d+)+)')
_plain_re = re.compile('[A-Za-z0-9-]+')


def ligand_symbols(cell: str) -> list[str]:
    """Symbols of a ligand cell.

    The cell ends in a bracket. Before it stands either a description, and
    then the first bracketed designation written as a symbol is taken
    (`KIT ligand (KITLG / stem cell factor)`), or symbols in capitals, and
    then the bracket holds their aliases (`CCL14 / CCL15 (MIP-5 / MIP-1delta)`).
    """
    head, bracket, designations = cell.rpartition('(')
    if not bracket:
        raise ValueError(f'Unknown Kirouac 2010 ligand notation {cell!r}')
    heads = [item.strip() for item in head.split(' / ')]
    if all(_capitals_re.fullmatch(item) for item in heads):
        return heads
    for item in designations.rstrip(') ').split('/'):
        if _symbol_re.fullmatch(item.strip()):
            return [item.strip()]
    raise ValueError(f'Unknown Kirouac 2010 ligand notation {cell!r}')


def receptor_designations(cell: str) -> list[tuple[str, str | None]]:
    """Symbols of a receptor cell, each with its second designation or None.

    The cell has alternatives separated by ` / `, each a symbol, a symbol with
    a second designation in brackets (`IL8RB (CXCR2)`, `Attractin (ATRN)`), a
    numbered series or a stem with several endings. Of two designations the one
    in capitals is the symbol (the bracketed one if neither is); the other is
    the synonym.
    """
    designations: list[tuple[str, str | None]] = []
    for item in cell.split(' / '):
        item = item.strip()
        if match := _aliased_re.fullmatch(item):
            first, second = match.groups()
            if _capitals_re.fullmatch(first):
                designations.append((first, second))
            else:
                designations.append((second, first))
        elif match := _series_re.fullmatch(item):
            stem, start, end = match.groups()
            designations.extend(
                (f'{stem}{number}', None)
                for number in range(int(start), int(end) + 1)
            )
        elif match := _endings_re.fullmatch(item):
            stem, endings = match.groups()
            designations.extend(
                (stem + ending, None) for ending in endings.split('/')
            )
        elif match := _numbers_re.fullmatch(item):
            stem, endings = match.groups()
            designations.extend(
                (stem + ending, None) for ending in endings.split('/')
            )
        elif _plain_re.fullmatch(item):
            designations.append((item, None))
        else:
            raise ValueError(f'Unknown Kirouac 2010 receptor notation {cell!r}')
    return designations


def receptor_symbols(cell: str) -> list[str]:
    """Symbols of a receptor cell (see `receptor_designations`)."""
    return [symbol for symbol, _ in receptor_designations(cell)]


def iter_pairs(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One row per ligand-receptor pair of Table S12, with the cells it comes from.

    A ligand without a receptor yields one row with no `ligand` and no `receptor`.
    """
    for row in iter_excel(opener, sheet='S12', header_row=2, **kwargs):
        ligands, receptors = row['LIGAND'], row['RECEPTOR(S)']
        cells = {'LIGAND': ligands, 'RECEPTOR(S)': receptors}
        if not ligands or not receptors or receptors in NO_RECEPTOR:
            yield {
                **cells,
                'ligand': None,
                'receptor': None,
                'receptor_synonym': None,
            }
            continue
        for ligand in ligand_symbols(ligands):
            for receptor, synonym in receptor_designations(receptors):
                yield {
                    **cells,
                    'ligand': ligand,
                    'receptor': receptor,
                    'receptor_synonym': synonym,
                }
