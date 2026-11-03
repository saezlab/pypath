"""Read the Cell Surface Protein Atlas (CSPA) files of Bausch-Fluck et al. 2015.

The atlas has cell surface glycoproteins of 41 human and 31 mouse cell types
(the paper's count), found by Cell Surface Capture and mass spectrometry. Two
workbooks of the web site are read.

``S2_File.xlsx`` lists the proteins of the atlas: ``Table A`` the 1,492 human
and ``Table B`` the 1,296 mouse proteins, each with a category (``1 - high
confidence``, ``2 - putative``, ``3 - unspecific``) and annotations:

    organism  ID_link  CSPA category  UP_Protein_name  UP_entry_name  CD
    ENTREZ_gene_ID  ENTREZ gene symbol  protein probability
    count detection in different cell types  Phobius TM predicted yes/no
    # predicted TM domains  GPI  UniProt Cell surface

``ID_link`` is the UniProt accession, ``CD`` the cluster of differentiation
name or ``no``. The mouse sheet is blank below its 1,296 rows, and its last
column is spelled ``UniProt cell surface``.

``S1_File.xlsx`` has the detections. Its ``Table_A`` lists every protein
identified in each cell type:

    organism  Experiment tag  ID_link  CD  ENTREZ ac  ENTREZ gene symbol
    UP_entry_name  count detection in different cell types
    protein probability  num unique peps

``Experiment tag`` is the name of the cell type or line. The other sheets of the
workbook hold the same detections as matrices (``Table_B``, ``Table_C``), counts
per cell type (``Table_D``) and the relative quantities of the proteins in log2
scale (``Table_E``, ``Table_F``). The quantity matrices name their cell types
with codes that differ from those of ``Table_A`` (``A431`` against ``A-431``,
``CD4pCD25n_Tcells`` against ``CD4+, CD25- Tcells``), and their detections are
not those of the other sheets, so they are not read.
"""

from __future__ import annotations

from collections.abc import Generator
from functools import partial
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

# Sheet and the spelling of the UniProt column, by organism.
ORGANISMS = {
    'Human': ('Table A', 'UniProt Cell surface'),
    'Mouse': ('Table B', 'UniProt cell surface'),
}
CATEGORIES = ('1 - high confidence', '2 - putative', '3 - unspecific')

PROTEIN_COLUMNS = (
    'organism',
    'ID_link',
    'CSPA category',
    'UP_Protein_name',
    'UP_entry_name',
    'CD',
    'ENTREZ_gene_ID',
    'ENTREZ gene symbol',
    'protein probability',
    'count detection in different cell types',
    'Phobius TM predicted yes/no',
    '# predicted TM domains',
    'GPI',
)
DETECTION_COLUMNS = (
    'organism',
    'Experiment tag',
    'ID_link',
    'ENTREZ ac',
    'ENTREZ gene symbol',
    'UP_entry_name',
    'protein probability',
    'num unique peps',
)

# The kinds of record of the terms: the list of the atlas, and a cell type.
LIST = 'list'
CELL_TYPE = 'cell_type'


def _flag(value: Any, field: str) -> bool:
    """A 0/1 column."""
    if value not in (0, 1):
        raise ValueError(f'{field} is not 0 or 1: {value!r}')
    return bool(value)


def _entrez(value: Any) -> str | None:
    """The Entrez Gene ID; the sheets write 0 (as for the symbol) for a protein
    that has no gene in Entrez Gene."""
    return None if value in (None, 0, '0') else str(int(value))


def _symbol(value: Any) -> str | None:
    """The gene symbol; the sheets write 0 (a number, or text) for a protein that
    has none (``P10404``, ``Q9D7S0`` of the mouse sheet)."""
    if value is None or value in (0, '0'):
        return None
    return str(value)


def _uniprot_cell_surface(value: Any) -> bool:
    """``yes`` where UniProt annotates the protein as cell surface; the sheet
    writes ``#N/A``, an Excel lookup that found nothing, otherwise."""
    if value == 'yes':
        return True
    if value == '#N/A':
        return False
    raise ValueError(f'Unknown value of the UniProt cell surface column: {value!r}')


def iter_surface_proteins(
    opener, organism: str = 'Human', **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """The proteins of the atlas for one organism, from ``S2_File.xlsx``.

    Fields: ``organism`` (the argument), ``uniprot``, ``entry_name``, ``name``, ``cd`` (None unless the
    protein has a CD name), ``entrez`` and ``gene`` (None where the sheet has none or writes 0),
    ``category``, ``probability`` (of the protein identification, 0 to 1),
    ``n_cell_types`` (number of cell types the protein was detected on),
    ``phobius_tm`` (yes/no), ``tm_domains``, ``gpi`` (booleans) and
    ``uniprot_cell_surface`` (see :func:`_uniprot_cell_surface`).
    """
    sheet, surface_column = ORGANISMS[organism]
    seen: set[str] = set()
    for row in iter_excel(
        opener,
        sheet=sheet,
        required=(*PROTEIN_COLUMNS, surface_column),
    ):
        if row['organism'] != organism:
            raise ValueError(f'{sheet} has a protein of organism {row["organism"]!r}')
        accession = row['ID_link']
        if accession in seen:
            raise ValueError(f'UniProt accession {accession!r} occurs twice in {sheet}')
        seen.add(accession)
        category = row['CSPA category']
        if category not in CATEGORIES:
            raise ValueError(f'{accession}: unknown CSPA category {category!r}')
        yield {
            'organism': organism,
            'uniprot': accession,
            'entry_name': row['UP_entry_name'],
            'name': row['UP_Protein_name'],
            'cd': None if row['CD'] == 'no' else row['CD'],
            'entrez': _entrez(row['ENTREZ_gene_ID']),
            'gene': _symbol(row['ENTREZ gene symbol']),
            'category': category,
            'probability': row['protein probability'],
            'n_cell_types': int(row['count detection in different cell types']),
            'phobius_tm': _flag(row['Phobius TM predicted yes/no'], 'Phobius TM'),
            'tm_domains': int(row['# predicted TM domains']),
            'gpi': _flag(row['GPI'], 'GPI'),
            'uniprot_cell_surface': _uniprot_cell_surface(row[surface_column]),
        }


iter_human_proteins = partial(iter_surface_proteins, organism='Human')
iter_mouse_proteins = partial(iter_surface_proteins, organism='Mouse')


def _detections(opener) -> Generator[dict[str, Any], None, None]:
    for row in iter_excel(opener, sheet='Table_A', required=DETECTION_COLUMNS):
        if row['organism'] not in ORGANISMS:
            raise ValueError(f'Table_A has a detection of organism {row["organism"]!r}')
        yield row


def iter_detections(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The identifications of proteins in cell types, from ``S1_File.xlsx``.

    Fields: ``organism`` (``Human`` or ``Mouse``), ``cell_type`` (the
    experiment tag), ``uniprot``, ``entry_name``, ``gene``, ``entrez``,
    ``probability`` (of the protein identification in this cell type, 0.77 to 1
    in the sheet) and ``peptides`` (the number of unique peptides).
    """
    seen: set[tuple[str, str, str]] = set()
    for row in _detections(opener):
        key = (row['organism'], row['Experiment tag'], row['ID_link'])
        if key in seen:
            raise ValueError(f'Detection {key} occurs twice')
        seen.add(key)
        yield {
            'organism': row['organism'],
            'cell_type': row['Experiment tag'],
            'uniprot': row['ID_link'],
            'entry_name': row['UP_entry_name'],
            'gene': _symbol(row['ENTREZ gene symbol']),
            'entrez': _entrez(row['ENTREZ ac']),
            'probability': row['protein probability'],
            'peptides': int(row['num unique peps']),
        }


def iter_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The list of the atlas, then the cell types of ``Table_A`` in order of
    appearance, each once. A cell type has ``organism`` and ``cell_type``.

    A cell type of the human and one of the mouse experiments never have the
    same name, so a name stands for one cell type.
    """
    yield {'kind': LIST}
    organism_of: dict[str, str] = {}
    for row in _detections(opener):
        name = row['Experiment tag']
        if organism_of.setdefault(name, row['organism']) != row['organism']:
            raise ValueError(f'Cell type {name!r} occurs in two organisms')
    for name, organism in organism_of.items():
        yield {'kind': CELL_TYPE, 'cell_type': name, 'organism': organism}
