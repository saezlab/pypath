"""
Parsers of the PanglaoDB cell type marker table.

``PanglaoDB_markers_27_Mar_2020.tsv.gz`` has a header row and one row for each
marker gene of a cell type:

    species                 'Mm Hs', 'Mm' or 'Hs': the species the marker is
                            stated for (see below for the other values)
    official gene symbol    the symbol, upper case also for a mouse gene
    cell type               PanglaoDB's own cell type label
    nicknames               aliases of the gene, ``|`` separated
    ubiquitousness index    of the gene, 0 to 1 (the paper: how often the gene is
                            expressed in cell clusters)
    product description     name of the gene product
    gene type               'protein-coding gene' and others
    canonical marker        1 or 0, or NA
    germ layer, organ       of the cell type
    sensitivity_human, sensitivity_mouse, specificity_human, specificity_mouse
                            one value for each species; NA where the marker is not
                            stated for the species

``NA`` is the table's missing value. The rows are read twice: as the statements
that a gene is a marker of a cell type, one for each species the row states
(:func:`iter_markers`), and as the vocabulary of cell types, organs and germ
layers that the statements use (:func:`iter_terms`).

Defects of the table that are repaired here, each counted in the audit of the
2020-03-27 file:

* two symbols have a trailing space (``ALCAM ``, ``PTGDS ``): stripped;
* the germ layer of the 17 rows of the cell type ``Myeloid-derived suppressor
  cells`` is ``mesoderm``, the other 4,810 mesoderm rows have ``Mesoderm``: the
  label is written like the others.

The table writes every symbol in upper case, also where the statement is about
a mouse gene, and mouse gene symbols are not upper case (``Try4``, ``Cd36``). A
symbol lookup is exact, so the upper-case form of a mouse gene would find
nothing or, where it happens to be an alias of another mouse gene (``AR`` is an
alias of ``Areg``, the gene is ``Ar``), the wrong gene. The statement of a
mouse gene therefore gets the symbol in the mouse form (:func:`mouse_symbol`),
the rule of the legacy input. Against the NCBI Gene symbols of mouse
(2026-10-13) this gives the official symbol for 7,751 of the 7,909 mouse
statements; of the rest, 149 are symbols mouse does not have in any case (old
names such as ``WISP2``, or genes without a mouse counterpart such as
``MNDA``) and 9 are symbols with a lower-case part the rule cannot know
(``H2-DMa``, ``Hbb-bs``). The human statement of the same row keeps the symbol
as written.

Rows that state no species are yielded with ``taxon`` None and left to the
mapper to exclude: 30 rows have the species ``None`` and 3 rows have ``4``.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = (
    'species',
    'official gene symbol',
    'cell type',
    'ubiquitousness index',
    'canonical marker',
    'germ layer',
    'organ',
    'sensitivity_human',
    'sensitivity_mouse',
    'specificity_human',
    'specificity_mouse',
)

# The species flags the table writes and the NCBI Taxonomy ID of each.
SPECIES = {'Hs': '9606', 'Mm': '10090'}
# Values of the species column that state no species: the literal ``None``
# (30 rows) and the number ``4`` (3 rows), both with the missing value in the
# columns of the product and, for ``4``, in all four sensitivity and
# specificity columns.
NO_SPECIES = frozenset({'None', '4'})
# The germ layer written in lower case by the source in one cell type.
GERM_LAYER_REPAIRS = {'mesoderm': 'Mesoderm'}


def _value(cell: str | None) -> str | None:
    """A cell of the table; ``NA`` and empty are missing."""
    if cell is None:
        return None
    cell = cell.strip()
    return None if cell in ('', 'NA') else cell


def _rows(opener) -> Generator[dict[str, Any], None, None]:
    for row in iter_tsv(opener, required=COLUMNS):
        symbol = _value(row['official gene symbol'])
        cell_type = _value(row['cell type'])
        if symbol is None or cell_type is None:
            raise ValueError(f'PanglaoDB row without gene or cell type: {row}')
        germ_layer = _value(row['germ layer'])
        yield {
            'species': row['species'].strip(),
            'gene': symbol,
            'cell_type': cell_type,
            'organ': _value(row['organ']),
            'germ_layer': GERM_LAYER_REPAIRS.get(germ_layer, germ_layer),
            'canonical_marker': _value(row['canonical marker']),
            'ubiquitousness_index': _value(row['ubiquitousness index']),
            'sensitivity': {
                '9606': _value(row['sensitivity_human']),
                '10090': _value(row['sensitivity_mouse']),
            },
            'specificity': {
                '9606': _value(row['specificity_human']),
                '10090': _value(row['specificity_mouse']),
            },
        }


MOUSE = '10090'


def mouse_symbol(symbol: str) -> str:
    """An upper-case symbol of the table in the form of a mouse gene symbol.

    Each part between hyphens gets a capital first letter and lower case after
    it: ``TRY4`` is ``Try4``, ``H2-M2`` stays ``H2-M2``.
    """
    return '-'.join(part.capitalize() for part in symbol.split('-'))


def iter_markers(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One statement for each species a row states.

    Fields: ``gene`` (the symbol, stripped: as written for human, in the mouse
    form for mouse, see :func:`mouse_symbol`), ``taxon`` (NCBI Taxonomy
    ID, ``9606`` or ``10090``; None where the row states no species),
    ``cell_type``, ``organ``, ``germ_layer`` (None where NA),
    ``canonical_marker`` (``1``, ``0`` or None), ``ubiquitousness_index``,
    ``sensitivity`` and ``specificity`` (of the species' columns; None where NA),
    ``species`` (the flags as written). An unknown species flag raises.
    """
    for row in _rows(opener):
        species = row.pop('species')
        sensitivity = row.pop('sensitivity')
        specificity = row.pop('specificity')
        if species in NO_SPECIES:
            yield {**row, 'species': species, 'taxon': None,
                   'sensitivity': None, 'specificity': None}
            continue
        flags = species.split()
        unknown = [flag for flag in flags if flag not in SPECIES]
        if unknown or not flags:
            raise ValueError(f'Unknown PanglaoDB species value: {species!r}')
        for flag in flags:
            taxon = SPECIES[flag]
            yield {
                **row,
                'gene': mouse_symbol(row['gene']) if taxon == MOUSE else row['gene'],
                'species': species,
                'taxon': taxon,
                'sensitivity': sensitivity[taxon],
                'specificity': specificity[taxon],
            }


def iter_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The cell types, organs and germ layers of the table, each once.

    Fields: ``kind`` (``cell_type``, ``organ`` or ``germ_layer``), ``id`` (the
    label), and for a cell type its ``organ`` and ``germ_layer`` (None where NA).
    A cell type with two organs or germ layers raises: the table gives one of
    each to every cell type.
    """
    cell_types: dict[str, tuple[str | None, str | None]] = {}
    organs: dict[str, None] = {}
    germ_layers: dict[str, None] = {}
    for row in _rows(opener):
        place = (row['organ'], row['germ_layer'])
        if cell_types.setdefault(row['cell_type'], place) != place:
            raise ValueError(
                f'PanglaoDB cell type {row["cell_type"]!r} has two organs or germ layers'
            )
        if row['organ']:
            organs[row['organ']] = None
        if row['germ_layer']:
            germ_layers[row['germ_layer']] = None
    if not cell_types:
        raise ValueError('PanglaoDB marker table has no rows')
    for germ_layer in germ_layers:
        yield {'kind': 'germ_layer', 'id': germ_layer}
    for organ in organs:
        yield {'kind': 'organ', 'id': organ}
    for cell_type, (organ, germ_layer) in cell_types.items():
        yield {'kind': 'cell_type', 'id': cell_type, 'organ': organ, 'germ_layer': germ_layer}
