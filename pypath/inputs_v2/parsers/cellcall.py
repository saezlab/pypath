"""
Read the CellCall ligand-receptor-transcription factor tables.

Source files (tab separated, with a header row; a core and an extended table
for human and for the mouse orthologs):

    new_ligand_receptor_TFs.txt
    new_ligand_receptor_TFs_extended.txt
    new_ligand_receptor_TFs_homology.txt
    new_ligand_receptor_TFs_homology_extended.txt

A source row is one axis: a ligand-receptor pair and one transcription factor
downstream of it in a branch of a KEGG pathway, so a pair repeats once for
every transcription factor. The reader folds the rows of a pair into one.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

# Columns of the four tables, in source order.
COLUMNS = (
    'Ligand_ID',
    'Receptor_ID',
    'TF_ID',
    'Pathway',
    'pathway_ID',
    'Ligand_Symbol',
    'Receptor_Symbol',
    'TF_Symbol',
)
PAIR_COLUMNS = ('Ligand_ID', 'Ligand_Symbol', 'Receptor_ID', 'Receptor_Symbol')
# ``Pathway`` lists the branches (``hsa04630_2``) of the KEGG pathways in
# ``pathway_ID`` that connect the pair to the transcription factor.
TF_COLUMNS = ('TF_ID', 'TF_Symbol', 'Pathway', 'pathway_ID')
# The subunits of a multi-subunit receptor are separated by commas, in the
# same order in the ID and in the symbol column (``3454,3455`` and
# ``IFNAR1,IFNAR2``).
SUBUNIT_DELIMITER = ','


def iter_pairs(
    opener,
    *,
    taxon_id: str,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield one row per ligand-receptor pair of a CellCall table.

    The table repeats a pair once for every downstream transcription factor.
    A yielded row has the four ligand and receptor columns under their own
    names, ``taxon_id`` (the NCBI Taxonomy ID of the file) and
    ``transcription_factors``: the distinct ``TF_ID``, ``TF_Symbol``,
    ``Pathway`` and ``pathway_ID`` cells of the pair's rows, in source order.
    Pairs are yielded in the order of their first row.
    """
    pairs: dict[tuple[str, ...], dict[str, Any]] = {}
    for row in iter_tsv(opener):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(
                f'CellCall table lacks expected columns: {", ".join(missing)}'
            )
        key = tuple((row[column] or '').strip() for column in PAIR_COLUMNS)
        pair = pairs.get(key)
        if pair is None:
            pair = dict(zip(PAIR_COLUMNS, key, strict=True))
            for role in ('Ligand', 'Receptor'):
                ids = pair[f'{role}_ID'].split(SUBUNIT_DELIMITER)
                symbols = pair[f'{role}_Symbol'].split(SUBUNIT_DELIMITER)
                if len(ids) != len(symbols):
                    raise ValueError(
                        f'CellCall {role.lower()} has unaligned subunits: '
                        f'{pair[f"{role}_ID"]!r} and {pair[f"{role}_Symbol"]!r}'
                    )
            pair['taxon_id'] = taxon_id
            pair['transcription_factors'] = []
            pairs[key] = pair
        factor = {column: (row[column] or '').strip() for column in TF_COLUMNS}
        # The mouse core table repeats many rows exactly.
        if factor not in pair['transcription_factors']:
            pair['transcription_factors'].append(factor)
    if not pairs:
        raise ValueError('CellCall table has no rows')
    yield from pairs.values()
