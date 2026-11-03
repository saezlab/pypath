"""
Reader for Table S3A of Havugimana et al. 2012 (Cell 150:1068), the 622 complexes.

The supplement ``mmc3.xls`` has four sheets; the first, ``Predicted_complexes``
(Table S3A), lists the complexes. The other three (the merged CORUM reference
complexes used for the comparison and the matches with yeast and fly complexes)
are not read.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

__all__ = ['SHEET', 'HEADER', 'iter_havugimana_complexes']

SHEET = 'Predicted_complexes'
# The header of the sheet is its third row. The spelling is the supplement's,
# typos included.
HEADER = (
    'Predicted complex',
    'Size (# of subunits)',
    'Co-complex membership_UniProtKB AC',
    'Co-complex membership_UniProtKB EntryName',
    '# of matched subunit in public DB',
    '"Novel" proteins',
    'Annotated proteins',
    'Best-matching CORUM complex',
    'Best-matchng REACTOME complex',
    'Best-mathng PINdb complex',
    'Best-matching HPRD complex',
    'Best-mathng Gene Ontology Cell Component complex',
    'Best-matching MEMO complex, Malovannaya et al',
    'Best-matching FLY complex (orthology), Guruharsha et al.',
    'Best-matching YEAST complex (orthology), Babu el al.,',
    'Taxonomic group',
)
# ``Best-matching`` columns -> key of the parsed row. The value is
# ``<id>(::|~|\|)<name> (<p-value>)`` or ``-``.
_MATCHES = {
    'Best-matching CORUM complex': 'best_corum',
    'Best-matchng REACTOME complex': 'best_reactome',
    'Best-mathng PINdb complex': 'best_pindb',
    'Best-matching HPRD complex': 'best_hprd',
    'Best-mathng Gene Ontology Cell Component complex': 'best_go',
    'Best-matching MEMO complex, Malovannaya et al': 'best_memo',
    'Best-matching FLY complex (orthology), Guruharsha et al.': 'best_fly',
    'Best-matching YEAST complex (orthology), Babu el al.,': 'best_yeast',
}

_complex_id = re.compile(r'C_\d+')
_accession = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
)
_entry_name = re.compile(r'[A-Z0-9]+_HUMAN')


def _list(value: Any) -> list[str]:
    """Items of a comma-separated cell, which may end with a comma; ``-`` is none."""
    items = [item.strip() for item in str(value or '').split(',') if item.strip()]
    return [] if items == ['-'] else items


def iter_havugimana_complexes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Complexes of Table S3A, one dict each.

    Keys:

    - ``complex_id``: ``C_1`` to ``C_622``, the number the table gives;
    - ``size``: the number of subunits (int);
    - ``accessions``, ``entry_names``: the UniProtKB accessions of the subunits
      and the UniProtKB entry names (all ``_HUMAN``) in the same order;
    - ``matched_subunits``: number of subunits found in a public complex
      database (int);
    - ``novel``, ``annotated``: the subunits the table calls novel, that is
      in no public complex, and annotated (both lists of accessions, sorted
      as the table has them; together the subunits);
    - ``best_corum``, ``best_reactome``, ``best_pindb``, ``best_hprd``,
      ``best_go``, ``best_memo``, ``best_fly``, ``best_yeast``: the best
      matching complex of each reference with its p-value, as written, or None;
    - ``taxonomic_group``: as written (``NA`` for one complex).

    Any complex whose size, names, id or accessions are inconsistent is an
    error, and so is a sheet with another header or without complexes.
    """
    rows = 0
    for row in iter_excel(opener, sheet=SHEET, header_row=3):
        if tuple(row) != HEADER:
            raise ValueError(f'Unexpected header of {SHEET}: {tuple(row)}')
        rows += 1
        complex_id = str(row[HEADER[0]] or '').strip()
        if not _complex_id.fullmatch(complex_id):
            raise ValueError(f'Unknown complex id: {complex_id!r}')
        accessions = _list(row[HEADER[2]])
        names = _list(row[HEADER[3]])
        size = row[HEADER[1]]
        if (
            not accessions
            or len(accessions) != len(names)
            or size != len(accessions)
            or len(set(accessions)) != len(accessions)
        ):
            raise ValueError(f'{complex_id}: size, accessions and entry names disagree')
        if not all(_accession.fullmatch(a) for a in accessions):
            raise ValueError(f'{complex_id}: malformed UniProt accession')
        if not all(_entry_name.fullmatch(n) for n in names):
            raise ValueError(f'{complex_id}: entry name of another species or malformed')
        novel = _list(row[HEADER[5]])
        annotated = _list(row[HEADER[6]])
        if sorted(novel + annotated) != sorted(accessions):
            raise ValueError(f'{complex_id}: novel and annotated proteins are not the subunits')
        parsed = {
            'complex_id': complex_id,
            'size': int(size),
            'accessions': accessions,
            'entry_names': names,
            'matched_subunits': int(row[HEADER[4]]),
            'novel': novel,
            'annotated': annotated,
            'taxonomic_group': row[HEADER[15]],
        }
        for column, key in _MATCHES.items():
            value = str(row[column] or '').strip()
            parsed[key] = None if value in ('', '-') else value
        yield parsed
    if not rows:
        raise ValueError(f'{SHEET} has no complexes')
