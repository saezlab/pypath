"""
Reader for the SIFTS chain-to-UniProt table (``pdb_chain_uniprot.tsv.gz``).

SIFTS (Structure Integration with Function, Taxonomy and Sequence, PDBe and
UniProt) maps the polymer chains of the PDB structures to UniProt. The file
has a comment line with the release (``# 2026/10/04 - 08:26 | PDB: 40.26 |
UniProt: 2026.03``), a header row, and one row per mapped segment of a chain:

    PDB  CHAIN  SP_PRIMARY  RES_BEG  RES_END  PDB_BEG  PDB_END  SP_BEG  SP_END

A chain that has segments from several UniProt entries (a fusion or a
chimera) has a row for each. The rows of an entry are contiguous.

:func:`iter_pdb_entries` groups the rows by entry. What the table states is
which UniProt proteins the chains of an entry map to; whether the proteins
form a complex in the cell is not in it.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['COLUMNS', 'iter_pdb_entries']

COLUMNS = (
    'PDB',
    'CHAIN',
    'SP_PRIMARY',
    'RES_BEG',
    'RES_END',
    'PDB_BEG',
    'PDB_END',
    'SP_BEG',
    'SP_END',
)

_release = re.compile(r'# (\d{4}/\d\d/\d\d) - \d\d:\d\d \| PDB: (\S+) \| UniProt: (\S+)')
_pdb_id = re.compile(r'[0-9][a-z0-9]{3}')
_accession = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}(?:-\d+)?'
)


def iter_pdb_entries(
    opener, complexes_only: bool = True, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """The PDB entries of the table, one dict each.

    Keys:

    - ``pdb_id``: the four-character PDB identifier, as written (lower case);
    - ``release``: the comment line of the file (date, PDB and UniProt releases);
    - ``accessions``: the distinct UniProt accessions of the entry, in the
      order of their first row;
    - ``chains``: the chain identifiers of the entry, in file order;
    - ``chimeric_chains``: the chains that map to more than one accession.

    With ``complexes_only`` (the default) an entry is yielded when its chains
    do not all map to the same proteins: at least two chains differ in the
    UniProt accessions they carry, so that the entry has at least two distinct
    UniProt proteins found on separate polypeptides. An entry with one protein,
    one chain made of segments of several proteins (a fusion), or chains that
    all carry the same fusion, is not. All entries are yielded with
    ``complexes_only=False``.

    A file without the comment line and the header, a row of another width, a
    malformed identifier or accession, an entry whose rows are not contiguous,
    and a file without entries are errors.
    """
    names = iter_tsv(opener, fieldnames=COLUMNS, quotechar=None)
    comment = next(names, None)
    header = next(names, None)
    if comment is None or header is None:
        raise ValueError('SIFTS table is empty')
    match = _release.fullmatch((comment['PDB'] or '').strip())
    if not match or comment['CHAIN'] is not None:
        raise ValueError(f'SIFTS table does not start with a release line: {comment}')
    if tuple(header.values()) != COLUMNS:
        raise ValueError(f'Unexpected SIFTS header: {tuple(header.values())}')
    release = (comment['PDB'] or '').strip().lstrip('# ')

    seen: set[str] = set()
    current: str | None = None
    chains: dict[str, list[str]] = {}
    entries = 0

    def entry() -> dict[str, Any] | None:
        accessions = list(dict.fromkeys(a for v in chains.values() for a in v))
        if complexes_only and len({frozenset(v) for v in chains.values()}) < 2:
            return None
        return {
            'pdb_id': current,
            'release': release,
            'accessions': accessions,
            'chains': list(chains),
            'chimeric_chains': [c for c, v in chains.items() if len(set(v)) > 1],
        }

    for row in names:
        if None in row or any(row[c] is None for c in COLUMNS):
            raise ValueError(f'SIFTS row with another number of fields than 9: {row}')
        pdb_id, chain, accession = row['PDB'], row['CHAIN'], row['SP_PRIMARY']
        if not _pdb_id.fullmatch(pdb_id) or not _accession.fullmatch(accession) or not chain:
            raise ValueError(f'Malformed SIFTS row: {row}')
        if pdb_id != current:
            if current is not None:
                entries += 1
                done = entry()
                if done:
                    yield done
            if pdb_id in seen:
                raise ValueError(f'Rows of PDB entry {pdb_id} are not contiguous')
            seen.add(pdb_id)
            current, chains = pdb_id, {}
        chains.setdefault(chain, []).append(accession)
    if current is None:
        raise ValueError('SIFTS table has no entries')
    done = entry()
    if done:
        yield done
