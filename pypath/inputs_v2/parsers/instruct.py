"""
Reader for the INstruct interaction table of human proteins (``sapiens.sin``).

INstruct is a database of protein-protein interactions with the interface
resolved to a pair of Pfam domains, inferred from co-crystal structures. A row
of the tab separated file, which has a header row and 14 columns, is one pair of
proteins and one domain pair:

    ProtA[uniprot] ProtB[uniprot] ProtA[EntrezID] ProtB[EntrezID]
    ProtA[Official Symbol] ProtB[Official Symbol] Pfam-domainA Pfam-domainB
    domainA domainB domainA_loc domainB_loc PDB-structs Publications

``domainA_loc`` is where the domain lies in protein A, as ``(start, end)``
residue numbers of the UniProt sequence; a domain that occurs several times in
the protein has all its locations, separated by a semicolon. ``PDB-structs``
and ``Publications`` list, separated by semicolons, the PDB entries that
support the domain pair and the publications of the interaction.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Generator
import logging
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = ['COLUMNS', 'iter_interactions', 'parse_locations']

_log = logging.getLogger(__name__)

# Header of the file and the field each column is read into, in file order.
COLUMNS = (
    ('ProtA[uniprot]', 'uniprot_a'),
    ('ProtB[uniprot]', 'uniprot_b'),
    ('ProtA[EntrezID]', 'entrez_a'),
    ('ProtB[EntrezID]', 'entrez_b'),
    ('ProtA[Official Symbol]', 'symbol_a'),
    ('ProtB[Official Symbol]', 'symbol_b'),
    ('Pfam-domainA', 'pfam_a'),
    ('Pfam-domainB', 'pfam_b'),
    ('domainA', 'domain_a'),
    ('domainB', 'domain_b'),
    ('domainA_loc', 'locations_a'),
    ('domainB_loc', 'locations_b'),
    ('PDB-structs', 'pdb_ids'),
    ('Publications', 'publications'),
)
_ACCESSION_RE = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
)
_PFAM_RE = re.compile(r'PF\d{5}')
_PDB_RE = re.compile(r'[0-9][A-Za-z0-9]{3}')
_LOCATION_RE = re.compile(r'\((\d+), (\d+)\)')


def parse_locations(cell: str) -> list[tuple[int, int]]:
    """Read ``(131, 183)`` or ``(810, 893);(605, 693)`` into (start, end) pairs.

    The pairs keep the order of the cell. A cell of another form, a start below
    1 or an end before the start is an error.
    """
    locations = []
    for part in cell.split(';'):
        match = _LOCATION_RE.fullmatch(part.strip())
        if match is None:
            raise ValueError(f'INstruct: unexpected domain location {cell!r}')
        start, end = int(match[1]), int(match[2])
        if start < 1 or end < start:
            raise ValueError(f'INstruct: unexpected domain location {cell!r}')
        locations.append((start, end))
    return locations


def iter_interactions(
    opener,
    stats: Counter | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per row (a protein pair with a domain pair).

    The fields are those named in :data:`COLUMNS`, with these differences:
    ``locations_a`` and ``locations_b`` are lists of (start, end) pairs;
    ``pdb_ids`` and ``publications`` are lists of the items as written (empty for
    an empty cell, and ``publications`` can hold text that is not a PubMed
    ID: the mapper decides what a publication is).

    A file without the expected columns, an accession, Pfam accession, Entrez
    Gene ID or PDB entry of another form, or no rows at all (the shared reader
    raises) is an error.
    ``stats``, if given, is filled with the count ``rows`` and the counts
    ``without_pdb_entries`` and ``multiple_locations`` (rows whose domain has more
    than one location in either protein), which are also logged.
    """
    stats = Counter() if stats is None else stats
    rows = iter_tsv(
        opener,
        quotechar=None,
        strict=True,
        required=[name for name, _ in COLUMNS],
    )
    for row in rows:
        stats['rows'] += 1
        record: dict[str, Any] = {field: row[name].strip() for name, field in COLUMNS}
        for side in 'ab':
            if not _ACCESSION_RE.fullmatch(record[f'uniprot_{side}']):
                raise ValueError(
                    f'INstruct: not a UniProt accession: {record[f"uniprot_{side}"]!r}'
                )
            if not record[f'entrez_{side}'].isdigit():
                raise ValueError(
                    f'INstruct: not an Entrez Gene ID: {record[f"entrez_{side}"]!r}'
                )
            if not _PFAM_RE.fullmatch(record[f'pfam_{side}']):
                raise ValueError(
                    f'INstruct: not a Pfam accession: {record[f"pfam_{side}"]!r}'
                )
            record[f'locations_{side}'] = parse_locations(record[f'locations_{side}'])
        record['pdb_ids'] = [i for i in record['pdb_ids'].split(';') if i]
        for pdb_id in record['pdb_ids']:
            if not _PDB_RE.fullmatch(pdb_id):
                raise ValueError(f'INstruct: not a PDB entry: {pdb_id!r}')
        record['publications'] = [i for i in record['publications'].split(';') if i]
        if not record['pdb_ids']:
            stats['without_pdb_entries'] += 1
        if len(record['locations_a']) > 1 or len(record['locations_b']) > 1:
            stats['multiple_locations'] += 1
        yield record
    _log.info('INstruct table: %s', dict(sorted(stats.items())))
