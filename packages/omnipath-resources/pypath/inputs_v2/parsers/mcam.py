"""Read the list of human cell adhesion molecules of MCAM.

The Mammalian Cell Adhesion Molecule database offers its human list as a
download named ``exel.cfm?seq=HumanUniprot``. Despite the name of the page it
is a tab-separated text file with carriage returns as line ends, a header and
one row per protein, with a tab after the last field:

    ACC  ENTNAME  FIELD2  HUNIPORTID  LENGTH  MATCHEDQUERYID  NF50  NF90  ORG  TAXGRP
    P78325  ADAM 8 precursor; CD156; ...  ADAM8  ADAM8_HUMAN  824  ADAM8_HUMAN,   UniRef50_P78325  UniRef90_P78325  Homo sapiens  Euk/mammal

``ACC`` holds the accessions of the UniProt entry as UniProt writes them, the
primary one first (``Q86WK6; Q8IW71; Q9ULQ7``). ``HUNIPORTID`` is the UniProt
entry name.
"""

from __future__ import annotations

from collections.abc import Generator
from functools import partial
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv, read_text

COLUMNS = (
    'ACC',
    'ENTNAME',
    'FIELD2',
    'HUNIPORTID',
    'LENGTH',
    'MATCHEDQUERYID',
    'NF50',
    'NF90',
    'ORG',
    'TAXGRP',
)

ORGANISM = 'Homo sapiens'

_table = partial(iter_tsv, quotechar=None, required=COLUMNS)


def _rows(opener) -> Generator[dict[str, Any], None, None]:
    # The lines end in a carriage return alone, which a file handle may or may
    # not translate; ``splitlines`` splits at it either way.
    lines = read_text(opener).splitlines()
    return _table(SimpleNamespace(result=lines))


def iter_proteins(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One record per protein: ``accessions`` (a list, the primary one first),
    ``entry_name``, ``organism`` (``Homo sapiens`` for every row)."""
    seen: set[str] = set()
    for row in _rows(opener):
        # The tab after the last field gives each row one empty field more.
        if row.pop(None, None) != ['']:
            raise ValueError(f'MCAM row with unexpected fields: {row!r}')
        accessions = [accession.strip() for accession in row['ACC'].split(';')]
        if not all(accessions):
            raise ValueError(f'Empty accession in {row["ACC"]!r}')
        if row['ORG'] != ORGANISM:
            raise ValueError(f'{accessions[0]}: organism {row["ORG"]!r}')
        if accessions[0] in seen:
            raise ValueError(f'MCAM accession {accessions[0]!r} occurs twice')
        seen.add(accessions[0])
        yield {
            'accessions': accessions,
            'entry_name': row['HUNIPORTID'],
            'organism': row['ORG'],
        }
    if not seen:
        raise ValueError('The MCAM list has no proteins')


def iter_list(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """A single empty record for the list term, once the file has been read and
    found to have proteins."""
    # The whole file is read, so that a malformed one fails here as well.
    if list(iter_proteins(opener, **kwargs)):
        yield {}
