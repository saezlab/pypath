"""Read the data of Membranome (single-pass membrane proteins and their membranes).

Membranome serves its data from a public API (``opm-back.cc.lehigh.edu/
membranome-backend``, documented in ``swagger-membranome-api.json`` of the
site), as JSON and, for the tables of its download page, as CSV. Three
tables are read.

``membranes`` (JSON, one page of the 24 membranes)::

    objects: id  name  topology_in  topology_out  abbreviation  proteins_count ...

``topology_in`` and ``topology_out`` name the two sides of the membrane in the
membrane's own words (``cytoplasmic side`` and ``extracellular side`` for the
plasma membrane).

``proteins`` (CSV, 16,697 proteins of 18 species)::

    id  name  uniprot_id  uniprotcode  genename  topology_show_in
    membrane_id  species_name_cache ...

``uniprotcode`` is the UniProt accession and ``uniprot_id`` the entry name
(``GNPTA_HUMAN_alpha`` for the two chains of a protein that the database
lists as two entries). ``membrane_id`` is the membrane the protein is in, and
``topology_show_in`` says which side of that membrane the protein's topology
is shown on: true for the ``topology_in`` side, false for ``topology_out``,
empty where the database gives none. Of the proteins whose first topological
domain in the UniProt comments is known, the flag is true for those that start
on the cytoplasmic side and false for those that start outside.

``protein_membranes`` (CSV, 18,549 rows; the download page calls it
``protein_localizations``)::

    id  membrane_name_cache  protein_uniprot_id_cache  pubmed  main
    membrane_id  protein_id

It lists every membrane a protein is in, with the publication of the evidence
(empty for most) and a flag ``main`` that does not always agree with
``membrane_id`` of the protein.
"""

from __future__ import annotations

from collections.abc import Generator
import json
from typing import Any

from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.base import iter_csv, read_text

MEMBRANE_COLUMNS = ('id', 'name', 'topology_in', 'topology_out', 'abbreviation')
PROTEIN_COLUMNS = (
    'id', 'name', 'uniprot_id', 'uniprotcode', 'genename', 'topology_show_in',
    'membrane_id', 'species_name_cache',
)
LOCALIZATION_COLUMNS = (
    'membrane_name_cache', 'protein_uniprot_id_cache', 'pubmed', 'main', 'membrane_id', 'protein_id',
)
# What the flags of the CSV files are written as.
_FLAGS = {'true': True, 'false': False}


def iter_membranes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The membranes of the ``/membranes`` response; raises if it lacks any.

    Fields: ``id`` (the database's number, as a string), ``name``,
    ``topology_in``, ``topology_out`` and ``abbreviation``.
    """
    response = json.loads(read_text(opener))
    objects = response.get('objects')
    if not objects or len(objects) != response.get('total_objects'):
        raise ValueError('The response of /membranes does not hold all membranes')
    for membrane in objects:
        missing = [column for column in MEMBRANE_COLUMNS if column not in membrane]
        if missing:
            raise ValueError(f'A membrane lacks the fields {missing}')
        yield {
            'id': str(membrane['id']),
            'name': membrane['name'],
            'topology_in': membrane['topology_in'],
            'topology_out': membrane['topology_out'],
            'abbreviation': membrane['abbreviation'],
        }


def _membranes(membrane_download: Download, force_refresh: bool) -> dict[str, dict[str, Any]]:
    opener = membrane_download.open(force_refresh=force_refresh)
    return {membrane['id']: membrane for membrane in iter_membranes(opener)}


def _proteins(protein_download: Download, force_refresh: bool) -> dict[str, dict[str, str]]:
    opener = protein_download.open(force_refresh=force_refresh)
    return {
        row['id']: row for row in iter_csv(opener, required=PROTEIN_COLUMNS)
    }


def iter_localizations(
    opener,
    *,
    membrane_download: Download,
    protein_download: Download,
    force_refresh: bool = False,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """One row for each membrane of each protein of ``protein_membranes``.

    Fields: ``accession`` (UniProt), ``entry_name`` (None where it is not a
    UniProt entry name, as for the two chains of GNPTA), ``gene``, ``name`` (the
    protein), ``species`` (the name the database gives), ``membrane`` (the
    membrane's name), ``pmid`` (None where the row has none), ``main`` (the
    ``main`` flag of the row), ``side`` (the topology side of the protein, only
    on the row of the membrane that the protein table gives, None otherwise
    and where the table has no topology) and ``topology_flag`` (the
    ``topology_show_in`` of the protein as written). A row whose protein or
    membrane is unknown, or whose names do not agree with the other tables,
    raises.
    """
    membranes = _membranes(membrane_download, force_refresh)
    proteins = _proteins(protein_download, force_refresh)
    for row in iter_csv(opener, required=LOCALIZATION_COLUMNS):
        protein = proteins.get(row['protein_id'])
        membrane = membranes.get(row['membrane_id'])
        if protein is None or membrane is None:
            raise ValueError(f'Unknown protein or membrane in the row {row["id"]}')
        if row['protein_uniprot_id_cache'] != protein['uniprot_id']:
            raise ValueError(f'The protein of the row {row["id"]} differs between the tables')
        if row['membrane_name_cache'] != membrane['name']:
            raise ValueError(f'The membrane of the row {row["id"]} differs between the tables')
        if row['main'] not in _FLAGS:
            raise ValueError(f'Unknown value of main: {row["main"]!r}')
        flag = protein['topology_show_in']
        if flag not in _FLAGS and flag != '':
            raise ValueError(f'Unknown value of topology_show_in: {flag!r}')
        side = None
        if flag and protein['membrane_id'] == row['membrane_id']:
            side = membrane['topology_in'] if _FLAGS[flag] else membrane['topology_out']
        entry_name = protein['uniprot_id']
        yield {
            'accession': protein['uniprotcode'],
            'entry_name': entry_name if entry_name.count('_') == 1 else None,
            'gene': protein['genename'] or None,
            'name': protein['name'],
            'species': protein['species_name_cache'],
            'membrane': membrane['name'],
            'pmid': row['pubmed'] or None,
            'main': _FLAGS[row['main']],
            'side': side,
            'topology_flag': flag or None,
        }
