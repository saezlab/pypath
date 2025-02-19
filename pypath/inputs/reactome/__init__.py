import csv
from typing import Literal
from collections.abc import Generator
from pypath.resources import urls
from pypath.share import curl

def _reactome_txt(url: str, colnames: list[str]) -> Generator[tuple]:

    c = curl.Curl(url, large=True)
    yield from csv.DictReader(c.fileobj, delimiter='\t', fieldnames=colnames)


def pathway_hierarchy() -> Generator[tuple]:

    url = urls.urls['reactome']['pathway_relations']

    yield from _reactome_txt(url = url, colnames = ['parent', 'child']) 


def reactome_raw(
        id_type: str,
        physical_entities: bool = False,
        event_type: Literal['all', 'lowest', 'reaction'] = 'all',
    ) -> Generator[tuple]:

    IDTYPES = {
        'uniprot': 'UniProt',
        'ensembl': 'Ensembl',
        'entrez': 'NCBI',
        'mirbase': 'miRBase',
        'gtop': 'GtoP',
        'chebi': 'ChEBI',
    }
    EVENT_TYPES = {
        'lowest': 'Pathway',
        'all': 'All_Levels',
        'reaction': 'Reactions',
        'reactions': 'Reactions',
        'pathway': 'Pathway',
    }

    fname = [f'{IDTYPES.get(id_type, id_type)}2Reactome']

    if physical_entities:

        fname.append('PE')

    if physical_entities or event_type != 'pathway':

        fname.append(EVENT_TYPES.get(event_type, event_type))

    fname = '_'.join(fname)

    if not physical_entities and event_type == 'reaction':

        fname = fname.replace('_', '')

    colnames = ['id']

    if physical_entities:

        colnames.extend(['pe_id', 'pe_name'])

    colnames.extend(
        [
            'pathway_id',
            'url',
            'pathway_name',
            'evidence_code',
            'species',
        ],
    )

    url = urls.urls['reactome']['txt'] % fname

    yield from _reactome_txt(
        url = url,
        colnames = colnames,
    )


