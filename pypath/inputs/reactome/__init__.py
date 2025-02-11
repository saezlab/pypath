# ChEBI2Reactome_PE_All_Levels: https://reactome.org/download/current/ChEBI2Reactome_PE_All_Levels.txt
# pathway hierarchy: https://reactome.org/download/current/ReactomePathwaysRelation.txt
import csv
from typing import Literal
from collections.abc import Generator
from pypath.resources import urls
from pypath.share import curl

def _reactome_txt(fname: str, colnames: list[str]) -> Generator[tuple]:

    url = urls.urls['reactome']['txt'] % fname

    c = curl.Curl(url, large=True)
    yield from csv.DictReader(c.fileobj, delimiter='\t', fieldnames=colnames)


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
        'lowest': '',
        'all': '_All_Levels',
        'reaction': 'Reactions',
    }

    colnames = ['id', 'pathway_id', 'url', 'pathway_name', 'evidence_code', 'species']
    id_type = IDTYPES.get(id_type, id_type)
    event_type = EVENT_TYPES.get(event_type, event_type)
    physical_entities = '_PE_' if physical_entities else ''

    fname = f'{id_type}2Reactome{physical_entities}{event_type}'

    yield from _reactome_txt(
        fname = fname,
        colnames = colnames,
    )


