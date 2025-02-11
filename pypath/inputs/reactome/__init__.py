# ChEBI2Reactome_PE_All_Levels: https://reactome.org/download/current/ChEBI2Reactome_PE_All_Levels.txt
# pathway hierarchy: https://reactome.org/download/current/ReactomePathwaysRelation.txt
import csv
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
        event_type: Literal['all', 'lowest', 'reaction'],
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
        'all': 'All_levels',
        'reaction': 'Reactions',
    }

    colnames = ['id', 'pathway_id', 'url', 'pathway_name', 'evidence_code', 'species']

    fname = f'{id_type}2Reactome{pe}{event_type}'

    yield from _reactome_txt(
        fname = 'ChEBI2Reactome_PE_All_Levels',
        colnames = colnames,
    )


