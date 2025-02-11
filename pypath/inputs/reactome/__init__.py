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


def reactome_raw():
    colnames = ['id', 'pathway_id', 'url', 'pathway_name', 'evidence_code', 'species']
    yield from _reactome_txt(fname='ChEBI2Reactome_PE_All_Levels', colnames=colnames)


