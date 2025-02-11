# ChEBI2Reactome_PE_All_Levels: https://reactome.org/download/current/ChEBI2Reactome_PE_All_Levels.txt
# pathway hierarchy: https://reactome.org/download/current/ReactomePathwaysRelation.txt
import csv
from collections.abc import Generator
from pypath.resources import urls
from pypath.share import curl

def reactome_raw(fname: str) -> Generator[tuple]:

    url = urls.urls['reactome']['txt'] % fname

    c = curl.Curl(url, large=True)
    yield from csv.DictReader(c.fileobj)
