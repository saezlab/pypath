# ChEBI2Reactome_PE_All_Levels: https://reactome.org/download/current/ChEBI2Reactome_PE_All_Levels.txt
# pathway hierarchy: https://reactome.org/download/current/ReactomePathwaysRelation.txt
from pypath.resources import urls

def reactome_raw(fname: str) -> Generator[tuple]:

    url = urls.urls['reactomr']['txt'] % fname
