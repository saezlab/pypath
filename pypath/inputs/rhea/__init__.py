import csv
from typing import Literal
from collections.abc import Generator
from pypath.resources import urls
from pypath.share import curl

RHEA_DATASETS = Literal[
    'rhea-reaction-smiles',
    'rhea-directions',
    'rhea-relationships',
    'rhea-chebi-smiles',
    'rhea2xrefs',
    'rhea2uniprot_sprot',
    #'rhea2uniprot_trembl.tsv.gz' # TODO: deal with gz file ending
]

def _rhea_txt(dataset: RHEA_DATASETS, colnames: list[str]) -> Generator[tuple]:

    c = curl.Curl(dataset, large=True)
    yield from csv.DictReader(c.fileobj, delimiter='\t', fieldnames=colnames)


