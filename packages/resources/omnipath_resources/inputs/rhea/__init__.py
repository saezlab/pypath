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

RHEA_COLNAMES = {
    'rhea-reaction-smiles': ['ID', 'SMILES'],
    'rhea-chebi-smiles': ['CHEBI', 'SMILES'],
}


def rhea_raw(dataset: RHEA_DATASETS) -> Generator[tuple]:
    """
    Download a table from the Rhea reaction databas.

    Args:
        dataset:
            The dataset to download. See ``RHEA_DATASETS`` for available
            datasets.
    """

    url = urls.urls['rhea']['url'] % dataset
    c = curl.Curl(url, large = True)
    names = RHEA_COLNAMES.get(dataset, None)

    yield from csv.DictReader(c.fileobj, delimiter = '\t', fieldnames = names)


# Helper to fetch the full Rhea TSV export from the web GUI API
def rhea_gui() -> Generator[dict]:
    """
    Download the full Rhea TSV export from the web GUI API.

    Yields:
        One dict per row, using the header row as fieldnames.
    """
    api_url = urls.urls['rhea']['api_url']
    c = curl.Curl(api_url, large=True)
    yield from csv.DictReader(c.fileobj, delimiter='\t')