from __future__ import annotations

from collections.abc import Generator
import csv

import pypath.resources.urls as urls
from pypath.share.downloads import download_and_open
from pypath_common import _misc

__all__ = [
    'table',
    'mapping',
]

def table(
        dataset: str = 'All',
        id_mapping: bool = False,
        max_lines: int | None = None,
    ) -> Generator[dict]:
    """
    Reads a tab-separated values file from the given URL and returns its
    contents in a generator of dictionaries, where each dictionary represents a
    table row.

    Args:
        dataset:
            BindingDB dataset: All, Articles, ChEMBL, PubChem, Patents, etc.
        id_mapping:
            If True, reads a mapping file instead of a table file.
        max_lines:
            Maximum number of lines to read from the file. If None, the entire
            file is read.

    Yields:
        dict
            A dictionary representing a table row.
    """

    url_key = 'mapping' if id_mapping else 'url'

    url = urls.urls['bindingdb'][url_key] % dataset

    # Download and open file using download manager
    # For zip files, opener.result is a dict like curl.Curl
    opener = download_and_open(
        url,
        filename=f'{dataset}_{"mapping" if id_mapping else "table"}.{"tsv" if id_mapping else "zip"}',
        subfolder='bindingdb',
        large=True,
        encoding='utf-8',
    )

    line_count = 0

    # mapping file is plain tsv, table file is zipped
    # For zip files, result is dict; for plain files, it's a file handle
    file_iterator = opener.result if id_mapping else _misc.first(opener.result.values())

    for line in csv.DictReader(file_iterator, delimiter = '\t'):
        yield line
        line_count += 1
        if max_lines is not None and line_count > max_lines:
            break


def mapping(dataset: str = 'UniProt') -> dict[str, list[str]]:
    """
    Reads a tab-separated values file from the given URL and returns its
    contents in a dictionary, where each key is a BindingDB protein name and each
    value is a list of UniProt IDs.

    Args:
        dataset:
            BindingDB mapping dataset: Uniprot, Pubmed, Drugbank, CID etc.

    Returns:
        dict
            A dictionary mapping BindingDB protein names to UniProt IDs.
        
    """

    return {
        i['BindingDB Name']: i[dataset].split()
        for i in table(dataset = dataset, id_mapping = True)
    }
