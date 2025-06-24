from collections.abc import Generator
import csv

import pypath.resources.urls as urls
import pypath.share.curl as curl
from pypath_common import _misc

__all__ = [
    'table'
]

def table(
        dataset: str = 'All',
        max_lines: int | None = None,
        mapping: bool = False
    ) -> Generator[dict]:
    """
    Reads a tab-separated values file from the given URL and returns its
    contents in a generator of dictionaries, where each dictionary represents a
    table row.

    Args:
        dataset:
            BindingDB dataset: All, Articles, ChEMBL, PubChem, Patents, etc.
        max_lines:
            Maximum number of lines to read from the file. If None, the entire
            file is read.

    Yields:
        dict
            A dictionary representing a table row.
    """

    url_key = 'mapping' if mapping else 'url'

    url = urls.urls['bindingdb'][url_key] % dataset
    c = curl.Curl(url, silent = False, large = True, slow = True)

    line_count = 0

    file_iterator = c.result if mapping else _misc.first(c.result.values())
    for line in csv.DictReader(file_iterator, delimiter = '\t'):

        yield line

        line_count += 1
        if max_lines is not None and line_count > max_lines:
            break
