from collections.abc import Generator
import re
import csv

import pypath.resources.urls as urls
import pypath.share.curl as curl

__all__ = [
    'tables'
]

def tables(url: str,
           max_lines: int | None = None
           ) -> Generator[dict]:
    """
    Reads a tab-separated values file from the given URL and returns its
    contents in a generator of dictionaries, where each dictionary represents a
    table row.

    Parameters
        url : str
            URL of the file to read.
        max_lines : int | None
            Maximum number of lines to read from the file. If None, the entire
            file is read.

    Yields
        dict
            A dictionary representing a table row.
    """

    c = curl.Curl(url, silent = False, large = True, slow = True)

    line_count = 0
    for line in csv.DictReader(c.result, delimiter = '\t'):
        
        yield line

        line_count += 1
        if max_lines is not None and line_count > max_lines:
            break
