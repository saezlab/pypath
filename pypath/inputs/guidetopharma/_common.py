from collections.abc import Generator

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls


def guide2pharma_table(name: str) -> Generator[dict]:
    """
    Downloads the table from Guide2Pharma.

    Args:
        name:
            The name of the table to download.

    Returns:
        Generator yielding the table as named tuples.
    """

    url = urls.urls['gtp']['url'] % name

    c = curl.Curl(url, silent = False, large = True, encoding = 'utf-8')

    return c

    # data = csv.DictReader(c.result)

