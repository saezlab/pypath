from __future__ import annotations

from collections.abc import Generator
import itertools
import re

import pypath.resources.urls as urls
from ._records import StitchLinks
from ._raw import tables


def _prepend(first, rest):
    return itertools.chain([first], rest)

__all__ = [
    'links',
]

RECHEMID = re.compile(r'CID([ms]?)0*(\d+)')          # e.g. CIDm01234, CIDs01234
REPROTID = re.compile(r'(\d+)\.?(ENS[A-Z]*P\d+)')  # e.g. 9606.ENSP00000123, 10090.ENSMUSP00000012

def links(max_lines: int | None = None,
          ncbi_tax_id: int = 9606) -> Generator[StitchLinks]:
    """
    Downloads the 'links' file from STITCH and formats the data.

    Tries the rescued OmniPath mirror first; falls back to the original
    STITCH server if the mirror returns no data (e.g. for non-human
    organisms that are not hosted on the mirror).

    Yields
        StitchLinks
            A named tuple containing information about each ligand-receptor
            interaction.
    """

    url = urls.urls['stitch']['links_rescued'] % ncbi_tax_id
    links = tables(url, max_lines)

    # Peek at the first record to detect a failed/empty download.
    # If the rescued mirror returned nothing, fall back to the original URL.
    try:
        first = next(links)
    except StopIteration:
        url = urls.urls['stitch']['links'] % ncbi_tax_id
        links = tables(url, max_lines)
        try:
            first = next(links)
        except StopIteration:
            return

    links = _prepend(first, links)


    for link in links:  # noqa: F402

        # extract ids
        stereo, cid_id = RECHEMID.match(link['chemical']).groups()
        tax, ens_id = REPROTID.match(link['protein']).groups()

        # yield the information as a named tuple
        yield StitchLinks(
            chemical_id = cid_id,
            protein_id = ens_id,
            experimental = link["experimental"],
            prediction = link["prediction"],
            database = link["prediction"],
            textmining = link["textmining"],
            combined_score = link["textmining"],
            ncbi_tax_id = int(tax),
            stereospecific = stereo == 's', # True if the ligand is stereospecific
        )
