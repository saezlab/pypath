from collections.abc import Generator
import re

import pypath.resources.urls as urls
from ._records import StitchLinks
from ._raw import tables

__all__ = [
    'links',
]

RECHEMID = re.compile(r'CID([ms]?)0*(\d+)') # seperates the stitch pubchemID
REPROTID = re.compile(r'(\d+)\.?(ENSP0*\d+)')  # seperates the stitch ensembl protein ID

def links(max_lines: int | None = None,
          ncbi_tax_id: int = 9606) -> Generator[StitchLinks]:
    """
    Downloads the 'links' file from STITCH and formats the data.

    Yields
        StitchLinks
            A named tuple containing information about each ligand-receptor
            interaction.
    """

    url = urls.urls['stitch']['links'] % ncbi_tax_id

    links = tables(url, max_lines)


    for link in links:

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