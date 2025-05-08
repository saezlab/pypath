from collections.abc import Generator
import re

import pypath.resources.urls as urls
from ._records import StitchLinks
from ._raw import tables

__all__ = [
    'links',
]

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

    sep = re.compile(r'[sm\.]')

    for link in links:
        # extract ids
        chemical_id = sep.split(link['chemical'])[-1]
        protein_id = sep.split(link['protein'])[-1]

        # yield the information as a named tuple
        yield StitchLinks(
            chemical_id = chemical_id,
            protein_id = protein_id,
            experimental = link["experimental"],
            prediction = link["prediction"],
            database = link["prediction"],
            textmining = link["textmining"],
            combined_score = link["textmining"],
        )