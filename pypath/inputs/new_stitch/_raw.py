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


def interactions(max_lines: int | None = None,
                 ncbi_tax_id: int = 9606):
    """
    Merges the actions and links tables into a single interactions table

    Args:
        max_lines: Maximum number of lines to read from the file.
        ncbi_tax_id: NCBI taxonomy ID of the organism to read actions for.

    Yields:
        StitchInteractions: A named tuple containing the interaction information
    """

    the_actions = actions(max_lines, ncbi_tax_id)

    # Create a dictionary to look up the link information for each action
    link_lookup = dict(((link.chemical_id, link.protein_id), link) for link in links(max_lines, ncbi_tax_id))

    # Iterate over the actions and look up the corresponding link information
    for action in the_actions:

        link = link_lookup.get((action.chemical_id, action.protein_id))
        
        if link:
            yield StitchInteractions(
                chemical_id = action.chemical_id,
                chemical_acting = action.chemical_acting,
                is_stereospecific = action.is_stereospecific,
                protein_id = action.protein_id,
                protein_acting = action.protein_acting,
                mode = action.mode,
                action = action.action,
                ncbi_taxa = action.ncbi_taxa,
                experimental = link.experimental,
                prediction = link.prediction,
                database = link.database,
                textmining = link.textmining,
                combined_score = link.combined_score,
            )
