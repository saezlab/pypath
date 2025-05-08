from collections.abc import Generator
import re

import pypath.resources.urls as urls
from ._records import StitchActions, Entity
from ._raw import tables

__all__ = [
    'actions',
]

REID = re.compile(r'((?:\d+)?)\.?(CID|ENSP)([ms]?)(\d+)')


def actions(max_lines: int | None = None,
            ncbi_tax_id: int = 9606
            ) -> Generator[StitchActions]:
    """
    Downloads the 'actions' file from STITCH and formats the data.

    Parameters
        max_lines : int | None
            Maximum number of lines to read from the file.
        ncbi_tax_id : int
            NCBI taxonomy ID of the organism to read actions for.

    Yields
        StitchActions
            A named tuple containing information about each action.
    """

    url = urls.urls['stitch']['actions'] % ncbi_tax_id

    actions = tables(url, max_lines)

    for action in actions:

        parsed_ids = convert_ids(action)

        yield StitchActions(
            chemical_id = parsed_ids.chemical_id,
            chemical_acting = parsed_ids.chemical_acting,
            is_stereospecific = parsed_ids.stereospecific,
            protein_id = parsed_ids.protein_id,
            protein_acting = parsed_ids.protein_acting,
            mode = action["mode"],
            action = action["action"],
            score = action["score"],
            ncbi_taxa = parsed_ids.ncbi_taxa
        )

def convert_ids(action: dict) -> Entity:
    """
    Converts the STITCH chemical and protein identifiers and
    flags into a ParsedIds named tuple. Also uses "a_is_acting"
    flag to determine if the chemical or protein is acting. Only
    'partner a' can act.

    Parameters
        action : dict
            A single action record from STITCH.
        sep : re.Pattern
            A regex pattern used to split the identifiers.

    Returns
        ParsedIds
    """
    partners = []

    for side in ('a', 'b'):

        tax, ens, stereo, _id = REID.match(action[f'item_id_{side}']).groups()
        id_prefix = ens if ens[:3] == 'ENS' else ''

        partners.append(
            Entity(
                id = f'{id_prefix}{_id},
                type = 'small_molecule' if ens == 'CID' else 'protein',
                ncbi_tax_id = tax,
                stereo = stereo == 's',
            )
        )

    if action['a_is_acting'].lower() == 't':

        partners = reversed(partners)

    return tuple(partners)


def cid_converter(cids: list) -> tuple:
    """
    Converts a STITCH CID identifier into a PubChem CID and a boolean indicating
    whether the compound is stereospecific.
    """
    cid = cids[-1]

    if cids[1] == "s":
        stereospecific = True
    else:
        stereospecific = False

    return cid, stereospecific

def ensembl_converter(ensembl_ids: list) -> tuple:
    """
    Converts the STITCH Ensembl identifier into an ensembl ID and an NCBI taxon ID
    """

    ensembl_id = ensembl_ids[-1]
    ncbi_taxa = ensembl_ids[0]

    return ensembl_id, ncbi_taxa
