from collections.abc import Generator
import re

import pypath.resources.urls as urls
from ._records import StitchAction, Entity
from ._raw import tables

__all__ = [
    'actions',
]

REID = re.compile(r'((?:\d+)?)\.?(CID|ENSP)([ms]?)(\d+)')


def actions(max_lines: int | None = None,
            ncbi_tax_id: int = 9606
            ) -> Generator[StitchAction]:
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

    parse_activation = lambda x: None if x == '' else x == 'activation'

    url = urls.urls['stitch']['actions'] % ncbi_tax_id

    actions = tables(url, max_lines)

    out = []

    for action in actions:

        a, b = parse_entities(action)

        out.append(
            StitchAction(
                source = a,
                target = b,
                directed = action["a_is_acting"].lower() == 't',
                mode = action["mode"],
                activation = parse_activation(action["action"]),
                score = int(action["score"]),
            )
        )

        if len(out) == 2:

            if set(out[0][:2]) == set(out[1][:2]):

                if not any(a.directed for a in out):

                    yield out[0]

                else:

                    for a in out:

                        if a.directed:

                            yield a

                out = []

            else:

                yield out.pop(0)


def parse_entities(action: dict) -> tuple[Entity, Entity]:
    """
    Converts the STITCH chemical and protein identifiers and
    flags into a ParsedIds named tuple. Also uses "a_is_acting"
    flag to determine if the chemical or protein is acting. Only
    'partner a' can act.

    Parameters
        action
            A single action record from STITCH.

    Returns
        pair of Entities
    """
    partners = []

    for side in ('a', 'b'):

        tax, ens, stereo, _id = REID.match(action[f'item_id_{side}']).groups()
        id_prefix = ens if ens[:3] == 'ENS' else ''

        partners.append(
            Entity(
                id = f'{id_prefix}{_id}',
                type = 'small_molecule' if ens == 'CID' else 'protein',
                ncbi_tax_id = int(tax) if tax else None,
                stereospecific = stereo == 's',
            )
        )

    if action['a_is_acting'].lower() == 't':

        partners = reversed(partners)

    return tuple(partners)
