from collections.abc import Generator

from ._records import StitchInteractions, StitchAction
from ._links import links
from ._actions import actions

__all__ = [
    'interactions',
]

def interactions(max_lines: int | None = None, 
                 ncbi_tax_id: int = 9606
                 ) -> Generator[StitchInteractions]:
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
    link_lookup = create_lookup(max_lines, ncbi_tax_id)

    # Iterate over the actions and look up the corresponding link information
    for action in the_actions:

        action_identifier = get_action_identifier(action)
        
        if action_identifier in link_lookup:
            link = link_lookup[action_identifier]
            yield StitchInteractions(
                source = action.source,
                target = action.target,
                directed = action.directed,
                mode = action.mode,
                activation = action.activation,
                inhibition = action.inhibition,
                experimental = link.experimental,
                prediction = link.prediction,
                database = link.database,
                textmining = link.textmining,
                combined_score = link.combined_score,
                final_score = action.score,
            )


def get_action_identifier(action: StitchAction) -> tuple:
    """
    Returns a tuple that can be used as a key to look up the link information for an action

    Args:
        action: A StitchAction object

    Returns:
        tuple: A tuple that can be used as a key to look up the link information for an action
    """

    # get the correct order to search the links lookup
    if action.source.type == 'small_molecule':
        return (action.source.id, action.target.id, action.target.ncbi_tax_id, action.source.stereospecific)
    else:
        return (action.target.id, action.source.id, action.source.ncbi_tax_id, action.target.stereospecific)
    
def create_lookup(max_lines: int | None = None,
                 ncbi_tax_id: int = 9606) -> dict:
    """
    Create a dictionary to look up the link information for each action

    Parameter:
        max_lines: Maximum number of lines to read from the file.
        ncbi_tax_id: NCBI taxonomy ID of the organism to read actions for.

    Returns:
        dict: A dictionary to look up the link information for each action.
    """

    all_links = links(max_lines, ncbi_tax_id)
    
    link_lookup = {(link.chemical_id, link.protein_id, link.ncbi_tax_id, link.stereospecific): link
                   for link in all_links}

    return link_lookup