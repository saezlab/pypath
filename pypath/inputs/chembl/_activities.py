from collections.abc import Generator

from ._records import ChemblActivity, ChemblAction
from . import _raw

__all__ = [
    "activity",
]


def activity(max_pages: int | None = None) -> Generator[ChemblActivity]:
    """
    Retrieves activity data from Chembl.
    This generator function retrieves the assay data from ChEMBL and
    yields the data as named tuples of the type `ChemblActivity.

    The function uses the `chembl_general` function to retrieve the
    data from ChEMBL.

    Args:
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        ChemblActivity: The named tuple of the retrieved data.
    """
    activities = _raw.json_pages(data_type="activity", max_pages=max_pages)

    for activ in activities:
        action_type = action_sorter(activ) 

        yield (ChemblActivity
                (
                activity_id = activ["activity_id"],
                action_type = action_type,
                standard_relation = activ["standard_relation"],
                standard_value = activ["standard_value"],
                standard_upper_value = activ['standard_upper_value'],
                standard_units = activ["standard_units"],
                standard_type = activ["standard_type"],
                ligand_efficiency = activ["ligand_efficiency"],
                validity_comment = activ["data_validity_comment"],
                potential_duplicate = activ["potential_duplicate"],
                pchembl_value = activ["pchembl_value"],
                source_id = activ["src_id"],
                molecule_chembl_id = activ["molecule_chembl_id"],
                target_chembl_id = activ["target_chembl_id"],
                target_taxa_id = activ["target_tax_id"],
                assay_id = activ["assay_chembl_id"],
                document_chembl_id = activ["document_chembl_id"],
                )
            )

def action_sorter(activ):
    """
    This function takes an activity data record from ChEMBL and returns an
    instance of the `ChemblAction` namedtuple with the action type
    information.

    args
        activ : dict
            An activity data record from ChEMBL.

    Returns
        ChemblAction: An instance of the `ChemblAction` namedtuple with the action type
        information.
    """
    action_type = activ["action_type"]

    if action_type:
        return ChemblAction(
            action_type=action_type["action_type"],
            description=action_type["description"],
            parent_type=action_type["parent_type"],
            is_stimulation=action_type["action_type"].lower() == "positive modulator",
            is_inhibition=action_type["action_type"].lower() == "negative modulator",
        )
