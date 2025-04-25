from collections.abc import Generator

from ._records import ChemblActivity
from . import _raw
def get_activity(max_pages: int | None = None) -> Generator[ChemblActivity]:
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

    yield from (ChemblActivity
        (
            activity_id = activity["activity_id"],
            action_type = activity["action_type"],
            standard_relation = activity["standard_relation"],
            standard_value = activity["standard_value"],
            standard_upper_value = activity['standard_upper_value'],
            standard_units = activity["standard_units"],
            standard_type = activity["standard_type"],
            ligand_efficiency = activity["ligand_efficiency"],
            validity_comment = activity["data_validity_comment"],
            potential_duplicate = activity["potential_duplicate"],
            pchembl_value = activity["pchembl_value"],
            source_id = activity["src_id"],
            molecule_chembl_id = activity["molecule_chembl_id"],
            target_chembl_id = activity["target_chembl_id"],
            target_taxa_id = activity["target_tax_id"],
            assay_id = activity["assay_chembl_id"],
            document_chembl_id = activity["document_chembl_id"],
        )
        for activity in activities
    )
