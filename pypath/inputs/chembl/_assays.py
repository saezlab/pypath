from collections.abc import Generator

from ._records import ChemblAssay, ChemblParam
from . import _raw

__all__ = [
    "assay",
]


def assay(max_pages: int | None = None) -> Generator[ChemblAssay]:
    """
    Retrieves assay data from ChEMBL.

    This generator function retrieves the assay data from ChEMBL and
    yields the data as named tuples of the type `ChemblAssay`.

    The function uses the `chembl_general` function to retrieve the
    data from ChEMBL.

    Args:
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        ChemblAssay: The named tuple of the retrieved data.
    """

    assays = _raw.json_pages(data_type="assay", max_pages=max_pages)

    for assay in assays:

        # Get the assay parameters
        parameters = tuple(_param_assay(assay['assay_parameters']))

        # Create the ChemblAssay named tuple
        yield ChemblAssay(
            assay_chembl_id = assay['assay_chembl_id'],
            assay_type = assay['assay_type'],
            assay_type_description = assay['assay_type_description'],
            assay_category = assay['assay_category'],
            target_chembl_id= assay['target_chembl_id'],
            organism = assay['assay_organism'],
            tax_id = assay['assay_tax_id'],
            tissue = assay['assay_tissue'],
            cell_type = assay['assay_cell_type'],
            subcellular_fraction = assay['assay_subcellular_fraction'],
            parameters = parameters,
            source_id = assay['src_id'],
            confidence_score = assay['confidence_score'],
            confidence_description = assay['confidence_description'],
            document_chembl_id = assay['document_chembl_id'],
            description = assay['description'],
        )

def _param_assay(parameters: dict) -> Generator[ChemblParam]:
    """
    Retrieves assay parameter data from ChEMBL.

    Args:
        parameters (dict): The dictionary of parameters.

    Yields:
        ChemblParam: The named tuple of the retrieved parameter data.
    """
    if parameters:
        yield from (ChemblParam
        (
            standard_type = parameter["standard_type"],
            standard_value = parameter["standard_value"],
            standard_units = parameter["standard_units"],
            standard_relation = parameter["standard_relation"]
        )
        for parameter in parameters
    )
