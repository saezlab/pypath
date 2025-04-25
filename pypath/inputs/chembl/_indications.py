
from collections.abc import Generator

from ._records import ChemblIndication
from . import _raw

__all__ = [
    "indication",
]

def indication(max_pages: int | None = None) -> Generator[ChemblIndication]:
    """
    Retrieves drug indications from ChEMBL.

    This function is a wrapper around the `chembl_general` function.
    It retrieves the indication information from ChEMBL and
    yields the data as named tuples of the type `ChemblIndication`.

    Args:
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        ChemblIndication: The named tuple of the retrieved data.
    """
    indications = _raw.json_pages(
        data_type="drug_indication",
        max_pages=max_pages,
    )

    yield from (ChemblIndication
        (
            chembl_id = indication["molecule_chembl_id"],
            efo_id = indication["efo_id"],
            efo_term = indication["efo_term"],
            mesh_id = indication["mesh_id"],
            mesh_term = indication["mesh_heading"],
            max_phase = indication["max_phase_for_ind"],
            refs = ";".join(ref["ref_id"] for ref in indication["indication_refs"])
        )
        for indication in indications
    )
