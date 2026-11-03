from collections.abc import Generator

from ._records import ChemblTarget, ChemblComponent
from . import _raw

__all__ = [
    'target',
]


def target(max_pages: int | None = None) -> Generator[ChemblTarget]:
    """
    Retrieves target data from ChEMBL.

    This generator function retrieves the target data from ChEMBL and
    yields the data as named tuples of the type `ChemblTarget`.

    The function uses the `chembl_general` function to retrieve the
    data from ChEMBL.

    Args:
        max_pages (int): The maximum number of pages to retrieve.
    Yields:
        ChemblTarget: The named tuple of the retrieved data.
    """

    targets= _raw.json_pages(data_type="target", max_pages=max_pages)

    # loop through the pages and yield the target data
    for target in targets:

        # components
        components = tuple(target_components(target))
        yield ChemblTarget(
            chembl_id = target['target_chembl_id'],
            target_type = target['target_type'],
            preferred_name = target["pref_name"],
            ncbi_taxa_id = target["tax_id"],
            organism = target["organism"],
            components = components,
            num_components = len(components),
        )


def target_components(target: dict) -> Generator[ChemblComponent]:
    """
    Retrieves target component data from ChEMBL.
    """
    comp_count = 0
    for component in target['target_components']:
        comp_count += 1
        yield ChemblComponent(
            uniprot_accession = component['accession'],
            component_type = component['component_type'],
            component_description = component['component_description'],
            component_id = component['component_id'],
            component_relationship = component['relationship'],
            component_number = comp_count,
        )
