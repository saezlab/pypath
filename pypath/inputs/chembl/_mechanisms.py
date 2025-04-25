from collections.abc import Generator

from ._records import ChemblMechanism
from . import _raw

__all__ = [
    'mechanism',
]

def mechanism(max_pages: int | None = None) -> Generator[ChemblMechanism]:
    """
    Retrieves mechanism data from ChEMBL.
    """

    mechanisms = _raw.json_pages(data_type="mechanism", max_pages=max_pages)

    for mech in mechanisms:

        # Get the mechanism references
        references = [ref["ref_url"] for ref in mech["mechanism_refs"]]

        yield ChemblMechanism(
            action_type = mech["action_type"],
            molecule_chembl_id = mech["molecule_chembl_id"],
            target_chembl_id = mech["target_chembl_id"],
            mechanism_id = mech["mec_id"],
            drug_rec_id = mech["record_id"],
            direct_interaction = mech["direct_interaction"],
            variant_sequence = mech["variant_sequence"],
            molecular_mechanism = mech["molecular_mechanism"],
            mechanism_refs = references,
        )
