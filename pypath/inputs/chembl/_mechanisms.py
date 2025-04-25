from collections.abc import Generator

from ._records import ChemblMechanism
from . import _raw
def get_mechanisms(max_pages: int | None = None) -> Generator[ChemblMechanism]:
    """
    Retrieves mechanism data from ChEMBL.
    """

    mechanisms = _raw.json_pages(data_type="mechanism", max_pages=max_pages)

    for mechanism in mechanisms:

        # Get the mechanism references
        references = [ref["ref_url"] for ref in mechanism["mechanism_refs"]]

        yield ChemblMechanism(
            action_type = mechanism["action_type"],
            molecule_chembl_id = mechanism["molecule_chembl_id"],
            target_chembl_id = mechanism["target_chembl_id"],
            mechanism_id = mechanism["mec_id"],
            drug_rec_id = mechanism["record_id"],
            direct_interaction = mechanism["direct_interaction"],
            variant_sequence = mechanism["variant_sequence"],
            molecular_mechanism = mechanism["molecular_mechanism"],
            mechanism_refs = references,
        )
