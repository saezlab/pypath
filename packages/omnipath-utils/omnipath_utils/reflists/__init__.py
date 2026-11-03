"""Reference lists — complete sets of identifiers for an organism.

Example::

    from omnipath_utils.reflists import all_swissprots, is_swissprot

    swissprots = all_swissprots()  # all human SwissProt IDs
    is_swissprot('P04637')        # True (TP53 is reviewed)
"""

from omnipath_utils.reflists._manager import ReferenceListManager


def get_reflist(list_name: str, ncbi_tax_id: int = 9606) -> set[str]:
    """Get a named reference list for an organism."""
    return ReferenceListManager.get().get_reflist(list_name, ncbi_tax_id)


def all_swissprots(ncbi_tax_id: int = 9606) -> set[str]:
    """All reviewed UniProt ACs for an organism."""
    return get_reflist('swissprot', ncbi_tax_id)


def all_trembls(ncbi_tax_id: int = 9606) -> set[str]:
    """All unreviewed UniProt ACs for an organism."""
    return get_reflist('trembl', ncbi_tax_id)


def all_swissprots_global() -> set[str]:
    """The complete reviewed (SwissProt) AC set across all organisms."""
    return ReferenceListManager.get().load_swissprot_global()


def all_uniprots(ncbi_tax_id: int = 9606) -> set[str]:
    """All UniProt ACs (SwissProt + TrEMBL) for an organism."""
    return get_reflist('uniprot', ncbi_tax_id)


def is_swissprot(uniprot_ac: str, ncbi_tax_id: int = 9606) -> bool:
    """Check if a UniProt AC is reviewed (SwissProt)."""
    return ReferenceListManager.get().is_swissprot(uniprot_ac, ncbi_tax_id)


def is_trembl(uniprot_ac: str, ncbi_tax_id: int = 9606) -> bool:
    """Check if a UniProt AC is unreviewed (TrEMBL)."""
    return ReferenceListManager.get().is_trembl(uniprot_ac, ncbi_tax_id)
