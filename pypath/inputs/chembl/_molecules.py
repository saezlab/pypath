from collections.abc import Generator

from ._records import ChemblMolecule, ChemblMolProps, ChemblMolStruct
from . import _raw

__all__ = [
    'molecule',
]


def molecule(max_pages: int | None = None) -> Generator[ChemblMolecule]:
    """
    Retrieves molecule information from ChEMBL
    """

    molecules = _raw.json_pages(data_type="molecule", max_pages=max_pages)

    for molecule in molecules:

        # Get the molecule properties
        molecule_properties = _molecule_props(molecule["molecule_properties"])

        # Get the molecule structures
        structures = _molecule_strucs(molecule["molecule_structures"])

        yield ChemblMolecule(
            molecule_chembl_id = molecule['molecule_chembl_id'],
            preferred_name = molecule['pref_name'],
            molecule_type = molecule['molecule_type'],
            structure_type = molecule['structure_type'],
            chirality = molecule['chirality'],
            biotherapeutic = molecule['biotherapeutic'],
            inorganic_flag = molecule['inorganic_flag'],
            natural_flag = molecule['natural_product'],
            polymer_flag = molecule['polymer_flag'],
            helm_notation = molecule['helm_notation'],
            molecule_properties = molecule_properties,
            structure = structures,
        )


def _molecule_props(properties: dict) -> ChemblMolProps | None:
    """
    Retrieves molecule properties from ChEMBL.

    Args:
        properties (dict): The dictionary of molecule properties.

    Returns:
        ChemblMolProps: The named tuple of the molecule properties.
    """
    return ChemblMolProps(
        mol_formula=properties.get("full_molformula"),
        full_mwt=properties.get("full_mwt"),
        monoisotopic_mwt=properties.get("mw_monoisotopic"),
        molecular_species=properties.get("molecular_species"),
        logp=properties.get("cx_logp"),
        logd=properties.get("cx_logd"),
        alogp=properties.get("alogp"),
    ) if properties else None


def _molecule_strucs(structure: dict) -> ChemblMolStruct | None:
    """
    Retrieves molecule structure from ChEMBL.

    Args:
        structure (dict): The dictionary of molecule structure.

    Returns:
        ChemblMolStruct: The named tuple of the molecule structure.
    """
    return ChemblMolStruct(
        canonical_smiles=structure.get("canonical_smiles"),
        inchi=structure.get("standard_inchi"),
        inchi_key=structure.get("standard_inchi_key"),
    ) if structure else None
