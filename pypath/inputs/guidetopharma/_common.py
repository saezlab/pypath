from typing import Literal
from collections.abc import Generator
import collections

import csv

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.utils.taxonomy as taxonomy
import pypath.resources.urls as urls

TABLES = Literal[
    "ligands",
    "interactions",
    "targets_and_families",
    "ligand_id_mapping",
    "ligand_physchem_properties",
    "endogenous_ligand_pairings_all",
    "endogenous_ligand_detailed",
    "approved_drug_primary_target_interactions",
    "approved_drug_detailed_interactions",
    "peptides",
    "GtP_to_HGNC_mapping",
    "GtP_to_UniProt_mapping",
]

POSITIVE_REGULATION = {
    "agonist",
    "activator",
    "potentiation",
    "partial agonist",
    "inverse antagonist",
    "full agonist",
    "activation",
    "irreversible agonist",
    "positive",
    "biased agonist",
    "slows inactivation",
}
NEGATIVE_REGULATION = {
    "inhibitor",
    "antagonist",
    "competitive",
    "feedback inhibition",
    "inhibition",
    "irreversible inhibition",
    "inverse agonist",
    "negative",
    "weak inhibition",
    "reversible inhibition",
    "voltage-dependent inhibition",
    "pore blocker",
}


G2PInteraction = collections.namedtuple(
    "G2PInteraction",
    [
        "is_stimulation",
        "is_inhibition",
        "ligand",
        "endogenous",
    ],
)

G2PLigand = collections.namedtuple(
    "G2PLigand",
    [
        "name",
        "pubchem",
        "chembl",
        "iupac",
        "smiles",
        "inchi",
    ],
)

G2PProtein = collections.namedtuple(
    "G2PProtein",
    [
        "uniprot",
    ],
)


def guide2pharma_table(name: TABLES) -> Generator[dict]:
    """
    Downloads the table from Guide2Pharma.

    Args:
        name
            The name of the table to download.

    Returns:
        Generator yielding the table as named tuples.
    """

    url = urls.urls["gtp"]["url"] % name

    c = curl.Curl(url, silent=False, large=True, encoding="utf-8")

    g2p_version = next(c.result).strip()

    return csv.DictReader(c.result)


def guide2pharma_interactions(
    organism: str | int | None = "human",
    endogenous: bool | None = None,
) -> Generator[tuple]:
    """
    Args:
        organism (str | int | None):
            Name of the organism, e.g. `human`. If None, all organisms will be
            included.]

        endogenous (bool | None):
            Whether to include only endogenous ligands interactions. If None,
            all ligands will be included.
    """

    organism_ = None
    ncbi_tax_id = None

    if isinstance(organism, str):
        ncbi_tax_id = get_taxid(organism)

        try:
            organism_ = taxonomy.ensure_common_name(ncbi_tax_id)
            organism_ = organism_.capitalize() if organism_ else None

        except KeyError:
            pass  # no organism specified

    ligands = guide2pharma_ligands()

    for row in guide2pharma_table("interactions"):

        _endogenous = row["Endogenous"].lower() == "true"

        if endogenous is not None and endogenous != _endogenous:

            continue

        yield G2PInteraction(
            is_stimulation = row["Action"].lower() in POSITIVE_REGULATION,
            is_inhibition = row["Action"].lower() in NEGATIVE_REGULATION,
            ligand = ligands.get(row["Ligand ID"]),
            endogenous = _endogenous,
        )

def get_taxid(organism_input: str) -> int:
    """
    Retrieves the NCBI Taxonomy ID for a given organism input.

    Args:
        organism_input (str): The name or identifier of the organism.

    Returns:
        str: The NCBI Taxonomy ID corresponding to the organism input.
    """
    if organism_input in {"", "None", None}:
        return _const.NOT_ORGANISM_SPECIFIC

    return taxonomy.ensure_ncbi_tax_id(organism_input)

def guide2pharma_ligands() -> dict[str, G2PLigand | G2PProtein]:
    """
    Downloads ligands from Guide2Pharma.

    Returns:
        Named tuples containing name, PubChem ID, ChEMBL ID, IUPAC name, SMILES,
        and InChI of each ligand.
    """
    return {
        row["Ligand ID"]: _entity_record(row)
        for row in guide2pharma_table("ligands")
    }

def _entity_record(row: dict) -> G2PLigand | G2PProtein:
    """
    Creates a G2PProtein or G2PLigand object from a row of the Guide2Pharma
    table.

    Args:
        row: A dictionary representing a row of the Guide2Pharma table.

    Returns:
        A G2PProtein or G2PLigand object.
    """
    if row['UniProt ID']:
        record = G2PProtein(
            uniprot = row["UniProt ID"],
        )
    else:
        record = G2PLigand(
            name=row["Name"],
            pubchem=row["PubChem CID"],
            chembl=row["ChEMBL ID"],
            iupac=row["IUPAC name"],
            smiles=row["SMILES"],
            inchi=row["InChI"],
        )

    return record
