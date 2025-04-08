# TODO: process cases where Target= None properly (only 257 cases)

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
        "endogenous",
        "affinity_high",
        "affinity_low",
        "affinity_median",
        "affinity_units",
        "primary_target",
        "pubmed_id",
        "ligand",
        "target",
    ],
)

G2PCompound = collections.namedtuple(
    "G2PCompound",
    [
        "name",
        "pubchem_cid",
        "chembl",
        "iupac",
        "smiles",
        "inchi",
        "compound_type",
        "species",
        "record_type",
    ],
)

G2PTargetProtein = collections.namedtuple(
    "G2PTargetProtein",
    [
        "family_name",
        "hgnc_name",
        "hgnc_symbol",
        "hgnc_id",
        "human_ensembl_gene",
        "human_entrez_gene",
        "human_swissprot",
        "human_nucleotide_refseq",
        "human_protein_refseq",
        "target_type",
        "synonyms",
        "record_type",
    ],
)

G2PLigandProtein = collections.namedtuple(
    "G2PLigandProtein",
    [
        "uniprot_id",
        "name",
        "pubchem_cid",
        "inchi",
        "smiles",
        "type",
        "species",
        "record_type",
    ],
)

def guide2pharma_interactions(
    organism: str | int | None = "human",
    endogenous: bool | None = None,
) -> Generator[tuple]:
    """
    Downloads interactions from Guide2Pharma.

    Args:
        organism (str | int | None):
            Name of the organism, e.g. `human`. If None, all organisms will be
            included. The organism name is translated to its NCBI Taxonomy ID
            using `pypath.utils.taxonomy.get_taxid`.

        endogenous (bool | None):
            Whether to include only endogenous ligands interactions. If None,
            all ligands will be included.

    Yields:
        Named tuples containing information about each interaction.
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

    compounds = g2p_compounds()
    protein_targets = g2p_protein_targets()

    for row in guide2pharma_table("interactions"):

        _endogenous = row["Endogenous"].lower() == "true"

        if endogenous is not None and endogenous != _endogenous:
            continue

        target = g2p_get_target(row, compounds, protein_targets)

        yield G2PInteraction(
            is_stimulation = row["Action"].lower() in POSITIVE_REGULATION,
            is_inhibition = row["Action"].lower() in NEGATIVE_REGULATION,
            endogenous = _endogenous,
            affinity_high = row["Affinity High"],
            affinity_low = row["Affinity Low"],
            affinity_median = row["Affinity Median"],
            affinity_units = row["Affinity Units"],
            primary_target = row["Primary Target"],
            pubmed_id = row["PubMed ID"],
            ligand = compounds.get(row["Ligand ID"]),
            target = target,
        )
def guide2pharma_table(name: TABLES) -> Generator[dict]:
    """
    Downloads a table from Guide2Pharma.

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


def get_taxid(organism_input: str) -> int:
    """
    Retrieves the NCBI Taxonomy ID for a given organism input.

    Args:
        organism_input: The name or identifier of the organism.

    Returns:
        int: The NCBI Taxonomy ID corresponding to the organism input.
        Returns a special constant if the input is empty or None.
    """
    
    # Check if the input is empty or None and return a constant value
    if organism_input in {"", "None", None}:
        return _const.NOT_ORGANISM_SPECIFIC

    # Use the taxonomy utility to ensure the input is converted to an NCBI Taxonomy ID
    return taxonomy.ensure_ncbi_tax_id(organism_input)

def g2p_get_target(
        row: dict,
        compounds: dict[str, G2PLigandProtein | G2PCompound],
        protein_targets: dict[str, G2PTargetProtein],
        ) -> G2PTargetProtein | G2PLigandProtein | G2PCompound | None:
    """
    Retrieves a target object from the given row in the Guide2Pharma
    interactions table. This object is either a G2PTargetProtein or a G2PCompound,
    depending on whether the target is a protein or a ligand.

    Args:
        row: A dictionary representing a row of the Guide2Pharma table.
        compounds: A dictionary of ligands with their IDs as keys and values as
        named tuples.
        protein_targets: A dictionary of protein targets with their IDs as keys and
        values as named tuples.

    Returns:
        G2PTargetProtein | G2PLigandProtein | G2PCompound | None: The target object
        corresponding to the given row, or None if the Target ID or Ligand Target ID is missing.
    """

    target : G2PTargetProtein | G2PLigandProtein | G2PCompound | None = None

    # Check if the target is a protein
    if row["Target ID"]:
        target = protein_targets.get(row["Target ID"]) # protein target

    # Check if the target is a ligand
    elif row["Target Ligand ID"]:
        target = compounds.get(row["Target Ligand ID"]) # ligand target

    # If neither Target ID nor Ligand Target ID is present, set target to None
    else:
        target = None

    return target


def g2p_compounds() -> dict[str, G2PLigandProtein | G2PCompound]:
    """
    Downloads ligands from Guide2Pharma.

    Ligands can be either proteins or small molecules. The function returns a
    dictionary of named tuples containing the name, PubChem ID, ChEMBL ID, IUPAC
    name, SMILES, and InChI of each ligand.

    Returns:
        dict[str, G2PLigandProtein | G2PCompound]: A dictionary of ligands with their
            IDs as keys and values as named tuples.
    """

    ligands : dict[str, G2PLigandProtein | G2PCompound] = {}
    for row in guide2pharma_table("ligands"):
        if row["UniProt ID"]: # check if ligand is a protein
            # If the ligand has a UniProt ID, it is a protein
            ligands[row["Ligand ID"]] = _source_record(row)
        else:
            # If the ligand does not have a UniProt ID, it is a small molecule
            ligands[row["Ligand ID"]] = _compound_record(row)

    return ligands

def g2p_protein_targets() -> dict[str, G2PTargetProtein]:
    """
    Downloads protein targets from Guide2Pharma.

    This function retrieves protein targets from the Guide2Pharma database.
    It returns a dictionary of named tuples, where the keys are the target IDs
    and the values are the G2PTargetProtein objects.

    Returns:
        dict[str, G2PTargetProtein]: A dictionary of protein targets.
    """

    targets = {}
    for row in guide2pharma_table("targets_and_families"):
        # Create a G2PTargetProtein object from the row
        targets[row["Target id"]] = _target_record(row)

    return targets
def _target_record(row: dict) -> G2PTargetProtein:
    """
    Creates a G2PProtein object from a row of the Guide2Pharma
    table.

    Args:
        row: A dictionary representing a row of the Guide2Pharma table.

    Returns:
        A G2PTargetProtein object.
    """
    
    record = G2PTargetProtein(
            family_name=row["Family name"],
            hgnc_name=row["HGNC name"],
            hgnc_symbol=row["HGNC symbol"],
            hgnc_id=row["HGNC id"],
            human_ensembl_gene=row["Human Ensembl Gene"],
            human_entrez_gene=row["Human Entrez Gene"],
            human_swissprot=row["Human SwissProt"],
            human_nucleotide_refseq=row["Human nucleotide RefSeq"],
            human_protein_refseq=row["Human protein RefSeq"],
            target_type=row["Type"],
            synonyms=row["synonyms"],
            record_type="protein",
        )

    return record

def _source_record(row: dict) -> G2PLigandProtein:
    """
    Creates a G2PLigandProtein object from a row of the Guide2Pharma
    table.

    Args:
        row: A dictionary representing a row of the Guide2Pharma table.

    Returns:
        A G2PLigandProtein object.
    """
    
    record = G2PLigandProtein(
            uniprot_id=row["UniProt ID"],
            name=row["Name"],
            pubchem_cid=row["PubChem CID"],
            inchi=row["InChI"],
            smiles=row["SMILES"],
            type=row["Type"],
            species=row["Species"],
            record_type="protein",
        )

    return record

def _compound_record(row: dict) -> G2PCompound:
    """
    Creates a G2PCompound object from a row of the Guide2Pharma
    table.

    Args:
        row: A dictionary representing a row of the Guide2Pharma table.

    Returns:
        A G2PCompound object.
    """
    
    record = G2PCompound(
            name=row["Name"],
            pubchem_cid=row["PubChem CID"],
            chembl=row["ChEMBL ID"],
            iupac=row["IUPAC name"],
            smiles=row["SMILES"],
            inchi=row["InChI"],
            compound_type=row["Type"],
            species=row["Species"],
            record_type="compound",
        )

    return record
