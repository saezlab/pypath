# TODO: process cases where Target= None properly (only 257 cases)

from typing import Literal
from collections.abc import Generator
import collections

import csv

import pypath_common._constants as _const
from pypath.share import curl
from pypath.utils import taxonomy
from pypath.resources import urls

__all__ = [
    'g2p_interactions',
    'g2p_table',
    'g2p_compounds',
    'g2p_protein_targets',
]

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
    "partial agonist",
    "inverse antagonist",
    "full agonist",
    "activation",
    "irreversible agonist",
    "positive",
    "biased agonist",
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
        "organism",
        "ncbi_taxa_id",
        "ligand_action",
        "ligand",
        "target",
    ],
)

G2PLigand = collections.namedtuple(
    "G2PLigand",
    [
        "name",
        "uniprot",
        "pubchem",
        "iupac",
        "chembl",
        "smiles",
        "inchi",
        "organism",
        "entity_type",
        "subtype",
        "role",
    ],
)

G2PLigand.__new__.__defaults__ = ('ligand',)


G2PTarget = collections.namedtuple(
    "G2PTarget",
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
        "entity_type",
        "role",
    ],
)

G2PTarget.__new__.__defaults__ = ('target',)


def g2p_interactions(
    organism: str | int | None = "human",
    endogenous: bool | None = None,
) -> Generator[tuple]:
    """
    Downloads interactions from Guide2Pharma and formats the data. Adds information
    from the Guide2Pharma 'ligands' and 'targets_and_families' tables.

    Args:
        organism (str | int | None):
            Name of the organism, e.g. `human`. If None, all organisms will be
            included.

        endogenous (bool | None):
            Whether to include only endogenous ligands interactions. If None,
            all ligands will be included.

    Yields:
        Named tuples containing information about each interaction.
    """
    # Retrieve the NCBI Taxonomy ID and common name for the given organism
    organism_, ncbi_tax_id = _organism_sorter(organism)

    # Download the ligands and protein targets from Guide2Pharma
    compounds = g2p_compounds()
    protein_targets = g2p_protein_targets()

    for row in g2p_table("interactions"):
        # Check if the interaction is endogenous
        _endogenous = row["Endogenous"].lower() == "true"

        # filters by users choice for endogenous ligands
        if endogenous is not None and endogenous != _endogenous:
            continue

        # Retrieve the NCBI Taxonomy ID and common name for the Target Species
        g2p_organism, g2p_tax_id = _organism_sorter(row["Target Species"])
        # filters by users choice for organism
        if organism_ is not None and g2p_tax_id != ncbi_tax_id:
            continue
        # Retrieve correct target object for the interaction
        target = _g2p_get_target(row, compounds, protein_targets)

        # Yield the interaction as a named tuple
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
            organism = g2p_organism,
            ncbi_taxa_id = g2p_tax_id,
            ligand_action=row["Type"],
            ligand = compounds.get(row["Ligand ID"]),
            target = target,
        )


def g2p_table(name: TABLES) -> Generator[dict]:
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

    yield from csv.DictReader(c.result)


def _organism_sorter(
        organism: str | int | None
        ) -> tuple[str | None, int | None]:
    """
    Takes an organism name or identifier from Guide to Pharmacology
    (e.g., "human") and converts it into both the common name and
    the NCBI Taxonomy ID.

    Args:
        organism (str | int | None): The name of the organism from Guide to Pharmacology.

    Returns:
        tuple: A tuple containing the common name and NCBI Taxonomy ID.
    """

    organism_ = None
    ncbi_tax_id = None

    if (isinstance(organism, str)
        and organism.lower() in taxonomy.taxids.values()):
        # If the organism is a string and it is a valid organism name
        ncbi_tax_id = _get_taxid(organism)

        try:
            # Get the common name for the organism
            organism_ = taxonomy.ensure_common_name(ncbi_tax_id)
            organism_ = organism_.capitalize() if organism_ else None
        except KeyError:
            pass

    return organism_, ncbi_tax_id


def _get_taxid(organism_input: str) -> int:
    """
    Retrieves the NCBI Taxonomy ID for a given organism input.

    Args:
        organism_input: The name or identifier of the organism.

    Returns:
        int: The NCBI Taxonomy ID corresponding to the organism input.
        Returns a special constant if the input is empty or None.
    """
    # ensure the input is not empty
    if organism_input:
        # Use the taxonomy utility to ensure the input is converted to an NCBI Taxonomy ID
        return taxonomy.ensure_ncbi_tax_id(organism_input)
    else:
        return _const.NOT_ORGANISM_SPECIFIC # TODO ask about _const.NOT_ORGANISM_SPECIFIC

def _g2p_get_target(
        row: dict,
        ligands: dict[str, G2PLigand],
        protein_targets: dict[str, G2PTargetProtein],
        ) -> G2PTargetProtein | G2PLigand | None:
    """
    Retrieves a target object from the compounds or protein_targets
    dictionaries based on the given row in the Guide2Pharma interactions table.
    This object is either a G2PTargetProtein or a G2PCompound, depending on
    whether the target is a protein or a ligand.

    Args:
        row: A dictionary representing a row of the Guide2Pharma table.
        compounds: A dictionary of ligands with their IDs as keys and values as
            named tuples.
        protein_targets: A dictionary of protein targets with their IDs as keys and
            values as named tuples.

    Returns:
        G2PTargetProtein | G2PLigandProtein | G2PCompound | None: The target object
        corresponding to the given row, or None if the Target ID or Ligand
        Target ID is missing.
    """

    target : G2PTarget | G2PLigand | None = None

    # Check if the target is a protein
    if row["Target ID"]:
        target = protein_targets.get(row["Target ID"])

    # Check if the target is a ligand
    elif row["Target Ligand ID"]:
        target = ligands.get(row["Target Ligand ID"])

    # If neither Target ID nor Ligand Target ID is present, set target to None
    else:
        target = None

    return target


def g2p_compounds() -> dict[str, G2PLigandProtein | G2PCompound]:
    """
    Downloads ligands from Guide2Pharma 'ligands' table.

    Ligands can be either proteins or small molecules. The function returns a
    dictionary of named tuples.

    Returns:
        dict[str, G2PLigandProtein | G2PCompound]: A dictionary of ligands with their
            IDs as keys and values as G2PLigandProtein or G2PCompund objects.
    """

    ligands : dict[str, G2PLigandProtein | G2PCompound] = {}
    for row in g2p_table("ligands"):
        ligands[row["Ligand ID"]] = _ligand(row)

    return ligands


def g2p_protein_targets() -> dict[str, G2PTargetProtein]:
    """
    Downloads protein targets from Guide2Pharma 'targets_and_families' table.

    This function retrieves protein targets from the Guide2Pharma database.
    It returns a dictionary of named tuples, where the keys are the target IDs
    and the values are the G2PTargetProtein objects.

    Returns:
        dict[str, G2PTargetProtein]: A dictionary of protein targets with their IDs
        as keys and G2PTargetProtein as objects.
    """

    targets = {}
    for row in g2p_table("targets_and_families"):
        # Create a G2PTargetProtein object from the row
        targets[row["Target id"]] = _target_record(row)

    return targets


def _target_record(row: dict) -> G2PTargetProtein:
    """
    Creates a G2PProtein object from a row of the Guide2Pharma 'targets_and_families'
    table.

    Args:
        row: A dictionary representing a row of the Guide2Pharma
            'targets_and_families' table.

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
            entity_type="protein",
        )

    return record

def _ligand(row: dict) -> G2PLigand:
    """
    Creates a G2PLigand object from a row of the Guide2Pharma ligands
    table.

    Args:
        row: A dictionary representing a row of the Guide2Pharma ligands table.

    Returns:
        A G2PLigand object.
    """

    record = G2PLigand(
        name=row["Name"],
        uniprot = row["UniProt ID"],
        pubchem=row["PubChem CID"],
        inchi=row["InChI"],
        smiles=row["SMILES"],
        iupac=row["IUPAC name"],
        chembl = row["ChEMBL ID"],
        organism=row["Species"],
        entity_type="protein" if row["UniProt ID"] else "compound",
        subtype=row["Type"],
    )

    return record
