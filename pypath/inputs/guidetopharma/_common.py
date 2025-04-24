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

NAME_AUTHORITIES = {
    'Human': 'HGNC',
    'Rat': 'RGD',
    'Mouse': 'MGI',
}

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
        "uniprot",
        "symbol",
        "entrez",
        "ensembl",
        "refseq",
        "refseqp",
        "family",
        "target_type",
        "organism",
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
    organism_, ncbi_tax_id = _organism_common_ncbi(organism)

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
        g2p_organism, g2p_tax_id = _organism_common_ncbi(row["Target Species"])
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


def _organism_common_ncbi(
        organism: str | int | None
    ) -> tuple[str | None, int | None]:
    """
    Common name and NCBI Taxonomy ID of an organism.

    Takes an organism name or identifier from Guide to Pharmacology
    (e.g., "human") and converts it into both the common name and
    the NCBI Taxonomy ID.

    Args:
        organism:
            The name of the organism from Guide to Pharmacology.

    Returns:
        A tuple of the common name and the NCBI Taxonomy ID. If organism is
        unknown, the former defaults to None, the latter to -1, which means not
        organism specific molecule (e.g. a compound).
    """

    return (
        taxonomy.ensure_common_name(organism),
        taxonomy.ensure_ncbi_tax_id(organism) or _const.NOT_ORGANISM_SPECIFIC
    )


def _g2p_get_target(
        row: dict,
        ligands: dict[str, G2PLigand],
        protein_targets: dict[str, G2PTarget],
        ) -> G2PTarget | G2PLigand | None:
    """
    Retrieves a target object from the compounds or protein_targets
    dictionaries based on the given row in the Guide2Pharma interactions table.
    This object is either a G2PTarget or a G2PLigand, depending on
    whether the target is a protein or a ligand.

    Args:
        row: A dictionary representing a row of the Guide2Pharma table.
        compounds: A dictionary of ligands with their IDs as keys and values as
            named tuples.
        protein_targets: A dictionary of protein targets with their IDs as keys and
            values as named tuples.

    Returns:
        The target object corresponding to the given row, or None if the Target
        ID or Ligand Target ID is missing.
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


def g2p_compounds() -> dict[str, G2PLigand | G2PLigand]:
    """
    Downloads ligands from Guide2Pharma 'ligands' table.

    Ligands can be either proteins or small molecules. The function returns a
    dictionary of named tuples.

    Returns:
        dict[str, G2PLigand | G2PLigand]: A dictionary of ligands with their
            IDs as keys and values as G2PLigand or G2PCompund objects.
    """

    ligands: dict[str, G2PLigand | G2PLigand] = {}

    for row in g2p_table("ligands"):
        ligands[row["Ligand ID"]] = _ligand(row)

    return ligands


def g2p_protein_targets() -> dict[str, G2PTarget]:
    """
    Downloads protein targets from Guide2Pharma 'targets_and_families' table.

    This function retrieves protein targets from the Guide2Pharma database.
    It returns a dictionary of named tuples, where the keys are the target IDs
    and the values are the G2PTarget objects.

    Returns:
        dict[str, G2PTarget]: A dictionary of protein targets with their IDs
        as keys and G2PTarget as objects.
    """

    targets = {}

    for row in g2p_table('targets_and_families'):

        targets[row['Target id']] = list(_parse_targets(row))

    return targets


def _parse_targets(row: dict) -> Generator[G2PTarget]:
    """
    Parse protein target data from a row of the "targets_and_families" table.

    Creates a G2PTarget object from a row of the Guide2Pharma
    'targets_and_families' table.

    Args:
        row: A dictionary representing a row of the Guide2Pharma
            'targets_and_families' table.

    Returns:
        A G2PTarget object.
    """

    for org in ('Human', 'Mouse', 'Rat'):

        record = G2PTarget(
            uniprot = _or_none(row, f'{org} SwissProt'),
            symbol = _or_none(row, f'{NAME_AUTHORITIES[org]} symbol'),
            entrez = _or_none(row, f'{org} Entrez Gene'),
            ensembl = _or_none(row, f'{org} Ensembl Gene'),
            refseq = _or_none(row, f'{org} nucleotide RefSeq'),
            refseqp = _or_none(row, f'{org} protein RefSeq'),
            family = row["Family name"],
            target_type = row["Type"],
            organism = taxonomy.ensure_ncbi_tax_id(org),
            entity_type = 'protein',
        )

        if any(x is not None for x in record[:6]):

            yield record


def _or_none(row, key):

    return row[key] or None


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
