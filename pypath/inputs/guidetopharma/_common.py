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


def _record_factory(
        name: str,
        fields: list[str],
        defaults: tuple = (),
    ) -> collections.namedtuple:

    RecordBase = collections.namedtuple(name, fields)
    RecordBase.__new__.__defaults__ = defaults

    class Record(RecordBase):

        def __new__(cls, *args, **kwargs):

            args = [arg if arg != '' else None for arg in args]
            kwargs = {k: (v if v != '' else None) for k, v in kwargs.items()}

            return super().__new__(cls, *args, **kwargs)

    return Record


G2PInteraction = _record_factory(
    "G2PInteraction",
    [
        "ligand",
        "target",
        "action",
        "action_type",
        "is_stimulation",
        "is_inhibition",
        "endogenous",
        "affinity_high",
        "affinity_low",
        "affinity_median",
        "affinity_units",
        "primary_target",
        "pubmed",
    ],
)

G2PLigand = _record_factory(
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
    ('ligand',),
)


G2PTarget = _record_factory(
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
    ('target',),
)


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
    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

    # Download the ligands and protein targets from Guide2Pharma
    ligands = g2p_ligands()
    targets = g2p_protein_targets()

    for row in g2p_table("interactions"):

        _endogenous = row["Endogenous"].lower() == "true"

        if endogenous is not None and endogenous != _endogenous:
            continue

        the_targets = (None,)

        if row["Target ID"]:

            the_targets = targets.get(row["Target ID"]) or (None,)

        elif row["Target Ligand ID"]:

            the_targets = (ligands.get(row["Target Ligand ID"]),)


        if the_targets is None:
            print(' ')
            print(row)
            print(the_targets)
            continue

        for target in the_targets:

            if target is None:

                target = G2PTarget(
                    *(None,) * 8,
                    organism = taxonomy.ensure_ncbi_tax_id(row['Target Species']),
                    entity_type = None,
                )

            if (
                ncbi_tax_id and
                target.organism != _const.NOT_ORGANISM_SPECIFIC and
                ncbi_tax_id != target.organism
            ):

                continue

            yield G2PInteraction(
                ligand = ligands.get(row['Ligand ID']),
                target = target,
                action = row['Action'],
                action_type = row['Type'],
                is_stimulation = row['Action'].lower() in POSITIVE_REGULATION,
                is_inhibition = row['Action'].lower() in NEGATIVE_REGULATION,
                endogenous = _endogenous,
                affinity_high = row['Affinity High'],
                affinity_low = row['Affinity Low'],
                affinity_median = row['Affinity Median'],
                affinity_units = row['Affinity Units'],
                primary_target = row['Primary Target'],
                pubmed = row['PubMed ID'],
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


def g2p_ligands() -> dict[str, G2PLigand | G2PLigand]:
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

        ligands[row["Ligand ID"]] = _parse_ligand(row)

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
            uniprot = row[f'{org} SwissProt'],
            symbol = row[f'{NAME_AUTHORITIES[org]} symbol'],
            entrez = row[f'{org} Entrez Gene'],
            ensembl = row[f'{org} Ensembl Gene'],
            refseq = row[f'{org} nucleotide RefSeq'],
            refseqp = row[f'{org} protein RefSeq'],
            family = row["Family name"],
            target_type = row["Type"],
            organism = taxonomy.ensure_ncbi_tax_id(org),
            entity_type = 'protein',
        )

        if any(x is not None for x in record[:6]):

            yield record


def _parse_ligand(row: dict) -> G2PLigand:
    """
    Parse ligand data from a row of the "ligands" table.

    Creates a G2PLigand object from a row of the Guide2Pharma ligands
    table.

    Args:
        row: A dictionary representing a row of the Guide2Pharma ligands table.

    Returns:
        A G2PLigand object.
    """

    record = G2PLigand(
        name = row['Name'],
        uniprot = row['UniProt ID'],
        pubchem = row['PubChem CID'],
        inchi = row['InChI'],
        smiles = row['SMILES'],
        iupac = row['IUPAC name'],
        chembl  =  row['ChEMBL ID'],
        organism = (
            taxonomy.ensure_ncbi_tax_id(row['Species']) or
            _const.NOT_ORGANISM_SPECIFIC
        ),
        entity_type = 'protein' if row['UniProt ID'] else 'compound',
        subtype = row['Type'].lower(),
    )

    return record
