from typing import Literal
from collections.abc import Generator

import csv

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls

TABLES = Literal[
    'ligands',
    'interactions',
    'targets_and_families',
    'ligand_id_mapping',
    'ligand_physchem_properties',
    'endogenous_ligand_pairings_all',
    'endogenous_ligand_detailed',
    'approved_drug_primary_target_interactions',
    'approved_drug_detailed_interactions',
    'peptides',
    'GtP_to_HGNC_mapping',
    'GtP_to_UniProt_mapping'
]

POSITIVE_REGULATION = {
    'agonist',
    'activator', 
    'potentiation', 
    'partial agonist',
    'inverse antagonist', 
    'full agonist', 
    'activation',
    'irreversible agonist',
    'positive',
}
NEGATIVE_REGULATION = {
    'inhibitor',
    'antagonist',
    'inhibition', 
    'irreversible inhibition',
    'inverse agonist', 
    'negative', 
    'weak inhibition',
    'reversible inhibition',
}

def guide2pharma_table(name: TABLES) -> Generator[dict]:
    """
    Downloads the table from Guide2Pharma.

    Args:
        name:
            The name of the table to download.

    Returns:
        Generator yielding the table as named tuples.
    """

    url = urls.urls['gtp']['url'] % name

    c = curl.Curl(url, silent = False, large = True, encoding = 'utf-8')

    g2p_version = next(c.result).strip()

    return csv.DictReader(c.result)


def guide2pharma_interactions(
        organism: str | int | None = 'human',
        endogenous: bool | None = None,
    ) -> Generator[tuple]:
    """
    Args:
        organism
            Name of the organism, e.g. `human`. If None, all organisms will be
            included.
        endogenous
            Whether to include only endogenous ligands interactions. If None,
            all ligands will be included.
    """

    get_taxid = lambda x: (
        _const.NOT_ORGANISM_SPECIFIC
            if x in {'', 'None', None} else
        taxonomy.ensure_ncbi_tax_id(x)
    )
    organism_ = None
    ncbi_tax_id = None

    if isinstance(organism, str):

        ncbi_tax_id = get_taxid(organism)

        try:

            organism_ = taxonomy.ensure_common_name(ncbi_tax_id)
            organism_ = organism_.capitalize() if organism_ else None

        except KeyError:

            pass  # no organism specified



