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