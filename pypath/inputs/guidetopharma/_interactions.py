# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

from collections.abc import Generator

import pypath_common._constants as _const
from pypath.utils import taxonomy
from . import _ligands, _targets, _raw
from ._records import G2PInteraction, G2PTarget
from ._constants import POSITIVE_REGULATION, NEGATIVE_REGULATION

__all__ = [
    'interactions',
]


def interactions(
        organism: str | int | None = "human",
        endogenous: bool | None = None,
    ) -> Generator[tuple]:
    """
    Interactions from Guide to Pharmacology.

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

    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

    ligands = _ligands.ligands()
    targets = _targets.protein_targets()

    for row in _raw.table('interactions'):

        _endogenous = row['Endogenous'].lower() == 'true'

        if endogenous is not None and endogenous != _endogenous:
            continue

        the_targets = (None,)

        if row['Target ID']:

            the_targets = targets.get(row['Target ID']) or (None,)

        elif row['Target Ligand ID']:

            the_targets = (ligands.get(row['Target Ligand ID']),)

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
