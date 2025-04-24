#!/usr/bin/env python
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

from __future__ import annotations

from collections.abc import Generator

from pypath.utils import taxonomy
from ._records import G2PTarget
from . import _raw
from . import _constants

__all__ = [
    'protein_targets',
]


def protein_targets() -> dict[str, G2PTarget]:
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

    for row in _raw.table('targets_and_families'):

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
            symbol = row[f'{_constants.NAME_AUTHORITIES[org]} symbol'],
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



