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

import pypath_common._constants as _const
from pypath.utils import taxonomy
from ._records import G2PLigand
from . import _raw

__all__ = [
    'ligands',
]


def ligands() -> dict[str, G2PLigand]:
    """
    Downloads ligands from Guide2Pharma 'ligands' table.

    Ligands can be either proteins or small molecules. The function returns a
    dictionary of named tuples.

    Returns:
        dict[str, G2PLigand | G2PLigand]: A dictionary of ligands with their
            IDs as keys and values as G2PLigand or G2PCompund objects.
    """

    ligands: dict[str, G2PLigand | G2PLigand] = {}

    for row in _raw.table("ligands"):

        ligands[row["Ligand ID"]] = _parse_ligand(row)

    return ligands


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
