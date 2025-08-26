#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2024
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

from . import _raw
from . import _records

__all__ = ['lipidmaps_lipids']


def lipidmaps_lipids():

    for rec in _raw.lipidmaps_raw():

        chebi = rec['name'].get('CHEBI_ID')
        chebi = f"CHEBI:{chebi}" if chebi else None

        yield _records.LipidmapsLipid(
            id = rec['id'],
            name = rec['annot'].get('NAME'),
            abbreviation = rec['name'].get('ABBREVIATION'),
            formula = rec['name']['FORMULA'],
            exact_mass = rec['annot']['EXACT_MASS'],
            category = rec['annot']['CATEGORY'],
            main_class = rec['annot']['MAIN_CLASS'],
            smiles = rec['name'].get('SMILES'),
            inchi = rec['name']['INCHI'],
            inchikey = rec['name']['INCHI_KEY'],
            synonyms = rec['name'].get('SYNONYMS'),
            iupac = rec['name'].get('SYSTEMATIC_NAME'),
            pubchem = rec['name'].get('PUBCHEM_CID'),
            chebi = chebi,
        )
