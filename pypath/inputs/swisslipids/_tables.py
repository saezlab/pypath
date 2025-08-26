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


from . import _raw, _records

__all__ = ['swisslipids_lipids']


def swisslipids_lipids():

    for rec in _raw.swisslipids_lipids_raw():

        yield _records.SwisslipidsLipid(
            id = rec['Lipid ID'],
            level = rec['Level'],
            name = rec['Name'],
            abbreviation = rec['Abbreviation*'],
            synonyms = rec['Synonyms*'],
            class = rec['Class'],
            parent = rec['Parent'],
            components = rec['Components*'],
            smiles = rec['SMILES (pH7.3)'],
            inchi = rec['InChI (pH7.3)'],
            inchikey = rec['InChIKey (pH7.3)'],
            formula = rec['Formula (pH7.3)'],
            charge = rec['Charge (pH7.3)'],
            exact_mass = rec['Exact mass (neutral form)'],
            chebi = f'CHEBI:{rec["CHEBI"]}',
            lipidmaps = rec['LIPID MAPS'],
            hmdb = rec['HMDB'],
            metanetx = rec['MetaNetX'],
            pmids = rec['PMID'].split(' | '),
        )
