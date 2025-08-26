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


import pypath.internals.intera as _intera
from . import _raw, _records

__all__ = [
    'swisslipids_evidences',
    'swisslipids_lipids',
    'swisslipids_tissues',
    'swisslipids_reactions',
]


def swisslipids_lipids():

    for rec in _raw.swisslipids_lipids_raw():

        yield _records.SwisslipidsLipid(
            id = rec['Lipid ID'],
            level = rec['Level'],
            name = rec['Name'],
            abbreviation = rec['Abbreviation*'],
            synonyms = rec['Synonyms*'],
            lipid_class = rec['Class'],
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


def swisslipids_evidences() -> dict[str, _records.SwisslipidsEvidence]:

    result = {}

    for rec in _raw.swisslipids_evidences_raw():

        result[rec['Evidence ID']] = _records.SwisslipidsEvidence(
            eco = rec['ECO ID'],
            eco_name = rec['ECO definition'],
            pmid = rec['PMID ID'],
        )

    return result


def swisslipids_tissues():

    evidences = swisslipids_evidences()

    for rec in _raw.swisslipids_tissues_raw():

        yield _records.SwisslipidsTissue(
            lipid_id = rec['Lipid ID'],
            lipid_name = rec['Lipid name'],
            ncbi_tax_id = int(rec['Taxon ID']),
            taxon_name = rec['Taxon scientific name'],
            tissue_uberon = rec['Tissue/Cell ID'],
            tissue_name = rec['Tissue/Cell name'],
            evidence = evidences.get(rec['Evidence tag ID']),
        )


def swisslipids_reactions():

    evidences = swisslipids_evidences()

    for rec in _raw.swisslipids_enzymes_raw():

        if '-' in (enzyme := rec['UniProtKB AC(s)']):

            enzyme = _intera.Complex(
                components = enzyme.split('-'),
                ncbi_tax_id = int(rec['Protein taxon']),
                sources = 'SwissLipids',
            )

        yield _records.SwisslipidsReaction(
            enzyme_uniprot = enzyme,
            reaction_rhea = rec['Rhea ID'],
            ncbi_tax_id = int(rec['Protein taxon']),
            evidence = evidences.get(rec['Evidence tag ID']),
        )
