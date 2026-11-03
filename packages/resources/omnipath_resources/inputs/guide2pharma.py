#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
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

import csv
import collections
import itertools

import pypath.share.curl as curl
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy
import pypath.internals.intera as intera


def guide2pharma_download(
        organism: str | int = 'human',
        endogenous: bool = True,
        process_interactions: bool = True,
        process_complexes: bool = True,
    ) -> tuple[list, dict]:
    """
    Downloads and processes Guide to Pharmacology data.
    Returns list of dicts.

    Args:
        organism
            Name of the organism, e.g. `human`.
        endogenous
            Whether to include only endogenous ligands interactions.
    """

    get_taxid = lambda x: (
        _const.NOT_ORGANISM_SPECIFIC
            if x in {'', 'None'} else
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

    positives = {
        'agonist', 'activator', 'potentiation', 'partial agonist',
        'inverse antagonist', 'full agonist', 'activation',
        'irreversible agonist', 'positive',
    }
    negatives = {
        'inhibitor', 'antagonist', 'inhibition', 'irreversible inhibition',
        'inverse agonist', 'negative', 'weak inhibition',
        'reversible inhibition',
    }


    GuideToPharmacologyInteraction = collections.namedtuple(
        'GuideToPharmacologyInteraction',
        [
            'ligand',
            'ligand_id_type',
            'target',
            'target_id_type',
            'target_is_ligand',
            'ligand_organism',
            'target_organism',
            'effect',
            'ligand_location',
            'target_type',
            'ligand_endogenous',
            'pubmed_ids',
        ]
    )

    def is_positive(term):
        return term.lower().strip() in positives

    def is_negative(term):
        return term.lower().strip() in negatives

    interactions = []
    complexes = {}

    url = urls.urls['gtp']['url'] % 'interactions'

    c = curl.Curl(url, silent = False, large = True, encoding = 'utf-8')

    line0 = next(c.result)

    if line0[:2] != '"#':

        c.fileobj.seek(0)

    data = csv.DictReader(c.result)

    if organism_ is not None:

        data = [
            d for d in data
            if (
                get_taxid(d['Target Species']) == ncbi_tax_id and
                ncbi_tax_id in set(
                    get_taxid(t)
                    for t in d['Ligand Species'].split('|')
                )
            )
        ]

    if endogenous:

        data = [d for d in data if d['Endogenous'].strip() == 'true']

    for d in data:

        if is_positive(d['Type']) or is_positive(d['Action']):
            effect = 1

        elif is_negative(d['Type']) or is_negative(d['Action']):
            effect = -1

        else:
            effect = 0

        ligands = d['Ligand Gene Symbol'] or d['Ligand PubChem SID']
        ligands = ligands.split('|')
        ligand_taxons = [get_taxid(l) for l in d['Ligand Species'].split('|')]

        for ligand_taxon in zip(ligands, ligand_taxons):

            targets = (
                d['Target UniProt ID'] or
                d['Target Ligand UniProt ID'] or
                d['Target Ligand PubChem SID']
            )
            targets = targets.split('|')
            references = d['PubMed ID'].split('|') if d['PubMed ID'] else []

            if process_interactions:

                for ligand, target in itertools.product(ligands, targets):

                    interactions.append(
                        GuideToPharmacologyInteraction(
                            ligand = ligand,
                            ligand_id_type = (
                                'genesymbol'
                                    if d['Ligand Gene Symbol'] else
                                'pubchem_sid'
                                    if d['Ligand PubChem SID'] else
                                None
                            ),
                            target = target,
                            target_id_type = (
                                'uniprot'
                                    if (
                                        d['Target UniProt ID'] or
                                        d['Target Ligand UniProt ID']
                                    ) else
                                'pubchem_sid'
                                    if d['Target Ligand PubChem SID'] else
                                None
                            ),
                            target_is_ligand = bool(d['Target Ligand']),
                            ligand_organism = ligand_taxon,
                            target_organism = get_taxid(d['Target Species']),
                            effect = effect,
                            ligand_location = (
                                d['Ligand Context'].strip().lower() or None
                            ),
                            target_type = (
                                d['Receptor Site'].strip().lower() or None
                            ),
                            ligand_endogenous = (
                                d['Endogenous'].strip() == 't'
                            ),
                            pubmed_ids = references,
                        )
                    )

            if process_complexes:
                if (
                    len(targets) > 1 and (
                        d['Target UniProt ID'] or
                        d['Target Ligand UniProt ID']
                    )
                ):
                    cplex = intera.Complex(
                        components = targets,
                        sources = 'Guide2Pharma',
                        references = references,
                    )
                    key = cplex.__str__()

                    if key in complexes:
                        complexes[key] += cplex

                    else:
                        complexes[key] = cplex

                if (
                    len(ligands) > 1 and
                    d['Ligand Gene Symbol']
                ):
                    ligand_uniprots = [
                        mapping.map_name0(ligand, 'genesymbol', 'uniprot')
                        for ligand in ligands
                    ]
                    ligand_uniprots = [u for u in ligand_uniprots if u]

                    if len(ligand_uniprots) > 1:
                        cplex = intera.Complex(
                            components = ligand_uniprots,
                            sources = 'Guide2Pharma',
                            references = references,
                        )
                        key = cplex.__str__()

                        if key in complexes:
                            complexes[key] += cplex

                        else:
                            complexes[key] = cplex

    return interactions, complexes


def guide2pharma_interactions(**kwargs):

    interactions, complexes = guide2pharma_download(
        process_complexes = False,
        **kwargs
    )

    return interactions


def guide2pharma_complexes(**kwargs):

    interactions, complexes = guide2pharma_download(
        process_interactions = False,
        **kwargs
    )

    return complexes
