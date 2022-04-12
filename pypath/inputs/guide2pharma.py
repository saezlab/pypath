#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import csv
import collections
import itertools

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.utils.taxonomy as taxonomy
import pypath.internals.intera as intera


def guide2pharma_download(
        organism = 'human',
        endogenous = True,
        process_interactions = True,
        process_complexes = True,
    ):
    """
    Downloads and processes Guide to Pharmacology data.
    Returns list of dicts.

    @organism : str
        Name of the organism, e.g. `human`.
    @endogenous : bool
        Whether to include only endogenous ligands interactions.
    """

    get_taxid = taxonomy.taxid_from_common_name

    if isinstance(organism, common.basestring):
        try:
            organism = taxonomy.taxid_from_common_name(organism)

        except KeyError:
            organism = None

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

    url = urls.urls['gtp']['url']

    c = curl.Curl(url, silent = False, large = True, encoding = 'utf-8')

    line0 = next(c.result)

    if line0[0] != '#':

        c.fileobj.seek(0)

    data = csv.DictReader(c.result)

    if organism is not None:

        data = [
            d for d in data
            if (
                get_taxid(d['target_species']) == organism and
                organism in set(
                    get_taxid(t)
                    for t in d['ligand_species'].split('|')
                )
            )
        ]

    if endogenous:
        data = [d for d in data if d['endogenous'].strip() == 't']

    for d in data:
        if is_positive(d['type']) or is_positive(d['action']):
            effect = 1

        elif is_negative(d['type']) or is_negative(d['action']):
            effect = -1

        else:
            effect = 0

        for ligand_taxon in d['ligand_species'].split('|'):
            ligand_taxid = get_taxid(ligand_taxon)

            ligands = d['ligand_gene_symbol'] or d['ligand_pubchem_sid']
            ligands = ligands.split('|')
            targets = (
                d['target_uniprot'] or
                d['target_ligand_uniprot'] or
                d['target_ligand_pubchem_sid']
            )
            targets = targets.split('|')
            references = d['pubmed_id'].split('|') if d['pubmed_id'] else []

            if process_interactions:
                for ligand, target in itertools.product(ligands, targets):
                    interactions.append(
                        GuideToPharmacologyInteraction(
                            ligand = ligand,
                            ligand_id_type = (
                                'genesymbol'
                                    if d['ligand_gene_symbol'] else
                                'pubchem_sid'
                                    if d['ligand_pubchem_sid'] else
                                None
                            ),
                            target = target,
                            target_id_type = (
                                'uniprot'
                                    if (
                                        d['target_uniprot'] or
                                        d['target_ligand_uniprot']
                                    ) else
                                'pubchem_sid'
                                    if d['target_ligand_pubchem_sid'] else
                                None
                            ),
                            target_is_ligand = bool(d['target_ligand']),
                            ligand_organism = ligand_taxid,
                            target_organism = get_taxid(d['target_species']),
                            effect = effect,
                            ligand_location = (
                                d['ligand_context'].strip().lower() or None
                            ),
                            target_type = (
                                d['receptor_site'].strip().lower() or None
                            ),
                            ligand_endogenous = (
                                d['endogenous'].strip() == 't'
                            ),
                            pubmed_ids = references,
                        )
                    )

            if process_complexes:
                if (
                    len(targets) > 1 and (
                        d['target_uniprot'] or
                        d['target_ligand_uniprot']
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
                    d['ligand_gene_symbol']
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
