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

import re
import itertools
import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera
import pypath.inputs.pubmed as pubmed
import pypath.inputs.rdata as rdata


def cellchatdb_download(organism = 9606, dataset = 'CellChatDB'):
    """
    Retrieves data from CellChatDB, extracts fhe R objects from the RDA
    format and returns the raw data.

    :param int,str organism:
        Human and mouse are available, for incomprehensive values falls back
        to human.
    :param str dataset:
        Either *CellChatDB* or *PPI*.
    """

    organism = taxonomy.ensure_common_name(organism)
    organism = (
        'mouse'
            if organism and organism.lower() == 'mouse' else
        'human'
    )
    dataset = 'PPI' if dataset.lower() == 'ppi' else 'CellChatDB'

    key = '%s.%s' % (dataset, organism)
    url = urls.urls['cellchatdb']['url'] % key

    c = curl.Curl(url, silent = False, large = True)
    rdata_path = c.fileobj.name
    c.fileobj.close()

    rdata_parsed = rdata.rdata.parser.parse_file(rdata_path)
    rdata_converted = rdata.rdata.conversion.convert(rdata_parsed)[key]

    if dataset == 'CellChatDB':

        df_names = rdata._rdata_list_get_names(rdata_parsed.object.value[0])

        rownames = [
            rdata._rdata_data_frame_get_rownames(df_obj)
            for df_obj in rdata_parsed.object.value[0].value
        ]

        rownames = dict(zip(df_names, rownames))

        for name, df in rdata_converted.items():

            df['rownames'] = rownames[name]

    return rdata_converted


def cellchatdb_complexes(organism = 9606):
    """
    Retrieves data from CellChatDB and processes protein complexes.

    :param int,str organism:
        Human and mouse are available, for incomprehensive values falls back
        to human.

    :return:
        Dictionary of complexes.
    """

    raw = cellchatdb_download(organism = organism)['complex']

    return _cellchatdb_process_complexes(raw, organism = organism)


def _cellchatdb_process_complexes(raw, organism = 9606):

    if isinstance(raw, dict):

        raw = raw['complex']

    ncbi_tax_id = _cellchatdb_organism(organism)

    complexes = {}

    for row in raw.itertuples():

        genesymbols = [c for c in row[1:-1] if c]

        uniprots = [
            mapping.map_name(
                gs,
                'genesymbol',
                'uniprot',
                ncbi_tax_id = ncbi_tax_id,
            )
            for gs in genesymbols
        ]

        uniprots = [up for up in uniprots if up]

        for components in itertools.product(*uniprots):

            cplex = intera.Complex(
                name = row.rownames,
                components = components,
                sources = 'CellChatDB',
                ncbi_tax_id = ncbi_tax_id,
            )
            complexes[cplex.__str__()] = cplex

    return complexes


def cellchatdb_cofactors(organism = 9606):

    raw = cellchatdb_download(organism = organism)

    return _cellchatdb_process_cofactors(raw, organism = organism)


def _cellchatdb_process_cofactors(raw, organism = 9606):

    if isinstance(raw, dict):

        raw = raw['cofactor']

    ncbi_tax_id = _cellchatdb_organism(organism)

    cofactors = {}

    for row in raw.itertuples():

        genesymbols = [c for c in row[1:-1] if c]

        uniprots = mapping.map_names(
            genesymbols,
            'genesymbol',
            'uniprot',
            ncbi_tax_id = ncbi_tax_id,
        )

        cofactors[row.rownames] = uniprots

    return cofactors


def cellchatdb_interactions(
    organism = 9606,
    ligand_receptor = True,
    cofactors = True,
):
    """
    Retrieves data from CellChatDB and processes interactions.

    :param int,str organism:
        Human and mouse are available, for incomprehensive values falls back
        to human.
    :param bool ligand_receptor:
        Include ligand-receptor interactions.
    :param bool cofactors:
        Include interactions between cofactors (agonists, antagonists),
        coreceptors and receptors.

    :return:
        Interactions as list of named tuples.
    """

    def process_name(name):

        return (
            (complexes_by_name[name],)
                if name in complexes_by_name else
            mapping.map_name(
                name,
                'genesymbol',
                'uniprot',
                ncbi_tax_id = ncbi_tax_id,
            )
        )


    CellChatDBInteraction = collections.namedtuple(
        'CellChatDBInteraction',
        [
            'id_a',
            'id_b',
            'role_a',
            'role_b',
            'effect',
            'pathway',
            'refs',
            'category',
        ]
    )

    repmid = re.compile('PMID:\s*?(\d+)')
    repmcid = re.compile('PMC\d+')

    ncbi_tax_id = _cellchatdb_organism(organism)

    raw = cellchatdb_download(organism = organism)

    _complexes = _cellchatdb_process_complexes(raw, organism = organism)
    _cofactors = _cellchatdb_process_cofactors(raw, organism = organism)

    complexes_by_name = dict(
        (cplex.name, cplex)
        for cplex in _complexes.values()
    )

    raw_ia = raw['interaction']

    result = []

    for row in raw_ia.itertuples():

        refs = (
            set(repmid.findall(row.evidence)) |
            set(pubmed.only_pmids(repmcid.findall(row.evidence)))
        )
        # PMIDs starting with 23209 are a mistake in CellChatDB
        refs = sorted(pmid for pmid in refs if not pmid.startswith('23209'))

        ligands = process_name(row.ligand)
        receptors = process_name(row.receptor)

        if ligand_receptor:

            i = 0

            for ligand, receptor in itertools.product(ligands, receptors):

                result.append(
                    CellChatDBInteraction(
                        id_a = ligand,
                        id_b = receptor,
                        role_a = 'ligand',
                        role_b = 'receptor',
                        effect = 'unknown',
                        pathway = row.pathway_name,
                        refs = refs,
                        category = row.annotation,
                    )
                )

        if cofactors:

            for cofactor_col, effect, targets in (
                ('agonist', 'stimulation', ligands),
                ('antagonist', 'inhibition', ligands),
                ('co_A_receptor', 'stimulation', receptors),
                ('co_I_receptor', 'inhibition', receptors),
            ):

                cofact_label = getattr(row, cofactor_col)

                if cofact_label not in _cofactors:

                    continue

                for cofactor, target in itertools.product(
                    _cofactors[cofact_label],
                    targets
                ):

                    result.append(
                        CellChatDBInteraction(
                            id_a = cofactor,
                            id_b = target,
                            role_a = cofactor_col,
                            role_b = (
                                'ligand'
                                    if cofactor_col.endswith('onist') else
                                'receptor'
                            ),
                            effect = effect,
                            pathway = row.pathway_name,
                            refs = set(),
                            category = row.annotation,
                        )
                    )

    return result


def cellchatdb_annotations(organism = 9606):

    CellChatDBAnnotation = collections.namedtuple(
        'CellChatDBAnnotation',
        [
            'role',
            'pathway',
            'category',
        ]
    )


    annotations = collections.defaultdict(set)
    interactions = cellchatdb_interactions(organism = organism)

    for ia in interactions:

        for side in ('a', 'b'):

            annotations[getattr(ia, 'id_%s' % side)].add(
                CellChatDBAnnotation(
                    role = getattr(ia, 'role_%s' % side),
                    pathway = ia.pathway,
                    category = ia.category,
                )
            )

    return dict(annotations)


def _cellchatdb_organism(organism = 9606):

    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)
    ncbi_tax_id = 10090 if ncbi_tax_id == 10090 else 9606

    return ncbi_tax_id
