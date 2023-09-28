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

import itertools
import collections

import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.utils.orthology as orthology
import pypath.share.common as common
import pypath.inputs.cell as cell_input


def embrace_raw():
    """
    Returns Supplementary Table S11 from 10.1016/j.isci.2019.10.026
    (Sheikh et al. 2019) as a list of tuples.
    """

    path = cell_input.cell_supplementary(
        supp_url = urls.urls['embrace']['url'],
        article_url = urls.urls['embrace']['article'],
    )

    content = inputs_common.read_xls(path)

    EmbraceRawRecord = collections.namedtuple(
        'EmbraceRawRecord',
        content[0]
    )

    return [
        EmbraceRawRecord(*(line[:2] + [int(float(n)) for n in line[2:]]))
        for line in
        content[1:]
    ]


def _embrace_id_translation(mouse_genesymbol, organism = 9606):

    uniprots = mapping.map_name(
        mouse_genesymbol,
        'genesymbol',
        'uniprot',
        ncbi_tax_id = 10090,
    )

    if organism != 10090:

        uniprots = orthology.translate(
            uniprots,
            target = organism,
            source = 10090,
        )

    return uniprots or [None]


def embrace_translated(organism = 9606):
    """
    Returns Supplementary Table S11 from 10.1016/j.isci.2019.10.026
    (Sheikh et al. 2019) translated to UniProt IDs of the requested organism.
    """

    raw = embrace_raw()
    record = raw[0].__class__
    result = []

    for row in raw:

        ligands = _embrace_id_translation(
            row.ligand_symbol,
            organism = organism,
        )
        receptors = _embrace_id_translation(
            row.receptor_symbol,
            organism = organism,
        )

        for ligand, receptor in itertools.product(ligands, receptors):

            result.append(
                record(*((ligand, receptor) + row[2:]))
            )

    return result


def embrace_interactions(organism = 9606):
    """
    Returns ligand-receptor interactions from Supplementary Table S11 of
    10.1016/j.isci.2019.10.026 (Sheikh et al. 2019) translated to UniProt IDs
    of the requested organism.
    """


    EmbraceInteraction = collections.namedtuple(
        'EmbraceInteraction',
        [
            'ligand',
            'receptor',
        ]
    )

    return [
        EmbraceInteraction(*rec[:2])
        for rec in embrace_translated(organism = organism)
        if rec[0] and rec[1]
    ]


def embrace_annotations(organism = 9606, ncbi_tax_id = None):
    """
    Returns protein annotations from Supplementary Table S11 of
    10.1016/j.isci.2019.10.026 (Sheikh et al. 2019) translated to UniProt IDs
    of the requested organism.
    Returns dict with UniProt IDs as keys and sets of annotation tuples as
    values.
    """

    def _get_value(rec, mainclass, celltype):

        return bool(getattr(rec, '%s_%s' % (celltype, mainclass)))


    EmbraceAnnotation = collections.namedtuple(
        'EmbraceAnnotation',
        [
            'mainclass',
            'neuron',
            'mural_cell',
            'microglia',
            'endothelial_cell',
        ]
    )


    organism = ncbi_tax_id or organism
    result = collections.defaultdict(set)

    for rec in embrace_translated(organism = organism):

        for mainclass in ('ligand', 'receptor'):

            identifier = getattr(rec, '%s_symbol' % mainclass)

            if not identifier:

                continue

            result[identifier].add(
                EmbraceAnnotation(
                    mainclass = mainclass,
                    neuron = _get_value(rec, mainclass, 'N'),
                    mural_cell = _get_value(rec, mainclass, 'P'),
                    microglia = _get_value(rec, mainclass, 'M'),
                    endothelial_cell = _get_value(rec, mainclass, 'EC'),
                )
            )

    return dict(result)
