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
import collections
import itertools

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.utils.orthology as orthology_mod
import pypath.internals.intera as intera
import pypath.inputs.common as inputs_common


def baccin2019_interactions(ncbi_tax_id = 9606):

    recamel = re.compile(r'(.+?)([A-Z][a-z])')
    recap = re.compile(r'(^[A-Z][a-z]|_[A-Z][a-z])(.+)')


    def camel_to_snake(value):

        return (
            recamel.sub(
                lambda m: m.group(1).lower() + '_' + m.group(2),
                value.strip()
            ).lower()
        )


    def id_translate(mouse_gs):

        uniprots = mapping.map_name(
            mouse_gs,
            'genesymbol',
            'uniprot',
            10090,
        )

        if ncbi_tax_id != 10090:

            uniprots = set(
                itertools.chain(*(
                    orthology_mod.translate(
                        uniprot,
                        target = ncbi_tax_id,
                        source = 10090,
                    )
                    for uniprot in uniprots
                ))
            )

        return uniprots


    def raw_to_uniprots(raw):

        components = raw.split('&')

        return set(
            itertools.product(
                *(id_translate(comp) for comp in components)
            )
        )


    def get_partners(components, sources, references):

        return {
            (
                comp[0]
                    if len(comp) == 1 else
                intera.Complex(
                    components = comp,
                    sources = sources,
                    references = references,
                )
            )
            for comp in components
        }


    Baccin2019Interaction = collections.namedtuple(
        'Baccin2019Interaction',
        [
            'ligand',
            'receptor',
            'correct',
            'ligand_location',
            'ligand_category',
            'resources',
            'references',
        ]
    )


    source_names = {
        'Baccin': 'Baccin2019',
        'Ramilowski': 'Ramilowski2015',
    }

    url = urls.urls['baccin2019']['url']
    c = curl.Curl(url, silent = False, large = True)
    data = inputs_common.read_xls(c.fileobj.name, sheet = 'SuppTable3')

    result = []

    for rec in data[3:]:

        if rec[4].strip().lower() == 'incorrect':

            continue

        ligand_components = raw_to_uniprots(rec[1])

        if not ligand_components:

            continue

        receptor_components = raw_to_uniprots(rec[2])

        if not receptor_components:

            continue

        sources = {'Baccin2019', rec[3].strip()}
        sources = {
            source_names[s] if s in source_names else s
            for s in sources
        }

        references = {
            _ref for _ref in
            (
                ref.strip().replace('.0', '')
                for ref in rec[7].split(',')
            )
            if _ref.isdigit()
        }

        ligands = get_partners(ligand_components, sources, references)
        receptors = get_partners(receptor_components, sources, references)

        for ligand, receptor in itertools.product(ligands, receptors):
            result.append(
                Baccin2019Interaction(
                    ligand = ligand,
                    receptor = receptor,
                    correct = rec[4].strip(),
                    ligand_location = camel_to_snake(rec[5]),
                    ligand_category = camel_to_snake(rec[6]),
                    resources = sources,
                    references = references,
                )
            )

    return result


def baccin2019_annotations(ncbi_tax_id = 9606):


    Baccin2019Annotation = collections.namedtuple(
        'Baccin2019Annotation',
        [
            'mainclass',
            'subclass',
            'location',
        ]
    )


    ia_all = baccin2019_interactions(ncbi_tax_id = ncbi_tax_id)

    result = collections.defaultdict(set)

    for ia in ia_all:

        result[ia.ligand].add(
            Baccin2019Annotation(
                mainclass = 'ligand',
                subclass = ia.ligand_category,
                location = ia.ligand_location,
            )
        )

        result[ia.receptor].add(
            Baccin2019Annotation(
                mainclass = 'receptor',
                subclass = (
                    '%s_receptor' % ia.ligand_category
                        if ia.ligand_category != 'other' else
                    None
                ),
                location = None,
            )
        )

    return dict(result)
