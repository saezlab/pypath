#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import re
import collections
import itertools

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.common as inputs_common
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera
import pypath.core.entity as entity


IcellnetRecord = collections.namedtuple(
    'IcellnetRecord',
    [
        'ligand',
        'receptor',
        'family',
        'subfamily',
        'classification',
        'resources',
        'references',
    ]
)


def icellnet_interactions():

    url = urls.urls['icellnet']['url']

    c = curl.Curl(url, silent = False, large = True)

    xls = c.fileobj
    xlsfile = xls.name
    xls.close()
    tbl = inputs_common.read_xls(xlsfile)

    for line in tbl[1:]:

        references = _icellnet_get_references(line)
        resources = _icellnet_get_resources(line)

        if resources:
            references.extend([r for r in resources if r.isdigit()])
            resources = [r for r in resources if not r.isdigit()]

        ligand_components = _icellnet_get_components(line, (0, 1))
        receptor_components = _icellnet_get_components(line, (2, 3, 4))

        ligand = _icellnet_get_entity(ligand_components, references)
        receptor = _icellnet_get_entity(receptor_components, references)

        yield IcellnetRecord(
            ligand = ligand,
            receptor = receptor,
            family = line[6].strip() or None,
            subfamily = line[7].strip() or None,
            classification = (
                [
                    cls.strip().replace('.', '').capitalize()
                    for cls in
                    line[8].split('/')
                ]
                    if line[8].strip() else
                None
            ),
            resources = resources,
            references = references,
        )


def icellnet_complexes():

    complexes = {}

    for ia in icellnet_interactions():

        for attr in ('ligand', 'receptor'):

            if hasattr(getattr(ia, attr), 'components'):

                cplex = getattr(ia, attr)
                cplex_str = cplex.__str__()

                if cplex_str in complexes:

                    complexes[cplex_str] += cplex

                else:

                    complexes[cplex_str] = cplex

    return complexes


def icellnet_annotations(complexes = None):

    IcellnetAnnotation = collections.namedtuple(
        'IcellnetAnnotation',
        [
            'role',
            'family',
            'subfamily',
            'classification',
        ]
    )


    def get_entities(ia, entity_attr, complexes):

        entities = getattr(ia, entity_attr)

        if not entities:

            return ()

        complex_entities = (
            (entities,)
                if entity.Entity._is_complex(entities) else
            ()
        )
        protein_entities = (
            (entities,)
                if entity.Entity._is_protein(entities) else
            tuple(entities.components.keys())
        )

        return (
            complex_entities
                if complexes == True else
            protein_entities
                if complexes == False else
            complex_entities + protein_entities
        )


    annotations = collections.defaultdict(set)

    for ia in icellnet_interactions():

        for role in ('ligand', 'receptor'):

            for en in get_entities(ia, role, complexes):

                annotations[en].add(
                    IcellnetAnnotation(
                        role = role,
                        family = ia.family,
                        subfamily = ia.subfamily,
                        classification = (
                            tuple(sorted(ia.classification))
                                if ia.classification else
                            None
                        ),
                    )
                )

    return annotations


def _icellnet_get_components(line, idx):

    return [
        uniprot
        for uniprot in
        (
            mapping.map_name0(line[i], 'genesymbol', 'uniprot')
            for i in idx
        )
        if uniprot
    ]


def _icellnet_get_entity(components, references):

    if len(components) > 1:

        return intera.Complex(
            components = components,
            sources = 'ICELLNET',
            references = references,
        )

    elif len(components) == 1:

        return components[0]


def _icellnet_get_references(line):

    return [
        str(int(float(ref)))
        for ref in
        (_ref.strip() for _ref in line[10].split(';'))
        if ref
    ]


def _icellnet_get_resources(line):

    rerami = re.compile(r'(Ramilowski)\d{4}')
    resource_synonyms = {
        'Signor': 'SIGNOR',
        'guidetopharmacology.org': 'Guide2Pharma',
        'IUPHAR': 'Guide2Pharma',
        'IUPHAR-DB': 'Guide2Pharma',
        'GO_lig_rec': 'GO-lig-rec',
        'CellPhone': 'CellPhoneDB',
        'SignaLink': 'SignaLink3',
        'Innate': 'InnateDB',
        'Kegg': 'KEGG',
    }

    resources = line[9].replace(
        'Dinarello et al.2013 (Immunity)',
        'Dinarello2013'
    )
    resources = {
        rerami.sub(
            r'\g<1>2015',
            resource_synonyms.get(res, res)
        )
        for res in
        (_res.strip() for _res in re.split(r'[/,; ]', resources))
    }

    resources.discard('')
    resources.discard('DB')

    return sorted(resources) or None
