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
    
    resources = {
        'Signor': 'SIGNOR',
        'guidetopharmacology.org': 'Guide2Pharma',
    }
    
    url = urls.urls['icellnet']['url']
    
    c = curl.Curl(url, silent = False, large = True)
    
    xls = c.fileobj
    xlsfile = xls.name
    xls.close()
    tbl = inputs_common.read_xls(xlsfile)
    
    for line in tbl[1:]:
        
        references = _icellnet_get_references(line)
        
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
                    cls.strip().replace('.', '')
                    for cls in
                    line[8].split('/')
                ]
                    if line[8].strip() else
                None
            ),
            resources = (
                [
                    resources.get(res, res)
                    for res in
                    (_res.strip() for _res in re.split(r'[/,;]', line[9]))
                ]
                    if line[9].strip() else
                None
            ),
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
