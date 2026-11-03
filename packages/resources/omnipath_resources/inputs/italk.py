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
import pyreadr

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping
import pypath.share.common as common


def italk_raw():
    """
    Returns a ``pandas.DataFrame`` with the iTalk database contents.
    """
    
    url = urls.urls['italk']['url']
    c = curl.Curl(url, silent = False, large = True)
    rdata_path = c.fileobj.name
    c.fileobj.close()
    
    rdata = pyreadr.read_r(rdata_path)['database']
    
    return rdata


def italk_interactions():
    
    
    ItalkInteraction = collections.namedtuple(
        'ItalkInteraction',
        [
            'ligand',
            'receptor',
            'classification',
        ]
    )
    
    
    rdata = italk_raw()
    result = []
    
    for row in rdata.itertuples():
        
        if (
            not isinstance(row[2], str) or
            not isinstance(row[4], str)
        ):
            
            continue
        
        ligands = mapping.map_name(row[2], 'genesymbol', 'uniprot')
        receptors = mapping.map_name(row[4], 'genesymbol', 'uniprot')
        cls = row[6]
        
        for ligand, receptor in itertools.product(ligands, receptors):
            
            result.append(
                ItalkInteraction(
                    ligand = ligand,
                    receptor = receptor,
                    classification = cls,
                )
            )
    
    return result


def italk_annotations():
    
    
    ItalkAnnotation = collections.namedtuple(
        'ItalkAnnotation',
        [
            'mainclass',
            'subclass',
        ]
    )
    
    
    rdata = italk_raw()
    result = collections.defaultdict(set)
    
    for row in rdata.itertuples():
        
        ligands = (
            mapping.map_name(row[2], 'genesymbol', 'uniprot')
                if isinstance(row[2], str) else
            ()
        )
        receptors = (
            mapping.map_name(row[4], 'genesymbol', 'uniprot')
                if isinstance(row[4], str) else
            ()
        )
        subclass = row[6]
        
        for mainclass, uniprot in itertools.chain(
            itertools.zip_longest((), ligands, fillvalue = 'ligand'),
            itertools.zip_longest((), receptors, fillvalue = 'receptor'),
        ):
            
            result[uniprot].add(
                ItalkAnnotation(
                    mainclass = mainclass,
                    subclass = subclass,
                )
            )
    
    return dict(result)
