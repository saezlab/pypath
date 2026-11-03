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

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl


def lincs_compounds():
    """
    The returned dict has names, brand names or company specific IDs of
    compounds as keys, and tuples of PubChem, ChEMBL, ChEBI, InChi,
    InChi Key, SMILES and LINCS as values.
    """

    LincsCompound = collections.namedtuple(
        'LincsCompound',
        (
            'lincs',
            'chembl',
            'chebi',
            'inchi',
            'inchi_key',
            'smiles',
            'alternatives',
        )
    )

    c = curl.Curl(urls.urls['lincs-compounds']['url'], silent = False)
    
    # Parse SDF file
    result = {}
    
    # SDF files contain molecules separated by $$$$
    molecules = c.result.split('$$$$')
    
    for molecule in molecules:
        if not molecule.strip():
            continue
            
        # Parse the molecule data
        lines = molecule.strip().split('\n')
        if len(lines) < 4:
            continue
            
        # First line is the molecule name
        name = lines[0].strip()
        
        # Parse the property block (after the connection table)
        properties = {}
        in_properties = False
        
        for line in lines:
            if line.strip().startswith('>'):
                # Property header: >  <PROPERTY_NAME>
                prop_name = line.strip()[3:-1]  # Remove >  < and >
                in_properties = True
                current_prop = prop_name
            elif in_properties and line.strip() and not line.strip().startswith('>'):
                # Property value
                properties[current_prop] = line.strip()
                in_properties = False
        
        # Extract names from various fields
        names = []
        if name:
            names.append(name)
            
        # Add alternative names from properties
        for prop in ['name', 'common_name', 'synonym', 'alternative_name']:
            if prop in properties and properties[prop]:
                names.append(properties[prop])
                
        # Add LINCS ID if available
        lincs_id = properties.get('lincs_id') or properties.get('LINCS_ID')
        if lincs_id:
            names.append(lincs_id)
            
        # Create the compound object
        compound = LincsCompound(
            lincs = lincs_id,
            chembl = 'CHEMBL%s' % properties.get('chembl_id', '') if properties.get('chembl_id') else None,
            chebi = 'CHEBI:%s' % properties.get('chebi_id', '') if properties.get('chebi_id') else None,
            inchi = properties.get('inchi'),
            inchi_key = properties.get('inchi_key'),
            smiles = properties.get('smiles'),
            alternatives = properties.get('synonym', ''),
        )
        
        # Add each name as a key pointing to the compound
        for name in names:
            if name:
                result[name] = compound
    
    return result