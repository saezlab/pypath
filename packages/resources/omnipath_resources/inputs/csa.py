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

import json

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.progress as progress
import pypath.internals.intera as intera


def get_csa(uniprots = None):
    """
    Downloads and preprocesses catalytic sites data from M-CSA API.
    This data tells which residues are involved in the catalytic
    activity of one protein.
    """

    url = urls.urls['catalytic_sites']['url']
    c = curl.Curl(url, silent = False)
    data = c.result

    if data is None:

        return None

    try:
        json_data = json.loads(data)
    except json.JSONDecodeError:
        return None

    if not json_data:
        return None

# No longer need PDB chain mapping since we use UniProt sequence directly
    css = {}
    prg = progress.Progress(len(json_data), 'Processing catalytic sites', 11)

    for residue_data in json_data:

        mcsa_id = residue_data.get('mcsa_id')
        
        # Process both PDB chain and UniProt sequence information
        for seq_data in residue_data.get('residue_sequences', []):
            
            uniprot = seq_data.get('uniprot_id')
            resname = seq_data.get('code')
            resnum = seq_data.get('resid')
            
            if uniprot and resname and resnum and mcsa_id:
                
                if uniprots is None or uniprot in uniprots:
                    
                    # Use the UniProt sequence position directly
                    if uniprot not in css:
                        css[uniprot] = {}
                    
                    # Use 'uniprot_seq' as the key instead of PDB ID
                    if 'uniprot_seq' not in css[uniprot]:
                        css[uniprot]['uniprot_seq'] = {}
                    
                    if mcsa_id not in css[uniprot]['uniprot_seq']:
                        css[uniprot]['uniprot_seq'][mcsa_id] = []
                    
                    css[uniprot]['uniprot_seq'][mcsa_id].append(
                        intera.Residue(
                            name = resname,
                            number = int(resnum),
                            protein = uniprot,
                        )
                    )

        prg.step()

    prg.terminate()

    return css
