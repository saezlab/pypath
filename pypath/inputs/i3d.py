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
import sys

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.pfam as pfam_input
import pypath.share.progress as progress


def get_i3d():
    """
    Interaction3D contains residue numbers in given chains in
    given PDB stuctures, so we need to add an offset to get the residue
    numbers valid for UniProt sequences. Offsets can be obtained from
    Instruct, or from the Pfam PDB-chain-UniProt mapping table.
    """

    dname_pfam, pfam_dname = pfam_input.pfam_names()

    if dname_pfam is None:
        sys.stdout.write('\n\t:: Could not get Pfam domain names\n\n')

    non_digit = re.compile(r'[^\d.-]+')
    c = curl.Curl(urls.urls['i3d_human']['url'], silent = False)
    data = c.result

    if data is None:
        return None

    data = data.replace('\r', '').split('\n')
    del data[0]
    i3d = []
    prg = progress.Progress(
        len(data), 'Processing domain-domain interactions', 11)

    for l in data:
        prg.step()
        l = l.split('\t')

        if len(l) > 20:
            domain1 = None if l[13] not in dname_pfam else dname_pfam[l[13]]
            domain2 = None if l[20] not in dname_pfam else dname_pfam[l[20]]
            pdb = l[5]
            uniprot1 = l[0]
            uniprot2 = l[1]
            chain1 = l[7]
            seq1 = [[
                int(non_digit.sub('', l[11])), int(non_digit.sub('', l[12]))
            ]]
            chain2 = l[14]
            seq2 = [[
                int(non_digit.sub('', l[18])), int(non_digit.sub('', l[19]))
            ]]
            i3d.append({
                uniprot1: {
                    'pfam': domain1,
                    'chain': chain1,
                    'seq': seq1
                },
                uniprot2: {
                    'pfam': domain2,
                    'chain': chain2,
                    'seq': seq2
                },
                'uniprots': [uniprot1, uniprot2],
                'source': 'I3D',
                'pdb': [pdb],
                'references': []
            })
    prg.terminate()

    return i3d
