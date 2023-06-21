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

try:
    from cStringIO import StringIO
except ModuleNotFoundError:
    from io import StringIO

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.progress as progress
import pypath.inputs.pdb as pdb_input
import pypath.inputs.common as inputs_common
import pypath.utils.pdb as pdb_utils
import pypath.internals.intera as intera


def get_csa(uniprots = None):
    """
    Downloads and preprocesses catalytic sites data.
    This data tells which residues are involved in the catalytic
    activity of one protein.
    """

    url = urls.urls['catalytic_sites']['url']
    c = curl.Curl(url, silent = False)
    data = c.result

    if data is None:

        return None

    u_pdb, pdb_u = pdb_input.pdb_chains()
    buff = StringIO()
    buff.write(data)
    cols = {
        'pdb': 0,
        'id': 1,
        'resname': 2,
        'chain': 3,
        'resnum': 4,
        'chem_fun': 5,
        'evidence': 6,
    }
    table = inputs_common.read_table(
        cols = cols,
        fileObject = buff,
        sep = ',',
        hdr = 1,
    )
    css = {}
    prg = progress.Progress(len(table), 'Processing catalytic sites', 11)

    for l in table:

        if l['pdb'] in pdb_u:

            if l['chain'] in pdb_u[l['pdb']]:

                uniprot = pdb_u[l['pdb']][l['chain']]['uniprot']

                if uniprots is None or uniprot in uniprots:

                    offset = pdb_u[l['pdb']][l['chain']]['offset']

                    if offset is not None:

                        l['resnum'] = int(l['resnum']) + offset

                    else:

                        this_res = pdb_utils.residue_pdb(
                            l['pdb'],
                            l['chain'],
                            l['resnum'],
                        )

                        if len(this_res) > 0:
                            l['resnum'] = int(this_res['UPCOUNT'])

                        else:
                            l['resnum'] = None

                    if l['resnum'] is not None:

                        if uniprot not in css:

                            css[uniprot] = {}

                        if l['pdb'] not in css[uniprot]:

                            css[uniprot][l['pdb']] = {}

                        if l['id'] not in css[uniprot][l['pdb']]:

                            css[uniprot][l['pdb']][l['id']] = []

                        css[uniprot][l['pdb']][l['id']].append(
                            intera.Residue(
                                name = l['resname'],
                                number = l['resnum'],
                                protein = uniprot,
                            )
                        )

        prg.step()

    prg.terminate()

    return css
