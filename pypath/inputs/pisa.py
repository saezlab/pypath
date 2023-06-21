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

from future.utils import iteritems
from past.builtins import xrange, range

import os
import sys
import re
import collections

try:
    import cPickle as pickle
except:
    import pickle

import bs4

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.cache as cache
import pypath.share.progress as progress
import pypath.utils.pdb as pdb_utils
import pypath.inputs.pdb as pdb_input
import pypath.internals.intera as intera


PisaBond = collections.namedtuple(
    'PisaBond',
    (
        'uniprot_1',
        'chain_1',
        'residue_1',
        'seqnum_1',
        'uniprot_2',
        'chain_2',
        'residue_2',
        'seqnum_2',
    ),
)


def pisa_bonds(bonds, chains):
    """
    To be refactored in the future. If you are interested in using this
    function, please contact the authors.
    """

    non_digit = re.compile(r'[^\d.-]+')
    result = []

    for bond in bonds.find_all('bond'):

        seqnum1 = int(non_digit.sub('', bond.find('seqnum-1').text))
        seqnum2 = int(non_digit.sub('', bond.find('seqnum-2').text))
        res1 = bond.find('res-1').text
        res1 = common.aaletters.get(res1, res1)
        res2 = bond.find('res-2').text
        res2 = common.aaletters.get(res2, res2)
        chain1 = bond.find('chain-1').text
        chain2 = bond.find('chain-2').text
        uniprot1 = chains.get(chain1, None)
        uniprot2 = chains.get(chain2, None)

        if uniprot1 and uniprot2:

            result.append(
                PisaBond(
                    uniprot_1 = uniprot1,
                    chain_1 = chain1,
                    residue_1 = res1,
                    seqnum_1 = seqnum1,
                    uniprot_2 = uniprot2,
                    chain_2 = chain2,
                    residue_2 = res2,
                    seqnum_2 = seqnum2,
                )
            )

    return result


def pisa_interfaces(pdbs, return_unmapped = False):
    """
    To be refactored in the future. If you are interested in using this
    function, please contact the authors.

    Args
        pdbs (set): A set of PDB IDs to query.
        return_unmapped (bool): Return also a list of unmapped residues.
        In this case, a tuple is returned, its first element is a dict of
        interfaces, while its second element is a list of unmapped
        residues (normally empty, if all residues could be mapped
        between PDB and UniProt sequences).

    Returns
        A dict of dicts with interfaces. The upper level keys are PDB
        structure IDs, the lower level keys are tuples of UniProt IDs,
        the values are ``pypath.internals.intera.Interface`` objects.
    """

    bond_types = {
        'hbonds': 'h-bonds',
        'sbridges': 'salt-bridges',
        'covbonds': 'cov-bonds',
        'ssbonds': 'ss-bonds'
    }
    interfaces = collections.defaultdict(dict)
    cachefile = os.path.join(cache.get_cachedir(), 'pisa.pickle')
    u_pdb, pdb_u = pdb_input.pdb_chains()

    if os.path.exists(cachefile):

        try:
            interfaces = pickle.load(open(cachefile, 'rb'))

        except:
            pass

    errors = []
    unmapped_residues = []
    p = 5
    pdbs = list(set(pdbs) - set(interfaces.keys()))
    prg = progress.Progress(
        len(pdbs) / p,
        'Downloading data from PDBe PISA',
        1,
    )

    pdbs = sorted(common.to_set(pdbs))

    for i in xrange(0, len(pdbs), p):

        to = i + p
        thisPart = pdbs[i:to]
        url = urls.urls['pisa_interfaces']['url'] + ','.join(thisPart)
        c = curl.Curl(url, cache = False)
        data = c.result

        if not data:

            msg = 'Could not download: \n\t\t%s' % url
            errors.append(msg)

            continue

        soup = bs4.BeautifulSoup(data, 'html.parser')

        for pdb in soup.find_all('pdb_entry'):

            pdb_id = pdb.find('pdb_code').text.lower()
            interfaces[pdb_id] = {}
            chains = {}
            resconv = pdb_utils.ResidueMapper()

            if pdb_id in pdb_u:

                for chain, chain_data in iteritems(pdb_u[pdb_id]):

                    chains[chain] = chain_data['uniprot']

                for interface in pdb.find_all('interface'):

                    for b, t in iteritems(bond_types):

                        bonds = interface.find(t)

                        if bonds:

                            bonds = pisa_bonds(bonds, chains)

                            for bond in bonds:

                                uniprots = (
                                    bond.uniprot_1,
                                    bond.uniprot_2,
                                )

                                if uniprots not in interfaces[pdb_id]:

                                    css = common.non_digit.sub(
                                        '', interface.find('css').text)
                                    css = (
                                        None if len(css) == 0 else float(css)
                                    )
                                    area = common.non_digit.sub(
                                        '', interface.find('int_area').text)
                                    area = None if len(area) == 0 else float(
                                        area)
                                    solv_en = common.non_digit.sub(
                                        '',
                                        interface.find('int_solv_en').text
                                    )
                                    solv_en = (
                                        None
                                            if len(solv_en) == 0 else
                                        float(solv_en)
                                    )
                                    stab_en = common.non_digit.sub(
                                        '',
                                        interface.find('stab_en').text
                                    )
                                    stab_en = (
                                        None
                                            if len(stab_en) == 0 else
                                        float(stab_en)
                                    )
                                    interfaces[pdb_id][uniprots] = (
                                        intera.Interface(
                                            uniprots[0],
                                            uniprots[1],
                                            source = 'PISA',
                                            pdb = pdb_id,
                                            css = css,
                                            solv_en = solv_en,
                                            area = area,
                                            stab_en = stab_en,
                                        )
                                    )

                                res1 = resconv.get_residue(
                                    pdb_id,
                                    bond.seqnum_1,
                                )
                                res2 = resconv.get_residue(
                                    pdb_id,
                                    bond.seqnum_2,
                                )

                                if (
                                    res1 and
                                    res2 and
                                    res1.uniprot == uniprots[0] and
                                    res2.uniprot == uniprots[1]
                                ):

                                    interfaces[pdb_id][uniprots].add_residues(
                                        (
                                            res1.resnum,
                                            bond.residue_1,
                                            uniprots[0],
                                        ),
                                        (
                                            res2.resnum,
                                            bond.residue_2,
                                            uniprots[1],
                                        ),
                                        typ = b,
                                    )

                                else:

                                    unmapped_residues.append(
                                        (
                                            pdb_id,
                                            bond.seqnum_1,
                                            bond.seqnum_2,
                                            uniprots[0],
                                            uniprots[1],
                                        )
                                    )

        prg.step()

    prg.terminate()
    pickle.dump(interfaces, open(cachefile, 'wb'), 2)
    interfaces = dict(interfaces)

    if len(errors) > 0:

        sys.stdout.write(
            '\t:: Failed to download %u files of total %u:\n\n' % (
                len(errors),
                len(pdbs)('P00968', 'P0A6F1')
            )
        )

        for e in errors:

            sys.stdout.write('\t' + e + '\n')

        sys.stdout.flush()

    if return_unmapped:
        return interfaces, unmapped_residues
    else:
        return interfaces
