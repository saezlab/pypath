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

import itertools
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.progress as progress
import pypath.share.common as common
import pypath.inputs.pdb as pdb_input
import pypath.inputs.pfam as pfam_input
import pypath.internals.intera as intera


def threedcomplex_complexes():
    """
    To be implemented later. Should return dictionary of
    pypath.internals.intera.Complex objects.
    """

    raise NotImplementedError


def threedcomplex_ddi(contacts = None):
    """
    Downloads and preprocesses data from the 3DComplex database
    (http://shmoo.weizmann.ac.il/elevy/3dcomplexV6/Home.cgi).

    Args
        contacts (set): A map of 3D structure based contacts, as provided by
            ``threedcomplex_contacts``. If None, it will be obtained
            automatically.

    Returns
        A list of domain-domain interaction objects.
    """

    contacts = contacts or threedcomplex_contacts()
    uniprots = (
        common.values(contacts, 'uniprot_1') |
        common.values(contacts, 'uniprot_2')
    )
    u_pfam = pfam_input.pfam_regions(uniprots, value = 'uniprot')

    ddi = []
    prg = progress.Progress(
        len(contacts),
        'Processing contact information',
        9,
    )

    for con in contacts:

        prg.step()

        pdb = common.prefix(con.pdb, '_')
        pfams1 = {common.prefix(x, '.') for x in con.domain_p1}
        pfams2 = {common.prefix(x, '.') for x in con.domain_p2}

        for pfam1, pfam2 in itertools.product(pfams1, pfams2):

            pfam1_details = [{
                'start': None,
                'end': None,
                'isoform': 1
            }]
            pfam2_details = [{
                'start': None,
                'end': None,
                'isoform': 1
            }]

            if con.uniprot_1 in u_pfam and pfam1 in u_pfam[con.uniprot_1]:

                pfam1_details = u_pfam[con.uniprot_1][pfam1]

            if con.uniprot_2 in u_pfam and pfam2 in u_pfam[con.uniprot_2]:

                pfam2_details = u_pfam[con.uniprot_2][pfam2]

            for pfam1_d, pfam2_d in itertools.product(
                pfam1_details,
                pfam2_details,
            ):

                dom1 = intera.Domain(
                    protein = con.uniprot_1,
                    domain = pfam1,
                    start = pfam1_d['start'],
                    end = pfam1_d['end'],
                    isoform = pfam1_d['isoform'],
                    chains = {pdb: con.chain_1},
                )
                dom2 = intera.Domain(
                    protein = con.uniprot_2,
                    domain = pfam2,
                    start = pfam2_d['start'],
                    end = pfam2_d['end'],
                    isoform = pfam2_d['isoform'],
                    chains = {pdb: con.chain_2},
                )
                dd = intera.DomainDomain(
                    dom1,
                    dom2,
                    pdbs = pdb,
                    sources = '3DComplex',
                    contact_residues = con.n_residues,
                )
                ddi.append(dd)

    prg.terminate()

    return ddi


def threedcomplex_chains():
    """
    Returns the chain correspondancy map of the 3D Complex database.
    """

    c = curl.Curl(
        urls.urls['3dcomplex_correspondancy']['url'],
        silent = False,
    )
    corresp = c.result

    corresp = corresp.split('\n')
    corr_dict = collections.defaultdict(
        lambda: collections.defaultdict(dict)
    )

    # chain correspondancy map in a dict
    for l in corresp:

        l = l.strip().split('\t')

        if len(l) > 2:

            pdb = common.prefix(l[0], '.')
            corr_dict[pdb][l[1]] = l[2]

    return dict((k, dict(v)) for k, v in iteritems(corr_dict))


def threedcomplex_contacts(chains = None, pdb_uniprot = None):
    """
    Returns the 3D structure based domain-domain contact map from the
    3D Complex database.

    Args
        chains (dict): A dict of chain correspondancy map, as created by
            ``threedcomplex_chains``. If None, it will be obtained
            automatically.
        pdb_uniprot (dict): A dict of PDB-UniProt mappings, as created by
            ``pypath.inputs.pdb.pdb_chains``. If None, it will be obtained
            automatically.
    """

    ThreedcomplexContact = collections.namedtuple(
        'ThreedcomplexContact',
        (
            'pdb',
            'uniprot_1',
            'uniprot_2',
            'chain_1',
            'chain_2',
            'n_residues',
            'length_1',
            'length_2',
            'domain_s1',
            'domain_p1',
            'domain_s2',
            'domain_p2',
            'ident',
            'homo',
        ),
    )


    chains = chains or threedcomplex_chains()
    pdb_u = pdb_uniprot or pdb_input.pdb_chains()[1]

    c = curl.Curl(
        urls.urls['3dcomplex_contact']['url'],
        silent = False,
        slow = True,
    )
    contact = c.result

    result = set()

    for l in contact.split('\n'):

        l = l.strip().split('\t')

        if len(l) > 11:

            compl = l[0]
            pdb = common.prefix(compl, '_')

            if (
                compl in chains and
                l[1] in chains[compl] and
                l[2] in chains[compl]
            ):

                ch1 = chains[compl][l[1]]
                ch2 = chains[compl][l[2]]

                if (
                    pdb in pdb_u and
                    ch1 in pdb_u[pdb] and
                    ch2 in pdb_u[pdb]
                ):

                    up1 = pdb_u[pdb][ch1]['uniprot']
                    up2 = pdb_u[pdb][ch2]['uniprot']

                    result.add(
                        ThreedcomplexContact(
                            pdb = compl,
                            uniprot_1 = up1,
                            uniprot_2 = up2,
                            chain_1 = ch1,
                            chain_2 = ch2,
                            n_residues = float(l[3]),
                            length_1 = int(l[4]),
                            length_2 = int(l[5]),
                            domain_s1 = tuple(l[6].split(';')),
                            domain_s2 = tuple(l[8].split(';')),
                            domain_p1 = tuple(l[7].split(';')),
                            domain_p2 = tuple(l[9].split(';')),
                            ident = bool(int(l[10])),
                            homo = bool(int(l[11])),
                        )
                    )

    return result


def threedcomplex_nresidues():
    """
    Downloads and preprocesses data from the 3DComplex database
    (http://shmoo.weizmann.ac.il/elevy/3dcomplexV6/Home.cgi).

    Returns dict of dicts where top level keys are PDB IDs, second level
    keys are pairs of tuples of UniProt IDs and values are list with the
    number of amino acids in contact.
    """

    nresidues = collections.defaultdict(
        lambda: collections.defaultdict(list)
    )

    for contact in threedcomplex_contacts():

        uniprot_key = tuple(sorted((contact.uniprot_1, contact.uniprot_2)))
        nresidues[contact.pdb][uniprot_key] = contact.n_residues

    return dict((k, dict(v)) for k, v in iteritems(nresidues))
