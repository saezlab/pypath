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

import re
import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.internals.intera as intera
import pypath.utils.reflists as reflists


def pdb_uniprot():
    """
    Mapping between UniProt and PDB identifiers.
    Returns two dictionaries:
    * UniProt to PDB mapping: keys are UniProt IDs, values are sets of
      tuples, each tuple with three values: the PDB structure ID, the
      structure analysis method and the structure resolution
    * PDB to UniProt mapping: keys are PDB IDs, values are sets of UniProt
      IDs
    """

    c = curl.Curl(urls.urls['uniprot_pdb']['url'], silent = False)
    data = c.result

    if data is None:
        return None, None

    data = data.split('\n')
    u_pdb = collections.defaultdict(set)
    pdb_u = collections.defaultdict(set)
    pdb = None
    pdb_re = re.compile(r'[0-9A-Z]{4}')

    for l in data:

        l = re.split(
            '[ ]{2,}',
            re.sub(
                '[ ]+,[ ]+', ',',
                re.sub(
                    r'[ ]*\(', '(',
                    l
                )
            )
        )

        if len(l[0]) == 4 and pdb_re.match(l[0]):

            pdb = l[0].lower()
            res = None if l[2] == '-' else float(l[2].replace(' A', ''))
            met = l[1]

        if pdb is not None and len(l) > 1:

            uniprots = l[1] if len(l) < 4 else l[3]
            uniprots = {
                u.split('(')[1].replace(')', '')
                for u in uniprots.split(',')
                if '(' in u
            }
            pdb_u[pdb].update(uniprots)

            for u in uniprots:

                u_pdb[u].add((pdb, met, res))

    return dict(u_pdb), dict(pdb_u)


def pdb_chains():
    """
    Amino acid chain level mapping between PDB and UniProt.
    Returns two dictionaries:
    * The first has UniProt IDs as keys and lists of dicts as values. Each
      of these dicts defines a mapping between UniProt and PDB amino acid
      chains with the chain identifier, PDB structure identifier and the
      start and end of the chain in the UniProt sequence and the PDB
      structure; the offset value is an integer if the PDB and the UniProt
      chain are the same length, otherwise None.
    * The second dict has PDB IDs as keys and dicts of chain mapping dicts
      as values, which are similar to the ones in the previous point, but
      here the chain identifiers are the keys.
    """

    def to_int(i):

        if i == 'None':

            return None

        return int(non_digit.sub('', i))


    c = curl.Curl(urls.urls['pdb_chains']['url'], silent = False)
    chains = c.result

    if chains is None:

        return None, None

    chains = chains.replace('\r', '').split('\n')
    del chains[0]
    del chains[0]
    pdb_u = {}
    u_pdb = {}
    non_digit = re.compile(r'[^\d.-]+')

    for l in chains:

        l = l.split('\t')

        if len(l) > 8:

            if l[0] not in pdb_u:

                pdb_u[l[0]] = {}

            pdb_u[l[0]][l[1]] = {
                'uniprot': l[2],
                'chain_beg': to_int(l[3]),
                'chain_end': to_int(l[4]),
                'pdb_beg': to_int(l[5]),
                'pdb_end': to_int(l[6]),
                'uniprot_beg': to_int(l[7]),
                'uniprot_end': to_int(l[8])
            }

            if (
                pdb_u[l[0]][l[1]]['pdb_end'] is not None and
                pdb_u[l[0]][l[1]]['pdb_beg'] is not None and
                pdb_u[l[0]][l[1]]['uniprot_beg'] is not None and
                pdb_u[l[0]][l[1]]['uniprot_end'] is not None and
                (
                    pdb_u[l[0]][l[1]]['pdb_end'] -
                    pdb_u[l[0]][l[1]]['pdb_beg'] ==
                    pdb_u[l[0]][l[1]]['uniprot_end'] -
                    pdb_u[l[0]][l[1]]['uniprot_beg']
                )
            ):

                pdb_u[l[0]][l[1]]['offset'] = (
                    pdb_u[l[0]][l[1]]['uniprot_beg'] -
                    pdb_u[l[0]][l[1]]['pdb_beg']
                )

            else:

                pdb_u[l[0]][l[1]]['offset'] = None

            if l[2] not in u_pdb:

                u_pdb[l[2]] = []

            u_pdb[l[2]].append({
                'pdb': l[0],
                'chain': l[1],
                'chain_beg': to_int(l[3]),
                'chain_end': to_int(l[4]),
                'pdb_beg': to_int(l[5]),
                'pdb_end': to_int(l[6]),
                'uniprot_beg': to_int(l[7]),
                'uniprot_end': to_int(l[8]),
                'offset': pdb_u[l[0]][l[1]]['offset']
            })

    return u_pdb, pdb_u


def pdb_complexes(organism = None):
    """
    Extracts protein complex data from PDB. The complexes are returned in
    a dict with string keys and ``pypath.internals.intera.Complex``
    objects as values. These latter carry their constitution, stoichiometry
    and the PDB identifiers.
    """

    complexes = {}

    uniprot_pdb, pdb_uniprot = pdb_chains()
    del uniprot_pdb

    for pdb_id, chains in iteritems(pdb_uniprot):

        uniprots = tuple(chain['uniprot'] for chain in chains.values())

        if len(uniprots) == 1:
            continue

        # if the organism set and any of the UniProt IDs does not
        # belong to this organism we drop the complex
        if organism and reflists.is_not(uniprots, 'uniprot', organism):

            continue

        cplex = intera.Complex(
            components = uniprots,
            sources = 'PDB',
            ids = pdb_id,
        )

        if cplex.__str__() in complexes:

            complexes[cplex.__str__()] += cplex

        else:

            complexes[cplex.__str__()] = cplex

    return complexes
