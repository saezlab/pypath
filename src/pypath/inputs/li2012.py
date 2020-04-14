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

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common
import pypath.inputs.common as inputs_common
import pypath.utils.mapping as mapping
import pypath.internals.intera as intera


def get_li2012():
    """
    Reads supplementary data of Li 2012 from local file.
    Returns table (list of lists).
    """

    url = urls.urls['li2012']['url']
    c = curl.Curl(url, silent = False, large = True)
    xls = c.fileobj
    xlsfile = xls.name
    xls.close()
    tbl = inputs_common.read_xls(xlsfile, sheet = 'File S1')
    return filter(lambda l: len(l[-1]) > 0, map(lambda l: l[:7], tbl[2:]))


def li2012_interactions():
    """
    Converts table read by ``pypath.inputs.li2012.get_li2012`` to
    list of interactions.
    """

    result = []
    data = get_li2012()

    for l in data:

        subs_protein = l[1].split('/')[0]
        tk_protein = l[2].split()[0]
        reader_protein = l[3].split()[0]
        route = l[4]
        result.append((
            tk_protein,
            subs_protein,
            route,
            'phosphorylation'
        ))
        result.append((
            subs_protein,
            reader_protein,
            route,
            'phosphomotif_binding'
        ))

    return [list(l) for l in common.uniq_list(result)]


def li2012_enzyme_substrate():
    """
    Converts table read by `pypath.dataio.get_li2012()` to
    list of dicts of kinase-substrate interactions.
    """

    result = []
    non_digit = re.compile(r'[^\d]+')
    data = get_li2012()

    for l in data:

        subs_protein = l[1].split('/')[0]
        tk_protein = l[2].split()[0]
        subs_resnum = int(non_digit.sub('', l[1].split('/')[1]))
        result.append(
            (
                subs_protein, # substrate
                tk_protein, # kinase
                None, # instance
                None, # start
                None, # end
                'Y',  # residue letter
                subs_resnum, # residue offset
            )
        )

    result = [
        dict(
            zip(
                [
                    'substrate', 'kinase', 'instance',
                    'start', 'end', 'resaa',
                    'resnum',
                ],
                list(l)
            )
        )
        for l in common.uniq_list(result)
    ]

    return result


def li2012_dmi():
    """
    Converts table read by ``pypath.inputs.li2012.get_li2012`` to
    list of ``pypath.internals.intera.DomainMotif`` objects.
    Translates GeneSymbols to UniProt IDs.
    """

    result = {}
    nondigit = re.compile(r'[^\d]+')
    se = uniprot_input.swissprot_seq(isoforms = True)
    data = get_li2012()

    for l in data:

        subs_protein = l[1].split('/')[0]
        tk_protein = l[2].split()[0]
        reader_protein = l[3].split()[0]
        subs_uniprots = mapping.map_name(
            subs_protein,
            'genesymbol',
            'uniprot',
        )
        tk_uniprots = mapping.map_name(tk_protein, 'genesymbol', 'uniprot')
        reader_uniprots = mapping.map_name(reader_protein, 'genesymbol',
                                          'uniprot')
        subs_resnum = int(non_digit.sub('', l[1].split('/')[1]))

        for su in subs_uniprots:
            if su in se:
                subs_iso = None
                for iso, s in iteritems(se[su].isof):
                    if se[su].get(subs_resnum, isoform = iso) == 'Y':
                        subs_iso = iso
                        break
                if subs_iso:
                    start = min(1, subs_resnum - 7)
                    end = max(subs_resnum + 7, len(se[su].isof[subs_iso]))
                    for ku in tk_uniprots:
                        res = intera.Residue(
                            subs_resnum, 'Y', su, isoform = subs_iso)
                        mot = intera.Motif(
                            su,
                            start,
                            end,
                            isoform = subs_iso,
                            instance = se[su].get(
                                start,
                                end,
                                isoform = subs_iso
                            )
                        )
                        ptm = intera.Ptm(
                            su,
                            motif = mot,
                            residue = res,
                            isoform = subs_iso,
                            source = 'Li2012'
                        )
                        dom = intera.Domain(ku)
                        dommot = intera.DomainMotif(
                            domain = dom, ptm = ptm, sources = ['Li2012'])
                        result = {}

    return result
