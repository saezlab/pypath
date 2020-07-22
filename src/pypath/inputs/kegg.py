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

from future.utils import iteritems

import re
import itertools
import collections
import bs4

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.progress as progress
import pypath.share.common as common
import pypath.utils.mapping as mapping


def kegg_interactions():
    """
    Downloads and processes KEGG Pathways.
    Returns list of interactions.
    """

    effect_terms = {'activation', 'inhibition'}
    mechanism_terms = {
        'phosphorylation',
        'binding/association',
        'expression',
        'dissociation',
        'ubiquitination',
        'dephosphorylation',
        'repression',
        'glycosylation',
        'state change',
        'methylation',
    }
    direct_terms = {'indirect effect'}

    KeggInteraction = collections.namedtuple(
        'KeggInteraction',
        [
            'id_a',
            'id_b',
            'effect',
            'pathway',
            'mechanism',
            'is_direct',
        ],
    )

    rehsa = re.compile(r'.*(hsa[0-9]+).*')
    req_hdrs = [
        'Referer: http://www.genome.jp/kegg-bin/show_pathway'
        '?map=hsa04710&show_description=show'
    ]
    hsa_list = []
    interactions = []

    c = curl.Curl(urls.urls['kegg_pws']['list_url'], silent = True)
    htmllst = c.result
    lstsoup = bs4.BeautifulSoup(htmllst, 'html.parser')

    for a in lstsoup.find_all('a', href = True):
        m = rehsa.match(a['href'])

        if m:
            hsa_list.append((m.groups(0)[0], a.text))

    prg = progress.Progress(
        len(hsa_list), 'Processing KEGG Pathways', 1, percent = False
    )

    for hsa, pw in hsa_list:

        prg.step()
        c = curl.Curl(
            urls.urls['kegg_pws']['kgml_url_2'] % hsa,
            silent = True,
            req_headers = req_hdrs
        )
        kgml = c.result
        kgmlsoup = bs4.BeautifulSoup(kgml, 'html.parser')
        entries = {}

        for ent in kgmlsoup.find_all('entry'):
            gr = ent.find('graphics')

            if gr and 'name' in gr.attrs:
                entries[ent.attrs['id']] = [
                    n.strip()
                    for n in gr.attrs['name'].replace('...', '').split(',')
                ]

        uentries = dict([(eid, common.uniq_list(
            common.flat_list([
                mapping.map_name(
                    gn, 'genesymbol', 'uniprot', strict = True) for gn in gns
            ]))) for eid, gns in iteritems(entries)])

        for rel in kgmlsoup.find_all('relation'):

            subtypes = {st.attrs['name'] for st in rel.find_all('subtype')}

            if (
                rel.attrs['entry1'] in uentries and
                rel.attrs['entry2'] in uentries and
                subtypes
            ):

                is_direct = 'indirect effect' not in subtypes
                effect = common.first(
                    effect_terms & subtypes,
                    default = 'unknown',
                )
                mechanism = ';'.join(mechanism_terms & subtypes)

                for u1 in uentries[rel.attrs['entry1']]:

                    for u2 in uentries[rel.attrs['entry2']]:

                        interactions.append(
                            KeggInteraction(
                                id_a = u1,
                                id_b = u2,
                                effect = effect,
                                pathway = pw,
                                mechanism = mechanism,
                                is_direct = is_direct,
                            )
                        )

    prg.terminate()

    return common.uniq_list(interactions)


def kegg_pathways():

    data = kegg_interactions()
    pws = common.uniq_list(map(lambda i: i[3], data))
    proteins_pws = dict(map(lambda pw: (pw, set([])), pws))
    interactions_pws = dict(map(lambda pw: (pw, set([])), pws))

    for u1, u2, eff, pw in data:

        proteins_pws[pw].add(u1)
        proteins_pws[pw].add(u2)
        interactions_pws[pw].add((u1, u2))

    return proteins_pws, interactions_pws


def kegg_pathway_annotations():

    KeggPathway = collections.namedtuple(
        'KeggPathway', ['pathway'],
    )


    result = collections.defaultdict(set)

    proteins, interactions = kegg_pathways()

    for pathway, uniprots in iteritems(proteins):
        record = KeggPathway(pathway = pathway)

        for uniprot in uniprots:
            result[uniprot].add(record)

    return result
