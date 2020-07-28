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

    positive_terms = {'activation', 'expression'}
    negative_terms = {'inhibition', 'repression'}
    transc_terms = {'expression', 'repression'}
    mechanism_terms = {
        'phosphorylation',
        'binding/association',
        'dissociation',
        'ubiquitination',
        'dephosphorylation',
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
            'transcriptional',
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
                effect = (
                    'inhibition'
                        if negative_terms & subtypes else
                    'activation'
                        if positive_terms & subtypes else
                    'unknown'
                )
                mechanism = ';'.join(mechanism_terms & subtypes)
                transcriptional = bool(transc_terms & subtypes)

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
                                transcriptional = transcriptional,
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


def kegg_medicus_interactions():


    reentity = re.compile(r'[,\+\(\)]|\w+')


    KeggMedicusInteraction = collections.namedtuple(
        'KeggMedicusInteraction',
        [
            'id_a',
            'id_b',
            'name_a',
            'name_b',
            'effect',
            'itype',
            'pw_type',
            'type_a',
            'type_b',
        ],
    )


    i_code = {
        '->': ('post_translational', 'stimulation'),
        '=>': ('transcriptional', 'stimulation'),
        '//': ('post_translational', 'missing'),
        '-|': ('post_translational', 'inhibition'),
        '=|': ('transcriptional', 'inhibition'),
        '--': ('post_translational', 'undirected'),
    }


    def process_entity(e):

        if isinstance(e, common.basestring):

            e = reentity.findall(e)

        sub = 0
        stack = []
        cplex = False

        for it in e:

            if it == ',':

                continue

            elif it == ')':

                sub -= 1

                if not sub:

                    stack.append(process_entity(this_stack))

                else:

                    this_stack.append(it)

            elif sub:

                this_stack.append(it)

                if it == '(':

                    sub += 1

            elif it == '(':

                if not sub:

                    this_stack = []

                sub += 1

            elif it == '+':

                cplex = True

            else:

                stack.append(it)

        if cplex:

            stack = tuple(stack)

        return stack


    def flatten_entity(e):

        flat = []

        if isinstance(e, common.basestring):

            flat.append(e)

        elif isinstance(e, tuple):

            flat.extend(
                itertools.product(*(
                    (c,) if isinstance(c, common.basestring) else c
                    for c in e
                ))
            )

        elif isinstance(e, list):

            flat.extend(itertools.chain(*(flatten_entity(c) for c in e)))

        return flat


    def get_interactions(connections, enames, pw_type):

        entities = dict(
            (
                i,
                flatten_entity(process_entity(connections[i]))
            )
            for i in range(0, len(connections), 2)
        )

        for i in range(0, len(connections) - 1, 2):

            itype, effect = i_code[connections[i + 1]]

            for id_a, id_b in itertools.product(entities[i], entities[i + 2]):

                name_a, type_a = get_name_type(id_a, enames)
                name_b, type_b = get_name_type(id_b, enames)

                yield KeggMedicusInteraction(
                    id_a = id_a,
                    id_b = id_b,
                    name_a = name_a,
                    name_b = name_b,
                    effect = effect,
                    itype = itype,
                    pw_type = pw_type,
                    type_a = type_a,
                    type_b = type_b,
                )


    def get_name_type(_id, enames):

        return (
            tuple(zip(*(enames[i] for i in _id)))
                if isinstance(_id, tuple) else
            enames[_id]
        )


    recollect = re.compile(r'^(GENE|PERTURBANT|VARIANT|METABOLITE)')
    result = set()
    url = urls.urls['kegg_pws']['medicus']
    c = curl.Curl(url, silent = False, large = True)
    enames = {}
    collecting = None

    for row in c.result:

        begin_coll = recollect.match(row)

        if begin_coll:

            collecting = begin_coll.group()
            row = row.split(maxsplit = 1)[-1]

        if collecting:

            if not begin_coll and row[0] != ' ':

                collecting = None
                continue

            if collecting == 'GENE':

                row = row.split(';')[0]

            try:
                _id, name = row.split(maxsplit = 1)
            except ValueError:
                print(row)
                return
            enames[_id] = (name, collecting.lower())

    c.fileobj.seek(0)

    for row in c.result:

        if row.startswith('ENTRY'):

            pw_type = None
            collecting = None
            print(row.split()[1])

        elif row.startswith('TYPE'):

            pw_type = row.strip().split()[-1].lower()

        elif row.startswith('  EXPANDED'):

            connections = row.split()[1:]

        elif row.startswith('///'):

            result.update(
                set(get_interactions(connections, enames, pw_type))
            )

    return result