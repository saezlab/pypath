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

import bs4
import xml.etree.cElementTree as ET

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping


def netpath_interactions():

    result = []
    repwnum = re.compile(r'NetPath_([0-9]+)_')
    mi = '{net:sf:psidev:mi}'
    url = urls.urls['netpath_psimi']['url']
    c = curl.Curl(url, silent = False)

    data = c.result
    data = dict([(k, v) for k, v in iteritems(data) if k.endswith('xml')])
    pwnames = netpath_names()

    for pwfile, rawxml in iteritems(data):

        try:
            pwnum = repwnum.findall(pwfile)[0]

        except:
            sys.stdout.write('Error at processing file:\n')
            sys.stdout.write(pwfile)
            sys.stdout.write('\n')
            sys.stdout.flush()

        pwname = pwnames[pwnum]
        root = ET.fromstring(rawxml)

        for e in root.findall(mi + 'entry'):

            thisInt = ()
            db = [
                pr.find(mi + 'primaryRef').attrib['db']
                for pr in e.find(mi + 'source').findall(mi + 'xref')
            ]
            refs = []
            mets = []

            for ex in e.find(mi + 'experimentList').findall(
                    mi + 'experimentDescription'):
                for pm in ex.find(mi + 'bibref').iter(mi + 'primaryRef'):
                    if pm.attrib['db'] == 'pubmed':
                        refs.append(pm.attrib['id'])

                for me in ex.find(mi + 'interactionDetectionMethod').\
                        iter(mi + 'shortLabel'):
                    mets.append(me.text)

            mols = {}

            for mo in e.find(mi + 'interactorList').findall(mi + 'interactor'):
                iid = mo.attrib['id']
                name = mo.find(mi + 'names').find(mi + 'shortLabel').text
                entrez = ''

                if mo.find(mi + 'xref') is not None:
                    entrez = ';'.join([
                        ac.attrib['id']
                        for ac in mo.find(mi + 'xref')
                        .findall(mi + 'secondaryRef')
                        if ac.attrib['db'] == 'Entrez gene'
                    ])

                mols[iid] = (name, entrez)

            theInt = e.find(mi + 'interactionList').find(mi + 'interaction')

            for p in theInt.find(mi + 'participantList').findall(
                    mi + 'participant'):
                pid = p.find(mi + 'interactorRef').text
                roles = ''

                if p.find(mi + 'experimentalRoleList') is not None:
                    roles = ';'.join([
                        rl.find(mi + 'names').find(mi + 'shortLabel').text
                        for rl in p.find(mi + 'experimentalRoleList')
                        .findall(mi + 'experimentalRole')
                    ])

                mols[pid] += (roles, )

            intTyp = (
                theInt.find(
                    mi + 'interactionType'
                ).find(
                    mi + 'names'
                ).find(
                    mi + 'shortLabel'
                ).text
            )
            molkeys = list(mols.keys())

            for i in range(0, len(mols) - 1):

                for j in range(i, len(mols)):

                    A = mols[molkeys[i]][0:2]
                    B = mols[molkeys[j]][0:2]
                    result.append(
                        list(A) +
                        list(B) +
                        [
                            ';'.join(refs),
                            ';'.join(mets),
                            intTyp,
                            pwname
                        ]
                    )

    return result



def netpath_names():

    repwnum = re.compile(r'_([0-9]+)$')
    result = {}
    url = urls.urls['netpath_names']['url']
    c = curl.Curl(url, silent = False)
    html = c.result
    soup = bs4.BeautifulSoup(html, 'html.parser')

    for a in soup.find_all('a'):

        if a.attrs['href'].startswith('pathways'):

            num = repwnum.findall(a.attrs['href'])[0]
            name = a.text
            result[num] = name.strip()

    return result


def netpath_pathway_annotations():

    NetpathPathway = collections.namedtuple(
        'NetpathPathway',
        ['pathway'],
    )

    result = collections.defaultdict(set)

    url_template = urls.urls['netpath_pw']['url']

    url_main = urls.urls['netpath_pw']['mainpage']
    c = curl.Curl(url_main, cache = False)
    cookie = [
        h.decode().split(':')[1].split(';')[0].strip()
        for h in c.resp_headers
        if h.startswith(b'Set-Cookie')
    ]
    cookie_hdr = ['Cookie: %s' % '; '.join(cookie)]

    pathway_ids = netpath_names()

    for _id, pathway in iteritems(pathway_ids):

        url = url_template % int(_id)
        c = curl.Curl(
            url,
            req_headers = cookie_hdr,
            silent = False,
            encoding = 'iso-8859-1',
        )

        soup = bs4.BeautifulSoup(c.result, 'html.parser')

        for tbl in soup.find_all('table'):
            hdr = tbl.find('td', {'class': 'barhead'})

            if not hdr or not hdr.text.strip().startswith('Molecules Invol'):
                continue

            for td in tbl.find_all('td'):
                genesymbol = td.text.strip()

                if not genesymbol:
                    continue

                uniprots = mapping.map_name(
                    genesymbol,
                    'genesymbol',
                    'uniprot',
                )

                for uniprot in uniprots:
                    result[uniprot].add(
                        NetpathPathway(
                            pathway = pathway
                        )
                    )

    return dict(result)
