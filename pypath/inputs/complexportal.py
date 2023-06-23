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

import collections

import bs4

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.internals.intera as intera


def complexportal_complexes(organism = 9606, return_details = False):
    """
    Complex dataset from IntAct.
    See more:
    http://www.ebi.ac.uk/intact/complex/
    http://nar.oxfordjournals.org/content/early/2014/10/13/nar.gku975.full.pdf
    """

    spec = {9606: 'Homo_sapiens'}

    zipurl = '%s/%s.zip' % (
        urls.urls['complex_portal']['url'],
        spec[organism],
    )
    c = curl.Curl(zipurl, large = True, silent = False)
    files = c.result

    errors = []
    complexes = {}
    details = []
    name_key = 'complex recommended name'

    for xmlname, xml in iteritems(c.result):

        soup = bs4.BeautifulSoup(xml, 'html.parser')
        interactors_xml = soup.find_all('interactor')
        interactors = {}
        interactions = {}

        for i in interactors_xml:

            if i.find('primaryref').attrs['db'] == 'uniprotkb':

                interactors[i.attrs['id']] = i.find('primaryref').attrs['id']

        interactions_xml = soup.find_all('interaction')

        for i in interactions_xml:

            description = ''
            pubmeds = []
            fullname = ''
            names = {}
            pdbs = []
            uniprots = []
            ids = collections.defaultdict(set)

            for a in i.find_all('attribute'):
                if a.attrs['name'] == 'curated-complex':
                    description = a.text

            for sr in i.find_all('secondaryref'):
                if sr.attrs['db'] == 'pubmed':
                    pubmeds.append(sr.attrs['id'])

                if sr.attrs['db'] == 'wwpdb':
                    pdbs.append(sr.attrs['id'])

            for pr in i.find_all('primaryref'):
                if pr.attrs['db'] in {'wwpdb', 'rcsb pdb', 'pdbe'}:
                    pdbs.append(pr.attrs['id'])

            for sr in i.find('xref').find_all('secondaryref'):

                if (
                    'reftype' in sr.attrs and
                    sr.attrs['db'] in {'intact', 'reactome'} and
                    sr.attrs['reftype'] == 'identity'
                ):

                    ids[sr.attrs['db']].add(sr.attrs['id'])

            pubmeds = list(set(pubmeds))
            pdbs = list(set(pdbs))
            fullname = (
                None
                    if i.find('fullname') is None else
                i.find('fullname').text
            )

            for a in i.find_all('alias'):
                names[a.attrs['type']] = a.text

            for intref in i.find_all('interactorref'):
                int_id = intref.text

                if int_id in interactors:
                    uniprot = interactors[int_id]

                    if uniprot.startswith('PRO'):
                        continue

                    uniprot = uniprot.split('-')[0]
                    uniprots.append(uniprot)

            if uniprots:

                if pdbs:
                    ids['PDB'].update(set(pdbs))

                cplex = intera.Complex(
                    components = uniprots,
                    name = names[name_key] if name_key in names else None,
                    references = set(pubmeds),
                    sources = 'ComplexPortal',
                    ids = ids,
                )

                if cplex.__str__() in complexes:
                    complexes[cplex.__str__()] += cplex

                else:
                    complexes[cplex.__str__()] = cplex

            details.append({
                'uniprots': uniprots,
                'pdbs': pdbs,
                'pubmeds': pubmeds,
                'fullname': fullname,
                'names': names,
                'description': description
            })

    if return_details:
        return complexes, details

    else:
        return complexes
