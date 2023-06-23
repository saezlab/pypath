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

import os
import re
import collections
import itertools
import shutil

try:
    import cPickle as pickle
except:
    import pickle

import bs4

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.share.progress as progress
import pypath.share.cache as cache
import pypath.share.settings as settings
import pypath.share.session as session_mod
import pypath.resources.urls as urls
import pypath.internals.intera as intera
import pypath.utils.mapping as mapping

_logger = session_mod.Logger(name = 'hpmr_input')
_log = _logger._log


HpmrInteraction = collections.namedtuple(
    'HpmrInteraction',
    (
        'receptor',
        'partner_role',
        'partner',
        'references',
        'unambiguous',
    ),
)


def get_hpmr(use_cache = None):
    """
    Downloads ligand-receptor and receptor-receptor interactions from the
    Human Plasma Membrane Receptome database.

    Args
        use_cache (bool): Use the intermediate cache (pickle file of
            processed data).

    Returns
        (dict): Two elements: "interactions" and "families".
    """

    def get_partner(interactors, typ, recname = None, references = None):
        """
        typ : str
            `Receptor` or `Ligand`.
        """

        components = [i[1] for i in interactors if i[0] == typ]

        if typ == 'Receptor' and recname:

            components.append(recname)

        if len(components) == 1:

            return components[0]

        elif len(components) > 1:

            return intera.Complex(
                components = components,
                sources = 'HPMR',
                references = references,
            )


    cachefile = cache.cache_item('hpmr_preprocessed')
    use_cache = (
        use_cache
            if isinstance(use_cache, bool) else
        settings.get('use_intermediate_cache')
    )

    if os.path.exists(cachefile) and use_cache:

        _log('Reading HPMR data from cache file `%s`.' % cachefile)

        return pickle.load(open(cachefile, 'rb'))

    rerecname = re.compile(r'Receptor ([A-z0-9]+) interacts with:')
    reint = re.compile(r'(Receptor|Ligand) ([A-z0-9]+) -')
    rerefid = re.compile(r'list_uids=([- \.:,0-9A-z]+)')
    refamid = re.compile(r'.*FamId=([0-9\.]+)')

    a_family_title = 'Open Family Page'
    a_receptor_title = 'Open Receptor Page'
    a_titles = {a_family_title, a_receptor_title}

    interactions = []
    families = {}
    complexes = set()
    recpages = []

    c = curl.Curl(urls.urls['hpmri']['browse_rescued'])
    soup = bs4.BeautifulSoup(c.result, 'html.parser')

    this_family = ('0', None)
    this_subfamily = ('0', None)
    this_subsubfamily = ('0', None)

    for a in soup.find_all('a'):

        a_title = a.attrs['title'] if 'title' in a.attrs else None

        if a_title not in a_titles:

            continue

        if a_title == a_family_title:

            family_id = refamid.match(a.attrs['href']).groups()[0]

            if family_id.startswith(this_subfamily[0]):

                this_subsubfamily = (family_id, a.text)

            elif family_id.startswith(this_family[0]):

                this_subfamily = (family_id, a.text)
                this_subsubfamily = ('0', None)

            else:

                this_family = (family_id, a.text)
                this_subfamily = ('0', None)
                this_subsubfamily = ('0', None)

        elif a_title == a_receptor_title:

            recpages.append((
                a.attrs['href'],
                this_family[1],
                this_subfamily[1],
                this_subsubfamily[1],
            ))

    prg = progress.Progress(len(recpages), 'Downloading HPMR data', 1)

    i_complex = 0

    regene = re.compile(r'Param=([^&]+)&ProtId=(\d+)&ProtType=(\w+)')

    genes_curl = curl.Curl(
        urls.urls['hpmri']['genes_rescued'],
        silent = False,
        large = True,
    )

    for url, family, subfamily, subsubfamily in recpages:

        protein, prot_id, prot_type = regene.search(url).groups()
        fname = 'gene_%s-%s-%s.html' % (protein, prot_id, prot_type)
        prg.step(status = 'Processing `%s`' % fname)

        _log(
            'Accessing `%s` from `%s` (%s).' % (
                fname,
                genes_curl.cache_file_name,
                genes_curl.url,
            )
        )

        if fname not in genes_curl.result:

            _log('File `%s` not found in the archive.' % fname)
            continue

        soup = bs4.BeautifulSoup(
            genes_curl.result[fname].read(),
            'html.parser',
        )
        ints = soup.find('div', {'id': 'GeneInts'})

        if not ints:

            _log('No interactions: `%s`' % url)
            continue

        recname = rerecname.search(
            ints.find_previous_sibling('span').text
        )
        recname = recname.groups()[0] if recname else 'Unknown'

        if recname == 'Unknown':

            _log('Could not find receptor name: `%s`' % url)
            continue

        recname_u = mapping.map_name0(recname, 'genesymbol', 'uniprot')

        if not recname_u:

            continue

        families[recname_u] = (
            family,
            subfamily,
            subsubfamily,
        )

        for td in ints.find_all('td'):

            interactors = []

            for span in td.find_all('span', {'class': 'IntRow'}):

                ints = reint.search(span.text)

                if ints:

                    interactors.append(ints.groups())

            references = []

            for ref in td.find_all(
                'a', {'title': 'click to open reference in new window'}
            ):

                references.append(
                    rerefid.search(ref.attrs['href']).groups()[0].strip()
                )

            interactors_u = []

            for role, genesymbol in interactors:

                uniprot = (
                    mapping.map_name0(genesymbol, 'genesymbol', 'uniprot')
                )

                if uniprot:

                    interactors_u.append((role, uniprot))

            partner_role = (
                'receptor'
                    if all(i[0] == 'Receptor' for i in interactors_u) else
                'ligand'
            )

            receptors = (
                recname_u
                    if partner_role == 'receptor' else
                get_partner(
                    interactors_u,
                    'Receptor',
                    recname = recname_u,
                    references = references,
                )
            )

            partners = (
                {u[1] for u in interactors_u} - {recname_u}
                    if partner_role == 'receptor' else
                get_partner(
                    interactors_u,
                    'Ligand',
                    references = references,
                )
            )

            receptors = common.to_list(receptors)
            partners = common.to_list(partners)

            unambiguous = (
                partner_role == 'ligand' or
                (
                    len(receptors) == 1 and
                    len(partners) == 1
                )
            )

            for receptor, partner in itertools.product(receptors, partners):

                interactions.append(
                    HpmrInteraction(
                        receptor = receptor,
                        partner = partner,
                        partner_role = partner_role,
                        references = ';'.join(references),
                        unambiguous = unambiguous,
                    )
                )

            for entity in itertools.chain(receptors, partners):

                if hasattr(entity, 'components'):

                    complexes.add(entity)

    prg.terminate()

    result = {
        'interactions': interactions,
        'families': families,
    }

    _log('Saving HPMR data to cache file `%s`.' % cachefile)
    pickle.dump(result, open(cachefile, 'wb'))

    return result


def hpmr_complexes(use_cache = None):
    """
    HPMR does not contain unambiguous protein complex data, and considering
    the resource is unmaintained, probably it never will. Hence this function
    always returns an empty dict.
    """

    hpmr_data = get_hpmr(use_cache = use_cache)

    complexes = dict(
        (
            cplex.__str__(),
            cplex,
        )
        for cplex in hpmr_data.get('complexes', ())
    )

    return complexes


def hpmr_interactions(use_cache = None):

    hpmr_data = get_hpmr(use_cache = use_cache)

    return hpmr_data['interactions']


def hpmr_annotations(use_cache = None):

    annot = collections.defaultdict(set)

    HPMRAnnotation = collections.namedtuple(
        'HPMRAnnotation',
        ('role', 'mainclass', 'subclass', 'subsubclass'),
    )

    hpmr_data = get_hpmr(use_cache = use_cache)

    for i in hpmr_data['interactions']:

        # first partner is always a receptor
        # (because ligand pages simply don't work on HPMR webpage)
        args1 = ('Receptor',) + (
            hpmr_data['families'][i[0]]
                if i[0] in hpmr_data['families'] else
            (None, None, None)
        )

        # the second is either a ligand or another receptor
        args2 = (i[1],) + (
            hpmr_data['families'][i[2]]
                if i[2] in hpmr_data['families'] else
            (None, None, None)
        )

        annot[i[0]].add(HPMRAnnotation(*args1))
        annot[i[2]].add(HPMRAnnotation(*args2))

    for uniprot, classes in iteritems(hpmr_data['families']):

        args = ('Receptor',) + classes

        annot[uniprot].add(HPMRAnnotation(*args))

    return dict(annot)
