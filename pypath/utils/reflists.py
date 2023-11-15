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
import json
import datetime
import time

try:
    import cPickle as pickle

except ImportError:
    import pickle

import timeloop

import pypath.inputs.uniprot_db as uniprot_input
import pypath.inputs.mirbase as mirbase_input
import pypath.share.common as common
import pypath.share.session as session_mod
import pypath.share.settings as settings
import pypath.share.cache as cache_mod


# method names for ID types
inputs = {
    'uniprot': 'all_uniprots',
    'swissprot': 'all_swissprots',
    'trembl': 'all_trembls',
    'mirbase': 'mirbase_mature_all',
    'mir-pre': 'mirbase_precursor_all',
}


_reflists_cleanup_timeloop = timeloop.Timeloop()
_reflists_cleanup_timeloop.logger.setLevel(9999)


class ReferenceListManager(session_mod.Logger):


    def __init__(self, cleanup_period = 10, lifetime = 300):

        session_mod.Logger.__init__(self, name = 'reflists')


        @_reflists_cleanup_timeloop.job(
            interval = datetime.timedelta(
                seconds = cleanup_period
            )
        )
        def _cleanup():

            self._remove_expired()


        _reflists_cleanup_timeloop.start(block = False)

        self.lifetime = lifetime
        self.lists = {}
        self.expiry = {}
        self.cachedir = cache_mod.get_cachedir()

        self._log('ReferenceListManager has been created.')


    def which_list(self, id_type, ncbi_tax_id = None):

        ncbi_tax_id = ncbi_tax_id or settings.get('default_organism')

        key = (id_type, ncbi_tax_id)

        self.expiry[key] = time.time()

        if key not in self.lists:

            self.load(key)

        if key in self.lists:

            return self.lists[key]


    def load(self, key):

        cachefile = 'reflist_%s_%u.pickle' % key
        cachefile = os.path.join(self.cachedir, cachefile)

        if os.path.exists(cachefile):

            self.lists[key] = pickle.load(open(cachefile, 'rb'))

            self._log(
                'Reference list for ID type `%s` for organism `%u` '
                'has been loaded from `%s`.' % (key + (cachefile,))
            )

        else:

            self.lists[key] = self._load(key)
            pickle.dump(self.lists[key], open(cachefile, 'wb'))
            self._log(
                'Reference list for ID type `%s` for organism `%u` '
                'has been saved to `%s`.' % (key + (cachefile,))
            )


    def _load(self, key):

        data = set()
        input_method = inputs[key[0]]

        if os.path.exists(input_method):

            with open(input_method, 'r') as fp:

                data = {l.strip() for l in fp.readlines()}

            self._log(
                'Reference list for ID type `%s` for organism `%u` has '
                'been loaded from `%s`.' % (key + (input_method,))
            )

        else:

            if hasattr(uniprot_input, input_method):

                input_func = getattr(uniprot_input, input_method)

            elif hasattr(mirbase_input, input_method):

                input_func = getattr(mirbase_input, input_method)

            ncbi_tax_id = key[1]
            data = set(input_func(organism = ncbi_tax_id))
            self._log(
                'Reference list for ID type `%s` for organism `%u` has '
                'been loaded by method `%s`.' % (key + (str(input_method),))
            )

        return data


    def check(self, name, id_type, ncbi_tax_id = None):
        """
        Checks if the identifier ``name`` is in the reference list with
        the provided ``id_type`` and organism.
        """

        lst = self.which_list(id_type = id_type, ncbi_tax_id = ncbi_tax_id)

        return name in lst


    def select(self, names, id_type, ncbi_tax_id = None):
        """
        Selects the identifiers in ``names`` which are in the reference list
        with the provided ``id_type`` and organism.
        """

        names = set(names)

        lst = self.which_list(id_type = id_type, ncbi_tax_id = ncbi_tax_id)

        return names & lst


    def is_not(self, names, id_type, ncbi_tax_id = None):
        """
        Returns the identifiers from ``names`` which are not instances of
        the provided ``id_type`` and from the given organism.
        """

        names = set(names)

        lst = self.which_list(id_type = id_type, ncbi_tax_id = ncbi_tax_id)

        return names - lst


    def _remove_expired(self):

        for key, last_used in list(self.expiry.items()):

            if time.time() - last_used > self.lifetime and key in self.lists:

                del self.lists[key]
                del self.expiry[key]


    def __del__(self):

        if hasattr(_reflists_cleanup_timeloop, 'stop'):

            _reflists_cleanup_timeloop.stop()


def init():

    globals()['manager'] = ReferenceListManager()


def get_manager():

    if 'manager' not in globals():

        init()

    return globals()['manager']


def check(name, id_type, ncbi_tax_id = None):
    """
    Checks if the identifier ``name`` is in the reference list with
    the provided ``id_type`` and organism.
    """

    manager = get_manager()

    return manager.check(
        name = name,
        id_type = id_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def select(names, id_type, ncbi_tax_id = None):
    """
    Selects the identifiers in ``names`` which are in the reference list
    with the provided ``id_type`` and organism.
    """

    manager = get_manager()

    return manager.select(
        names = names,
        id_type = id_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def is_not(names, id_type, ncbi_tax_id = None):
    """
    Returns the identifiers from ``names`` which are not instances of
    the provided ``id_type`` and from the given organism.
    """

    manager = get_manager()

    return manager.is_not(
        names = names,
        id_type = id_type,
        ncbi_tax_id = ncbi_tax_id,
    )


def get_reflist(id_type, ncbi_tax_id = None):

    manager = get_manager()

    return manager.which_list(id_type = id_type, ncbi_tax_id = ncbi_tax_id)
