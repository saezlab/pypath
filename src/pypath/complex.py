#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
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

import sys
import imp
import traceback
import collections

try:
    import cPickle as pickle
except:
    import pickle

import numpy as np
import pandas as pd

import pypath.dataio as dataio
import pypath.intera as intera
import pypath.resource as resource
import pypath.settings as settings
import pypath.session_mod as session_mod


complex_resources = (
    'Signor',
    'Corum',
    'CellPhoneDB',
    'Havugimana',
    'Compleat',
    'ComplexPortal',
    'Pdb',
    'Hpmr',
    'GuideToPharmacology',
    'Humap',
)


class AbstractComplexResource(resource.AbstractResource):
    """
    A resource which provides information about molecular complexes.
    """


    def __init__(
            self,
            name,
            ncbi_tax_id = 9606,
            input_method = None,
            input_args = None,
            dump = None,
            **kwargs
        ):
        """
        name : str
            Custom name for the resource.
        input_method : callable
            Method providing the input data.
        process_method : callable
            Method processing the data and yielding ``intera.Complex``
            instances.
        """

        session_mod.Logger.__init__(self, name = 'complex')

        self.complexes = {}

        resource.AbstractResource.__init__(
            self,
            name = name,
            ncbi_tax_id = ncbi_tax_id,
            input_method = input_method,
            input_args = input_args,
            dump = dump,
            data_attr_name = 'complexes',
        )

        self.load()


    def load(self):

        resource.AbstractResource.load(self)
        self.update_index()


    def _process_method(self):

        self.complexes = self.data

        delattr(self, 'data')


    def __iter__(self):

        for cplex in self.complexes.values():

            yield cplex


    def update_index(self):

        self.proteins = collections.defaultdict(set)
        self.resources = collections.defaultdict(set)
        self.ids = {}

        for cplex in self:

            for protein in cplex:

                self.proteins[protein].add(cplex)

            for db in cplex.sources:

                self.resources[protein].add(cplex)

            for db, ids in iteritems(cplex.ids):

                for _id in ids:

                    self.ids[(db, _id)] = cplex


    def __contains__(self, other):

        # a Complex instance
        if isinstance(other, intera.Complex):

            other = other.__str__()

        # either a UniProt ID or
        # a complex string representation
        if isinstance(other, common.basestring):

            if len(other) <= 10:

                return other in self.proteins

            else:

                return other in self.complexes

        return False


    def __len__(self):

        return len(self.complexes)


    def make_df(self):

        colnames = [
            'name',
            'components',
            'components_genesymbols',
            'stoichiometry',
            'sources',
            'references',
            'identifiers',
        ]

        records = []

        for cplex in self.complexes.values():

            records.append([
                cplex.name if cplex.name else None,
                cplex.__str__()[8:],
                cplex.genesymbol_str,
                cplex.stoichiometry,
                ';'.join(cplex.sources),
                ';'.join(cplex.references),
                ';'.join(
                    '%s:%s' % (db, _id)
                    for db, ids in iteritems(cplex.ids)
                    for _id in ids
                ),
            ])

        self.df = pd.DataFrame(
            records,
            columns = colnames,
        )


    def _from_dump_callback(self):

        if hasattr(self, '_from_dump'):

            self.complexes = self._from_dump
            delattr(self, '_from_dump')
            delattr(self, 'dump')


class CellPhoneDB(AbstractComplexResource):


    def __init__(self, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'CellPhoneDB',
            input_method = 'cellphonedb_complexes',
        )


class Corum(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'CORUM',
            input_method = 'corum_complexes',
            input_args = input_args or {},
        )


class Havugimana(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'Havugimana2012',
            input_method = 'havugimana_complexes',
            input_args = input_args or {},
        )


class Compleat(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'Compleat',
            input_method = 'compleat_complexes',
            input_args = input_args or {},
        )


class ComplexPortal(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'ComplexPortal',
            input_method = 'complexportal_complexes',
            input_args = input_args or {},
        )


class Pdb(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        if 'organism' not in input_args:

            input_args['organism'] = settings.get('default_organism')

        AbstractComplexResource.__init__(
            self,
            name = 'PDB',
            input_method = 'pdb_complexes',
            input_args = input_args or {},
        )


class Signor(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        if 'organism' not in input_args:

            input_args['organism'] = settings.get('default_organism')

        AbstractComplexResource.__init__(
            self,
            name = 'Signor',
            input_method = 'signor_complexes',
            input_args = input_args or {},
        )


class Hpmr(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        AbstractComplexResource.__init__(
            self,
            name = 'HPMR',
            input_method = 'hpmr_complexes',
            input_args = input_args or {},
        )


class Humap(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'hu.MAP',
            input_method = 'humap_complexes',
        )


class GuideToPharmacology(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        AbstractComplexResource.__init__(
            self,
            name = 'Guide2Pharma',
            input_method = 'guide2pharma_complexes',
            input_args = input_args or {},
        )


class ComplexAggregator(AbstractComplexResource):


    def __init__(
            self,
            resources = None,
            pickle_file = None,
        ):
        """
        Combines complexes from multiple resources.

        :arg list resources:
            List of resources. Names of complex resource classes in this
            module or custom
        """

        self.pickle_file = pickle_file
        self.resources = resources or complex_resources

        AbstractComplexResource.__init__(
            self,
            name = 'OmniPath',
        )


    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def load(self):

        if self.pickle_file:

            self.load_from_pickle(self.pickle_file)
            return

        self.data = {}

        for res in self.resources:

            try:

                if not callable(res):

                    if res in globals():

                        res = globals()[res]

                if callable(res):

                    processor = res()

                elif hasattr(res, 'complexes'):

                    processor = res

                for key, cplex in iteritems(processor.complexes):

                    if key in self.data:

                        self.data[key] += cplex

                    else:

                        self.data[key] = cplex

            except Exception:

                self._log(
                    'Failed to load resouce `%s`: %s' % (
                        str(res),
                        traceback.format_exception(*sys.exc_info()),
                    ))

        resource.AbstractResource.load(self)
        self.update_index()


    def load_from_pickle(self, pickle_file):

        with open(pickle_file, 'rb') as fp:

            self.complexes = pickle.load(fp)


    def save_to_pickle(self, pickle_file):

        with open(pickle_file, 'wb') as fp:

            pickle.dump(
                obj = self.complexes,
                file = fp,
            )


def init_db(**kwargs):
    """
    Initializes or reloads the complex database.
    The database will be assigned to the ``db`` attribute of this module.
    """

    globals()['db'] = ComplexAggregator(**kwargs)


def get_db(**kwargs):
    """
    Retrieves the current database instance and initializes it if does
    not exist yet.
    """

    if 'db' not in globals():

        init_db(**kwargs)

    return globals()['db']


def all_complexes():
    """
    Returns a set of all complexes in the database which serves as a
    reference set for many methods, just like ``uniprot_input.all_uniprots``
    represents the proteome.
    """

    db = get_db()

    return set(db.complexes.values())
