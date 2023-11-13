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

import sys
import importlib as imp
import traceback
import collections

try:
    import cPickle as pickle
except:
    import pickle

import numpy as np
import pandas as pd

import pypath.inputs as inputs
import pypath.internals.intera as intera
import pypath.internals.resource as resource
import pypath.share.settings as settings
import pypath.share.session as session_mod
import pypath.share.common as common


complex_resources = (
    'Signor',
    'Corum',
    'CellPhoneDB',
    'Havugimana',
    'Compleat',
    'ComplexPortal',
    'Pdb',
    'GuideToPharmacology',
    'Humap',
    'Humap2',
    'Icellnet',
    'Kegg',
    'Cellchatdb',
    'Cellinker',
    'Spike',
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

        self._log(
            'Loaded resource `%s`: %u proteins, %u complexes.' % (
                self.name,
                len(self.proteins),
                len(self.complexes),
            )
        )


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
        if isinstance(other, str):

            if len(other) <= 10:

                return other in self.proteins

            else:

                return other in self.complexes

        return False


    def __len__(self):

        return len(self.complexes)


    def __repr__(self):

        return '<Complex database: %u complexes>' % len(self)


    @property
    def numof_references(self):

        return len(
            set.union(*(
                cplex.references for cplex in self.complexes.values()
            ))
        )


    @property
    def curation_effort(self):

        return len(
            set.union(*(
                {(key, ref) for ref in cplex.references}
                for key, cplex in iteritems(self.complexes)
            ))
        )


    @property
    def has_stoichiometry(self):

        return any(
            cnt > 1
            for cplex in self.complexes.values()
            for cnt in cplex.components.values()
        )


    @property
    def all_sources(self):

        return set.union(*(
            cplex.sources
            for cplex in self.complexes.values()
        ))


    @property
    def homomers(self):

        return sum(
            1 for cplex in self.complexes.values()
            if len(cplex.components) == 1
        )


    @property
    def heteromers(self):

        return sum(
            1 for cplex in self.complexes.values()
            if len(cplex.components) > 1
        )


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

        self._log('Creating a data frame of complexes.')

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

        self._log(
            'Created data frame of complexes. '
            'Memory usage: %s.' % common.df_memory_usage(self.df)
        )


    def _from_dump_callback(self):

        if hasattr(self, '_from_dump'):

            self.complexes = self._from_dump
            delattr(self, '_from_dump')
            delattr(self, 'dump')


    @property
    def summary(self):

        return {
            'n_complexes': self.__len__(),
            'n_references': self.numof_references,
            'curation_effort': self.curation_effort,
            'has_stoichiometry': self.has_stoichiometry,
            'name': self.name,
            'sources': self.all_sources,
            'homomers': self.homomers,
            'heteromers': self.heteromers,
        }


    @property
    def summary_str(self):

        s = self.summary
        bar = '=' * 70

        return (
            '\n%s\n'
            'Complex resource `%s`\n'
            '%s\n'
            '\tNumber of complexes: %u\n'
            '\tHomomers: %u\n'
            '\tHeteromers: %u\n'
            '\tNumber of literature references: %u\n'
            '\tCuration effort (reference-entity pairs): %u\n'
            '\tHas stoichiometry: %s\n'
            '\tSources: %s\n'
            '%s\n\n'
        ) % (
            bar,
            self.name,
            bar,
            s['n_complexes'],
            s['homomers'],
            s['heteromers'],
            s['n_references'],
            s['curation_effort'],
            str(s['has_stoichiometry']),
            ', '.join(s['sources']),
            bar
        )


class CellPhoneDB(AbstractComplexResource):


    def __init__(self, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'CellPhoneDB',
            input_method = 'cellphonedb.cellphonedb_complexes',
        )


class Corum(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'CORUM',
            input_method = 'corum.corum_complexes',
            input_args = input_args or {},
        )


class Havugimana(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'Havugimana2012',
            input_method = 'havugimana.havugimana_complexes',
            input_args = input_args or {},
        )


class Compleat(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'Compleat',
            input_method = 'compleat.compleat_complexes',
            input_args = input_args or {},
        )


class ComplexPortal(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'ComplexPortal',
            input_method = 'complexportal.complexportal_complexes',
            input_args = input_args or {},
        )


class Kegg(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'KEGG',
            input_method = 'kegg.kegg_medicus_complexes',
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
            input_method = 'pdb.pdb_complexes',
            input_args = input_args or {},
        )


class Signor(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        if 'organism' not in input_args:

            input_args['organism'] = settings.get('default_organism')

        AbstractComplexResource.__init__(
            self,
            name = 'SIGNOR',
            input_method = 'signor.signor_complexes',
            input_args = input_args or {},
        )


class Hpmr(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        AbstractComplexResource.__init__(
            self,
            name = 'HPMR',
            input_method = 'hpmr.hpmr_complexes',
            input_args = input_args or {},
        )


class Humap(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'hu.MAP',
            input_method = 'humap.humap_complexes',
        )


class Humap2(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        AbstractComplexResource.__init__(
            self,
            name = 'hu.MAP2',
            input_method = 'humap.humap2_complexes',
        )


class GuideToPharmacology(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        AbstractComplexResource.__init__(
            self,
            name = 'Guide2Pharma',
            input_method = 'guide2pharma.guide2pharma_complexes',
            input_args = input_args or {},
        )


class Icellnet(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        AbstractComplexResource.__init__(
            self,
            name = 'ICELLNET',
            input_method = 'icellnet.icellnet_complexes',
            input_args = input_args or {},
        )


class Spike(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        AbstractComplexResource.__init__(
            self,
            name = 'SPIKE',
            input_method = 'spike.spike_complexes',
            input_args = input_args or {},
        )


class Cellchatdb(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        if 'organism' not in input_args:

            input_args['organism'] = settings.get('default_organism')

        AbstractComplexResource.__init__(
            self,
            name = 'CellChatDB',
            input_method = 'cellchatdb.cellchatdb_complexes',
            input_args = input_args or {},
        )


class Cellinker(AbstractComplexResource):


    def __init__(self, input_args = None, **kwargs):

        input_args = input_args or {}

        if 'organism' not in input_args:

            input_args['organism'] = settings.get('default_organism')

        AbstractComplexResource.__init__(
            self,
            name = 'Cellinker',
            input_method = 'cellinker.cellinker_complexes',
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

            self._log('Loading database from pickle `%s`.' % self.pickle_file)
            self.load_from_pickle(self.pickle_file)
            return

        self.data = {}
        self.summaries = {}

        for res in self.resources:

            total_attempts = settings.get('complex_load_resource_attempts')

            for attempt in range(total_attempts):

                try:

                    self._log(
                        f'Loading resource `{str(res)}`; '
                        f'attempt {attempt + 1}/{total_attempts}.'
                    )

                    if not callable(res):

                        if res in globals():

                            res = globals()[res]

                    if callable(res):

                        processor = res()

                    elif hasattr(res, 'complexes'):

                        processor = res

                    if hasattr(processor, 'summary'):

                        self.summaries[processor.name] = processor.summary

                    for key, cplex in iteritems(processor.complexes):

                        if key in self.data:

                            self.data[key] += cplex

                        else:

                            self.data[key] = cplex

                    self._log(f'Successfully loaded resource `{str(res)}`.')
                    break

                except Exception:

                    exc = sys.exc_info()
                    self._log('Failed to load resource `%s`:' % str(res))
                    self._log_traceback()

        resource.AbstractResource.load(self)
        self.update_index()
        self.update_summaries()


    def load_from_pickle(self, pickle_file):

        self._log('Loading from pickle `%s`.' % pickle_file)

        with open(pickle_file, 'rb') as fp:

            self.complexes, self.summaries = pickle.load(fp)

        self._log('Loaded from pickle `%s`.' % pickle_file)


    def update_summaries(self):

        for src in self.summaries.keys():

            self.summaries[src]['unique_complexes'] = sum(
                1 for cplex in self.complexes.values()
                if len(cplex.sources) == 1 and src in cplex.sources
            )

            self.summaries[src]['shared_complexes'] = sum(
                1 for cplex in self.complexes.values()
                if len(cplex.sources) > 1 and src in cplex.sources
            )


    def summaries_tab(self, outfile = None, return_table = False):

        columns = (
            ('name', 'Resource'),
            ('n_complexes', 'All complexes'),
            ('homomers', 'Homomers'),
            ('heteromers', 'Heteromers'),
            ('has_stoichiometry', 'Stoichiometry'),
            ('unique_complexes', 'Unique complexes'),
            ('shared_complexes', 'Shared complexes'),
            ('n_references', 'References'),
            ('curation_effort', 'Curation effort'),
        )

        tab = []
        tab.append([f[1] for f in columns])

        tab.extend([
            [
                str(self.summaries[src][f[0]])
                for f in columns
            ]
            for src in sorted(self.summaries.keys())
        ])

        if outfile:

            with open(outfile, 'w') as fp:

                fp.write('\n'.join('\t'.join(row) for row in tab))

        if return_table:

            return tab


    def _update_complex_attribute_classes(self):

        self._update_complex_attribute_classes_static(self.complexes)


    @staticmethod
    def _update_complex_attribute_classes_static(cplexes, mod = None):

        mod = mod or sys.modules[__name__]

        for key in cplexes:

            if hasattr(key, 'attrs'):

                for attr, val in iteritems(key.attrs):

                    cls = val.__class__.__name__

                    if hasattr(mod, cls):

                        val.__class__ = getattr(mod, cls)


    def save_to_pickle(self, pickle_file):

        self._log('Saving to pickle `%s`.' % pickle_file)

        self._update_complex_attribute_classes()

        with open(pickle_file, 'wb') as fp:

            pickle.dump(
                obj = (self.complexes, self.summaries),
                file = fp,
            )

        self._log('Saved to pickle `%s`.' % pickle_file)


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
    reference set for many methods, just like
    ``inputs.uniprot_db.all_uniprots`` represents the proteome.
    """

    db = get_db()

    return set(db.complexes.values())
