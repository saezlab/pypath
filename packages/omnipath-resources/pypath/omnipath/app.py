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

from __future__ import annotations

from typing import Literal

import os
import sys
import shutil
import importlib as imp
import time
import pprint
import copy
import collections
import itertools

import pypath.resources.network as netres
from pypath.core import annot
from pypath.core import intercell
from pypath.core import complex
from pypath.core import enz_sub
from pypath.core import network
from pypath.share import session as session_mod

import pypath.share.settings as settings
import pypath.share.common as common


class DatabaseManager(session_mod.Logger):
    """
    Builds and serves the databases in OmniPath such as various networks,
    enzyme-substrate interactions, protein complexes, annotations and
    inter-cellular communication roles. Saves the databases to and loads
    them from pickle dumps on demand.
    """


    def __init__(self, rebuild = False, **kwargs):

        session_mod.Logger.__init__(self, name = 'omnipath.dbmanager')

        self.timestamp = time.strftime(settings.get('timestamp_format'))
        self.param = kwargs
        self.rebuild = rebuild
        self.datasets = self.get_param('datasets')
        self.ensure_dirs()
        self.network_dfs = {}

        self._log('The OmniPath database manager has been initialized.')


    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
        self.foreach_dataset(method = self.reload_module)


    def build(self):
        """
        Builds all built-in datasets.
        """

        self._log(
            'Building databases. Rebuild forced: %s.' % str(self.rebuild)
        )

        self.foreach_dataset(method = self.ensure_dataset)


    def ensure_dataset(
            self,
            dataset,
            force_reload = False,
            force_rebuild = False,
            ncbi_tax_id = 9606,
        ):
        """
        Makes sure a dataset is loaded. It loads only if it's not loaded
        yet or :py:arg:`force_reload` is ``True``. It only builds if it's
        not availabe as a pickle dump or :py:arg:`force_rebuild` is ``True``.

        :arg str dataset:
            The name of the dataset.
        :arg int ncbi_tax_id:
            NCBI Taxonomy ID of the organism. Considered only if the dataset
            builds for one organism and saved to organism specific pickle
            files.
        """

        for dep_dataset in self.dataset_dependencies(dataset):

            self.ensure_dataset(dep_dataset)

        rebuild_dataset = self.get_param('rebuild_%s' % dataset)

        _dataset = self._dataset_taxid(dataset, ncbi_tax_id = ncbi_tax_id)

        if (
            force_reload or
            force_rebuild or
            not hasattr(self, _dataset)
        ):

            if (
                force_rebuild or
                self.rebuild or
                rebuild_dataset or
                not self.pickle_exists(dataset, ncbi_tax_id = ncbi_tax_id)
            ):

                self.remove_db(dataset, ncbi_tax_id = ncbi_tax_id)
                self.build_dataset(dataset, ncbi_tax_id = ncbi_tax_id)

            elif (
                not hasattr(self, _dataset) or
                force_reload
            ):

                self.load_dataset(dataset, ncbi_tax_id = ncbi_tax_id)


    def dataset_dependencies(self, dataset):
        """
        Returns the dependencies of a dataset. E.g. to build `annotations`
        `complexes` must be loaded hence the former is dependent on the
        latter.
        """

        deps = self.get_param('dependencies')

        return deps[dataset] if dataset in deps else ()


    def ensure_dirs(self):
        """
        Checks if the directories for tables, figures and pickles exist and
        creates them if necessary.
        """

        if self.get_param('timestamp_dirs'):

            self.tables_dir = os.path.join(
                self.get_param('tables_dir'),
                self.timestamp
            )
            self.figures_dir = os.path.join(
                self.get_param('figures_dir'),
                self.timestamp,
            )
            settings.setup(
                tables_dir = self.tables_dir,
                figures_dir = self.figures_dir,
            )

        os.makedirs(self.get_param('pickle_dir'), exist_ok = True)

        for _dir in ('pickle', 'tables', 'figures'):

            path = self.get_param('%s_dir' % _dir)
            self._log(
                '%s directory: `%s` (exists: %s).' % (
                    _dir.capitalize(),
                    path,
                    'yes' if os.path.exists(path) else 'no',
                )
            )


    def pickle_path(self, dataset, ncbi_tax_id = 9606):
        """
        Returns the path of the pickle dump for a dataset according to
        the current settings.
        """

        pickle_fname = (
            self.get_param('%s_pickle' % dataset) or
            '%s.pickle' % dataset
        )

        if dataset == 'enz_sub':

            pickle_fname = pickle_fname % ncbi_tax_id

        return os.path.join(
            self.get_param('pickle_dir'),
            pickle_fname,
        )


    def pickle_exists(self, dataset, ncbi_tax_id = 9606):
        """
        Tells if a pickle dump of a particular dataset exists.
        """

        return os.path.exists(
            self.pickle_path(dataset, ncbi_tax_id = ncbi_tax_id)
        )


    def table_path(self, dataset):
        """
        Returns the full path for a table (to be exported or imported).
        """

        return os.path.join(
            self.get_param('tables_dir'),
            self.get_param('%s_tsv' % dataset),
        )


    def build_dataset(self, dataset, ncbi_tax_id = 9606):
        """
        Builds a dataset.
        """

        self._log('Building dataset `%s`.' % dataset)

        args = self.get_build_args(dataset)

        self._log('Build param: [%s].' % common.dict_str(args))

        mod = self.ensure_module(dataset)

        if dataset == 'enz_sub':

            args['ncbi_tax_id'] = ncbi_tax_id

        if hasattr(mod, 'db'):

            delattr(mod, 'db')

        db = mod.get_db(**args)

        pickle_path = self.pickle_path(dataset, ncbi_tax_id = ncbi_tax_id)

        old_pickle_path = '%s.old' % pickle_path

        if os.path.exists(pickle_path):

            shutil.move(pickle_path, old_pickle_path)

        self._log('Saving dataset `%s` to `%s`.' % (dataset, pickle_path))

        try:
            db.save_to_pickle(pickle_file = pickle_path)

            if os.path.exists(old_pickle_path):

                os.remove(old_pickle_path)

            self._log(
                'Saved dataset `%s` to `%s`.' % (
                    dataset,
                    pickle_path
                )
            )

        except Exception as e:

            exc = sys.exc_info()
            self._log_traceback()
            os.remove(pickle_path)

            self._log(
                'Failed to save dataset `%s` to `%s`. '
                'The dataset is currently loaded. '
                'Try restart Python and re-build the dataset. '
                'If the issue persists please report it.' % (
                    dataset,
                    pickle_path,
                )
            )

            if os.path.exists(old_pickle_path):

                self._log('Restoring the old version of `%s`.' % pickle_path)
                shutil.move(old_pickle_path, pickle_path)

        self._log('Successfully built dataset `%s`.' % dataset)

        _dataset = self._dataset_taxid(dataset, ncbi_tax_id = ncbi_tax_id)

        setattr(self, _dataset, db)

        self._add_network_df(dataset, ncbi_tax_id = ncbi_tax_id)


    def ensure_module(self, dataset, reset = True):
        """
        Makes sure the module providing a particular dataset is available
        and has no default database loaded yet (:py:attr:`db` attribute
        of the module).
        """

        mod_str = self.get_param('%s_mod' % dataset)
        mod = sys.modules['pypath.core.%s' % mod_str]

        if reset and hasattr(mod, 'db'):

            delattr(mod, 'db')

        return mod


    def reload_module(self, dataset):
        """
        Reloads the module of the database object of a particular dataset.
        E.g. in case of network datasets the ``pypath.network`` module
        will be reloaded.
        """

        mod = self.ensure_module(dataset, reset = False)
        imp.reload(mod)

        if hasattr(mod, 'db'):

            mod.db.reload()


    def get_build_args(self, dataset):
        """
        Retrieves the default database build parameters for a dataset.
        """

        args = self.get_param('%s_args' % dataset) or {}

        if hasattr(self, 'get_args_%s' % dataset):

            args.update(getattr(self, 'get_args_%s' % dataset)())

        return args


    def load_dataset(self, dataset, ncbi_tax_id = 9606):
        """
        Loads a dataset, builds it if no pickle dump is available.
        """

        pickle_path = self.pickle_path(dataset, ncbi_tax_id = ncbi_tax_id)

        self._log('Loading dataset `%s` from `%s`.' % (dataset, pickle_path))

        mod = self.ensure_module(dataset)

        _dataset = self._dataset_taxid(dataset, ncbi_tax_id = ncbi_tax_id)

        setattr(self, _dataset, mod.get_db(pickle_file = pickle_path))

        self._log('Loaded dataset `%s` from `%s`.' % (dataset, pickle_path))

        self._add_network_df(dataset, ncbi_tax_id = ncbi_tax_id)


    def _dataset_taxid(self, dataset, ncbi_tax_id = 9606):

        return '%s%s' % (
            dataset,
            '_%u' % ncbi_tax_id if dataset == 'enz_sub' else '',
        )


    # TODO
    # the get_args_* methods below will be replaced by the
    # pypath.omnipath.databases module

    def get_args_curated(self):
        """
        Returns the arguments for building the curated PPI network dataset.
        """

        resources = copy.deepcopy(netres.pathway)
        resources.update(copy.deepcopy(netres.enzyme_substrate))

        return {'resources': resources}


    def get_args_tf_target(self):
        """
        Returns the arguments for building the TF-target network dataset.
        """

        transcription = (
            netres.dorothea_expand_levels(
                resources = netres.transcription,
                levels = self.get_param('tfregulons_levels'),
            )
                if self.get_param('dorothea_expand_levels') else
            netres.transcription
        )

        return {'resources': transcription}


    def get_args_tf_mirna(self):
        """
        Returns the arguments for building the TF-miRNA network dataset.
        """

        return {'resources': netres.tf_mirna}


    def get_args_mirna_mrna(self):
        """
        Returns the arguments for building the miRNA-mRNA network dataset.
        """

        return {'resources': netres.mirnatarget}


    def get_args_lncrna_mrna(self):
        """
        Returns the arguments for building the lncRNA-mRNA network dataset.
        """

        return {'resources': netres.lncrna_mrna}


    def get_args_small_molecule(self):
        """
        Returns the arguments for building the small molecule-protein
        network dataset.
        """

        return {'resources': netres.small_molecule}


    def compile_tables(self):
        """
        Compiles the `summaries` table for all datasets. These tables contain
        various quantitative descriptions of the data contents.
        """

        self.foreach_dataset(method = self.compile_table)


    def compile_table(self, dataset):
        """
        Compiles the `summaries` table for a dataset. These tables contain
        various quantitative descriptions of the data contents.
        """

        table_path = self.table_path(dataset)
        db = self.get_db(dataset)
        db.update_summaries()
        db.summaries_tab(outfile = table_path)


    def foreach_dataset(self, method):
        """
        Applies a method for each dataset.
        """

        for dataset in self.datasets:

            _ = method(dataset)


    def get_db(self, dataset, ncbi_tax_id = 9606):
        """
        Returns a dataset object. Loads and builds the dataset if necessary.

        :arg int ncbi_tax_id:
            NCBI Taxonomy ID of the organism. Considered only if the dataset
            builds for one organism and saved to organism specific pickle
            files.
        """

        self.ensure_dataset(dataset, ncbi_tax_id = ncbi_tax_id)

        _dataset = self._dataset_taxid(dataset, ncbi_tax_id = ncbi_tax_id)

        return getattr(self, _dataset)


    def remove_db(self, dataset, ncbi_tax_id = 9606):
        """
        Removes a dataset. Deletes the references to the object
        in the module, however if you have references elsewhere in your
        code it remains in the memory.
        """

        _dataset = self._dataset_taxid(dataset, ncbi_tax_id = ncbi_tax_id)

        if hasattr(self, _dataset):

            delattr(self, _dataset)


    def remove_all(self):
        """
        Removes all loaded datasets. Deletes the references to the objects
        in the module, however if you have references elsewhere in your
        code they remain in the memory.
        """

        self.foreach_dataset(method = self.ensure_module)
        self.foreach_dataset(method = self.remove_db)


    def get_param(self, key):
        """
        Retrieves a parameter from the :py:attr:`param` dict of the current
        object or from the module settings.
        """

        return self.param.get(key, settings.get(key))


    def _create_network_df(self, dataset = 'omnipath', **kwargs):

        graph = self.get_db(dataset)

        return self._network_df(graph, **kwargs)


    def network_df(self, dataset, by_source = False):
        """
        Creates a data frame of a network dataset where rows aggregate
        information from all resources describing an interaction.
        """

        self.ensure_dataset(dataset)

        by_source_str = 'by_source' if by_source else 'plain'

        return self.network_dfs[dataset][by_source_str]


    def network_df_by_source(self, dataset = 'omnipath'):
        """
        Creates a data frame of a network dataset where each row contains
        information from one resource.
        """

        self.ensure_dataset(dataset)

        return self.network_dfs[dataset]['by_source']


    def _network_df(self, obj, **kwargs):

        if not isinstance(obj, network.Network):

            obj = network.Network.from_igraph(obj)

        obj.make_df(**kwargs)

        return obj.df


    def _add_network_df(self, dataset, ncbi_tax_id = 9606):

        _dataset = self._dataset_taxid(dataset, ncbi_tax_id = ncbi_tax_id)

        obj = getattr(self, _dataset)

        if (
            (
                hasattr(obj, 'graph') and
                hasattr(obj.graph, 'es')
            ) or
            isinstance(obj, network.Network)
        ):

            network_df = self._network_df(obj, by_source = False)
            network_df_by_source = self._network_df(obj, by_source = True)

            self.network_dfs[dataset] = {}
            self.network_dfs[dataset]['plain'] = network_df
            self.network_dfs[dataset]['by_source'] = network_df_by_source

            self._log('Created network data frames for `%s`.' % dataset)


    def set_network(self, dataset, by_source = False, **kwargs):
        """
        Sets dataset as the default
        """

        network_df = self.network_df(dataset, by_source = by_source, **kwargs)

        self.ensure_dataset('intercell')

        self.intercell.register_network(network_df)


    def define_dataset(
            self,
            name: str,
            module: Literal[
                'annot',
                'complex',
                'enz_sub',
                'intercell',
                'network',
            ],
            args: dict | None = None,
            pickle: str | None = None,
            **param,
        ):
        """
        Add a new dataset definition.

        Args
            name:
                Arbitrary name for the dataset.
            module:
                A database builder module: this determines the type of the
                dataset.
            args:
                Arguments for the database provider method (typically
                called ``get_db``) of the above module.
            pickle:
                A name for the pickle file, if not provided it will be
                named as "<name>_<module>.pickle".
            param:
                Further parameters, saved directly into the :attr:``param``
                dict of this object, however the three arguments above
                override values provided this way.
        """

        settings.setup(datasets = setting.get('datasets') + [name])

        param[f'{name}_pickle'] = pickle or f'{name}_{module}.pickle'
        param[f'{name}_mod'] = module
        param[f'{name}_args'] = args

        self.param.update(param)
