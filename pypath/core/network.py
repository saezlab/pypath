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

from future.utils import iteritems

from typing import Mapping

import importlib as imp
import re
import os
import sys
import collections
import itertools
import functools
import copy as copy_mod
import pickle
import random
import traceback
from typing_extensions import Literal

import numpy as np
import pandas as pd

import pypath.share.session as session_mod
import pypath.share.progress as progress
import pypath.core.interaction as interaction_mod
import pypath.core.evidence as evidence
import pypath.core.entity as entity_mod
import pypath.core.common as core_common
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.share.settings as settings
import pypath.share.cache as cache_mod
import pypath.utils.mapping as mapping
import pypath.inputs.pubmed as pubmed_input
import pypath.share.curl as curl
import pypath.internals.refs as refs_mod
import pypath.utils.reflists as reflists
import pypath.resources.network as network_resources
import pypath.internals.input_formats as input_formats
import pypath.internals.resource as resource_formats
import pypath.inputs as inputs

# Py 2/3
try:
    input = raw_input
except NameError:
    pass

NetworkEntityCollection = collections.namedtuple(
    'NetworkEntityCollection',
    [
        'total',
        'by_resource',
        'by_category',
        'shared',
        'unique',
        'shared_res_cat',
        'unique_res_cat',
        'shared_cat',
        'unique_cat',
        'resource_cat',
        'cat_resource',
        'method',
        'label',
    ],
)
NetworkEntityCollection.__new__.__defaults__ = (None,) * 8


class NetworkEntityCollection(object):

    __slots__ = [
        'collection',
        '_collection',
        'label',

        'shared_within_data_model',
        'unique_within_data_model',
        'shared_within_interaction_type',
        'unique_within_interaction_type',

        'n_collection',
        'n_shared_within_data_model',
        'n_unique_within_data_model',
        'n_shared_within_interaction_type',
        'n_unique_within_interaction_type',

        'pct_collection',
        'pct_within_data_model',
        'pct_within_interaction_type',
        'pct_shared_within_data_model',
        'pct_unique_within_data_model',
        'pct_shared_within_interaction_type',
        'pct_unique_within_interaction_type',

        'by_data_model',
        'by_interaction_type',
        'unique_by_data_model',
        'shared_by_data_model',
        'unique_by_interaction_type',
        'shared_by_interaction_type',

        'n_by_data_model',
        'n_by_interaction_type',
        'n_unique_by_data_model',
        'n_shared_by_data_model',
        'n_unique_by_interaction_type',
        'n_shared_by_interaction_type',

        'pct_by_data_model',
        'pct_by_interaction_type',
        'pct_unique_by_data_model',
        'pct_shared_by_data_model',
        'pct_unique_by_interaction_type',
        'pct_shared_by_interaction_type',

    ]


    def __init__(self, collection, label = None):

        self.collection = collection.copy()
        # we need a copy where we don't add the totals
        # so these don't bother the shared and unique methods
        self._collection = collection.copy()
        self.label = label

        self.main()


    def main(self):

        self.setup()


    def setup(self):

        self.update()
        self.collection_add_total()
        self.update_collection_counts()


    def update_collection_counts(self):

        self.n_collection = common.dict_counts(self.collection)
        self.pct_collection = common.dict_set_percent(self.collection)


    def collection_add_total(self):

        self.collection = self._add_total(
            self.collection,
            key = ('all', 'all', 'Total')
        )


    def update(self):

        for level in ('interaction_type', 'data_model'):

            self._update(level = level)
            self._update(level = level, summarize_groups = True)


    def _update(self, level, summarize_groups = False):

        midpart = '_by_' if summarize_groups else '_within_'

        if summarize_groups:

            collection = common.dict_subtotals(
                self._expand_keys(level = level)
            )

            by = 'by_%s' % level

            setattr(
                self,
                by,
                collection
            )
            setattr(
                self,
                'n%s%s' % (midpart, level),
                common.dict_counts(collection)
            )

            for k, v in iteritems(getattr(self, by)):

                k = k if isinstance(k, tuple) else (k, 'all')

                k += ('Total',)

                self.collection[k] = v

        else:

            collection = self._expand_keys(level = level)

        setattr(
            self,
            'pct%s%s' % (midpart, level),
            (
                common.dict_set_percent(collection)
                    if summarize_groups else
                self._percent_and_collapse(collection)
            )
        )

        for method in ('shared', 'unique'):

            shared_unique = (
                self._add_total(
                    common.shared_unique_foreach(collection, op = method),
                    key = (
                        'all'
                            if level == 'interaction_type' else
                        ('all', 'all')
                    )
                )
                    if summarize_groups else
                self._shared_unique(
                    dct = collection,
                    method = method,
                    total_key = (
                        ('all', 'Total')
                            if level == 'interaction_type' else
                        None
                    ),
                )
            )

            if not summarize_groups:

                shared_unique_flat = common.dict_collapse_keys(shared_unique)

            attr = '%s%s%s' % (method, midpart, level)
            n_attr = 'n_%s' % attr
            pct_attr = 'pct_%s' % attr

            setattr(
                self,
                attr,
                shared_unique
            )
            setattr(
                self,
                n_attr,
                common.dict_collapse_keys(
                    common.dict_counts(shared_unique)
                )
            )
            setattr(
                self,
                pct_attr,
                common.dict_collapse_keys(
                    common.dict_set_percent(shared_unique)
                        if summarize_groups else
                    self._percent_and_collapse(shared_unique)
                )
            )


    def _expand_keys(self, level):

        return common.dict_expand_keys(
            self._collection,
            depth = 1,
            front = level == 'interaction_type',
        )


    @classmethod
    def _shared_unique(cls, dct, method, total_key = None):

        return dict(
            (
                key,
                cls._add_total(
                    common.shared_unique_foreach(val, op = method),
                    key = total_key
                )
            )
            for key, val in iteritems(dct)
        )


    @staticmethod
    def _add_total(dct, key = None):

        if isinstance(key, (str, tuple)):

            _key = key

        else:

            first_key = next(dct.keys().__iter__())

            if callable(key):

                _key = key(first_key)

            else:

                _key = (
                    'Total'
                        if isinstance(first_key, str) else
                    first_key[:-1] + ('Total',)
                )

        dct[_key] = common.dict_union(dct)

        return dct


    @classmethod
    def _percent_and_collapse(cls, dct):

        return (
            common.dict_collapse_keys(
                dict(
                    (
                        key,
                        common.dict_set_percent(val)
                    )
                    for key, val in iteritems(dct)
                )
            )
        )


NetworkStatsRecord = collections.namedtuple(
    'NetworkStatsRecord',
    [
        'total',
        'by_resource',
        'by_category',
        'shared',
        'unique',
        'percent',
        'shared_res_cat',
        'unique_res_cat',
        'percent_res_cat',
        'shared_cat',
        'unique_cat',
        'percent_cat',
        'resource_cat',
        'cat_resource',
        'method',
        'label',
    ],
)
NetworkStatsRecord.__new__.__defaults__ = (None,) * 11


class Network(session_mod.Logger):
    """
    Represents a molecular interaction network. Provides various methods to
    query the network and its components. Optionally converts the network
    to a ``pandas.DataFrame`` of interactions.

    :arg list,dict resources:
        One or more lists or dictionaries containing
        ``pypath.internals.resource.NetworkResource`` objects.
    :arg bool make_df:
        Create a ``pandas.DataFrame`` already when creating the instance.
        If no network data loaded no data frame will be created.
    :arg int ncbi_tax_id:
        Restrict the network only to this organism. If ``None`` identifiers
        from any organism will be allowed.
    :arg bool allow_loops:
        Allow interactions with the their two endpoints being the same entity.
    """

    _partners_methods = (
        {
            '': {},
            'transcriptionally_': {
                'interaction_type': {
                    'transcriptional',
                    'mirna_transcriptional',
                },
            },
            'post_transcriptionally_': {
                'interaction_type': {
                    'post_transcriptional',
                    'lncrna_post_transcriptional',
                },
            },
            'post_translationally_': {
                'interaction_type': 'post_translational',
            },
        },
        {
            'regulat': {
                'direction': True,
            },
            'activat': {
                'effect': 'positive',
            },
            'suppress': {
                'effect': 'negative',
            },
        },
        {
            'es': {
                'mode': 'IN',
            },
            'ed_by': {
                'mode': 'OUT',
            }
        },
    )


    def __init__(
            self,
            resources = None,
            make_df = False,
            df_by_source = False,
            df_with_references = False,
            df_columns = None,
            df_dtype = None,
            pickle_file = None,
            ncbi_tax_id = 9606,
            allow_loops = None,
            **kwargs
        ):

        session_mod.Logger.__init__(self, name = 'network')

        self._log('Creating network object.')

        self.df_by_source = df_by_source
        self.df_with_references = df_with_references
        self.df_columns = df_columns
        self.df_dtype = df_dtype
        self.ncbi_tax_id = ncbi_tax_id
        self.allow_loops = allow_loops
        self.cache_dir = cache_mod.get_cachedir()
        self.keep_original_names = settings.get('network_keep_original_names')
        self.default_name_types = settings.get('default_name_types')

        self.reset()

        if pickle_file and os.path.exists(pickle_file):

            self.load_from_pickle(pickle_file = pickle_file)
            return

        self.load(resources = resources, make_df = make_df, **kwargs)


    def reload(self, recursive: bool = False):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

        if recursive:

            imp.reload(entity_mod)
            imp.reload(interaction_mod)

            for entity in self.nodes.values():

                entity.__class__ = entity_mod.Entity

            for interaction in self.interactions.values():

                interaction.__class__ = interaction_mod.Interaction
                interaction.a.__class__ = entity_mod.Entity
                interaction.b.__class__ = entity_mod.Entity


    def __len__(self):

        return len(self.interactions)


    def __bool__(self):

        return bool(self.interactions)


    def __iter__(self):

        for ia in self.interactions.values():

            yield ia


    def __contains__(self, other):

        return any(other in ia for ia in self.interactions.values())


    def reset(self):
        """
        Removes network data i.e. creates empty interaction and node
        dictionaries.
        """

        self.raw_data = {}
        self.interactions = {}
        self.nodes = {}
        self.nodes_by_label = {}
        self.interactions_by_nodes = collections.defaultdict(set)


    def load(
            self,
            resources = None,
            make_df = False,
            exclude = None,
            reread = False,
            redownload = False,
            keep_raw = False,
            top_call = True,
            cache_files = None,
            only_directions = False,
            pickle_file = None,
            allow_loops = None,
            first_n = None,
        ):
        """
        Loads data from a network resource or a collection of resources.

        :arg str,dict,list,resource.NetworkResource resources:
            An object defining one or more network resources. If *str* it
            will be looked up among the collections in the
            ``pypath.resources.network`` module (e.g. ``'pathway'`` will load
            all resources in the `pathway` collection). If *dict* or *list*
            it will be processed recursively i.e. the ``load`` method will be
            called for each element. If it is a
            ``pypath.internals.resource.NetworkResource`` object it will be
            processed and added to the network.
        :arg bool make_df:
            Whether to create a ``pandas.DataFrame`` after loading all
            resources.
        :arg NoneType,set exclude:
            A *set* of resource names to be ignored. It is useful if you want
            to load a collection with the exception of a few resources.
        """

        if pickle_file:

            self.load_from_pickle(pickle_file = pickle_file)
            return

        kwargs = {
            'reread': reread,
            'redownload': redownload,
            'keep_raw': keep_raw,
            'top_call': False,
            'only_directions': only_directions,
            'allow_loops': allow_loops,
            'first_n': first_n,
        }

        exclude = common.to_set(exclude)

        resources = (
            (resources,)
                if not isinstance(resources, (list, Mapping, tuple, set)) else
            resources.values()
                if isinstance(resources, Mapping) else
            resources
        )

        for resource in resources:

            if (
                isinstance(resource, str) and
                hasattr(network_resources, resource)
            ):

                self.load(
                    resources = getattr(network_resources, resource),
                    **kwargs
                )

            elif isinstance(resource, (list, dict, tuple, set)):

                self.load(
                    resources = resource,
                    **kwargs
                )

            elif (
                isinstance(
                    resource,
                    (
                        input_formats.NetworkInput,
                        resource_formats.NetworkResource,
                    )
                ) and resource.name not in exclude
            ):

                self.load_resource(resource, **kwargs)

            elif resource is not None:

                self._log(
                    'Could not recognize network input '
                    'definition: `%s`.' % str(resource)
                )

        if make_df and top_call:

            self.make_df()


    # synonyms (old method names of PyPath)
    load_resources = load
    init_network = load


    def load_resource(
            self,
            resource,
            clean = True,
            reread = None,
            redownload = None,
            keep_raw = False,
            only_directions = False,
            allow_loops = None,
            first_n = None,
            **kwargs
        ):
        """
        Loads the data from a single resource and attaches it to the
        network

        :arg pypath.input_formats.NetworkInput resource:
            :py:class:`pypath.input_formats.NetworkInput` instance
            containing the detailed definition of the input format to
            the downloaded file.
        :arg bool clean:
            Legacy parameter, has no effect at the moment.
            Optional, ``True`` by default. Whether to clean the graph
            after importing the data or not. See
            :py:meth:`pypath.main.PyPath.clean_graph` for more
            information.
        :arg dict cache_files:
            Legacy parameter, has no effect at the moment.
            Optional, ``{}`` by default. Contains the resource name(s)
            [str] (keys) and the corresponding cached file name [str].
            If provided (and file exists) bypasses the download of the
            data for that resource and uses the cache file instead.
        :arg bool reread:
            Optional, ``False`` by default. Specifies whether to reread
            the data files from the cache or omit them (similar to
            *redownload*).
        :arg bool redownload:
            Optional, ``False`` by default. Specifies whether to
            re-download the data and ignore the cache.
        :arg bool only_directions:
            If ``True``, no new interactions will be created but direction
            and effect sign evidences will be added to existing interactions.
        :arg int first_n:
            Load only the first n interactions.
        """

        total_attempts = settings.get('network_load_resource_attempts')

        for attempt in range(total_attempts):

            try:

                self._log(
                    f'Loading network data from resource `{resource.name}`'
                    f' (dataset: {resource.dataset}); '
                    f'attempt {attempt + 1} of {total_attempts}.'
                )

                self._read_resource(
                    resource,
                    reread = reread,
                    redownload = redownload,
                    keep_raw = keep_raw,
                    first_n = first_n,
                )

                self._log(
                    'Successfully read interactions '
                    f'from resource `{resource.name}`.'
                )
                break

            except Exception as e:

                exc = sys.exc_info()
                self._log(
                    'Failed to read interactions '
                    f'from resource `{resource.name}`:'
                )
                self._log_traceback(console = True)

                if attempt == total_attempts - 1:

                    self._log(
                        f'Not loading `{resource.name}`: giving up after '
                        f'{total_attempts} attempts.'
                    )
                    return

        allow_loops = self._allow_loops(
            allow_loops = allow_loops,
            resource = resource,
        )

        self._log('Loops allowed for resource `%s`: %s' % (
            resource.name,
            allow_loops,
        ))

        self._add_edge_list(
            only_directions = only_directions,
            allow_loops = allow_loops,
        )

        self.organisms_check()
        self.remove_zero_degree()

        self._log(
            'Completed: loading network data from '
            'resource `%s`.' % resource.name
        )


    def _read_resource(
            self,
            resource,
            reread = False,
            redownload = False,
            keep_raw = False,
            cache_files = None,
            first_n = None,
        ):
        """
        Reads interaction data file containing node and edge attributes
        that can be read from simple text based files and adds it to the
        networkdata. This function works not only with files, but with
        lists as well. Any other function can be written to download and
        preprocess data, and then give it to this function to finally
        attach to the network.

        :arg pypath.input_formats.NetworkInput resource:
            :py:class:`pypath.input_formats.NetworkInput` instance
            containing the detailed definition of the input format of
            the file. Instead of the file name (on the
            :py:attr:`pypath.input_formats.NetworkInput.input`
            attribute) you can give a custom function name, which will
            be executed, and the returned data will be used instead.
        :arg bool keep_raw:
            Optional, ``False`` by default. Whether to keep the raw data
            read by this function, in order for debugging purposes, or
            further use.
        :arg dict cache_files:
            Optional, ``{}`` by default. Contains the resource name(s)
            [str] (keys) and the corresponding cached file name [str].
            If provided (and file exists) bypasses the download of the
            data for that resource and uses the cache file instead.
        :arg bool reread:
            Optional, ``False`` by default. Specifies whether to reread
            the data files from the cache or omit them (similar to
            *redownload*).
        :arg bool redownload:
            Optional, ``False`` by default. Specifies whether to
            re-download the data and ignore the cache.
        :arg int first_n:
            Load only the first n interactions.
        """

        self._log('Reading network data from `%s`.' % resource.name)

        SMOL_TYPES = settings.get('small_molecule_entity_types')

        # workaround in order to make it work with both NetworkInput
        # and NetworkResource type param
        _resource = (
            resource
                if isinstance(resource, resource_formats.NetworkResource) else
            resource_formats.NetworkResource(
                name = resource.name,
                interaction_type = resource.interaction_type,
                networkinput = resource,
                data_model = resource.data_model or 'unknown',
                resource_attrs = resource.resource_attrs,
            )
        )

        networkinput = _resource.networkinput

        _resources_secondary = ()

        expand_complexes = (
            networkinput.expand_complexes
                if isinstance(networkinput.expand_complexes, bool) else
            settings.get('network_expand_complexes')
        )
        reread = (
            reread
                if isinstance(reread, bool) else
            not settings.get('network_pickle_cache')
        )

        self._log('Expanding complexes for `%s`: %s' % (
            networkinput.name, str(expand_complexes),
        ))

        edge_list = []
        edge_list_mapped = []
        self.edge_list_mapped = []
        infile = None
        _name = networkinput.name.lower()

        edges_cache = os.path.join(
            self.cache_dir,
            '%s_%s_%s.edges.pickle' % (
                _name,
                _resource.data_model,
                _resource.interaction_type,
            )
        )

        interaction_cache = os.path.join(
            self.cache_dir,
            '%s_%s_%s.interactions.pickle' % (
                _name,
                _resource.data_model,
                _resource.interaction_type,
            )
        )

        if not reread and not redownload:

            infile, edge_list_mapped = self._lookup_cache(
                _name,
                cache_files,
                interaction_cache,
                edges_cache,
            )

        if not len(edge_list_mapped):

            if infile is None:

                if not isinstance(
                    resource,
                    (
                        input_formats.NetworkInput,
                        resource_formats.NetworkResource,
                    )
                ):

                    self._log(
                        '_read_network_data: No proper input file '
                        'definition. `param` should be either '
                        'a `pypath.internals.input_formats.NetworkInput` or a '
                        '`pypath.internals.resource.NetworkResource` instance.',
                        -5,
                    )

                    return None

                if networkinput.huge:

                    sys.stdout.write(
                        '\n\tProcessing %s requires huge memory.\n'
                        '\tPlease hit `y` if you have at '
                        'least 2G free memory,\n'
                        '\tor `n` to omit %s.\n'
                        '\tAfter processing once, it will be saved in \n'
                        '\t%s, so next time can be loaded quickly.\n\n'
                        '\tProcess %s now? [y/n]\n' %
                        (
                            networkinput.name,
                            networkinput.name,
                            edges_cache,
                            networkinput.name
                        )
                    )
                    sys.stdout.flush()

                    while True:
                        answer = input().lower()

                        if answer == 'n':
                            return None

                        elif answer == 'y':
                            break

                        else:
                            sys.stdout.write(
                                '\n\tPlease answer `y` or `n`:\n\t')
                            sys.stdout.flush()

                # if no method available it gonna be None
                input_func = inputs.get_method(networkinput.input)

                # reading from remote or local file, or executing import
                # function:
                if (
                    isinstance(networkinput.input, str) and (
                        networkinput.input.startswith('http') or
                        networkinput.input.startswith('ftp')
                    )
                ):

                    curl_use_cache = not redownload
                    c = curl.Curl(
                        networkinput.input,
                        silent = False,
                        large = True,
                        cache = curl_use_cache
                    )
                    infile = c.fileobj.read()

                    if type(infile) is bytes:

                        try:
                            infile = infile.decode('utf-8')

                        except UnicodeDecodeError as e:

                            try:
                                infile = infile.decode('iso-8859-1')

                            except UnicodeDecodeError:

                                raise e

                    infile = [
                        x for x in infile.replace('\r', '').split('\n')
                        if len(x) > 0
                    ]
                    self._log(
                        "Retrieving data from `%s` ..." % networkinput.input
                    )

                elif input_func is not None:

                    self._log(
                        'Retrieving data by method `%s` of the '
                        'pypath.inputs module...' % input_func.__name__
                    )

                    _store_cache = curl.CACHE

                    if isinstance(redownload, bool):

                        curl.CACHE = not redownload

                    try:

                        infile = input_func(**networkinput.input_args)

                    except Exception as e:

                        self._log(
                            f'Error in method `{input_func.__name__}` of the '
                            'pypath.inputs module. '
                        )

                        raise e

                    finally:

                        curl.CACHE = _store_cache

                elif os.path.isfile(networkinput.input):

                    infile = curl.Curl(
                        networkinput.input,
                        large = True,
                        silent = False,
                    ).result

                    self._log('%s opened...' % networkinput.input)

                if infile is None:

                    self._log(
                        '`%s`: Could not find file or input function '
                        'or failed preprocessing.' %
                        networkinput.input,
                        -5,
                    )
                    return None

            is_directed = networkinput.is_directed
            sign = networkinput.sign
            ref_col = (
                networkinput.refs[0]
                    if isinstance(networkinput.refs, tuple) else
                networkinput.refs
                    if isinstance(networkinput.refs, int) else
                None
            )
            ref_sep = (
                networkinput.refs[1]
                    if isinstance(networkinput.refs, tuple) else
                ';'
            )
            # column index of the sign
            sig_col = None if not isinstance(sign, tuple) else sign[0]
            # column index and value(s) for the direction
            dir_col = None
            dir_val = None
            dir_sep = None

            if isinstance(is_directed, tuple):

                dir_col = is_directed[0]
                dir_val = is_directed[1]
                dir_sep = is_directed[2] if len(is_directed) > 2 else None

            elif isinstance(sign, tuple):

                dir_col = sign[0]
                dir_val = sign[1:3]
                dir_val = (
                    dir_val
                        if type(dir_val[0]) in _const.SIMPLE_TYPES else
                    common.flat_list(dir_val)
                )
                dir_sep = sign[3] if len(sign) > 3 else None

            dir_val = common.to_set(dir_val)

            must_have_references = (
                settings.get('keep_noref') or
                networkinput.must_have_references
            )
            self._log(
                'Resource `%s` %s have literature references '
                'for all interactions. Interactions without references '
                'will be %s. You can alter this condition globally by '
                '`pypath.settings.keep_noref` or for individual resources '
                'by the `must_have_references` attribute of their '
                '`NetworkInput` object.' % (
                    networkinput.name,
                    'must' if must_have_references else 'does not need to',
                    'dropped' if must_have_references else 'included',
                ),
                1,
            )
            self._log(
                '`%s` must have references: %s' % (
                    networkinput.name,
                    str(must_have_references)
                )
            )

            # iterating lines from input file
            input_filtered = 0
            ref_filtered = 0
            taxon_filtered = 0
            read_error = False
            lnum = 0 # we need to define it here to avoid errors if the
                     # loop below runs zero cycles

            prg = progress.Progress(
                iterable = infile,
                name = 'Reading network data - %s' % networkinput.name,
            )

            try:

                for lnum, line in enumerate(prg):

                    if len(line) <= 1 or (lnum == 1 and networkinput.header):
                        # empty lines
                        # or header row
                        continue

                    if not isinstance(line, (list, tuple)):

                        if hasattr(line, 'decode'):
                            line = line.decode('utf-8')

                        line = line.strip('\n\r').split(networkinput.separator)

                    else:
                        line = [
                            x.replace('\n', '').replace('\r', '')
                                if hasattr(x, 'replace') else
                            x
                            for x in line
                        ]

                    # 1) filters
                    if self._filters(
                        line,
                        networkinput.positive_filters,
                        networkinput.negative_filters
                    ):

                        input_filtered += 1
                        continue

                    # 2) direction
                    # reading names and attributes:
                    if is_directed and not isinstance(is_directed, tuple):

                        this_edge_dir = True

                    else:

                        this_edge_dir = self._process_direction(
                            line,
                            dir_col,
                            dir_val,
                            dir_sep,
                        )

                    # 3) references
                    refs = []

                    if ref_col is not None:

                        if line[ref_col] is None:

                            refs = ()

                        elif isinstance(line[ref_col], (list, set, tuple)):

                            refs = line[ref_col]

                        elif isinstance(line[ref_col], int):

                            refs = (line[ref_col],)

                        else:

                            refs = line[ref_col].split(ref_sep)

                        refs = common.del_empty(list(set(refs)))

                    refs = pubmed_input.only_pmids(
                        [str(r).strip() for r in refs]
                    )

                    if len(refs) == 0 and must_have_references:

                        ref_filtered += 1
                        continue

                    # 4) entity types
                    entity_type_a = self._process_field(
                        networkinput.entity_type_a,
                        line,
                    )
                    entity_type_b = self._process_field(
                        networkinput.entity_type_b,
                        line,
                    )

                    # 5) ID types
                    id_type_a = self._process_field(networkinput.id_type_a, line)
                    id_type_b = self._process_field(networkinput.id_type_b, line)

                    # 6) organisms
                    # to give an easy way for input definition:
                    if isinstance(networkinput.ncbi_tax_id, int):

                        taxon_a = (
                            _const.NOT_ORGANISM_SPECIFIC
                                if entity_type_a in SMOL_TYPES else
                            networkinput.ncbi_tax_id
                        )
                        taxon_b = (
                            _const.NOT_ORGANISM_SPECIFIC
                                if entity_type_b in SMOL_TYPES else
                            networkinput.ncbi_tax_id
                        )

                    # to enable more sophisticated inputs:
                    elif isinstance(networkinput.ncbi_tax_id, dict):

                        taxx = self._process_taxon(
                            networkinput.ncbi_tax_id,
                            line,
                        )

                        if isinstance(taxx, tuple):

                            taxon_a, taxon_b = taxx

                        else:

                            taxon_a = taxon_b = taxx

                        taxd_a = (
                            networkinput.ncbi_tax_id['A']
                                if 'A' in networkinput.ncbi_tax_id else
                            _const.NOT_ORGANISM_SPECIFIC
                                if entity_type_a in SMOL_TYPES else
                            networkinput.ncbi_tax_id
                        )
                        taxd_b = (
                            networkinput.ncbi_tax_id['B']
                                if 'B' in networkinput.ncbi_tax_id else
                            _const.NOT_ORGANISM_SPECIFIC
                                if entity_type_b in SMOL_TYPES else
                            networkinput.ncbi_tax_id
                        )

                        only_default = networkinput.only_default_organism

                        if not (
                            self._match_taxon(taxd_a, taxon_a, only_default) and
                            self._match_taxon(taxd_b, taxon_b, only_default)
                        ):

                            taxon_filtered += 1
                            continue

                    # assuming by default the default organism
                    else:

                        taxon_a = taxon_b = self.ncbi_tax_id

                    if taxon_a is None or taxon_b is None:

                        taxon_filtered += 1
                        continue

                    # 7) effect (sign)
                    positive = False
                    negative = False

                    if isinstance(sign, tuple):

                        positive, negative = (
                            self._process_sign(line[sign[0]], sign)
                        )

                    # 8) resources (source databases)
                    resource = (
                        line[networkinput.resource]
                            if isinstance(networkinput.resource, int) else
                        line[networkinput.resource[0]].split(
                            networkinput.resource[1]
                        )
                            if (
                                isinstance(networkinput.resource, tuple) and
                                hasattr(line[networkinput.resource[0]], 'split')
                            ) else
                        []
                            if isinstance(networkinput.resource, tuple) else
                        networkinput.resource
                    )

                    resource = common.to_set(resource)

                    _resources_secondary = tuple(
                        resource_formats.NetworkResource(
                            name = sec_res,
                            interaction_type = _resource.interaction_type,
                            data_model = _resource.data_model,
                            via = _resource.name,
                            dataset = _resource.dataset,
                        )
                        for sec_res in resource
                        if sec_res != _resource.name
                    )

                    resource.add(networkinput.name)

                    # 9) interacting partners
                    id_a = self._process_partner(networkinput.id_col_a, line)
                    id_b = self._process_partner(networkinput.id_col_b, line)

                    # 10) further attributes
                    # getting additional edge and node attributes
                    attrs_edge = self._process_attrs(
                        line,
                        networkinput.extra_edge_attrs,
                        lnum,
                    )
                    attrs_node_a = self._process_attrs(
                        line,
                        networkinput.extra_node_attrs_a,
                        lnum,
                    )
                    attrs_node_b = self._process_attrs(
                        line,
                        networkinput.extra_node_attrs_b,
                        lnum,
                    )

                    # 11) creating the Evidence object
                    evidences = evidence.Evidences(
                        evidences = (
                            evidence.Evidence(
                                resource = _res,
                                references = None if _res.via else refs,
                                attrs = attrs_edge,
                            )
                            for _res in
                            _resources_secondary + (_resource,)
                        )
                    )

                    # 12) node attributes that
                    #     depend on the interaction direction
                    if networkinput.mark_source:

                        attrs_node_a[networkinput.mark_source] = this_edge_dir

                    if networkinput.mark_target:

                        attrs_node_b[networkinput.mark_target] = this_edge_dir

                    # 13) all interaction data goes into a dict
                    new_edge = {
                        'id_a': id_a,
                        'id_b': id_b,
                        'id_type_a': id_type_a,
                        'id_type_b': id_type_b,
                        'entity_type_a': entity_type_a,
                        'entity_type_b': entity_type_b,
                        'source': resource,
                        'is_directed': this_edge_dir,
                        'references': refs,
                        'positive': positive,
                        'negative': negative,
                        'taxon_a': taxon_a,
                        'taxon_b': taxon_b,
                        'interaction_type': networkinput.interaction_type,
                        'evidences': evidences,
                        'attrs_node_a': attrs_node_a,
                        'attrs_node_b': attrs_node_b,
                        'attrs_edge': attrs_edge,
                    }

                    if read_error:

                        self._log(
                            'Errors occured, certain lines skipped.'
                            'Trying to read the remaining.\n',
                            5,
                        )

                    edge_list.append(new_edge)

                    if first_n and len(edge_list) >= first_n:

                        break

            except Exception as e:

                self._log(
                    'Error at loading resource `%s`.' % networkinput.name
                )

                raise e

            if hasattr(infile, 'close'):

                infile.close()

            # 14) ID translation of edges
            edge_list_mapped = self._map_list(
                edge_list,
                expand_complexes = expand_complexes,
            )

            self._log(
                '%u lines have been read from %s, '
                '%u links after mapping; '
                '%u lines filtered by filters; '
                '%u lines filtered because lack of references; '
                '%u lines filtered by taxon filters.' %
                (
                    lnum - 1,
                    networkinput.input,
                    len(edge_list_mapped),
                    input_filtered,
                    ref_filtered,
                    taxon_filtered,
                )
            )

            if reread or redownload:

                pickle.dump(edge_list_mapped, open(edges_cache, 'wb'), -1)
                self._log('ID translated edge list saved to %s' % edges_cache)

        else:

            self._log(
                'Previously ID translated edge list '
                'has been loaded from `%s`.' % edges_cache
            )

        if keep_raw:

            self.raw_data[networkinput.name] = edge_list_mapped

        self.edge_list_mapped = edge_list_mapped


    def _lookup_cache(self, name, cache_files, int_cache, edges_cache):
        """
        Checks up the cache folder for the files of a given resource.
        First checks if *name* is on the *cache_files* dictionary.
        If so, loads either the interactions or edges otherwise. If
        not, checks *edges_cache* or *int_cache* otherwise.

        :arg str name:
            Name of the resource (lower-case).
        :arg dict cache_files:
            Contains the resource name(s) [str] (keys) and the
            corresponding cached file name [str] (values).
        :arg str int_cache:
            Path to the interactions cache file of the resource.
        :arg str edges_cache:
            Path to the edges cache file of the resource.

        :return:
            * (*file*) -- The loaded pickle file from the cache if the
              file is contains the interactions. ``None`` otherwise.
            * (*list*) -- List of mapped edges if the file contains the
              information from the edges. ``[]`` otherwise.
        """

        cache_files = cache_files or {}
        infile = None
        edge_list_mapped = []
        cache_file = cache_files[name] if name in cache_files else None

        if cache_file is not None and os.path.exists(cache_file):
            cache_type = cache_file.split('.')[-2]

            if cache_type == 'interactions':
                infile = self.read_from_cache(int_cache)

            elif cache_type == 'edges':
                edge_list_mapped = self.read_from_cache(edges_cache)

        elif os.path.exists(edges_cache):
            edge_list_mapped = self.read_from_cache(edges_cache)

        elif os.path.exists(int_cache):
            infile = self.read_from_cache(int_cache)

        return infile, edge_list_mapped


    @classmethod
    def _filters(
            cls,
            line,
            positive_filters = None,
            negative_filters = None,
        ):
        """
        Applies negative and positive filters on a line (record from an
        interaction database). If returns ``True`` the interaction will be
        discarded, if ``False`` the interaction will be further processed
        and if all other criteria fit then will be added to the network
        after identifier translation.

        Return
            (bool): True if the line should be filtered (removed), False
                if all filters passed, the record can be further processed.
        """

        return (
            cls._process_filters(line, negative_filters, False) or
            cls._process_filters(line, positive_filters, True)
        )


    @classmethod
    def _process_filters(cls, line, filters = None, negate = False):
        """
        Args
            negate (bool): Whether to negate the filter matches. Sorry for
                the confusion, but it should be True for positive filters
                and False for negatives.


        Return
            (bool): True if the line should be filtered (removed), False
                if all filters passed, the record can be further processed.
        """

        _negate = (lambda x: not x) if negate else (lambda x: x)

        filters = filters or ()

        for filtr in filters:

            if _negate(cls._process_filter(line, filtr)):

                return True

        return False


    @classmethod
    def _process_filter(cls, line, filtr):
        """
        Return
            (bool): True if the filter matches.
        """

        if callable(filtr):

            if filtr(line):

                return True

        else:

            if len(filtr) > 2:

                sep = filtr[2]
                thisVal = set(line[filtr[0]].split(sep))

            else:

                thisVal = common.to_set(line[filtr[0]])

            filtrVal = common.to_set(filtr[1])

            return bool(thisVal & filtrVal)


    def _process_sign(self, sign_data, sign_def):
        """
        Processes the sign of an interaction, used when processing an
        input file.

        :arg str sign_data:
            Data regarding the sign to be processed.
        :arg tuple sign_def:
            Contains information about how to process *sign_data*. This
            is defined in :py:mod:`pypath.data_formats`. First element
            determines the position on the direction information of each
            line on the data file [int], second element is either [str]
            or [list] and defines the terms for which an interaction is
            defined as stimulation, third element is similar but for the
            inhibition and third (optional) element determines the
            separator for *sign_data* if contains more than one element.

        :return:
            * (*bool*) -- Determines whether the processed interaction
              is considered stimulation (positive) or not.
            * (*bool*) -- Determines whether the processed interaction
              is considered inhibition (negative) or not.
        """

        positive = False
        negative = False
        sign_sep = sign_def[3] if len(sign_def) > 3 else None
        sign_data = sign_data.split(sign_sep) if sign_sep else sign_data
        sign_data = common.to_set(sign_data)
        pos = common.to_set(sign_def[1])
        neg = common.to_set(sign_def[2])

        if bool(sign_data & pos):

            positive = True

        if bool(sign_data & neg):

            negative = True

        return positive, negative


    def _process_direction(self, line, dir_col, dir_val, dir_sep):
        """
        Processes the direction information of an interaction according
        to a data file from a source.

        :arg list line:
            The stripped and separated line from the resource data file
            containing the information of an interaction.
        :arg int dir_col:
            The column/position number where the information about the
            direction is to be found (on *line*).
        :arg list dir_val:
            Contains the terms [str] for which that interaction is to be
            considered directed.
        :arg str dir_sep:
            Separator for the field in *line* containing the direction
            information (if any).

        :return:
            (*bool*) -- Determines whether the given interaction is
            directed or not.
        """

        if isinstance(dir_col, bool):

            return dic_col

        if (
            dir_val is None and
            isinstance(dir_col, int) and
            isinstance(line[dir_col], bool)
        ):

            return line[dir_col]

        if dir_col is None or dir_val is None:

            return False

        else:

            value = line[dir_col].split(dir_sep) if dir_sep else line[dir_col]
            value = common.to_set(value)

            return bool(value & dir_val)


    def _process_field(self, fmt, line):
        """
        Extract a value from a line describing an interaction.

        Args
            fmt (str, tuple, callable): The value, or a definition how to
                process it.
            line (list): The raw interaction record.

        Return
            (str): The extracted value.
        """

        if common.is_str(fmt) or isinstance(fmt, list):

            return fmt

        elif callable(fmt):

            return fmt(line)

        if isinstance(fmt, int):

            idx, dct = fmt, {}

        elif isinstance(fmt, tuple):

            idx, dct = fmt

        val = line[idx]
        val = dct.get(val, val)

        return val


    @staticmethod
    def _process_partner(fmt, line):

        if isinstance(fmt, int):

            partner = line[fmt]

        elif isinstance(fmt, tuple):

            idx, proc = fmt
            obj = line if idx is None else line[idx]

            partner = proc(obj)

        return partner.strip() if hasattr(partner, 'strip') else partner


    def _map_list(
            self,
            lst,
            single_list = False,
            expand_complexes = True,
        ):
        """
        Maps the names from a list of edges or items (molecules).

        :arg list lst:
            List of items or edge dictionaries whose names have to be
            mapped.
        :arg bool single_list:
            Optional, ``False`` by default. Determines whether the
            provided elements are items or edges. This is, either calls
            :py:meth:`pypath.main.PyPath.map_edge` or
            :py:meth:`pypath.main.PyPath.map_item` to map the item
            names.
        :arg bool expand_complexes:
            Expand complexes, i.e. create links between each member of
            the complex and the interacting partner.

        :return:
            (*list*) -- Copy of *lst* with their elements' names mapped.
        """

        list_mapped = []

        if single_list:

            for item in lst:

                list_mapped += self._map_item(
                    item,
                    expand_complexes = expand_complexes,
                )

        else:

            for edge in lst:

                list_mapped += self._map_edge(
                    edge,
                    expand_complexes = expand_complexes,
                )

        return list_mapped


    def _map_item(self, item, expand_complexes = True):
        """
        Translates the name in *item* representing a molecule. Default
        name types are defined in
        :py:attr:`pypath.main.PyPath.default_name_type` If the mapping
        is unsuccessful, the item will be added to
        :py:attr:`pypath.main.PyPath.unmapped` list.

        :arg dict item:
            Item whose name is to be mapped to a default name type.
        :arg bool expand_complexes:
            Expand complexes, i.e. create links between each member of
            the complex and the interacting partner.

        :return:
            (*list*) -- The default mapped name(s) [str] of *item*.
        """

        # TODO: include
        default_id = mapping.map_name(
            item['name'], item['id_type'],
            self.default_name_types[item['type']],
            expand_complexes = expand_complexes,
        )

        if len(default_id) == 0:

            self.unmapped.append(item['name'])

        return default_id


    def _map_edge(self, edge, expand_complexes = True):
        """
        Translates the identifiers in *edge* representing an edge. Default
        name types are defined in
        :py:attr:`pypath.main.PyPath.default_name_type` If the mapping
        is unsuccessful, the item will be added to
        :py:attr:`pypath.main.PyPath.unmapped` list.

        :arg dict edge:
            Item whose name is to be mapped to a default name type.
        :arg bool expand_complexes:
            Expand complexes, i.e. create links between each member of
            the complex and the interacting partner.

        :return:
            (*list*) -- Contains the edge(s) [dict] with default mapped
            names.
        """

        edge_stack = []

        defnt = self.default_name_types
        def_name_type_a = defnt.get(edge['entity_type_a'], edge['id_type_a'])
        def_name_type_b = defnt.get(edge['entity_type_b'], edge['id_type_b'])

        default_id_a = mapping.map_name(
            edge['id_a'],
            edge['id_type_a'],
            def_name_type_a,
            ncbi_tax_id = edge['taxon_a'],
            expand_complexes = expand_complexes,
        )

        default_id_b = mapping.map_name(
            edge['id_b'],
            edge['id_type_b'],
            def_name_type_b,
            ncbi_tax_id = edge['taxon_b'],
            expand_complexes = expand_complexes,
        )

        # this is needed because the possibility ambigous mapping
        # and expansion of complexes
        # one name can be mapped to multiple ones
        # this multiplies the nodes and edges
        # in case of proteins this does not happen too often
        for id_a, id_b in itertools.product(default_id_a, default_id_b):

            this_edge = copy_mod.copy(edge)
            this_edge['default_name_a'] = id_a
            this_edge['default_name_type_a'] = def_name_type_a

            this_edge['default_name_b'] = id_b
            this_edge['default_name_type_b'] = def_name_type_b

            edge_stack.append(this_edge)

        return edge_stack


    def _process_attrs(self, line, spec, lnum):
        """
        Extracts the extra (custom, resource specific) attributes from a
        line of the input based on the given specification (defined in the
        network input definition).
        """

        attrs = {}

        for col in spec.keys():
            # extra_edge_attrs and extra_node_attrs are dicts
            # of additional parameters assigned to edges and nodes
            # respectively;
            # key is the name of the parameter, value is the col number,
            # or a tuple of col number and the separator,
            # if the column contains additional subfields e.g. (5, ";")

            try:

                if spec[col].__class__ is tuple:

                    if hasattr(spec[col][1], '__call__'):
                        field_value = spec[col][1](line[spec[col][0]])

                    else:
                        field_value = line[spec[col][0]].split(spec[col][1])

                else:
                    field_value = line[spec[col]]

            except:
                self._log(
                    'Wrong column index (%s) in extra attributes? '
                    'Line #%u' % (str(col), lnum),
                    -5,
                )

            field_name = col
            attrs[field_name] = field_value

        return attrs


    def _process_taxon(self, tax_dict, fields): # TODO
        """
        """

        if isinstance(tax_dict, int):

            return tax_dict

        elif 'A' in tax_dict and 'B' in tax_dict:

            return (
                self._process_taxon(tax_dict['A'], fields),
                self._process_taxon(tax_dict['B'], fields),
            )

        else:

            if 'dict' not in tax_dict:
                return int(fields[tax_dict['col']])

            elif fields[tax_dict['col']] in tax_dict['dict']:
                return tax_dict['dict'][fields[tax_dict['col']]]

            else:
                return None


    def _match_taxon(self, tax_dict, taxon, only_default_organism = False):

        has_dict = isinstance(tax_dict, dict)
        has_include = has_dict and 'include' in tax_dict
        has_exclude = has_dict and 'exclude' in tax_dict

        return (
            (
                taxon == _const.NOT_ORGANISM_SPECIFIC
            ) or (
                has_include and
                taxon in tax_dict['include']
            ) or (
                has_exclude and
                taxon not in tax_dict['exclude']
            ) or (
                not has_include and
                not has_exclude and
                (
                    not only_default_organism or
                    taxon == self.ncbi_tax_id
                )
            )
        )


    def _add_edge_list(
            self,
            edge_list = False,
            regulator = False,
            only_directions = False,
            allow_loops = None,
        ):
        """
        Adds edges to the network from *edge_list* obtained from file or
        other input method. If none is passed, checks for such data in
        :py:attr:`pypath.network.Network.edge_list_mapped`.

        :arg str edge_list:
            Optional, ``False`` by default. The source name of the list
            of edges to be added. This must have been loaded previously
            (e.g.: with :py:meth:`pypath.main.PyPath.read_data_file`).
            If none is passed, loads the data directly from
            :py:attr:`pypath.main.PyPath.raw_data`.
        :arg bool regulator:
            Optional, ``False`` by default. If set to ``True``, non
            previously existing nodes, will not be added (and hence, the
            edges involved).
        """

        self._log('Adding preprocessed edge list to existing network.')

        allow_loops = self._allow_loops(allow_loops = allow_loops)

        if not edge_list:

            if (
                hasattr(self, 'edge_list_mapped') and
                self.edge_list_mapped is not None
            ):

                edge_list = self.edge_list_mapped

            else:

                self._log('_add_edge_list(): No data, nothing to do.')
                return True

        if isinstance(edge_list, str):

            if edge_list in self.raw_data:

                edge_list = self.raw_data[edge_list]

            else:

                self._log(
                    '`%s` looks like a source name, but no data '
                    'available under this name.' % edge_list
                )

                return False

        self._filtered_loops = 0

        prg = progress.Progress(
            iterable = edge_list,
            name = 'Processing interactions',
        )

        for e in prg:

            self._add_update_edge(
                e,
                allow_loops = allow_loops,
                only_directions = only_directions,
            )

        self._log(
            'New network resource added, current number '
            'of nodes: %u, edges: %u.' % (
                self.vcount,
                self.ecount
            )
        )

        if not allow_loops:

            self._log('Loop edges discarded: %u' % self._filtered_loops)

        delattr(self, '_filtered_loops')

        self.raw_data = None


    def _add_update_edge(
            self,
            edge,
            allow_loops = None,
            only_directions = False,
        ):
        """
        Adds a new interaction (edge) or updates the attributes of the edge
        if it already exists.

        :arg dict edge:
            A dictionary describing an edge (interaction) with the following
            items:
            :item str id_a:
                Name of the source node of the edge to be added/updated.
            :item str id_b:
                Name of the source node of the edge to be added/updated.
            :item set source:
                Or [list], contains the names [str] of the resources
                supporting that edge.
            :item pypath.evidence.Evidence evidence:
                A ``pypath.evidence.Evidence`` object.
            :item bool is_directed:
                Whether if the edge is directed or not.
            :item set refs:
                Or [list], contains the instances of the references
                :py:class:`pypath.refs.Reference` for that edge.
            :item bool stim:
                Whether the edge is stimulatory or not.
            :item bool inh:
                Whether the edge is inhibitory or note
            :item int taxon_a:
                NCBI Taxonomic identifier of the source molecule.
            :item int taxon_b:
                NCBI Taxonomic identifier of the target molecule.
            :item str typ:
                The type of interaction (e.g.: ``'trascriptional'``)
            :item dict extra_attrs:
                Optional, ``{}`` by default. Contains any extra attributes
                for the edge to be updated.

        :arg bool only_directions:
            Optional, ``False`` by default. If set to ``True`` and the
            edge is not in the network, it won't be created. If it already
            exists the attributes of the new edge will be added to the
            existing one.
        """

        (
            id_a,
            id_b,
            id_type_a,
            id_type_b,
            entity_type_a,
            entity_type_b,
            source,
            evidences,
            is_directed,
            refs,
            positive,
            negative,
            taxon_a,
            taxon_b,
            interaction_type,
            extra_attrs,
            extra_attrs_a,
            extra_attrs_b,
        ) = (
            edge['default_name_a'],
            edge['default_name_b'],
            edge['default_name_type_a'],
            edge['default_name_type_b'],
            edge['entity_type_a'],
            edge['entity_type_b'],
            edge['source'],
            edge['evidences'],
            edge['is_directed'],
            edge['references'],
            edge['positive'],
            edge['negative'],
            edge['taxon_a'],
            edge['taxon_b'],
            edge['interaction_type'],
            edge['attrs_edge'],
            edge['attrs_node_a'],
            edge['attrs_node_b'],
        )

        allow_loops = allow_loops or self.allow_loops

        refs = {refs_mod.Reference(pmid) for pmid in refs}

        entity_a = entity_mod.Entity(
            identifier = id_a,
            id_type = id_type_a,
            entity_type = entity_type_a,
            taxon = taxon_a,
            attrs = extra_attrs_a,
        )
        entity_b = entity_mod.Entity(
            identifier = id_b,
            id_type = id_type_b,
            entity_type = entity_type_b,
            taxon = taxon_b,
            attrs = extra_attrs_b,
        )

        interaction = interaction_mod.Interaction(
            a = entity_a,
            b = entity_b,
            attrs = extra_attrs,
        )

        if not allow_loops and interaction.is_loop():

            self._filtered_loops += 1
            return

        if is_directed:

            interaction.add_evidence(
                evidence = evidences,
                direction = (entity_a, entity_b),
            )

        else:

            interaction.add_evidence(
                evidence = evidences,
                direction = 'undirected',
            )

        # setting signs:
        if positive:

            interaction.add_evidence(
                evidence = evidences,
                direction = (entity_a, entity_b),
                effect = 1,
            )

        if negative:

            interaction.add_evidence(
                evidence = evidences,
                direction = (entity_a, entity_b),
                effect = -1,
            )

        if is_directed and not positive and not negative:

            interaction.add_evidence(
                evidence = evidences,
                direction = (entity_a, entity_b),
                effect = 0,
            )

        self.add_interaction(
            interaction,
            attrs = extra_attrs,
            only_directions = only_directions,
        )


    def organisms_check(
            self,
            organisms = None,
            remove_mismatches = True,
            remove_nonspecific = False,
        ):
        """
        Scans the network for one or more organisms and removes the nodes
        and interactions which belong to any other organism.

        :arg int,set,NoneType organisms:
            One or more NCBI Taxonomy IDs. If ``None`` the value in
            :py:attr:`ncbi_tax_id` will be used. If that's too is ``None``
            then only the entities with discrepancy between their stated
            organism and their identifier.
        :arg bool remove_mismatches:
            Remove the entities where their ``identifier`` can not be found
            in the reference list from the database for their ``taxon``.
        :arg bool remove_nonspecific:
            Remove the entities with taxonomy ID zero, which is used to
            represent the non taxon specific entities such as metabolites
            or drug compounds.
        """

        self._log(
            'Checking organisms. %u nodes and %u interactions before.' % (
                self.vcount,
                self.ecount,
            )
        )

        organisms = common.to_set(organisms or self.ncbi_tax_id)

        to_remove = set()

        for node in self.nodes.values():

            if (
                organisms and
                node.taxon != _const.NOT_ORGANISM_SPECIFIC and
                node.taxon not in organisms
            ):

                to_remove.add(node)

            if (
                (
                    remove_mismatches and
                    not node.entity_type in {
                        'complex',
                        'lncrna',
                        'drug',
                        'small_molecule'
                    } and
                    not reflists.check(
                        name = node.identifier,
                        id_type = node.id_type,
                        ncbi_tax_id = node.taxon,
                    )
                ) or (
                    remove_nonspecific and
                    not node.taxon
                )
            ):

                to_remove.add(node)

        for node in to_remove:

            self.remove_node(node)

        self._log(
            'Finished checking organisms. '
            '%u nodes have been removed, '
            '%u nodes and %u interactions remained.' % (
                len(to_remove),
                self.vcount,
                self.ecount,
            )
        )


    def get_organisms(self):
        """
        Returns the set of all NCBI Taxonomy IDs occurring in the network.
        """

        return {n.taxon for n in self.nodes.values()}


    @property
    def vcount(self):

        return len(self.nodes)


    @property
    def ecount(self):

        return len(self.interactions)


    def make_df(
            self,
            records = None,
            by_source = None,
            with_references = None,
            columns = None,
            dtype = None,
        ):
        """
        Creates a ``pandas.DataFrame`` from the interactions.
        """

        self._log('Creating interactions data frame.')

        by_source = by_source if by_source is not None else self.df_by_source
        with_references = (
            with_references
                if with_references is not None else
            self.df_with_references
        )
        columns = columns or self.df_columns
        dtype = dtype or self.df_dtype

        if not dtype:

            dtype = {
                'id_a': 'category',
                'id_b': 'category',
                'type_a': 'category',
                'type_b': 'category',
                'effect': 'int8',
                'type': 'category',
                'dmodel': 'category' if by_source else 'object',
                'sources': 'category' if by_source else 'object',
                'references': 'object' if with_references else 'category',
            }

        if not records:

            records = self.generate_df_records(
                by_source = by_source,
                with_references = with_references,
            )

        if not isinstance(records, (list, tuple, np.ndarray)):

            records = list(records)

        if not columns and hasattr(records[0], '_fields'):

            columns = records[0]._fields

        self.records = records
        self.dtype = dtype

        self.df = pd.DataFrame(
            records,
            columns = columns,
        )

        ### why?
        if dtype:

            self.df = self.df.astype(dtype)

        self._log(
            'Interaction data frame ready. '
            'Memory usage: %s ' % common.df_memory_usage(self.df)
        )


    def get_df(self):

        if not hasattr(self, 'df'):

            self.make_df()

        return self.df


    def filtered(
            self,
            resource = None,
            entity_type = None,
            data_model = None,
            interaction_type = None,
            only_directed = None,
            only_undirected = None,
            only_signed = None,
            only_proteins = None,
            effect = None,
            entities = None,
            source_entities = None,
            target_entities = None,
            swap_undirected = True,
            **kwargs
        ):

        return self.filter_df(
            df = self.get_df(),
            resource = resource,
            entity_type = entity_type,
            data_model = data_model,
            interaction_type = interaction_type,
            only_directed = only_directed,
            only_undirected = only_undirected,
            only_signed = only_signed,
            only_proteins = only_proteins,
            effect = effect,
            entities = entities,
            source_entities = source_entities,
            target_entities = target_entities,
            swap_undirected = swap_undirected,
            **kwargs
        )


    @staticmethod
    def filter_df(*args, **kwargs):

        return core_common.filter_network_df(*args, **kwargs)


    def generate_df_records(self, by_source = False, with_references = False):

        for ia in self.interactions.values():

            for rec in ia.generate_df_records(
                by_source = by_source,
                with_references = with_references,
            ):

                yield rec


    @classmethod
    def from_igraph(cls, pa, **kwargs):
        """
        Creates an instance from an ``igraph.Graph`` based
        ``pypath.main.PyPath`` object.

        :arg pypath.main.PyPath pa:
            A ``pypath.main.PyPath`` object with network data loaded.
        """

        obj = cls(**kwargs)

        for ia in pa.graph.es['attrs']:

            obj.add_interaction(ia)

        return obj


    def add_interaction(
            self,
            interaction,
            attrs = None,
            only_directions = False,
        ):
        """
        Adds a ready ``pypath.interaction.Interaction`` object to the network.
        If an interaction between the two endpoints already exists, the
        interactions will be merged: this stands for the directions, signs,
        evidences and other attributes.

        :arg interaction.Interaction interaction:
            A ``pypath.interaction.Interaction`` object.
        :arg NoneType,dict attrs:
            Optional, a dictionary of extra (usually resource specific)
            attributes.
        :arg bool only_directions:
            If the interaction between the two endpoints does not exist it
            won't be added to the network. Otherwise all attributes
            (direction, effect sign, evidences, etc) will be merged to the
            existing interaction. Apart from the endpoints also the
            ``interaction_type`` of the existing interaction has to match the
            interaction added here.
        """

        attrs = attrs or {}

        key = (interaction.a, interaction.b)

        if key not in self.interactions:

            if only_directions:

                return

            else:

                self.interactions[key] = interaction

        else:

            if only_directions:

                if (
                    self.interactions[key].get_interaction_types() &
                    interaction.get_interaction_types()
                ):

                    for itype_to_remove in (
                        interaction.get_interaction_types() -
                        self.interactions[key].get_interaction_types()
                    ):

                        interaction.unset_interaction_type(itype_to_remove)

                else:

                    return

            self.interactions[key] += interaction

        self.interactions[key].update_attrs(**attrs)

        self.add_node(interaction.a, add = not only_directions)
        self.add_node(interaction.b, add = not only_directions)

        self.interactions_by_nodes[interaction.a].add(key)
        self.interactions_by_nodes[interaction.b].add(key)


    def add_node(self, entity, attrs = None, add = True):
        """
        Adds a molecular entity to the py:attr:``nodes`` and
        py:attr:``nodes_by_label`` dictionaries.

        :arg entity.Entity entity:
            An object representing a molecular entity.
        :arg NoneType,dict attrs:
            Optional extra attributes to be assigned to the entity.
        :arg bool add:
            Whether to add a new molecular entity to the network if it does
            not exist yet. If ``False`` will only update attributes for
            existing entities otherwise will do nothing.
        """

        if attrs:

            entity.update_attrs(**attrs)

        if entity.identifier in self.nodes:

            self.nodes[entity.identifier] += entity

        elif add:

            self.nodes[entity.identifier] = entity
            self.nodes_by_label[entity.label or entity.identifier] = entity


    def remove_node(self, entity):
        """
        Removes a node with all its interactions.
        If the removal of the interactions leaves any of the partner nodes
        without interactions it will be removed too.

        :arg str,Entity entity:
            A molecular entity identifier, label or ``Entity`` object.
        """

        entity = self.entity(entity)

        if not entity:

            return

        _ = self.nodes.pop(entity.identifier, None)
        _ = self.nodes_by_label.pop(entity.label, None)

        if entity in self.interactions_by_nodes:

            partners = set()

            for i_key in self.interactions_by_nodes[entity].copy():

                self.remove_interaction(*i_key)

            _ = self.interactions_by_nodes.pop(entity, None)


    def remove_interaction(self, entity_a, entity_b):
        """
        Removes the interaction between two nodes if exists.

        :arg str,Entity entity_a,entity_b:
            A pair of molecular entity identifiers, labels or ``Entity``
            objects.
        """

        entity_a = self.entity(entity_a)
        entity_b = self.entity(entity_b)

        key_ab = (entity_a, entity_b)
        key_ba = (entity_b, entity_a)

        _ = self.interactions.pop(key_ab, None)
        _ = self.interactions.pop(key_ba, None)

        keys = {key_ab, key_ba}
        self.interactions_by_nodes[entity_a] -= keys
        self.interactions_by_nodes[entity_b] -= keys

        if (
            entity_a in self.interactions_by_nodes and
            not self.interactions_by_nodes[entity_a]
        ):

            self.remove_node(entity_a)

        if (
            entity_b in self.interactions_by_nodes and
            not self.interactions_by_nodes[entity_b]
        ):

            self.remove_node(entity_b)


    def remove_zero_degree(self):
        """
        Removes all nodes with no interaction.
        """

        self._log(
            'Removing zero degree nodes. '
            '%u nodes and %u interactions before.' % (
                self.vcount,
                self.ecount,
            )
        )

        to_remove = set()

        for node, interactions in iteritems(self.interactions_by_nodes):

            if not interactions:

                to_remove.add(node)

        for node in to_remove:

            self.remove_node(node)

        self._log(
            'Finished removing zero degree nodes. '
            '%u nodes have been removed, '
            '%u nodes and %u interactions remained.' % (
                len(to_remove),
                self.vcount,
                self.ecount,
            )
        )


    def remove_loops(self):
        """
        Removes the loop interactions from the network i.e. the ones with
        their two endpoints being the same entity.
        """

        self._log(
            'Removing loop edges. Number of edges before: %u.' % len(self)
        )

        for ia in list(self):

            if ia.is_loop():

                self.remove_interaction(ia.a, ia.b)

        self._log(
            'Removed loop edges. Number of edges after: %u.' % len(self)
        )


    @property
    def resources(self):
        """
        Returns a set of all resources.
        """

        return set.union(*(ia.get_resources() for ia in self))


    @property
    def resource_names(self):
        """
        Returns a set of all resource names.
        """

        return set.union(*(ia.get_resource_names() for ia in self))


    def entities_by_resource(self):
        """
        Returns a dict of sets with resources as keys and sets of entity IDs
        as values.
        """

        return dict(
            (
                resource,
                set(
                    itertools.chain(
                        *self.df[
                            [
                                resource in resources
                                for resources in self.df.sources
                            ]
                        ][['id_a', 'id_b']].values
                    )
                )
            )
            for resource in self.resources
        )


    def entity_by_id(self, identifier):
        """
        Returns a ``pypath.entity.Entity`` object representing a molecular
        entity by looking it up by its identifier. If the molecule does not
        present in the current network ``None`` will be returned.

        :arg str identifier:
            The identifier of a molecular entity. Unless it's been set
            otherwise for genes/proteins it is the UniProt ID.
            E.g. ``'P00533'``.
        """

        if identifier in self.nodes:

            return self.nodes[identifier]


    def entity_by_label(self, label):
        """
        Returns a ``pypath.entity.Entity`` object representing a molecular
        entity by looking it up by its label. If the molecule does not
        present in the current network ``None`` will be returned.

        :arg str label:
            The label of a molecular entity. Unless it's been set otherwise
            for genes/proteins it is the Gene Symbol. E.g. ``'EGFR'``.
        """

        if label in self.nodes_by_label:

            return self.nodes_by_label[label]


    def interaction(self, a, b):
        """
        Retrieves the interaction `a --> b` if it exists in the network,
        otherwise `b --> a`. If no interaction exist between `a` and `b`
        returns `None`.
        """

        entity_a = self.entity(a)
        entity_b = self.entity(b)

        key_ab = (entity_a, entity_b)
        key_ba = (entity_b, entity_a)

        if key_ab in self.interactions:

            return self.interactions[key_ab]

        elif key_ba in self.interactions:

            return self.interactions[key_ba]


    def random_interaction(self, **kwargs):
        """
        Picks a random interaction from the network.

        Returns
            An Interaction object, or None if the network is empty.
        """

        key = None

        keys = (
            self.get_interactions(**kwargs)
                if kwargs else
            self.interactions.keys()
        )

        for _, key in zip(range(random.randint(0, len(self)) + 1), keys):

            pass

        if key:

            key = tuple(sorted(key, key = lambda e: e.identifier))

        return self.interactions[key] if key else None


    def _get_interaction(self, id_a, id_b, name_type = 'id'):

        method = 'entity_by_%s' % name_type

        entity_a = getattr(self, method)(id_a)
        entity_b = getattr(self, method)(id_b)

        a_b = (entity_a, entity_b)
        b_a = (entity_b, entity_a)

        if a_b in self.interactions:

            return self.interactions[a_b]

        elif b_a in self.interactions:

            return self.interactions[b_a]


    def entity(self, entity):

        if not isinstance(entity, entity_mod.Entity):

            entity = self.entity_by_id(entity) or self.entity_by_label(entity)

        return entity


    def interaction_by_id(self, id_a, id_b):
        """
        Returns a ``pypath.interaction.Interaction`` object by looking it up
        based on a pair of identifiers. If the interaction does not exist
        in the network ``None`` will be returned.

        :arg str id_a:
            The identifier of one of the partners in the interaction. Unless
            it's been set otherwise for genes/proteins it is the UniProt ID.
            E.g. ``'P00533'``.
        :arg str id_b:
            The other partner, similarly to ``id_a``. The order of the
            partners does not matter here.
        """

        return self._get_interaction(id_a, id_b)


    def interaction_by_label(self, label_a, label_b):
        """
        Returns a ``pypath.interaction.Interaction`` object by looking it up
        based on a pair of labels. If the interaction does not exist
        in the network ``None`` will be returned.

        :arg str label_a:
            The label of one of the partners in the interaction. Unless
            it's been set otherwise for genes/proteins it is the Gene Symbol.
            E.g. ``'EGFR'``.
        :arg str label_b:
            The other partner, similarly to ``label_a``. The order of the
            partners does not matter here.
        """

        return self._get_interaction(label_a, label_b, name_type = 'label')


    def to_igraph(self):
        """
        Converts the network to the legacy ``igraph.Graph`` based ``PyPath``
        object.
        """

        raise NotImplementedError


    def __repr__(self):

        return '<Network: %u nodes, %u interactions>' % (
            self.vcount,
            self.ecount,
        )


    def save_to_pickle(self, pickle_file):
        """
        Saves the network to a pickle file.

        :arg str pickle_file:
            Path to the pickle file.
        """

        self._log('Saving to pickle `%s`.' % pickle_file)

        with open(pickle_file, 'wb') as fp:

            pickle.dump(
                obj = (
                    self.interactions,
                    self.nodes,
                    self.nodes_by_label,
                ),
                file = fp,
            )

        self._update_interactions_by_nodes()

        self._log('Saved to pickle `%s`.' % pickle_file)


    def _update_interactions_by_nodes(self):

        self.interactions_by_nodes = collections.defaultdict(set)

        for key, ia in iteritems(self.interactions):

            self.interactions_by_nodes[ia.a].add(key)
            self.interactions_by_nodes[ia.b].add(key)


    def load_from_pickle(self, pickle_file):
        """
        Loads the network to a pickle file.

        :arg str pickle_file:
            Path to the pickle file.
        """

        self._log('Loading from pickle `%s`.' % pickle_file)

        with open(pickle_file, 'rb') as fp:

            (
                self.interactions,
                self.nodes,
                self.nodes_by_label,
            ) = pickle.load(fp)

        self._update_interactions_by_nodes()

        self._log('Loaded from pickle `%s`.' % pickle_file)


    @classmethod
    def from_pickle(cls, pickle_file: str, **kwargs):
        """
        Initializes a new ``Network`` object by loading it from a pickle
        file. Returns a ``Network`` object.

        Args
            pickle_file:
                Path to a pickle file.
            kwargs:
                Passed to ``Network.__init__``.
        """

        new = cls(
            pickle_file = pickle_file,
            **kwargs
        )

        return new


    def extra_directions(
            self,
            resources = 'extra_directions',
            use_laudanna = False,
            use_string = False,
            dataset = 'directionextra',
        ):
        """
        Adds additional direction & effect information from resources having
        no literature curated references, but giving sufficient evidence
        about the directionality for interactions already supported by
        literature evidences from other sources.
        """

        resources = (
            getattr(network_resources, resources)
                if isinstance(resources, str) else
            list(resources)
        )

        if use_laudanna:

            resources.append(
                network_resources.pathway_bad['laudanna_effects']
            )
            resources.append(
                network_resources.pathway_bad['laudanna_directions']
            )

        if use_string:

            pass

        resources = resource_formats.NetworkDataset(
            name = dataset,
            resources = resources,
        )

        self.load(resources = resources, only_directions = True)


    @staticmethod
    def omnipath_resources(
            omnipath = None,
            kinase_substrate_extra = False,
            ligand_receptor_extra = False,
            pathway_extra = False,
            old_omnipath_resources = False,
            exclude = None,
        ) -> list[resource_formats.NetworkResource]:


        def reference_constraints(resources, data_model, relax = True):

            result = []

            resources = (
                resources.values()
                    if isinstance(resources, dict) else
                resources
            )

            resources = copy_mod.deepcopy(resources)

            for res in resources:

                if res.data_model == data_model:

                    res.networkinput.must_have_references = not relax
                    result.append(res)

            return result


        omnipath = omnipath or copy_mod.deepcopy(network_resources.omnipath)
        exclude = common.to_set(exclude)

        if old_omnipath_resources:

            interaction_resources = (
                copy_mod.deepcopy(network_resources.interaction)
            )

            omnipath = copy_mod.deepcopy(omnipath)
            omnipath['biogrid'] = interaction_resources['biogrid']
            omnipath['alz'] = interaction_resources['alz']
            omnipath['netpath'] = interaction_resources['netpath']
            exclude.update({'IntAct', 'HPRD'})

        else:

            omnipath['huri'] = copy_mod.deepcopy(
                network_resources.interaction_misc['huri']
            )

        omnipath = list(omnipath.without(exclude))

        for dataset, data_model, enabled in (
            ('pathwayextra', 'activity_flow', pathway_extra),
            ('ligrecextra', 'ligand_receptor', ligand_receptor_extra),
            ('kinaseextra', 'enzyme_substrate', kinase_substrate_extra),
        ):

            if enabled:

                extra = list(
                    resource_formats.NetworkDataset(
                        name = dataset,
                        resources = reference_constraints(
                            omnipath,
                            data_model,
                        ),
                    )
                )

                omnipath.extend(extra)

        return omnipath


    def load_omnipath(
            self,
            omnipath = None,
            kinase_substrate_extra = False,
            ligand_receptor_extra = False,
            pathway_extra = False,
            extra_directions = True,
            remove_htp = False,
            htp_threshold = 1,
            keep_directed = True,
            remove_undirected = True,
            min_refs_undirected = None,
            min_resources_undirected = 2,
            old_omnipath_resources = False,
            exclude = None,
            pickle_file = None,
            allow_loops = None,
        ):

        self._log('Loading the `OmniPath` network.')

        if pickle_file:

            self.load(pickle_file = pickle_file)
            return

        omnipath = self.omnipath_resources(
            omnipath = omnipath,
            kinase_substrate_extra = kinase_substrate_extra,
            ligand_receptor_extra = ligand_receptor_extra,
            pathway_extra = pathway_extra,
            old_omnipath_resources = old_omnipath_resources,
            exclude = exclude,
        )

        self.load(omnipath, exclude = exclude, allow_loops = allow_loops)


        for dataset, label, enabled in (
            ('pathwayextra', 'activity flow', pathway_extra),
            ('ligrecextra', 'ligand-receptor', ligand_receptor_extra),
            ('kinaseextra', 'enzyme-PTM', kinase_substrate_extra),
        ):

            if enabled:

                self._log(f'Loading extra {label} interactions.')

                self.load(
                    getattr(network_resources, dataset).rename(dataset),
                    exclude = exclude,
                )

        if extra_directions:

            self.extra_directions()

        if remove_htp:

            self.remove_htp(
                threshold = htp_threshold,
                keep_directed = keep_directed,
            )

        if remove_undirected:

            self.remove_undirected(
                min_refs = min_refs_undirected,
                min_resources = min_resources_undirected,
            )

        self._log('Finished loading the `OmniPath` network.')


    def remove_htp(self, threshold = 50, keep_directed = False):

        self._log(
            'Removing high-throughput interactions above threshold %u'
            ' interactions per reference. Directed interactions %s.' % (
                threshold,
                'will be kept' if keep_directed else 'also will be removed'
            )
        )

        to_remove = self.htp_interactions(
            threshold = threshold,
            ignore_directed = keep_directed,
        )

        ecount_before = self.ecount
        vcount_before = self.vcount

        for key in to_remove:

            self.remove_interaction(*key)

        self._log(
            'Interactions with only high-throughput references '
            'have been removed. %u interactions removed. '
            'Number of edges decreased from %u to %u, '
            'number of nodes from %u to %u.' % (
                len(to_remove),
                ecount_before,
                self.ecount,
                vcount_before,
                self.vcount,
            )
        )


    def htp_references(self, threshold = 50):
        """
        Collects the high-throughput references i.e. the ones cited at a
        higher number of interactions than ``threshold``.
        """

        interactions_per_reference = self.numof_interactions_per_reference()

        htp_refs = {
            ref
            for ref, cnt in iteritems(interactions_per_reference)
            if cnt > threshold
        }

        self._log('High-throughput references collected: %u' % len(htp_refs))

        return htp_refs


    def htp_interactions(self, threshold = 50, ignore_directed = False):
        """
        Collects the interactions only from high-throughput studies.

        :returns:
            Set of interaction keys (tuples of entities).
        """

        htp_refs = self.htp_references(threshold = threshold)

        htp_int = set()

        for key, ia in iteritems(self.interactions):

            if (
                (
                    not ignore_directed or
                    not ia.is_directed()
                ) and
                not ia.get_references() - htp_refs
            ):

                htp_int.add(key)

        self._log('High-throughput interactions collected: %u' % len(htp_int))

        return htp_int


    def remove_undirected(self, min_refs = None, min_resources = None):

        self._log(
            'Removing undirected interactions%s%s%s.' % (
                (
                    ' with less than %u references' % min_refs
                )
                if min_refs else '',
                ' and' if min_refs and min_resources else '',
                (
                    ' with less than %u resources ' % min_resources
                ),
            )
        )

        ecount_before = self.ecount
        vcount_before = self.vcount

        to_remove = set()

        for key, ia in iteritems(self.interactions):

            if (
                not ia.is_directed() and
                (
                    not min_refs or
                    ia.count_references() < min_refs
                ) and
                (
                    not min_resources or
                    ia.count_resource_names() < min_resources
                )
            ):

                to_remove.add(key)

        for key in to_remove:

            self.remove_interaction(*key)

        self._log(
            'Undirected interactions %s have been removed. '
            '%u interactions removed. Number of edges '
            'decreased from %u to %u, number of vertices '
            'from %u to %u.' % (
                ''
                    if min_refs is None else
                'with less than %u references' % min_refs,
                len(to_remove),
                ecount_before,
                self.ecount,
                vcount_before,
                self.vcount,
            )
        )


    def numof_interactions_per_reference(self):
        """
        Counts the number of interactions for each literature reference.
        Returns a ``collections.Counter`` object (similar to ``dict``).
        """

        return collections.Counter(
            itertools.chain(
                *(
                    ia.get_references()
                    for ia in self
                )
            )
        )


    def interactions_by_reference(self):
        """
        Creates a ``dict`` with literature references as keys and interactions
        described by each reference as values.
        """

        interactions_by_reference = collections.defaultdict(set)

        for i_key, ia in iteritems(self.interactions):

            for ref in ia.get_references():

                interactions_by_reference[ref].add(i_key)

        return dict(interactions_by_reference)

    #
    # Methods for loading specific datasets or initializing the object
    # with loading datasets
    #

    @classmethod
    def omnipath(
            cls,
            omnipath = None,
            kinase_substrate_extra = False,
            ligand_receptor_extra = False,
            pathway_extra = False,
            extra_directions = True,
            remove_htp = False,
            htp_threshold = 1,
            keep_directed = True,
            min_refs_undirected = 2,
            old_omnipath_resources = False,
            exclude = None,
            ncbi_tax_id = 9606,
            **kwargs
        ):

        make_df = kwargs.pop('make_df', None)

        new = cls(ncbi_tax_id = ncbi_tax_id, **kwargs)

        new.load_omnipath(
            omnipath = omnipath,
            kinase_substrate_extra = kinase_substrate_extra,
            ligand_receptor_extra = ligand_receptor_extra,
            pathway_extra = pathway_extra,
            extra_directions = extra_directions,
            remove_htp = remove_htp,
            htp_threshold = htp_threshold,
            keep_directed = keep_directed,
            min_refs_undirected = min_refs_undirected,
            old_omnipath_resources = old_omnipath_resources,
            exclude = exclude,
        )

        if make_df:

            cls.make_df()

        return new


    @staticmethod
    def dorothea_resources(levels = None, expand_levels = None):

        expand_levels = (
            expand_levels
                if isinstance(expand_levels, bool) else
            settings.get('dorothea_expand_levels')
        )

        dorothea = copy_mod.deepcopy(network_resources.transcription_dorothea)

        if levels:

            dorothea['dorothea'].networkinput.input_args['levels'] = levels

        dorothea = (
            network_resources.dorothea_expand_levels(dorothea, levels = levels)
                if expand_levels else
            dorothea
        )

        dorothea = dorothea.rename('dorothea')

        return dorothea


    def load_dorothea(self, levels = None, expand_levels = None, **kwargs):

        dorothea = self.dorothea_resources(
            levels = levels,
            expand_levels = expand_levels,
        )

        self.load(dorothea, **kwargs)


    @classmethod
    def dorothea(cls, levels = None, ncbi_tax_id = 9606, **kwargs):
        """
        Initializes a new ``Network`` object with loading the transcriptional
        regulation network from DoRothEA.

        :arg NontType,set levels:
            The confidence levels to include.
        """

        make_df = kwargs.pop('make_df', False)

        new = cls(ncbi_tax_id = ncbi_tax_id, **kwargs)

        new.load_dorothea(levels = levels, make_df = make_df)

        return new


    def load_collectri(self, **kwargs):

        self.load(network_resources.collectri, **kwargs)


    @classmethod
    def collectri(cls, ncbi_tax_id = 9606, **kwargs):
        """
        Initializes a new ``Network`` object with loading the transcriptional
        regulation network from CollecTRI.
        """

        make_df = kwargs.pop('make_df', False)

        new = cls(ncbi_tax_id = ncbi_tax_id, **kwargs)

        new.load_collectri(make_df = make_df)

        return new


    def load_transcription(
            self,
            collectri = True,
            dorothea = True,
            original_resources = True,
            dorothea_levels = None,
            exclude = None,
            reread = False,
            redownload = False,
            allow_loops = None,
            **kwargs
        ):

        make_df = kwargs.pop('make_df', None)

        if collectri:

            self.load_collectri(
                reread = reread,
                redownload = redownload,
                allow_loops = allow_loops,
            )

        if dorothea:

            self.load_dorothea(
                levels = dorothea_levels,
                reread = reread,
                redownload = redownload,
                allow_loops = allow_loops,
            )

        if original_resources:

            transcription = (
                original_resources
                    if not isinstance(original_resources, bool) else
                network_resources.transcription_onebyone.rename('tf_target')
            )

            self.load(
                resources = transcription,
                reread = reread,
                redownload = redownload,
                exclude = exclude,
                allow_loops = allow_loops,
            )

        if make_df:

            self.make_df()


    @classmethod
    def transcription(
            cls,
            dorothea = True,
            original_resources = True,
            dorothea_levels = None,
            exclude = None,
            reread = False,
            redownload = False,
            make_df = False,
            ncbi_tax_id = 9606,
            allow_loops = None,
            **kwargs
        ):
        """
        Initializes a new ``Network`` object with loading a transcriptional
        regulation network from all databases by default.

        Args
            kwargs:
                Passed to ``Network.__init__``.
        """

        load_args = locals()
        kwargs = load_args.pop('kwargs')
        ncbi_tax_id = load_args.pop('ncbi_tax_id')
        kwargs['ncbi_tax_id'] = ncbi_tax_id
        cls = load_args.pop('cls')

        new = cls(**kwargs)

        new.load_transcription(**load_args)

        return new


    def load_mirna_target(self, **kwargs):

        if 'resources' not in kwargs:

            kwargs['resources'] = (
                network_resources.mirna_target.rename('mirnatarget')
            )

        self.load(**kwargs)


    @classmethod
    def mirna_target(
            cls,
            resources = None,
            make_df = None,
            reread = False,
            redownload = False,
            exclude = None,
            ncbi_tax_id = 9606,
            **kwargs
        ):
        """
        Initializes a new ``Network`` object with loading a miRNA-mRNA
        regulation network from all databases by default.

        Args
            kwargs:
                Passed to ``Network.__init__``.
        """

        new = cls(ncbi_tax_id = ncbi_tax_id, **kwargs)

        new.load_mirna_target(
            exclude = exclude,
            make_df = make_df,
            reread = reread,
            redownload = redownload,
        )

        return new

    #
    # Methods for querying partners by node
    #

    def partners(
            self,
            entity,
            mode = 'ALL',
            direction: bool | tuple | None = None,
            effect: bool | str | None = None,
            resources: str | set[str] | None = None,
            interaction_type: str | set[str] | None = None,
            data_model: str | set[str] | None = None,
            via: bool | str | set[str] | None = None,
            references: bool | str | set[str] | None = None,
            return_interactions: bool = False,
        ):
        """
        :arg str,Entity,list,set,tuple,EntityList entity:
            An identifier or label of a molecular entity or an
            :py:class:`Entity` object. Alternatively an iterator with the
            elements of any of the types valid for a single entity argument,
            e.g. a list of gene symbols.
        :arg str mode:
            Mode of counting the interactions: `IN`, `OUT` or `ALL` , whether
            to consider incoming, outgoing or all edges, respectively,
            respective to the `node defined in `entity``.

        :returns:
            :py:class:`EntityList` object containing the partners having
            interactions to the queried node(s) matching all the criteria.
            If ``entity`` doesn't present in the network the returned
            ``EntityList`` will be empty just like if no interaction matches
            the criteria.
        """

        if (
            not common.is_str(entity) and
            not hasattr(entity, 'identifier') and
            hasattr(entity, '__iter__')
        ):

            kwargs = locals()
            _ = kwargs.pop('self')
            _ = kwargs.pop('entity')
            _ = kwargs.pop('return_interactions')

            return entity_mod.EntityList(
                set(itertools.chain(*(
                    self.partners(_entity, **kwargs)
                    for _entity in entity
                )))
            )

        entity = self.entity(entity)

        # we need to swap it to make it work relative to the queried entity
        _mode = (
            'IN'
                if mode == 'OUT' else
            'OUT'
                if mode == 'IN' else
            'ALL'
        )

        return (
            entity_mod.EntityList(
                {
                    partner
                    for ia in self.interactions_by_nodes[entity]
                    for partner in self.interactions[ia].get_degrees(
                        mode = _mode,
                        direction = direction,
                        effect = effect,
                        resources = resources,
                        interaction_type = interaction_type,
                        data_model = data_model,
                        via = via,
                        references = references,
                    )
                    if partner != entity or self.interactions[ia].is_loop()
                }
                if entity in self.interactions_by_nodes else
                ()
            )
        )


    def count_partners(self, entity, **kwargs):
        """
        Returns the count of the interacting partners for one or more
        entities according to the specified criteria.
        Please refer to the docs of the ``partners`` method.
        """

        return len(self.partners(entity = entity, **kwargs))


    @classmethod
    def _generate_partners_methods(cls):

        def _create_partners_method(method_args):

            count = method_args.pop('count')
            method = 'count_partners' if count else 'partners'

            @functools.wraps(method_args)
            def _partners_method(*args, **kwargs):

                self = args[0]
                kwargs.update(method_args)

                return getattr(self, method)(*args[1:], **kwargs)

            _partners_method.__doc__ = getattr(cls, method).__doc__

            return _partners_method

        for name_parts, arg_parts in (
            zip(*param)
            for param in
            itertools.product(
                *(iteritems(variety) for variety in cls._partners_methods)
            )
        ):

            for count in (False, True):

                method_args = dict(
                    itertools.chain(
                        *(iteritems(part) for part in arg_parts)
                    )
                )
                method_name = ''.join(name_parts)
                method_name = (
                    'count_%s' % method_name if count else method_name
                )
                method_args['count'] = count
                method = _create_partners_method(method_args)
                method.__name__ = method_name

                setattr(
                    cls,
                    method_name,
                    method,
                )

    #
    # Methods for selecting paths and motives in the network
    #

    def find_paths(
            self,
            start: (
                str | entity.Entity | entity.EntityList |
                Iterable[str | entity.Entity]
            ),
            end: (
                str | entity.Entity | entity.EntityList |
                Iterable[str | entity.Entity] |
                None
            ) = None,
            loops: bool = False,
            mode: Literal['OUT', 'IN', 'ALL'] = 'OUT',
            maxlen: int = 2,
            minlen: int = 1,
            direction: bool | tuple | None = None,
            effect: bool | str | None = None,
            resources: str | set[str] | None = None,
            interaction_type: str | set[str] | None = None,
            data_model: str | set[str] | None = None,
            via: bool | str | set[str] | None = None,
            references: bool | str | set[str] | None = None,
            silent: bool = False,
        ):
        """
        Find paths or motifs in a network.

        Finds all paths up to length ``maxlen`` between groups of nodes.
        In addition is able to search for motifs or select the nodes of a
        subnetwork around certain nodes.

        Args
            start:
                Starting node(s) of the paths.
            end:
                Target node(s) of the paths. If ``None`` any target node will
                be accepted and all paths from the starting nodes with length
                ``maxlen`` will be returned.
            loops:
                Search for loops, i.e. the start and end nodes of each path
                should be the same.
            mode:
                Direction of the paths. ``'OUT'`` means from ``start`` to ``end``,
                ``'IN'`` the opposite direction while ``'ALL'`` both directions.
            maxlen:
                Maximum length of paths in steps, i.e. if maxlen = 3, then
                the longest path may consist of 3 edges and 4 nodes.
            minlen:
                Minimum length of the path.
            silent:
                Indicate progress by showing a progress bar.

        Details
            The arguments: ``direction``, ``effect``, ``resources``,
            ``interaction_type``, ``data_model``, ``via`` and ``references``
            will be passed to the ``partners`` method of this object and from
            there to the relevant methods of the ``Interaction`` and
            ``Evidence`` objects. By these arguments it is possible to filter
            the interactions in the paths according to custom criteria. If any
            of these arguments is a ``tuple`` or ``list``, its first value will
            be used to match the first interaction in the path, the second for
            the second one and so on. If the list or tuple is shorter then
            ``maxlen``, its last element will be used for all interactions.
            If it's longer than ``maxlen``, the remaining elements will be
            discarded. This way the method is able to search for custom
            motives. For example, let's say you want to find the motives
            where the estrogen receptor transcription factor *ESR1*
            transcriptionally regulates a gene encoding a protein which
            then has some effect post-translationally on *ESR1*:

        Examples

            n.find_paths(
                'ESR1',
                loops = True,
                minlen = 2,
                interaction_type = ('transcriptional', 'post_translational'),
            )

            # Or if you are interested only in the -/+ feedback loops i.e.
            # *ESR1 --(-)--> X --(+)--> ESR1*:

            n.find_paths(
                'ESR1',
                loops = True,
                minlen = 2,
                interaction_type = ('transcriptional', 'post_translational'),
                effect = ('negative', 'positive'),
            )
        """

        def list_of_entities(entities):

            entities = (
                (entities,)
                    if isinstance(
                        entities,
                        (str, entity_mod.Entity)
                    ) else
                entities
            )

            entities = [self.entity(en) for en in entities]

            return entities


        def interaction_arg(value):

            value = (
                tuple(value)
                    if isinstance(value, (tuple, list)) else
                (value,)
            )

            value = value + (value[-1],) * (maxlen - len(value))
            value = value[:maxlen]

            return value


        def find_all_paths_aux(start, end, path, maxlen = None):

            path = path + [start]

            if (
                len(path) >= minlen + 1 and
                (
                    start == end or
                    (
                        end is None and
                        not loops and
                        len(path) == maxlen + 1
                    ) or
                    (
                        loops and
                        path[0] == path[-1]
                    )
                )
            ):

                return [path]

            paths = []

            if len(path) <= maxlen:

                next_steps = set(
                    self.partners(
                        entity = start,
                        **interaction_args[len(path) - 1]
                    )
                )

                next_steps = next_steps if loops else next_steps - set(path)

                for node in next_steps:

                    paths.extend(
                        find_all_paths_aux(
                            node,
                            end,
                            path, maxlen
                        )
                    )

            return paths


        minlen = max(1, minlen)
        start = list_of_entities(start)
        end = list_of_entities(end) if end else (None,)

        interaction_args = {
            'mode': interaction_arg(mode),
            'direction': interaction_arg(direction),
            'effect': interaction_arg(effect),
            'resources': interaction_arg(resources),
            'interaction_type': interaction_arg(interaction_type),
            'data_model': interaction_arg(data_model),
            'via': interaction_arg(via),
            'references': interaction_arg(references),
        }
        interaction_args = tuple(
            dict(
                (key, interaction_args[key][i])
                for key in interaction_args.keys()
            )
            for i in range(maxlen)
        )

        all_paths = []

        if not silent:
            prg = progress.Progress(
                len(start) * len(end),
                'Looking up all paths up to length %u' % maxlen, 1)

        for s in start:

            for e in end:

                if not silent:
                    prg.step()

                all_paths.extend(find_all_paths_aux(s, e, [], maxlen))

        if not silent:
            prg.terminate()

        return all_paths

    #
    # Methods for collecting interaction attributes across the network
    #

    def _collect(
            self,
            what,
            by = None,
            add_total = False,
            **kwargs
        ):
        """
        Collects the values of an attribute over all interactions in the
        network.

        Args
            kwargs:
                Passed to methods of
                :py:class:`pypath.interaction.Interaction`.
        """

        result = set() if not by else collections.defaultdict(set)

        method = self._get_by_method_name(what, by)

        if not hasattr(interaction_mod.Interaction, method):

            self._log('Collecting attributes: no such method: `%s`.' % method)

        else:

            for ia in self:

                ia_attrs = getattr(ia, method)(**kwargs)

                if by:

                    for grp, val in iteritems(ia_attrs):

                        result[grp].update(val)

                else:

                    result.update(ia_attrs)

        if by and add_total:

            result['total'] = set.union(*result.values())

        return dict(result) if by else result


    @classmethod
    def _generate_collect_methods(cls):

        def _create_collect_method(what):

            @functools.wraps(what)
            def _collect_method(self, **kwargs):

                kwargs['what'] = what

                self._log('Collecting `%s`.' % what)

                collection = self._collect(
                    by = 'interaction_type_and_data_model_and_resource',
                    **kwargs
                )

                return (
                    NetworkEntityCollection(
                        collection = collection,
                        label = what,
                    )
                )

            return _collect_method


        for _get in interaction_mod.Interaction._get_methods:

            method = _create_collect_method(_get)
            method_name = 'collect_%s' % _get
            doc = (
                'Builds a comprehensive collection of `%s` entities '
                'across the network, counts unique and shared objects '
                'by resource, data model and interaction types.' % _get
            )
            signature = interaction_mod.Interaction._get_method_signature

            if 'degree' in _get:

                signature = [('mode',)] + signature

            cls._add_method(
                method_name,
                method,
                signature = signature,
                doc = doc,
            )


    def update_summaries(self, collect_args = None):


        def get_labels(lab, key, segments):

            return tuple(
                (
                    '%s%s%s%s' % (
                        key,
                        '_' if seg else '',
                        seg.replace(' ', '_'),
                        '_pct' if pct else '_n',
                    ),
                    '%s%s%s%s' % (lab, ' ' if seg else '', seg, pct)
                )
                for seg in segments
                for pct in ('', r' [%]')
            )


        def add_resource_segments(rec, res, key, lab, segments, coll):

            get = coll[key].__getattribute__

            values = tuple(itertools.chain(*zip(*(
                (
                    get('%s_collection' % n_pct).get(res, 0),
                    get('%s_shared_within_data_model' % n_pct).get(res, 0),
                    get('%s_unique_within_data_model' % n_pct).get(res, 0),
                    get(
                        '%s_shared_within_interaction_type' % n_pct
                    ).get(res, 0),
                    get(
                        '%s_unique_within_interaction_type' % n_pct
                    ).get(res, 0),
                )
                for n_pct in ('n', 'pct')
            ))))

            labels = get_labels(lab, key, segments)

            rec.extend(list(zip(labels, values)))

            return rec


        def add_dmodel_segments(rec, itype, dmodel, key, lab, segments, coll):

            it_dm_key = (itype, dmodel)
            total_key = it_dm_key + ('Total',)

            get = coll[key].__getattribute__

            values = tuple(itertools.chain(*zip(*(
                (
                    get('%s_by_data_model' % n_pct).get(it_dm_key, 0),
                    get(
                        '%s_shared_within_data_model' % n_pct
                    ).get(total_key, 0),
                    get(
                        '%s_unique_within_data_model' % n_pct
                    ).get(total_key, 0),
                    get('%s_shared_by_data_model' % n_pct).get(it_dm_key, 0),
                    get('%s_unique_by_data_model' % n_pct).get(it_dm_key, 0),
                )
                for n_pct in ('n', 'pct')
            ))))

            labels = get_labels(lab, key, segments)

            rec.extend(list(zip(labels, values)))

            return rec


        def add_itype_segments(rec, itype, key, lab, segments, coll):

            get = coll[key].__getattribute__
            total_key = (itype, 'all', 'Total')

            values = tuple(itertools.chain(*zip(*(
                (
                    get('%s_by_interaction_type' % n_pct).get(itype, 0),
                    get(
                        '%s_shared_within_interaction_type' % n_pct
                    ).get(total_key, 0),
                    get(
                        '%s_unique_within_interaction_type' % n_pct
                    ).get(total_key, 0),
                    get('%s_shared_by_data_model' % n_pct).get(total_key, 0),
                    get('%s_unique_by_data_model' % n_pct).get(total_key, 0),
                )
                for n_pct in ('n', 'pct')
            ))))

            labels = get_labels(lab, key, segments)

            rec.extend(list(zip(labels, values)))

            return rec


        collect_args = collect_args or {'via': False}


        required = collections.OrderedDict(
            entities = 'Entities',
            proteins = 'Proteins',
            mirnas = 'miRNAs',
            interactions_0 = 'Edges',
            references = 'References',
            curation_effort = 'Curation effort',
            interactions_non_directed_0 = 'Undirected interactions',
            interactions_directed = 'Directed interactions',
            interactions_positive = 'Stimulatory interactions',
            interactions_negative = 'Inhibitory interactions',
            interactions_mutual = 'Mutual interactions',
        )

        segments = (
            '',
            'shared within database category',
            'unique within database category',
            'shared within interaction type',
            'unique within interaction type',
        )

        self.summaries = []

        coll = {}

        self._log('Updating summaries.')

        for method in required.keys():

            coll[method] = getattr(self, 'collect_%s' % method)(
                **collect_args
            )

        for itype in self.get_interaction_types():

            for dmodel in self.get_data_models(interaction_type = itype):

                for res in sorted(
                    self.get_resource_names(
                        interaction_type = itype,
                        data_model = dmodel,
                        **collect_args
                    ),
                    key = lambda r: r.lower()
                ):

                    # compiling a record for each resource
                    # within the data model

                    rec = [(('resource', 'Resource'), res)]

                    _res = (itype, dmodel, res)

                    for key, lab in iteritems(required):

                        rec = add_resource_segments(
                            rec, _res, key, lab, segments, coll,
                        )

                    self.summaries.append(rec)

                # compiling a summary record for the data model

                rec = [(
                    ('resource', 'Resource'),
                    '%s total' % dmodel.replace('_', ' ').capitalize()
                )]

                for key, lab in iteritems(required):

                    rec = add_dmodel_segments(
                        rec, itype, dmodel, key, lab, segments, coll,
                    )

                self.summaries.append(rec)

            # compiling a summary record for the interaction type

            rec = [(
                ('resource', 'Resource'),
                '%s total' % itype.replace('_', ' ').capitalize()
            )]

            for key, lab in iteritems(required):

                rec = add_itype_segments(rec, itype, key, lab, segments, coll)

            self.summaries.append(rec)

        # maybe we could compile a summary record for the entire network

        self.summaries = [
            collections.OrderedDict(rec)
            for rec in self.summaries
        ]

        self._log('Finished updating summaries.')


    def summaries_tab(
            self,
            outfile = None,
            return_table = False,
            label_type = 1,
        ):
        """
        Creates a table from resource vs. entity counts and optionally
        writes it to ``outfile`` and returns it.
        """

        tab = []

        tab.append(key[label_type] for key in self.summaries[0].keys())

        for rec in self.summaries:

            tab.append([str(val) for val in rec.values()])

        if outfile:

            with open(outfile, 'w') as fp:

                fp.write('\n'.join('\t'.join(row) for row in tab))

        if return_table:

            return tab


    def homology_translate(self, taxon, exclude = None):

        self._log(
            'Translating network by homology from organism `%u` to `%u`.' % (
                self.ncbi_tax_id,
                taxon,
            )
        )

        new = Network(ncbi_tax_id = taxon)

        n_ia_translated = 0
        entities_translated = set()

        for ia in self:

            ia_translated = False

            for new_ia in ia.homology_translate(
                taxon = taxon,
                exclude = exclude,
            ):

                new.add_interaction(new_ia)
                ia_translated = True
                entities_translated.update(ia.get_entities())

            n_ia_translated += ia_translated

        self._log(
            'Orthology translation ready. '
            '%u out of %u interactions (%.02f%%), '
            '%u out of %u entities (%.02f%%) '
            'have been translated.' % (
                n_ia_translated,
                len(self),
                n_ia_translated / len(self) * 100,
                len(entities_translated),
                len(self.nodes),
                len(entities_translated) / len(self.nodes) * 100,
            )
        )

        return new


    @staticmethod
    def _get_by_method_name(get, by):

        return (
            ''.join(
                (
                    'get_' if not by else '',
                    get,
                    '_by_' if by else '',
                    by or '',
                )
            )
        )


    @staticmethod
    def _iter_get_by_methods():

        return (
            itertools.product(
                interaction_mod.Interaction._get_methods | {'entities'},
                interaction_mod.Interaction._by_methods + (None,),
            )
        )

    @classmethod
    def _generate_get_methods(cls):

        def _create_get_method(what, by):

            wrap_args = (what, by)

            @functools.wraps(wrap_args)
            def _get_by_method(*args, **kwargs):

                what, by = wrap_args

                self = args[0]
                kwargs['what'] = what
                kwargs['by'] = by

                return self._collect(**kwargs)

            return _get_by_method


        for _get, _by in cls._iter_get_by_methods():

            method_name = cls._get_by_method_name(_get, _by)

            setattr(
                cls,
                method_name,
                _create_get_method(what = _get, by = _by),
            )


    @classmethod
    def _generate_count_methods(cls):

        def _create_count_method(what, by):

            method_name = cls._get_by_method_name(what, by)

            @functools.wraps(method_name)
            def _count_method(*args, **kwargs):

                self = args[0]

                collection = getattr(self, method_name)(**kwargs)

                return (
                    len(collection)
                        if isinstance(collection, set) else
                    common.dict_counts(collection)
                )

            return _count_method


        for _get, _by in cls._iter_get_by_methods():

            method_name = (
                'count_%s' % (
                    cls._get_by_method_name(_get, _by).replace('get_', '')
                )
            )

            setattr(
                cls,
                method_name,
                _create_count_method(what = _get, by = _by)
            )


    @classmethod
    def _add_method(cls, method_name, method, signature = None, doc = None):

        common.add_method(
            cls,
            method_name,
            method,
            signature = signature,
            doc = doc,
        )


    def _allow_loops(self, allow_loops = None, resource = None):
        """
        Integrates settings for the `allow_loops` parameter from the
        method, instance and module level settings.
        """

        default = settings.get('network_allow_loops')

        return (
            # from the arguments of the actual `load` call
            allow_loops
                if isinstance(allow_loops, bool) else
            # from the current instance
            self.allow_loops
                if isinstance(self.allow_loops, bool) else
            # resource specific settings
            resource.networkinput.allow_loops
                if (
                    hasattr(resource, 'networkinput') and
                    isinstance(resource.networkinput.allow_loops, bool)
                ) else
            # interaction type specific settings from the module level
            resource.networkinput.interaction_type in default
                if (
                    isinstance(default, _const.LIST_LIKE) and
                    hasattr(resource, 'networkinput')
                ) else
            # general settings from the module level
            bool(default)
        )


    def count_loops(self):

        return sum(ia.is_loop() for ia in self)


    def direction_consistency(self):
        """
        Collects statistics about the consistency of interaction
        directions between resources.
        * total_directed: number of directed edges
        * shared_directed: number of directed edges in overlap with other
          resources
        * consistent_edges: number of edges consistent with other resources
        * inconsistent_edges: number of edges inconsistent with other
          resources
        * total_consistency: sum of consistencies (for all edges and all
          resources)
        * total_inconsistency: sum of inconsistencies (for all edges and all
          resources)
        """

        def dd_matrix(dd):

            names = list(dd.keys())

            return pd.DataFrame(
                [
                    [key] + list(val.values())
                    for key, val in dd.items()
                ],
                columns = ['resource'] + names,
            )


        DirectionConsistency = collections.namedtuple(
            'DirectionConsistency',
            [
                'total_directed',
                'shared_directed',
                'consistent_edges',
                'inconsistent_edges',
                'total_consistency',
                'total_inconsistency',
                'total_signed',
                'shared_signed',
                'consistent_signed_edges',
                'inconsistent_signed_edges',
                'total_sign_consistency',
                'total_sign_inconsistency',
            ]
        )

        summary = {}

        resources = sorted(self.get_resource_names(via = False))
        consistencies = collections.OrderedDict(
            (
                resource1,
                collections.OrderedDict(
                    (resource2, 0)
                    for resource2 in resources
                )
            )
            for resource1 in resources
        )
        inconsistencies = copy_mod.deepcopy(consistencies)
        sign_consistencies = copy_mod.deepcopy(consistencies)
        sign_inconsistencies = copy_mod.deepcopy(consistencies)

        for resource in resources:

            total_directed = 0
            shared_directed = 0
            consistent_edges = 0
            inconsistent_edges = 0
            total_consistency = 0
            total_inconsistency = 0
            total_signed = 0
            shared_signed = 0
            consistent_signed_edges = 0
            inconsistent_signed_edges = 0
            total_sign_consistency = 0
            total_sign_inconsistency = 0

            for ia in self:

                if not ia.is_directed():

                    continue

                res_a_b = ia.direction[ia.a_b].get_resource_names(via = False)
                res_b_a = ia.direction[ia.b_a].get_resource_names(via = False)
                res_a_b_pos = ia.positive[ia.a_b].get_resource_names(
                    via = False
                )
                res_a_b_neg = ia.negative[ia.a_b].get_resource_names(
                    via = False
                )
                res_b_a_pos = ia.positive[ia.b_a].get_resource_names(
                    via = False
                )
                res_b_a_neg = ia.negative[ia.b_a].get_resource_names(
                    via = False
                )

                if resource in res_a_b or resource in res_b_a:

                    total_directed += 1

                else:

                    continue

                if resource in res_a_b_pos or resource in res_a_b_neg:

                    total_signed += 1

                if resource in res_b_a_pos or resource in res_b_a_neg:

                    total_signed += 1

                if len(res_a_b | res_b_a) > 1:

                    shared_directed += 1

                if len(res_a_b_pos | res_a_b_neg) > 1:

                    shared_signed += 1

                if len(res_b_a_pos | res_b_a_neg) > 1:

                    shared_signed += 1

                if (
                    (resource in res_a_b and len(res_a_b) > 1) or
                    (resource in res_b_a and len(res_b_a) > 1)
                ):

                    consistent_edges += 1

                if (
                    (resource in res_a_b_pos and len(res_a_b_pos) > 1) or
                    (resource in res_a_b_neg and len(res_a_b_neg) > 1)
                ):

                    consistent_signed_edges += 1

                if (
                    (resource in res_b_a_pos and len(res_b_a_pos) > 1) or
                    (resource in res_b_a_neg and len(res_b_a_neg) > 1)
                ):

                    consistent_signed_edges += 1

                if (
                    (
                        resource in res_a_b and
                        resource not in res_b_a and
                        res_b_a
                    ) or
                    (
                        resource in res_b_a and
                        resource not in res_a_b and
                        res_a_b
                    )
                ):

                    inconsistent_edges += 1

                if (
                    (
                        resource in res_a_b_pos and
                        resource not in res_a_b_neg and
                        res_a_b_neg
                    ) or
                    (
                        resource in res_a_b_neg and
                        resource not in res_a_b_pos and
                        res_a_b_pos
                    )
                ):

                    inconsistent_signed_edges += 1

                if (
                    (
                        resource in res_b_a_pos and
                        resource not in res_b_a_neg and
                        res_b_a_neg
                    ) or
                    (
                        resource in res_b_a_neg and
                        resource not in res_b_a_pos and
                        res_b_a_pos
                    )
                ):

                    inconsistent_signed_edges += 1

                if resource in res_a_b:

                    total_consistency += len(res_a_b) - 1

                else:

                    total_inconsistency += len(res_a_b)

                if resource in res_a_b_pos:

                    total_sign_consistency += len(res_a_b_pos) - 1

                if resource in res_a_b_neg:

                    total_sign_consistency += len(res_a_b_neg) - 1

                if resource in res_b_a_pos:

                    total_sign_consistency += len(res_b_a_pos) - 1

                if resource in res_b_a_neg:

                    total_sign_consistency += len(res_b_a_neg) - 1

                if resource not in res_a_b_pos:

                    total_sign_inconsistency += len(res_a_b_pos)

                if resource not in res_a_b_neg:

                    total_sign_inconsistency += len(res_a_b_neg)

                if resource not in res_b_a_pos:

                    total_sign_inconsistency += len(res_b_a_pos)

                if resource not in res_b_a_neg:

                    total_sign_inconsistency += len(res_b_a_neg)

                if resource in res_b_a:

                    total_consistency += len(res_b_a) - 1

                else:

                    total_inconsistency += len(res_b_a)

                for dir_resources in (res_a_b, res_b_a):

                    for res_other in dir_resources:

                        if resource in dir_resources:

                            consistencies[resource][res_other] += 1

                        else:

                            inconsistencies[resource][res_other] += 1

                for sign_resources in (
                    res_a_b_pos,
                    res_a_b_neg,
                    res_b_a_pos,
                    res_a_b_neg,
                ):

                    for res_other in sign_resources:

                        if resource in sign_resources:

                            sign_consistencies[resource][res_other] += 1

                        else:

                            sign_inconsistencies[resource][res_other] += 1

            summary[resource] = DirectionConsistency(
                total_directed = total_directed,
                shared_directed = shared_directed,
                consistent_edges = consistent_edges,
                inconsistent_edges = inconsistent_edges,
                total_consistency = total_consistency,
                total_inconsistency = total_inconsistency,
                total_signed = total_signed,
                shared_signed = shared_signed,
                consistent_signed_edges = consistent_signed_edges,
                inconsistent_signed_edges = inconsistent_signed_edges,
                total_sign_consistency = total_sign_consistency,
                total_sign_inconsistency = total_sign_inconsistency,
            )

        consistencies = dd_matrix(consistencies)
        inconsistencies = dd_matrix(inconsistencies)
        sign_consistencies = dd_matrix(sign_consistencies)
        sign_inconsistencies = dd_matrix(sign_inconsistencies)

        summary = pd.DataFrame(
            [
                [resource] + list(values)
                for resource, values in summary.items()
            ],
            columns = ['resource'] + list(DirectionConsistency._fields),
        )

        return {
            'summary': summary,
            'consistencies': consistencies,
            'inconsistencies': inconsistencies,
            'sign_consistencies': sign_consistencies,
            'sign_inconsistencies': sign_inconsistencies,
        }


Network._generate_get_methods()
Network._generate_partners_methods()
Network._generate_count_methods()
Network._generate_collect_methods()


def init_db(use_omnipath = False, method = None, **kwargs):

    method_name = (
        'load_omnipath'
            if use_omnipath else
        (method or 'load')
    )

    new_network = Network()
    maybe_network = getattr(new_network, method_name)(**kwargs)

    globals()['db'] = maybe_network or new_network


def get_db(**kwargs):

    if 'db' not in globals():

        init_db(**kwargs)

    return globals()['db']
