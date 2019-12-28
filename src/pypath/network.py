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

import importlib as imp
import os
import collections
import itertools
import functools
import copy as copy_mod
import pickle

import pandas as pd

import pypath.session_mod as session_mod
import pypath.interaction as interaction_mod
import pypath.evidence as evidence
import pypath.entity as entity_mod
import pypath.common as common
import pypath.settings as settings
import pypath.cache as cache_mod
import pypath.mapping as mapping
import pypath.dataio as dataio
import pypath.curl as curl
import pypath.refs as refs_mod
import pypath.resources.network as network_resources


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


class Network(session_mod.Logger):
    """
    Represents a molecular interaction network. Provides various methods to
    query the network and its components. Optionally converts the network
    to a ``pandas.DataFrame`` of interactions.
    
    :arg list,dict resources:
        One or more lists or dictionaries containing
        ``pypath.resource.NetworkResource`` objects.
    :arg bool make_df:
        Create a ``pandas.DataFrame`` already when creating the instance.
        If no network data loaded no data frame will be created.
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
            **kwargs
        ):
        
        session_mod.Logger.__init__(self, name = 'network')
        
        self._log('Creating network object.')
        
        self.df_by_source = df_by_source
        self.df_with_references = df_with_references
        self.df_columns = df_columns
        self.df_dtype = df_dtype
        
        self.cache_dir = cache_mod.get_cachedir()
        self.keep_original_names = settings.get('network_keep_original_names')
        self.default_name_types = settings.get('default_name_types')
        
        self.reset()
        
        if pickle_file and os.path.exists(pickle_file):
            
            self.load_from_pickle(pickle_file = pickle_file)
            return
        
        self.load(resources = resources, make_df = make_df, **kwargs)
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def __len__(self):
        
        return len(self.interactions)
    
    
    def __iter__(self):
        
        for ia in self.interactions.values():
            
            yield ia
    
    
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
            ``pypath.resource.NetworkResource`` object it will be processed
            and added to the network.
        :arg bool make_df:
            Whether to create a ``pandas.DataFrame`` after loading all
            resources.
        :arg NoneType,set exclude:
            A *set* of resource names to be ignored. It is useful if you want
            to load a collection with the exception of a few resources.
        """
        
        kwargs = {
            'reread': reread,
            'redownload': redownload,
            'keep_raw': keep_raw,
            'top_call': False,
        }
        
        exclude = common.to_set(exclude)
        
        resources = (
            (resources,)
                if not isinstance(resources, (list, dict, tuple, set)) else
            resources.values()
                if isinstance(resources, dict) else
            resources
        )
        
        for resource in resources:
            
            if (
                isinstance(resource, common.basestring) and
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
                        network_resources.data_formats.\
                            input_formats.NetworkInput,
                        network_resources.resource.NetworkResource,
                    )
                ) and resource.name not in exclude
            ):
                
                self.load_resource(resource, **kwargs)
                
            else:
                
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
        """

        self._log('Loading network data from resource `%s`.' % resource.name)

        self._read_resource(
            resource,
            reread = reread,
            redownload = redownload,
            keep_raw = keep_raw,
        )
        self._add_edge_list(only_directions = only_directions)
        
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
        """
        
        self._log('Reading network data from `%s`.' % resource.name)

        # workaround in order to make it work with both NetworkInput
        # and NetworkResource type param
        _resource = (
            resource
                if isinstance(
                    resource,
                    network_resources.resource.NetworkResource
                ) else
            network_resources.resource.NetworkResource(
                name = resource.name,
                interaction_type = resource.interaction_type,
                networkinput = resource,
                data_model = resource.data_model or 'unknown',
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
                        network_resources.data_formats.\
                            input_formats.NetworkInput,
                        network_resources.resource.NetworkResource,
                    )
                ):
                    
                    self._log(
                        '_read_network_data: No proper input file '
                        'definition. `param` should be either '
                        'a `pypath.input_formats.NetworkInput` or a '
                        '`pypath.resource.NetworkResource` instance.',
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
                        answer = raw_input().lower()

                        if answer == 'n':
                            return None

                        elif answer == 'y':
                            break

                        else:
                            sys.stdout.write(
                                '\n\tPlease answer `y` or `n`:\n\t')
                            sys.stdout.flush()

                input_func = (
                    getattr(dataio, networkinput.input)
                        if hasattr(dataio, networkinput.input) else
                    None
                )

                # reading from remote or local file, or executing import
                # function:
                if (
                    isinstance(networkinput.input, common.basestring) and (
                        networkinput.input.startswith('http') or
                        networkinput.input.startswith('ftp')
                    )
                ):

                    curl_use_cache = not redownload
                    c = curl.Curl(
                        networkinput.input,
                        silent=False,
                        large=True,
                        cache=curl_use_cache
                    )
                    infile = c.fileobj.read()

                    if type(infile) is bytes:

                        try:
                            infile = infile.decode('utf-8')

                        except:

                            try:
                                infile = infile.decode('iso-8859-1')

                            except:
                                pass

                    infile = [
                        x for x in infile.replace('\r', '').split('\n')
                        if len(x) > 0
                    ]
                    self._log(
                        "Retrieving data from%s ..." % networkinput.input
                    )

                # elif hasattr(dataio, networkinput.input):
                elif input_func is not None:

                    self._log("Retrieving data by dataio.%s() ..." %
                                    input_func.__name__)

                    _store_cache = curl.CACHE

                    if isinstance(redownload, bool):

                        curl.CACHE = not redownload

                    # this try-except needs to be removed
                    # once correct exception handling will
                    # be implemented in every input function
                    try:
                        infile = input_func(**networkinput.input_args)

                    except Exception as e:
                        self._log(
                            'Error in `pypath.dataio.%s()`. '
                            'Skipping to next resource. '
                            'See below the traceback.' % input_func.__name__
                        )
                        self._log(str(e.args))

                        try:
                            traceback.print_tb(
                                e.__traceback__, file = sys.stdout)

                        except Exception as e:
                            self._log('Failed handling exception.')
                            self._log(str(e.args))

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
                        '`%s`: Could not find file or dataio function '
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
            sig_col = None if not isinstance(sign, tuple) else sign[0]
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
                dir_val = dir_val if type(dir_val[
                    0]) in common.simple_types else common.flat_list(dir_val)
                dir_sep = sign[3] if len(sign) > 3 else None

            dir_val = common.to_set(dir_val)

            must_have_references = (
                settings.get('keep_noref') or
                networkinput.must_have_references
            )
            self._log(
                'Resource `%s` %s have literature references '
                'for all interactions. Interactions without references '
                'will be dropped. You can alter this condition globally by '
                '`pypath.settings.keep_noref` or for individual resources '
                'by the `must_have_references` attribute of their '
                '`NetworkInput` object.' % (
                    networkinput.name,
                    'must' if must_have_references else 'does not need to'
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

            for lnum, line in enumerate(infile):

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

                # applying filters:
                if self._filters(
                    line,
                    networkinput.positive_filters,
                    networkinput.negative_filters
                ):

                    input_filtered += 1
                    continue

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

                refs = []
                if ref_col is not None:

                    if isinstance(line[ref_col], (list, set, tuple)):

                        refs = line[ref_col]

                    elif isinstance(line[ref_col], int):

                        refs = (line[ref_col],)

                    else:

                        refs = line[ref_col].split(ref_sep)

                    refs = common.del_empty(list(set(refs)))

                refs = dataio.only_pmids([str(r).strip() for r in refs])

                if len(refs) == 0 and must_have_references:
                    ref_filtered += 1
                    continue

                # to give an easy way for input definition:
                if isinstance(networkinput.ncbi_tax_id, int):
                    taxon_a = networkinput.ncbi_tax_id
                    taxon_b = networkinput.ncbi_tax_id

                # to enable more sophisticated inputs:
                elif isinstance(networkinput.ncbi_tax_id, dict):

                    taxx = self._process_taxon(
                        networkinput.ncbi_tax_id,
                        line,
                    )

                    if isinstance(taxx, tuple):
                        taxon_a = taxx[0]
                        taxon_b = taxx[1]

                    else:
                        taxon_a = taxon_b = taxx

                    taxdA = (
                        networkinput.ncbi_tax_id['A']
                        if 'A' in networkinput.ncbi_tax_id else
                        networkinput.ncbi_tax_id
                    )
                    taxdB = (
                        networkinput.ncbi_tax_id['B']
                        if 'B' in networkinput.ncbi_tax_id else
                        networkinput.ncbi_tax_id
                    )

                    if (('include' in taxdA and
                        taxon_a not in taxdA['include']) or
                        ('include' in taxdB and
                        taxon_b not in taxdB['include']) or
                        ('exclude' in taxdA and
                        taxon_a in taxdA['exclude']) or
                        ('exclude' in taxdB and
                        taxon_b in taxdB['exclude'])):

                        taxon_filtered += 1
                        continue

                else:
                    taxon_a = taxon_b = self.ncbi_tax_id

                if taxon_a is None or taxon_b is None:
                    taxon_filtered += 1
                    continue

                positive = False
                negative = False

                if isinstance(sign, tuple):
                    positive, negative = (
                        self._process_sign(line[sign[0]], sign)
                    )

                resource = (
                    line[networkinput.resource]
                        if isinstance(networkinput.resource, int) else
                    line[networkinput.resource[0]].split(
                        networkinput.resource[1]
                    )
                        if isinstance(networkinput.resource, tuple) else
                    networkinput.resource
                )

                resource = common.to_set(resource)

                _resources_secondary = tuple(
                    network_resources.resource.NetworkResource(
                        name = sec_res,
                        interaction_type = _resource.interaction_type,
                        data_model = _resource.data_model,
                        via = _resource.name,
                    )
                    for sec_res in resource
                    if sec_res != _resource.name
                )

                resource.add(networkinput.name)

                id_a = line[networkinput.id_col_a]
                id_b = line[networkinput.id_col_b]
                id_a = id_a.strip() if hasattr(id_a, 'strip') else id_a
                id_b = id_b.strip() if hasattr(id_b, 'strip') else id_b

                evidences = evidence.Evidences(
                    evidences = (
                        evidence.Evidence(
                            resource = _res,
                            references = refs,
                        )
                        for _res in
                        _resources_secondary + (_resource,)
                    )
                )


                new_edge = {
                    'id_a': id_a,
                    'id_b': id_b,
                    'id_type_a': networkinput.id_type_a,
                    'id_type_b': networkinput.id_type_b,
                    'entity_type_a': networkinput.entity_type_a,
                    'entity_type_b': networkinput.entity_type_b,
                    'source': resource,
                    'is_directed': this_edge_dir,
                    'references': refs,
                    'positive': positive,
                    'negative': negative,
                    'taxon_a': taxon_a,
                    'taxon_b': taxon_b,
                    'interaction_type': networkinput.interaction_type,
                    'evidences': evidences,
                }

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

                if networkinput.mark_source:

                    attrs_node_a[networkinput.mark_source] = this_edge_dir

                if networkinput.mark_target:

                    attrs_node_b[networkinput.mark_target] = this_edge_dir

                # merging dictionaries
                node_attrs = {
                    'attrs_node_a': attrs_node_a,
                    'attrs_node_b': attrs_node_b,
                    'attrs_edge': attrs_edge,
                }
                new_edge.update(node_attrs)

                if read_error:

                    self._log(
                        'Errors occured, certain lines skipped.'
                        'Trying to read the remaining.\n',
                        5,
                    )

                edge_list.append(new_edge)

            if hasattr(infile, 'close'):

                infile.close()

            # ID translation of edges
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
    
    
    def _filters(
            self,
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
        """

        negative_filters = negative_filters or ()

        for filtr in negative_filters:

            if len(filtr) > 2:
                sep = filtr[2]
                thisVal = set(line[filtr[0]].split(sep))

            else:
                thisVal = set([line[filtr[0]]])

            filtrVal = common.to_set(filtr[1])

            if thisVal & filtrVal:
                return True

        positive_filters = positive_filters or ()

        for filtr in positive_filters:

            if len(filtr) > 2:
                sep = filtr[2]
                thisVal = set(line[filtr[0]].split(sep))

            else:
                thisVal = {line[filtr[0]]}

            filtrVal = common.to_set(filtr[1])

            if not thisVal & filtrVal:
                return True

        return False
    
    
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
        sign_data = set(str(sign_data).split(sign_sep))
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

        if dir_col is None or dir_val is None:
            
            return False

        else:
            
            this_directed = set(line[dir_col].split(dir_sep))
            return bool(this_directed & dir_val)
    
    
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

        default_id_a = mapping.map_name(
            edge['id_a'],
            edge['id_type_a'],
            self.default_name_types[edge['entity_type_a']],
            ncbi_tax_id = edge['taxon_a'],
            expand_complexes = expand_complexes,
        )

        default_id_b = mapping.map_name(
            edge['id_b'],
            edge['id_type_b'],
            self.default_name_types[edge['entity_type_b']],
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
            this_edge['default_name_type_a'] = (
                self.default_name_types[edge['entity_type_a']]
            )

            this_edge['default_name_b'] = id_b
            this_edge['default_name_type_b'] = (
                self.default_name_types[edge['entity_type_b']]
            )
            edge_stack.append(this_edge)

        return edge_stack
    
    
    def _process_attrs(self, line, spec, lnum): # TODO
        """
        """

        attrs = {}

        for col in spec.keys():
            # extra_edge_attrs and extraNodeAttrs are dicts
            # of additional parameters assigned to edges and nodes respectively;
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

        if 'A' in tax_dict and 'B' in tax_dict:

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
    
    
    def _add_edge_list(
            self,
            edge_list = False,
            regulator = False,
            only_directions = False,
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

        if not edge_list:

            if self.edge_list_mapped is not None:
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
        
        edges = []

        for e in edge_list:
            
            self._add_update_edge(e, only_directions = only_directions)

        self._log(
            'New network resource added, current number '
            'of nodes: %u, edges: %u.' % (
                self.vcount,
                self.ecount
            )
        )
        
        self.raw_data = None
    
    
    def _add_update_edge(
            self,
            edge,
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
        )

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
        
        self.add_interaction(
            interaction,
            attrs = extra_attrs,
            only_directions = only_directions,
        )
    
    
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
                'sources': 'category' if by_source else 'object',
                'references': 'object' if with_references else 'category',
            },
        
        if not records:
            
            records = self.generate_df_records(
                by_source = by_source,
                with_references = with_references,
            )
        
        if not isinstance(records, (list, tuple, pd.np.ndarray)):
            
            records = list(records)
        
        if not columns and hasattr(records[0], '_fields'):
            
            columns = records[0]._fields
        
        self.df = pd.DataFrame(
            records,
            columns = columns,
            dtype = dtype,
        )
        
        ### why?
        if dtype:
            
            self.df = self.df.astype(dtype)
        
        self._log(
            'Interaction data frame ready. '
            'Memory usage: %s ' % common.df_memory_usage(self.df)
        )
    
    
    def generate_df_records(self, by_source = False, with_references = False):
        
        for ia in self.interactions:
            
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
        
        key = (interaction.a, interaction.b)
        
        if key not in self.interactions:
            
            if only_directions:
                
                return
                
            else:
                
                self.interactions[key] = interaction
            
        else:
            
            if only_directions:
                
                if (
                    self.interaction[key].get_interaction_types() &
                    interaction.get_interaction_types()
                ):
                    
                    for itype_to_remove in (
                        interaction.get_interaction_types() -
                        self.interaction[key].get_interaction_types()
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
            self.nodes_by_label[entity.label] = entity
    
    
    def remove_node(self, entity):
        
        entity = self.entity(entity)
        
        if not entity:
            
            return
        
        _ = self.nodes.pop(entity.identifier, None)
        _ = self.nodes_by_label.pop(entity.label, None)
        
        if entity in self.interactions_by_nodes:
            
            for i_key in self.interactions_by_nodes[entity]:
                
                _ = self.interactions.pop(i_key, None)
            
            del self.interactions_by_nodes[entity]
    
    
    def remove_interaction(self, entity_a, entity_b):
        
        entity_a = self.entity(entity_a)
        entity_b = self.entity(entity_b)
        
        key_ab = (entity_a, entity_b)
        key_ba = (entity_b, entity_a)
        
        _ = self.interactions.pop(key_ab, None)
        _ = self.interactions.pop(key_ba, None)
        
        keys = {key_ab, key_ba}
        self.interactions_by_nodes[entity_a] -= keys
        self.interactions_by_nodes[entity_b] -= keys
        
        if not self.interactions_by_nodes[entity_a]:
            
            self.remove_node(entity_a)
        
        if not self.interactions_by_nodes[entity_b]:
            
            self.remove_node(entity_b)
    
    
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
        
        self._log('Saved to pickle `%s`.' % pickle_file)
    
    
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
        
        self._log('Loaded from pickle `%s`.' % pickle_file)
    
    
    @classmethod
    def from_pickle(cls, pickle_file, **kwargs):
        """
        Initializes a new ``Network`` object by loading it from a pickle
        file. Returns a ``Network`` object.
        
        :arg str pickle_file:
            Path to a pickle file.
        **kwargs:
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
        ):
        """
        Adds additional direction & effect information from resources having
        no literature curated references, but giving sufficient evidence
        about the directionality for interactions already supported by
        literature evidences from other sources.
        """
        
        resources = (
            getattr(network_resources, resources)
                if isinstance(resources, common.basestring) else
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
        
        self.load(resources = resources, only_directions = True)
    
    
    def load_omnipath(
            self,
            omnipath = None,
            kinase_substrate_extra = False,
            ligand_receptor_extra = False,
            pathway_extra = False,
            extra_directions = True,
            remove_htp = True,
            htp_threshold = 1,
            keep_directed = True,
            min_refs_undirected = 2,
            old_omnipath_resources = False,
            exclude = None,
        ):
        
        
        def reference_constraints(resources, interaction_type, release):
            
            resources = (
                resources.values()
                    if isinstance(resources, dict) else
                resources
            )
            
            for res in resources:
                
                if res.resource.interaction_type == interaction_type:
                    
                    res.networkinput.must_have_references = not release
        
        
        self._log('Loading the `OmniPath` network.')
        
        omnipath = omnipath or copy_mod.deepcopy(network_resources.omnipath)
        
        if old_omnipath_resources:
            
            omnipath = copy_mod.deepcopy(omnipath)
            omnipath['biogrid'] = network_resources.interaction['biogrid']
            omnipath['alz'] = network_resources.interaction['alz']
            omnipath['netpath'] = network_resources.interaction['netpath']
            exclude = exclude or []
            exclude.extend(['IntAct', 'HPRD'])
        
        reference_constraints(
            omnipath,
            'pathway',
            pathway_extra,
        )
        reference_constraints(
            omnipath,
            'ligand_receptor',
            ligand_receptor_extra,
        )
        reference_constraints(
            omnipath,
            'enzyme_substrate',
            kinase_substrate_extra,
        )
        
        self.load(omnipath, exclude = exclude)
        
        if kinase_substrate_extra:
            
            self._log('Loading extra enzyme-substrate interactions.')
            
            self.load(network_resources.ptm_misc, exclude = exclude)
        
        if ligand_receptor_extra:
            
            self._log('Loading extra ligand-receptor interactions.')
            
            self.load(network_resources.ligand_receptor, exclude = exclude)
        
        if pathway_extra:
            
            self._log('Loading extra activity flow interactions.')
            
            self.load(network_resources.pathway_noref, exclude = exclude)
        
        if extra_directions:
            
            self.extra_directions()
        
        if remove_htp:
            
            self.remove_htp(
                threshold = htp_threshold,
                keep_directed = keep_directed,
            )
        
        if not keep_directed:
            
            self.remove_undirected(min_refs = min_refs_undirected)
        
        self._log('Finished loading the `OmniPath` network.')
    
    
    def remove_htp(self, threshold = 50, keep_directed = False):
        
        self._log(
            'Removing high-throughput interactions above threshold %u'
            ' interactions per reference. Directed interactions %s.' % (
                threshold,
                'will be kept' if keep_directed else 'also will be removed'
            )
        )
        
        interactions_per_reference = self.numof_interactions_per_reference()
        interactions_by_reference = self.interactions_by_reference()
        
        htp_refs = {
            ref
            for ref, cnt in iteritems(interactions_per_reference)
            if cnt > threshold
        }
        
        to_remove = set()
        
        ecount_before = self.ecount
        vcount_before = self.vcount
        
        for key, ia in iteritems(self.interactions):
            
            if (
                not ia.get_references() - htp_refs and (
                    not keep_directed or
                    not ia.is_directed()
                )
            ):
                
                to_remove.add(key)
        
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
    
    
    def remove_undirected(self, min_refs = None):
        
        self._log(
            'Removing undirected interactions%s.' % (
                (
                    'with less than %u references' % min_refs
                )
                if min_refs else ''
            )
        )
        
        ecount_before = self.ecount
        vcount_before = self.vcount
        
        removed = 0
        
        for key, ia in iteritems(self.interactions):
            
            if (
                ia.is_directed() and (
                    not min_refs or
                    len(ia.get_references()) < min_refs
                )
            ):
                
                self.remove_interaction(*key)
                removed += 1
        
        self._log(
            'Undirected interactions %s have been removed. '
            '%u interactions removed. Number of edges '
            'decreased from %u to %u, number of vertices '
            'from %u to %u.' % (
                ''
                    if min_refs is None else
                'with less than %u references' % min_refs,
                removed,
                ecount_before,
                self.ecount,
                vcount_before,
                self.vcount,
            )
        )
    
    
    def numof_interactions_per_reference(self):
        
        return collections.Counter(
            itertools.chain(
                ia.get_references()
                for ia in self.interactions
            )
        )
    
    
    def interactions_by_reference(self):
        
        interactions_by_reference = collections.defaultdict(set)
        
        for i_key, ia in iteritems(self.interactions):
            
            for ref in ia.get_references():
                
                interactions_by_reference[ref].add(i_key)
        
        return interactions_by_reference
    
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
            remove_htp = True,
            htp_threshold = 1,
            keep_directed = True,
            min_refs_undirected = 2,
            old_omnipath_resources = False,
            exclude = None,
            **kwargs
        ):
        
        make_df = kwargs.pop('make_df', None)
        
        new = cls(**kwargs)
        
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
    
    
    def load_dorothea(self, levels = None, **kwargs):
        
        dorothea = copy_mod.deepcopy(network_resources.dorothea['dorothea'])
        
        if levels:
            dorothea.networkinput.input_args['levels'] = levels
        
        self.load(dorothea, **kwargs)
    
    
    @classmethod
    def dorothea(cls, levels = None, **kwargs):
        """
        Initializes a new ``Network`` object with loading the transcriptional
        regulation network from DoRothEA.
        
        :arg NontType,set levels:
            The confidence levels to include.
        """
        
        make_df = kwargs.pop('make_df', False)
        
        new = cls(**kwargs)
        
        new.load_dorothea(levels = levels, make_df = make_df)
        
        return new
    
    
    def load_transcription(
            self,
            dorothea = True,
            original_resources = True,
            dorothea_levels = None,
            exclude = None,
            reread = False,
            redownload = False,
            **kwargs,
        ):
        
        make_df = kwargs.pop('make_df', None)
        
        if dorothea:
            
            self.load_dorothea(
                levels = dorothea_levels,
                reread = reread,
                redownload = redownload,
            )
        
        if original_resources:
            
            transcription = (
                original_resources
                    if not isinstance(original_resources, bool) else
                network_resources.transcription_onebyone
            )
            
            self.load(
                resources = transcription,
                reread = reread,
                redownload = redownload,
                exclude = exclude,
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
            **kwargs,
        ):
        """
        Initializes a new ``Network`` object with loading a transcriptional
        regulation network from all databases by default.
        
        **kwargs: passed to ``Network.__init__``.
        """
        
        load_args = locals()
        kwargs = load_args.pop('kwargs')
        cls = load_args.pop('cls')
        
        new = cls(**kwargs)
        
        new.load_transcription(**load_args)
        
        return new
    
    
    def load_mirna_target(self, **kwargs):
        
        if 'resources' not in kwargs:
            
            kwargs['resources'] = network_resources.mirna_target
        
        self.load(**kwargs)
    
    
    @classmethod
    def mirna_target(
            cls,
            resources = None,
            make_df = None,
            reread = False,
            redownload = False,
            exclude = None,
            **kwargs,
        ):
        """
        Initializes a new ``Network`` object with loading a miRNA-mRNA
        regulation network from all databases by default.
        
        **kwargs: passed to ``Network.__init__``.
        """
        
        new = cls(**kwargs)
        
        new.mirna_target(
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
            direction = None,
            effect = None,
            resources = None,
            interaction_type = None,
            data_model = None,
            via = None,
            references = None,
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
            not isinstance(entity, common.basestring) and
            hasattr(entity, '__iter__')
        ):
            
            kwargs = locals()
            _ = kwargs.pop('self')
            _ = kwargs.pop('entity')
            
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
            
            @functools.wraps(method_args)
            def _partners_method(*args, **kwargs):
                
                self = args[0]
                kwargs.update(method_args)
                
                method = 'count_partners' if count else 'partners'
                
                return getattr(self, method)(*args[1:], **kwargs)
            
            return _partners_method
        
        for name_parts, arg_parts in (
            zip(*param)
            for param in
            itertools.product(
                *(iteritems(variety) for variety in  cls._partners_methods)
            )
        ):
            
            for count in (False, True):
                
                method_args = dict(
                    itertools.chain(
                        *(iteritems(part) for part in arg_parts)
                    )
                )
                method_name = ''.join(name_parts)
                method_args['count'] = count
                
                setattr(
                    cls,
                    'count_%s' % method_name if count else method_name,
                    _create_partners_method(method_args),
                )
    
    
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
    
    
    @staticmethod
    def _get_by_method_name(get, by):
        
        return (
            ''.join(
                (
                    'get_' if not by else '',
                    get,
                    '_by_' if by else '',
                    by if by else '',
                )
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
        
        
        for _get, _by in itertools.product(
            interaction_mod.Interaction._get_methods,
            interaction_mod.Interaction._by_methods + (None,),
        ):
            
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
        
        
        for _get, _by in itertools.product(
            interaction_mod.Interaction._get_methods,
            interaction_mod.Interaction._by_methods + (None,),
        ):
            
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


Network._generate_get_methods()
Network._generate_partners_methods()
Network._generate_count_methods()
