#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Exports tables for webservice.
#
#  Copyright (c) 2014-2019 - EMBL
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
from past.builtins import xrange, range

import os
import sys
import importlib as imp
import re
import json
import copy
import pandas as pd
import itertools

import pypath.progress as progress
import pypath.urls as urls
import pypath.data_formats as data_formats
import pypath.settings as settings
import pypath.entity as entity

strip_json = re.compile(r'[\[\]{}\"]')
simple_types = {bool, int, float, type(None)}


class Export(object):

    default_header_uniquepairs = ['UniProt_A', 'GeneSymbol_A', 'UniProt_B',
                                  'GeneSymbol_B', 'Databases', 'PubMed_IDs',
                                  'Undirected', 'Direction_A-B',
                                  'Direction_B-A', 'Stimulatory_A-B',
                                  'Inhibitory_A-B', 'Stimulatory_B-A',
                                  'Inhibitory_B-A', 'Category']
    
    default_dtypes_uniquepairs = {
        'UniProt_A': 'category',
        'GeneSymbol_A': 'category',
        'UniProt_B': 'category',
        'GeneSymbol_B': 'category',
        'Undirected': 'category',
        'Direction_A-B': 'category',
        'Direction_B-A': 'category',
        'Stimulatory_A-B': 'category',
        'Inhibitory_A-B': 'category',
        'Stimulatory_B-A': 'category',
        'Inhibitory_B-A': 'category',
        'Category': 'category',
    }

    default_header_bydirs = ['source', 'target', 'source_genesymbol',
                             'target_genesymbol', 'is_directed',
                             'is_stimulation', 'is_inhibition',
                             'consensus_direction',
                             'consensus_stimulation',
                             'consensus_inhibition',
                             'sources',
                             'references', 'dip_url']
    
    default_dtypes_bydirs = {
        'source': 'category',
        'target': 'category',
        'source_genesymbol': 'category',
        'target_genesymbol': 'category',
        'is_directed': 'int8',
        'is_stimulation': 'int8',
        'is_inhibition': 'int8',
        'consensus_direction': 'int8',
        'consensus_stimulation': 'int8',
        'consensus_inhibition': 'int8',
        'sources': 'category',
        'references': 'category',
    }

    def __init__(
            self,
            pa,
            only_sources = None,
            extra_node_attrs = None,
            extra_edge_attrs = None,
            outfile = None,
            default_vertex_attr_processor = None,
            default_edge_attr_processor = None,
        ):

        self.extra_node_attrs = extra_node_attrs or {}
        self.extra_edge_attrs = extra_edge_attrs or {}
        self.outfile = outfile
        self.pa    = pa
        self.graph = pa._get_undirected()
        self.only_sources = only_sources

        if isinstance(self.only_sources, list):

            self.only_sources = set(self.only_sources)

        self.default_vertex_attr_processor = (
            default_vertex_attr_processor or
            self.default_vertex_attr_processor
        )
        self.default_edge_attr_processor = (
            default_edge_attr_processor or
            self.default_edge_attr_processor
        )

    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def make_df(
            self,
            unique_pairs = True,
            extra_node_attrs = None,
            extra_edge_attrs = None,
        ):
        """
        Creates a data frame from the network.

        By default UniProt IDs, Gene Symbols, source databases, literature
        references, directionality and sign information and interaction type
        are included.

        Args:
        -----
        :param bool unique_pairs:
            If `True` each line corresponds to a unique pair of molecules,
            all directionality and sign information are covered in other
            columns. If `False`, order of `A` and `B` IDs corresponds to
            the direction while sign covered in further columns.
        :param dict extra_node_attrs:
            Additional node attributes to be included in the exported table.
            Keys are column ames used in the header while values are names
            of vertex attributes. Values also might be methods which then
            will be called then on each vertex. These should return strings
            or their result will be converted to string.
            In the header `_A` and `_B` suffixes will be appended to the
            column names so the values can be assigned to A and B side
            interaction partners.
        :param dict extra_edge_attrs:
            Additional edge attributes to be included in the exported table.
            Keys are column ames used in the header while values are names
            of edge attributes or callables accepting an edge as single
            argument.
        :param str outfile:
            Name of the output file. If `None` a file name
            "netrowk-<session id>.tab" is used.
        """

        result = []

        self.pa.genesymbol_labels()

        self.extra_node_attrs = extra_node_attrs or self.extra_node_attrs
        self.extra_edge_attrs = extra_edge_attrs or self.extra_edge_attrs

        suffix_a = 'A' if unique_pairs else 'source'
        suffix_b = 'B' if unique_pairs else 'target'

        dtypes = (
            self.default_dtypes_uniquepairs
                if unique_pairs else
            self.default_dtypes_bydirs
        )
        header = copy.copy(
            self.default_header_uniquepairs
            if unique_pairs else
            self.default_header_bydirs
        )
        header += self.extra_edge_attrs.keys()
        header += [
            '%s_%s' % (x, suffix_a)
            for x in self.extra_node_attrs.keys()
        ]
        header += [
            '%s_%s' % (x, suffix_b)
            for x in self.extra_node_attrs.keys()
        ]

        prg = progress.Progress(
            total = self.graph.ecount(),
            name = 'Creating table',
            interval=31
        )

        for e in self.graph.es:

            # adding default fields
            lines = (
                self.process_edge_uniquepairs(e)
                    if unique_pairs else
                self.process_edge_bydirection(e)
            )

            result.extend(lines)

            prg.step()

        prg.terminate()

        self.df = pd.DataFrame(result, columns = header)
        self.df = self.df.astype(dtypes)

    def process_edge_uniquepairs(self, e):
        """
        Returns a table row representing a network edge with covering all
        annotations in a single row: directionality represented by fields
        like `Direction_A-B` and `Direction_B-A` effect sign a similar way.

        Args:
        -----
        :param igraph.Edge e:
            An edge from a pypath igraph object.
        """

        if self.only_sources and not e['sources'] & self.only_sources:

            return []

        vertex_a = self.graph.vs[e.source]
        vertex_b = self.graph.vs[e.target]
        
        name_a, label_a = entity.Entity.igraph_vertex_name_label(vertex_a)
        name_b, label_b = entity.Entity.igraph_vertex_name_label(vertex_b)

        return [
            list(itertools.chain(
                (
                    # uniprots, genesymbols
                        name_a,
                        label_a,
                        name_b,
                        label_b,
                    # sources, references
                        ';'.join(list(e['sources'])),
                        ';'.join(map(lambda r: r.pmid, e['references'])),
                    # directions
                        ';'.join(e['dirs'].get_dir(
                            'undirected', sources=True)),
                        ';'.join(e['dirs'].get_dir(
                            (name_a, name_b), sources=True)),
                        ';'.join(e['dirs'].get_dir(
                            (name_b, name_a), sources=True))
                ),
                (
                    # signs
                    ';'.join(a) for a in
                    itertools.chain(
                        e['dirs'].get_sign((name_a, name_b), sources=True),
                        e['dirs'].get_sign((name_b, name_a), sources=True),
                    )
                ),
                (
                    # category
                    ';'.join(e['type']),
                )
            ))
        ]

    def process_edge_bydirection(self, e):
        """
        Returns one or more table rows representing a network edge a way
        that opposite direction connections contained in separate rows.
        Directionality and sign information covered in 3 columns:
        `is_directed`, `is_inhibition` and `is_stimulation`.
        This is the row format used in the webservice.

        Args:
        -----
        :param igraph.Edge e:
            An edge from a pypath igraph object.
        """

        lines = []
        
        consensus_edges = set(map(tuple, e['dirs'].consensus_edges()))
        consensus_dir = set(c[:3] for c in consensus_edges)

        for d in ['straight', 'reverse']:

            uniprots = getattr(e['dirs'], d)

            if e['dirs'].dirs[uniprots]:
                
                is_stimulation = int(e['dirs'].is_stimulation(uniprots))
                is_inhibition = int(e['dirs'].is_inhibition(uniprots))

                this_edge = [
                    entity.Entity.entity_name_str(uniprots[0]),
                    entity.Entity.entity_name_str(uniprots[1]),
                    self.pa.name_to_label(uniprots[0]),
                    self.pa.name_to_label(uniprots[1]),
                    1, # is_directed
                    is_stimulation,
                    is_inhibition,
                    int(
                        (uniprots[0], uniprots[1], 'directed')
                        in consensus_dir
                    ),
                    int(
                        (uniprots[0], uniprots[1], 'directed', 'positive')
                        in consensus_edges
                    ),
                    int(
                        (uniprots[0], uniprots[1], 'directed', 'negative')
                        in consensus_edges
                    ),
                ]

                dsources = (
                    e['dirs'].get_dir(uniprots, sources=True) |
                    e['dirs'].get_dir('undirected', sources=True)
                )

                if self.only_sources and not dsources & self.only_sources:

                    continue

                this_edge.extend([
                    ';'.join(sorted(dsources)),
                    ';'.join(
                        r.pmid
                        for r in itertools.chain(*(
                            rs for s, rs in iteritems(e['refs_by_source'])
                            if s in dsources
                        ))
                    )
                ])

                this_edge.append(self._dip_urls(e))

                this_edge = self.add_extra_fields(e, this_edge, uniprots)

                lines.append(this_edge)

        if not e['dirs'].is_directed():
            
            this_edge = [
                entity.Entity.entity_name_str(e['dirs'].nodes[0]),
                entity.Entity.entity_name_str(e['dirs'].nodes[1]),
                self.pa.name_to_label(e['dirs'].nodes[0]),
                self.pa.name_to_label(e['dirs'].nodes[1]),
                0,
                0,
                0,
                0,
                0,
                0,
                ';'.join(sorted(e['sources'])),
                ';'.join([r.pmid for r in e['references']]),
                self._dip_urls(e)
            ]

            this_edge = self.add_extra_fields(e, this_edge, 'undirected')

            lines.append(this_edge)

        return lines


    def add_extra_fields(self, e, line, dr = None):
        """
        Takes one table row and using the `igraph.Edge` object and the
        direction provided adds the extra node and edge attribute fields
        as they are defined in `extra_node_attrs` and `extra_edge_attrs`.

        Returns the row with extra fields added.

        Args:
        -----
        :param igraph.Edge e:
            One edge.
        :param list line:
            A table row.
        :param tuple,str dr:
            Direction key. A tuple of names (most often UniProt IDs) or
            `undirected`.
        """

        # extra edge attributes on demand of the user
        for k, v in iteritems(self.extra_edge_attrs):

            line.append(
                self.generic_attr_processor(v, e, dr)
                if hasattr(v, '__call__') else
                self.default_edge_attr_processor(
                    e[v]
                    if v in self.graph.es.attributes() else
                    None
                )
            )

        # extra vertex attributes
        for vertex in (self.graph.vs[e.source], self.graph.vs[e.target]):

            for k, v in iteritems(self.extra_node_attrs):

                line.append(
                    self.generic_attr_processor(v, vertex, dr)
                    if hasattr(v, '__call__') else
                    self.default_vertex_attr_processor(
                        vertex[v]
                        if v in self.graph.vs.attributes() else
                        None
                    )
                )

        return line


    @staticmethod
    def default_vertex_attr_processor(vattr):

            return (
                vattr if type(vattr) in simple_types else
                ';'.join([
                    x.strip()
                    for x in strip_json.sub('',
                        json.dumps(
                            list(vattr)
                            if type(vattr) is set
                            else vattr
                        )).split(',')
                ])
            )

    @staticmethod
    def default_edge_attr_processor(eattr):

        return (
            eattr if type(eattr) in simple_types else
            ';'.join([
                x.strip()
                for x in strip_json.sub('', json.dumps(eattr)).split(',')
            ])
        )

    @staticmethod
    def generic_attr_processor(proc, obj, dr = None):
        """
        Wraps the attribute processor to handle unknown number of arguments.

        Not knowing if the attribute processor expects one or two arguments,
        have no better way than try: if calling with 2 arguments fails with
        `TypeError` we call with one argument.
        """
        
        try:

            return proc(obj, dr)

        except TypeError:

            return proc(obj)

    def write_tab(self, outfile = None, **kwargs):
        """
        Writes the data frame into a tab separated file.

        Args:
        -----
        :param **kwargs:
            Forwarded to `make_df()`.
        """

        if not hasattr(self, 'df'):

            self.make_df(**kwargs)

        self.write(outfile = outfile)

    def write(self, outfile = None):

        outfile = outfile or self.outfile or os.path.join(
            self.pa.outdir, 'network-%s.tab' % self.pa.session
        )

        self.df.to_csv(outfile, sep = '\t', index = False)

    def _dip_urls(self, e):

        result = []
        if 'dip_id' in e.attributes():
            for dip_id in e['dip_id']:
                try:
                    result.append(urls.urls['dip']['ik'] %
                        int(dip_id.split('-')[1][:-1]))
                except:

                    sys.stdout.write('Could not find DIP ID: %s\n' % dip_id)

        return ';'.join(result)

    def webservice_interactions_df(self):

        sources_omnipath = set(
            f.name for f in data_formats.omnipath.values()
        )
        sources_extra_directions = settings.get('network_extra_directions')
        sources_kinase_extra = set(
            f.name for f in data_formats.ptm_misc.values()
        )
        sources_ligrec_extra = set(
            f.name for f in data_formats.ligand_receptor.values()
        )
        sources_pathway_extra = set(
            f.name for f in data_formats.pathway_noref.values()
        )
        sources_mirna = set(
            f.name for f in data_formats.mirna_target.values()
        )

        self.make_df(
            unique_pairs = False,
            extra_node_attrs = {
                'ncbi_tax_id': 'ncbi_tax_id'
            },
            extra_edge_attrs = {
                'omnipath': lambda e, d: (
                    (
                        bool(e['dirs'].sources[d] & sources_omnipath) or
                        (
                            bool(
                                e['dirs'].sources['undirected'] &
                                sources_omnipath
                            ) and
                            bool(
                                e['dirs'].sources[d] &
                                sources_extra_directions
                            )
                        )
                    ) and
                    'PPI' in e['type']
                ),
                'kinaseextra': lambda e, d: (
                    bool(e['dirs'].sources[d] & sources_kinase_extra) and
                    'PPI' in e['type']
                ),
                'ligrecextra': lambda e, d: (
                    bool(e['dirs'].sources[d] & sources_ligrec_extra) and
                    'PPI' in e['type']
                ),
                'pathwayextra': lambda e, d: (
                    bool(e['dirs'].sources[d] & sources_pathway_extra) and
                    'PPI' in e['type']
                ),
                'mirnatarget': lambda e, d: (
                    bool(e['dirs'].sources[d] & sources_mirna) and
                    'MTI' in e['type']
                ),
                'tfregulons': lambda e, d: (
                    'TF' in e['sources_by_type'] and bool(
                        e['sources_by_type']['TF'] &
                        e['dirs'].sources[d]
                    )
                ),
                'tfregulons_curated': 'tfregulons_curated',
                'tfregulons_chipseq': 'tfregulons_chipseq',
                'tfregulons_tfbs':    'tfregulons_tfbs',
                'tfregulons_coexp':   'tfregulons_coexp',
                'tfregulons_level': lambda e, d: (
                    ';'.join(sorted(e['tfregulons_level'])) if
                    'tfregulons_level' in e.attributes() and
                    'TF' in e['sources_by_type'] and bool(
                        e['sources_by_type']['TF'] &
                        e['dirs'].sources[d]
                    ) else ''),
                'type': lambda e: e['type'][0]
            }
        )

    @classmethod
    def sources_table(
            cls,
            pa,
            only_sources = None,
            unique_pairs = True,
            extra_edge_attrs = None,
            extra_node_attrs = None,
            outfile = None,
            default_vertex_attr_processor = None,
            default_edge_attr_processor   = None,
        ):
        """
        Creates a data frame which contains a column for each source database
        with binary values showing presence-absence of interactions across
        resources.
        """

        new = cls(
            pa = pa,
            only_sources = only_sources,
            extra_edge_attrs = extra_edge_attrs or {},
            extra_node_attrs = extra_node_attrs or {},
            outfile = outfile,
            default_vertex_attr_processor = default_vertex_attr_processor,
            default_edge_attr_processor   = default_edge_attr_processor,
        )

        new.make_df(unique_pairs = unique_pairs)

        src_attr = 'Databases' if unique_pairs else 'sources'

        src_all = sorted(only_sources or pa.sources)

        src_cols = dict((src, []) for src in src_all)

        for i, row in new.df.iterrows():

            this_row_src = set(row[src_attr].split(';'))

            for src in src_all:

                src_cols[src].append(int(src in this_row_src))

        for src in src_all:

            new.df.insert(
                loc = new.df.columns.get_loc(src_attr),
                column = src,
                value = src_cols[src]
            )

        return new
