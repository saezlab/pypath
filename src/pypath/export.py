#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import print_function

from future.utils import iteritems
from past.builtins import xrange, range

import os
import sys
import imp
import re
import json
import copy
import pandas as pd
import itertools

import pypath.progress as progress
import pypath.urls as urls
import pypath.data_formats as data_formats

strip_json = re.compile(r'[\[\]{}\"]')
simple_types = {bool, int, float, type(None)}


class Export(object):
    
    default_header_uniquepairs = [
        'UniProt_A', 'GeneSymbol_A', 'UniProt_B', 'GeneSymbol_B',
        'Databases', 'PubMed_IDs', 'Undirected', 'Direction_A-B',
        'Direction_B-A', 'Stimulatory_A-B', 'Inhibitory_A-B',
        'Stimulatory_B-A', 'Inhibitory_B-A', 'Category'
    ]
    
    default_header_bydirs = [
        'source', 'target', 'source_genesymbol', 'target_genesymbol',
        'is_directed', 'is_stimulation', 'is_inhibition', 'sources',
        'references', 'dip_url'
    ]
    
    def __init__(
            self,
            pa,
            extra_node_attrs = {},
            extra_edge_attrs = {},
            outfile = None,
            default_vertex_attr_processor = None,
            default_edge_attr_processor   = None
        ):
        
        self.extra_node_attrs = extra_node_attrs
        self.extra_edge_attrs = extra_edge_attrs
        self.outfile = outfile
        self.pa    = pa
        self.graph = pa._get_undirected()
        
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
            extra_node_attrs={},
            extra_edge_attrs={},
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
        
        header = copy.copy(
            self.default_header_uniquepairs
            if unique_pairs else
            self.default_header_bydirs
        )
        header += extra_edge_attrs.keys()
        header += ['%s_%s' % (x, suffix_a) for x in extra_node_attrs.keys()]
        header += ['%s_%s' % (x, suffix_b) for x in extra_node_attrs.keys()]
        
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
        
        name_a = self.graph.vs[e.source]['name']
        name_b = self.graph.vs[e.target]['name']
        
        return [
            list(itertools.chain(
                (
                    # uniprots, genesymbols
                        name_a.replace(' ', ''),
                        self.graph.vs[e.source]['label'].replace(' ', ''),
                        name_b.replace(' ', ''),
                        self.graph.vs[e.target]['label'].replace(' ', ''),
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
                            e['dirs'].get_sign((name_b, name_a), sources=True)
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
        
        for d in ['straight', 'reverse']:
            
            uniprots = getattr(e['dirs'], d)
            
            if e['dirs'].dirs[uniprots]:
                
                this_edge = [
                    uniprots[0],
                    uniprots[1],
                    self.pa.nodLab[self.pa.nodDct[uniprots[0]]],
                    self.pa.nodLab[self.pa.nodDct[uniprots[1]]],
                    1, # is_directed
                    int(e['dirs'].is_stimulation(uniprots)),
                    int(e['dirs'].is_inhibition(uniprots))
                ]
                
                dsources = (
                    e['dirs'].get_dir(uniprots, sources=True) |
                    e['dirs'].get_dir('undirected', sources=True)
                )
                
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
                e['dirs'].nodes[0],
                e['dirs'].nodes[1],
                self.pa.nodLab[self.pa.nodDct[e['dirs'].nodes[0]]],
                self.pa.nodLab[self.pa.nodDct[e['dirs'].nodes[1]]],
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
        
        self.make_df(**kwargs)
        
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
        sources_kinase_extra = set(
            f.name for f in data_formats.ptm_misc.values()
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
                    bool(e['dirs'].sources[d] & sources_omnipath) and
                    'PPI' in e['type']
                ),
                'kinaseextra': lambda e, d: (
                    bool(e['dirs'].sources[d] & sources_kinase_extra) and
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
