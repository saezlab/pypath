#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2018 - EMBL
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
# this module makes possible
# dynamic data integration, downloads
# files from various resources, in standard
# or non-standard text based and xml formats,
# processes them, sometimes parses html
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

import pypath.progress as progress
import pypath.urls as urls


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
            outfile = None
        ):
        
        self.extra_node_attrs = extra_node_attrs
        self.extra_edge_attrs = extra_edge_attrs
        self.outfile = outfile
        self.pa    = pa
        self.graph = pa._get_undirected()
    
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
        
        self.genesymbol_labels()
        
        extra_node_attrs = extra_node_attrs or self.extra_node_attrs
        extra_edge_attrs = extra_edge_attrs or self.extra_edge_attrs
        
        suffix_a = 'A' if unique_pairs else 'source'
        suffix_b = 'B' if unique_pairs else 'target'
        
        header = copy.copy(
            self.default_header_uniquepairs
            if unique_pairs else
            self.default_header_bydirs
        )
        header += extra_edge_attrs.keys()
        header += ['%s_%s' (x, suffix_a) for x in extra_node_attrs.keys()]
        header += ['%s_%s' (x, suffix_b) for x in extra_node_attrs.keys()]
        
        stripJson = re.compile(r'[\[\]{}\"]')
        
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
            
            # extra edge attributes on demand of the user
            for k, v in iteritems(extra_edge_attrs):
                
                value = (
                    v(e)
                    if hasattr(v, '__call__') else
                    ';'.join([
                        x.strip()
                        for x in stripJson.sub('', json.dumps(e[v])).split(',')
                    ])
                )
                
                for l in lines:
                    
                    l.append(value)
            
            # extra vertex attributes
            for vertex in (self.graph.vs[e.source], self.graph.vs[e.target]):
                
                for k, v in iteritems(extra_node_attrs):
                    
                    value = (
                        v(vertex)
                        if hasattr(v, '__call__') else
                        ';'.join([
                            x.strip()
                            for x in stripJson.sub('',
                                json.dumps(
                                    list(vertex[v])
                                    if type(vertex[v]) is set
                                    else vertex[v]
                                )).split(',')
                        ])
                    )
                    
                    for l in lines:
                        
                        l.append(value)
            
            result.extend(lines)
            
            prg.step()
        
        prg.terminate()
        
        self.df = pd.DataFrame(result, columns = header, index = False)
    
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
            [
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
                    (name_b, name_a), sources=True)),
                ';'.join(a)
                for a in e['dirs'].get_sign(
                    (name_a, name_b), sources=True) + e['dirs'].get_sign(
                        (name_b, name_a), sources=True),
                ';'.join(e['type'])
            ]
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
                    self.g.vs[self.p.nodDct[uniprots[0]]]['label'],
                    self.g.vs[self.p.nodDct[uniprots[1]]]['label'],
                    1, # is_directed
                    int(e['dirs'].is_stimulation(uniprots)),
                    int(e['dirs'].is_inhibition(uniprots))
                ]
                
                dsources = (
                    e['dirs'].get_dir(uniprots, sources=True) |
                    e['dirs'].get_dir('undirected', sources=True)
                )
                
                this_edge.extend([
                    sorted(dsources),
                    [
                        r.pmid
                        for r in flatList([
                            rs for s, rs in iteritems(e['refs_by_source'])
                            if s in dsources
                        ])
                    ]
                ])
                
                this_edge.append(self._dip_urls(e))
                
                lines.append(this_edge)
        
        if not e['dirs'].is_directed():
            
            this_edge = [
                e['dirs'].nodes[0],
                e['dirs'].nodes[1],
                self.g.vs[self.p.nodDct[e['dirs'].nodes[0]]]['label'],
                self.g.vs[self.p.nodDct[e['dirs'].nodes[1]]]['label'],
                0,
                0,
                0,
                this_edge.append(list(e['sources'])),
                [r.pmid for r in e['references']],
                self._dip_urls(e)
            ]
            
            lines.append(this_edge)
        
        return lines
    
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
            self.outdir, 'network-%s.tab' % self.pa.session
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
