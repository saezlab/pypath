#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
    #After satisfying identification of ChEMBL IDs for 
    #compounds in GDSC screening, this script looks up
    #in the network their nominal tarets and other known 
    #targets form ChEMBL binding assays.
#

# generic modules #

import sys
import os
import re
import cPickle as pickle
import copy
import igraph
import louvain
import cairo

# stats and plotting modules #

import math
import ranking
import numpy as np
from numpy.random import randn
import pandas as pd
from scipy import stats
from scipy.misc import comb
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as hc
import hcluster as hc2
import matplotlib.patches as mpatches

# from bioigraph #

import bioigraph
from bioigraph import chembl
from bioigraph.common import *
from bioigraph.data_formats import best, good, ugly, transcription
import _sensitivity as sens
from bioigraph import progress
from bioigraph.ig_drawing import DefaultGraphDrawerFFsupport

net = bioigraph.BioGraph(9606)

#net.init_network(pfile = 'cache/default_network.pickle')
#net.init_network(pfile = 'cache/plus_phospho.pickle')
#net.load_resources(lst={'mimp': good['mimp']})
#net.load_resources(lst={'pnetworks': good['pnetworks']})
#net.load_resources(lst={'psite_noref': good['psite_noref']})
#net.save_network(pfile = 'cache/plus_phospho.pickle')
#net.get_directed(conv_edges = True, mutual = True)
#net.graph = net.dgraph
#net.save_network(pfile = 'cache/plus_phospho_directed.pickle')
#net.init_network(pfile = 'cache/default_network.pickle')
#net.load_resources(lst={'acsn': ugly['acsn']})
#net.save_network(pfile = 'cache/default_plus_acsn.pickle')
#net.init_network(pfile = 'cache/default_plus_acsn.pickle')

net.genesymbol_labels()

sim = net.databases_similarity()

edges = [e for e in [(it[0], iit[0], iit[1], sim['nodes'][it[0]][iit[0]]) \
    for it in sim['edges'].items() for iit in it[1].items()] \
    if (e[2] > 0.15 or (e[0] == 'MatrixDB' and e[2] > 0.02) or \
    (e[0] == 'NRF2ome' and e[2] > 0.07) or (e[0] == 'ACSN' and e[2] > 0.10)) \
    and e[0] != e[1]]

g = igraph.Graph.TupleList(edges, edge_attrs = ['weight', 'nodes'])
g.simplify(combine_edges = 'mean')


lo = g.layout_fruchterman_reingold(weights = 'weight', repulserad = g.vcount() ** 2.8, 
    maxiter = 1000, area = g.vcount() ** 2.3)

g.vs['ncount'] = [len([e for e in net.graph.es if v['name'] in e['sources']])**0.55 for v in g.vs]

sf = cairo.PDFSurface('resource_graph_edge_simpson.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)

plot = igraph.plot(g, vertex_label = g.vs['name'],
    layout = lo, 
    bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_size = g.vs['ncount'], 
    vertex_frame_width = 0, vertex_color = '#6EA945', 
    vertex_label_color = '#777777FF',  vertex_label_family = 'Sentinel Book',
    vertex_label_size = 24,  vertex_label_dist = 1.4, 
    edge_width = map(lambda x: (x * 10.0)**1.8, g.es['weight']), 
    edge_color = '#007B7F55',
    edge_curved = False)

plot.redraw()
plot.save()

redges = [(s1, s2, bioigraph.common.simpson_index(
    [r.pmid for r in uniqList(flatList([[] if s1 not in e['refs_by_source'] \
        else e['refs_by_source'][s1] \
        for e in net.graph.es]))], 
    [r.pmid for r in uniqList(flatList([[] if s2 not in e['refs_by_source'] \
        else e['refs_by_source'][s2] \
        for e in net.graph.es]))]
    )) for s1 in list(set(g.vs['name']) - set(['ACSN'])) \
        for s2 in list(set(g.vs['name']) - set(['ACSN']))]


f = igraph.Graph.TupleList([e for e in redges if e[2] > 0.0545 and e[0] != e[1]], 
    edge_attrs = ['weight'])
f.simplify(combine_edges = 'mean')


flo = f.layout_fruchterman_reingold(weights = 'weight', repulserad = f.vcount() ** 2.8, 
    maxiter = 1000, area = f.vcount() ** 2.3)

f.vs['ncount'] = [len(uniqList(flatList([e['refs_by_source'][v['name']] \
    for e in net.graph.es \
    if v['name'] in e['refs_by_source']])))**0.48 for v in f.vs]

sf = cairo.PDFSurface('resource_graph_refs_simpson.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)

plot = igraph.plot(f, vertex_label = f.vs['name'],
    layout = flo, 
    bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_size = f.vs['ncount'], 
    vertex_frame_width = 0, vertex_color = '#6EA945', 
    vertex_label_color = '#777777FF',  vertex_label_family = 'Sentinel Book',
    vertex_label_size = 24,  vertex_label_dist = 1.4, 
    edge_width = map(lambda x: (x * 10.0)**1.8, f.es['weight']), 
    edge_color = '#007B7F55',
    edge_curved = False)

plot.redraw()
plot.save()

# ### #

cedges = [(s1, s2, sum([0.0 if s1 not in e['refs_by_source'] or \
    s2 not in e['refs_by_source'] \
    else len(set([r1.pmid for r1 in e['refs_by_source'][s1]]).\
    symmetric_difference(set([r2.pmid for r2 in e['refs_by_source'][s2]]))) \
    for e in net.graph.es]) / \
    float(len([e for e in net.graph.es \
        if s1 in e['sources'] and s2 in e['sources']]) + 0.001)) 
    for s1 in net.sources for s2 in net.sources]

len([c for c in cedges if c[2] > 4])

c = igraph.Graph.TupleList([e for e in cedges if e[2] > 5.9 and e[0] != e[1]], 
    edge_attrs = ['weight'])
c.simplify(combine_edges = 'mean')


citeffort = []
for v in c.vs:
    allrefs = len(uniqList(flatList([[r.pmid for r in e1['refs_by_source'][v['name']]] \
        for e1 in net.graph.es \
        if v['name'] in e1['refs_by_source']])))
    alledges = float(len([e2.index for e2 in net.graph.es if v['name'] in e2['sources']]))
    uniqcits = sum([len([rr.pmid for rr in e3['refs_by_source'][v['name']]]) \
    for e3 in net.graph.es \
    if v['name'] in e3['refs_by_source']])
    citeffort.append(allrefs/alledges*uniqcits)

c.vs['ncount'] = citeffort


clo = c.layout_fruchterman_reingold(weights = 'weight', repulserad = c.vcount() ** 2.98, 
    maxiter = 1000, area = c.vcount() ** 2.3)

sf = cairo.PDFSurface('resource_graph_curation.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)

plot = igraph.plot(c, vertex_label = c.vs['name'],
    layout = clo, 
    bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_size = map(lambda x: (x*3)**0.4, c.vs['ncount']), 
    vertex_frame_width = 0, vertex_color = '#6EA945', 
    vertex_label_color = '#777777FF',  vertex_label_family = 'Sentinel Book',
    vertex_label_size = 24,  vertex_label_dist = 1.4, 
    edge_width = map(lambda x: (x * 0.35)**1.3, c.es['weight']), 
    edge_color = '#007B7F55',
    edge_curved = False)

plot.redraw()
plot.save()

## TikZ summary figure ##


tikzfname = 'resource_tree_tikz.tex'

tikz = r'''\documentclass[a4paper,10pt]{article}
    \usepackage{fontspec}
    \usepackage{xunicode}
    \usepackage{polyglossia}
    \setdefaultlanguage{english}
    \usepackage{xltxtra}
    \usepackage{microtype}
    \usepackage{fullpage}
    \usepackage[usenames,dvipsnames,svgnames,table]{xcolor}
    \usepackage{color}
    \setmainfont{HelveticaNeueLTStd-Lt}
    \usepackage{tikz}
    \definecolor{zircon}{RGB}{228, 236, 236}
    \definecolor{teal}{RGB}{0, 123, 127}
    \definecolor{twilightblue}{RGB}{239, 244, 233}
    \definecolor{mantis}{RGB}{110, 169, 69}
    \begin{document}
    \thispagestyle{empty}
    \begin{tikzpicture}
'''

for i in xrange(11):
    tikz += r'''        \fill[anchor = south west, fill = twilightblue] '''\
        r'''(0.0, %f) rectangle (17.0, %f);
        \node[anchor = west] at (0.33, %f) {\LARGE{\color{teal}%u}};
''' % (float(i), i + 0.96, i + 0.50, 2005 + i)

tikz += r'''\end{tikzpicture}
    \end{document}
'''

with open(tikzfname, 'w') as f:
    f.write(tikz)

