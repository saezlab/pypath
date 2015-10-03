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
import subprocess
import cPickle as pickle
import copy
import igraph
import louvain
import cairo
from datetime import date

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
from bioigraph import descriptions
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
net.init_network(pfile = 'cache/default_network_wo-intact_ltp-only.pickle')

net.genesymbol_labels()

## ## ##

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

scale = [50, 100, 500, 1000, 5000]
escale = [0.05, 0.1, 0.2, 0.5]

g.add_vertices([str(i) for i in scale])
g.add_vertices(['%.2f_%u' % (i, a) for i in escale for a in [0, 1]])
g.add_edges([('%.2f_%u' % (i, 0), '%.2f_%u' % (i, 1)) for i in escale])

xmax = max([c[0] for c in lo._coords])
ymin = min([c[1] for c in lo._coords])
xrng = xmax - min([c[0] for c in lo._coords])
yrng = max([c[1] for c in lo._coords]) - ymin
xleg = xmax + xrng * 0.2

for i, s in enumerate(scale):
    v = g.vs[g.vs['name'].index(str(s))]
    v['ncount'] = s**0.55
    lo._coords.append([xleg, ymin + i*1.25 + sum(scale[:i + 1])**0.55 * 0.065])

g.es['label'] = ['' for _ in g.es]

for i, s in enumerate(escale):
    v1 = g.vs[g.vs['name'].index('%.2f_%u' % (s, 0))]
    v2 = g.vs[g.vs['name'].index('%.2f_%u' % (s, 1))]
    e = g.es[g.get_eid(v1.index, v2.index)]
    e['weight'] = s
    e['label'] = '%.2f' % s
    ycoo =  ymin + yrng * 0.7 + i * 1.8
    lo._coords.append([xleg - xrng * 0.07, ycoo])
    lo._coords.append([xleg + xrng * 0.07, ycoo])
    v1['ncount'] = 0.0
    v2['ncount'] = 0.0
    v1['name'] = ''
    v2['name'] = ''

sf = cairo.PDFSurface('resource_graph_edge_simpson.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)

plot = igraph.plot(g, vertex_label = g.vs['name'],
    layout = lo, 
    bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_size = g.vs['ncount'], 
    vertex_frame_width = 0, vertex_color = '#6EA945', 
    vertex_label_color = '#777777FF', vertex_label_family = 'Sentinel Book',
    edge_label_color = '#777777FF', edge_label_family = 'Sentinel Book',
    vertex_label_size = 24,  vertex_label_dist = 1.4, 
    edge_label_size = 24, 
    edge_label = g.es['label'], 
    edge_width = map(lambda x: (x * 10.0)**1.8, g.es['weight']), 
    edge_color = '#007B7F55',
    edge_curved = False)

plot.redraw()
plot.save()

## ## ##

redges = [(s1, s2, bioigraph.common.simpson_index(
    [r.pmid for r in uniqList(flatList([[] if s1 not in e['refs_by_source'] \
        else e['refs_by_source'][s1] \
        for e in net.graph.es]))], 
    [r.pmid for r in uniqList(flatList([[] if s2 not in e['refs_by_source'] \
        else e['refs_by_source'][s2] \
        for e in net.graph.es]))]
    )) for s1 in list(set(net.sources) - set(['ACSN'])) \
        for s2 in list(set(net.sources) - set(['ACSN']))]


g = igraph.Graph.TupleList([e for e in redges if e[2] > 0.0545 and e[0] != e[1]], 
    edge_attrs = ['weight'])
g.simplify(combine_edges = 'mean')


lo = g.layout_fruchterman_reingold(weights = 'weight', repulserad = g.vcount() ** 2.8, 
    maxiter = 1000, area = g.vcount() ** 2.3)

g.vs['ncount'] = [len(uniqList(flatList([e['refs_by_source'][v['name']] \
    for e in net.graph.es \
    if v['name'] in e['refs_by_source']])))**0.48 for v in g.vs]

scale = [50, 100, 500, 1000, 2000]
escale = [0.05, 0.1, 0.2, 0.5]

g.add_vertices([str(i) for i in scale])
g.add_vertices(['%.2f_%u' % (i, a) for i in escale for a in [0, 1]])
g.add_edges([('%.2f_%u' % (i, 0), '%.2f_%u' % (i, 1)) for i in escale])

xmax = max([c[0] for c in lo._coords])
ymin = min([c[1] for c in lo._coords])
xrng = xmax - min([c[0] for c in lo._coords])
yrng = max([c[1] for c in lo._coords]) - ymin
xleg = xmax + xrng * 0.2

for i, s in enumerate(scale):
    v = g.vs[g.vs['name'].index(str(s))]
    v['ncount'] = s**0.48
    lo._coords.append([xleg, ymin + i*1.25 + sum(scale[:i + 1])**0.48 * 0.065])

g.es['label'] = ['' for _ in g.es]

for i, s in enumerate(escale):
    v1 = g.vs[g.vs['name'].index('%.2f_%u' % (s, 0))]
    v2 = g.vs[g.vs['name'].index('%.2f_%u' % (s, 1))]
    e = g.es[g.get_eid(v1.index, v2.index)]
    e['weight'] = s
    e['label'] = '%.2f' % s
    ycoo =  ymin + yrng * 0.7 + i * 1.8
    lo._coords.append([xleg - xrng * 0.07, ycoo])
    lo._coords.append([xleg + xrng * 0.07, ycoo])
    v1['ncount'] = 0.0
    v2['ncount'] = 0.0
    v1['name'] = ''
    v2['name'] = ''

sf = cairo.PDFSurface('resource_graph_refs_simpson-2.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)

plot = igraph.plot(g, vertex_label = g.vs['name'],
    layout = lo, 
    bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_size = g.vs['ncount'], 
    vertex_frame_width = 0, vertex_color = '#6EA945', 
    vertex_label_color = '#777777FF',  vertex_label_family = 'Sentinel Book',
    edge_label_color = '#777777FF', edge_label_family = 'Sentinel Book',
    vertex_label_size = 24,  vertex_label_dist = 1.4, 
    edge_label_size = 24, 
    edge_label = g.es['label'], 
    edge_width = map(lambda x: (x * 10.0)**1.8, g.es['weight']), 
    edge_color = '#007B7F55',
    edge_curved = False)

plot.redraw()
plot.save()

## ## ##

cedges = [(s1, s2, sum([0.0 if s1 not in e['refs_by_source'] or \
    s2 not in e['refs_by_source'] \
    else len(set([r1.pmid for r1 in e['refs_by_source'][s1]]).\
    symmetric_difference(set([r2.pmid for r2 in e['refs_by_source'][s2]]))) \
    for e in net.graph.es]) / \
    float(len([e for e in net.graph.es \
        if s1 in e['sources'] and s2 in e['sources']]) + 0.001)) 
    for s1 in net.sources for s2 in net.sources]

len([c for c in cedges if c[2] > 4])

g = igraph.Graph.TupleList([e for e in cedges if e[2] > 5.9 and e[0] != e[1]], 
    edge_attrs = ['weight'])
g.simplify(combine_edges = 'mean')


citeffort = []
for v in g.vs:
    allrefs = len(uniqList(flatList([[r.pmid for r in e1['refs_by_source'][v['name']]] \
        for e1 in net.graph.es \
        if v['name'] in e1['refs_by_source']])))
    alledges = float(len([e2.index for e2 in net.graph.es if v['name'] in e2['sources']]))
    uniqcits = sum([len([rr.pmid for rr in e3['refs_by_source'][v['name']]]) \
    for e3 in net.graph.es \
    if v['name'] in e3['refs_by_source']])
    citeffort.append(allrefs/alledges*uniqcits)

g.vs['ncount'] = citeffort

lo = g.layout_fruchterman_reingold(weights = 'weight', repulserad = g.vcount() ** 2.98, 
    maxiter = 1000, area = g.vcount() ** 2.3)

scale = [100, 1000, 5000, 10000, 20000]
escale = [5.0, 7.5, 10.0, 15.0, 30.0]

g.add_vertices([str(i) for i in scale])
g.add_vertices(['%.2f_%u' % (i, a) for i in escale for a in [0, 1]])
g.add_edges([('%.2f_%u' % (i, 0), '%.2f_%u' % (i, 1)) for i in escale])

xmax = max([c[0] for c in lo._coords])
ymin = min([c[1] for c in lo._coords])
xrng = xmax - min([c[0] for c in lo._coords])
yrng = max([c[1] for c in lo._coords]) - ymin
xleg = xmax + xrng * 0.2

for i, s in enumerate(scale):
    v = g.vs[g.vs['name'].index(str(s))]
    v['ncount'] = s
    lo._coords.append([xleg, ymin + i*0.22 + (sum(scale[:i + 1])*3)**0.55 * 0.016])

g.es['label'] = ['' for _ in g.es]

for i, s in enumerate(escale):
    v1 = g.vs[g.vs['name'].index('%.2f_%u' % (s, 0))]
    v2 = g.vs[g.vs['name'].index('%.2f_%u' % (s, 1))]
    e = g.es[g.get_eid(v1.index, v2.index)]
    e['weight'] = s
    e['label'] = '%.2f' % s
    ycoo =  ymin + yrng * 0.7 + i * 1.8
    lo._coords.append([xleg - xrng * 0.07, ycoo])
    lo._coords.append([xleg + xrng * 0.07, ycoo])
    v1['ncount'] = 0.0
    v2['ncount'] = 0.0
    v1['name'] = ''
    v2['name'] = ''

sf = cairo.PDFSurface('resource_graph_curation-2.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)

plot = igraph.plot(g, vertex_label = g.vs['name'],
    layout = lo, 
    bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_size = map(lambda x: (x*3)**0.4, g.vs['ncount']), 
    vertex_frame_width = 0, vertex_color = '#6EA945', 
    vertex_label_color = '#777777FF',  vertex_label_family = 'Sentinel Book',
    edge_label_color = '#777777FF', edge_label_family = 'Sentinel Book',
    vertex_label_size = 24,  vertex_label_dist = 1.4, 
    edge_label_size = 24, 
    edge_label = g.es['label'], 
    edge_width = map(lambda x: (x * 0.35)**1.3, g.es['weight']), 
    edge_color = '#007B7F55',
    edge_curved = False)

plot.redraw()
plot.save()

## TikZ summary figure ##

# the figure is based on the bioigraph.descriptions dict:
import subprocess
from datetime import date
from bioigraph.descriptions import *
from bioigraph.common import *
d = descriptions

# parameters of the figure
tikzfname = 'resource_tree_tikz.tex'
latex = '/usr/bin/xelatex'
firstyear = min(flatList([[r['year']] for r in d.values() if 'year' in r] + \
    [r['releases'] for r in d.values() if 'releases' in r]))

lastyear = date.today().year
years = range(lastyear - firstyear + 1)
yearbarwidth = 0.4
sepwidth = 0.04
lineheight = [1.0] * len(years)
labelbg = 'twilightblue'
labelfg = 'teal'
nodelabbg = 'teal'
nodelabfg = 'white'
dotcol = 'teal'
linecol = 'teal'
dataimportcol = 'mantis'
dotsize = 3.0
linewidth = 1.0
rowbg = 'twilightblue'
width = 17.0
xoffset = 0.5
dotlineopacity = 0.7
horizontal = True # whether the timeline should be the horizontal axis

# TikZ styles
tikzstyles = r'''
    \tikzstyle{omnipath}=[rectangle, anchor = center, inner sep = 2pt, fill = %s, 
        rotate = 90, text = %s, draw = %s]
    \tikzstyle{others}=[rectangle, anchor = center, inner sep = 2pt, fill = %s, 
        rotate = 90, text = %s, draw = %s]
''' % (
        nodelabfg, 
        nodelabbg,
        nodelabbg,
        nodelabbg,
        nodelabfg,
        nodelabbg
    )

# LaTeX preamble for XeLaTeX
tikz = r'''\documentclass[a4paper,10pt]{article}
    \usepackage{fontspec}
    \usepackage{xunicode}
    \usepackage{polyglossia}
    \setdefaultlanguage{english}
    \usepackage{xltxtra}
    \usepackage{microtype}
    \usepackage[cm]{fullpage}
    \usepackage{rotating}
    \usepackage[usenames,dvipsnames,svgnames,table]{xcolor}
    \usepackage{color}
    \setmainfont{HelveticaNeueLTStd-Lt}
    \usepackage{tikz}
    \definecolor{zircon}{RGB}{228, 236, 236}
    \definecolor{teal}{RGB}{0, 123, 127}
    \definecolor{twilightblue}{RGB}{239, 244, 233}
    \definecolor{mantis}{RGB}{110, 169, 69}%s
    \begin{document}
    \thispagestyle{empty}
    \pgfdeclarelayer{background}
    \pgfdeclarelayer{nodes}
    \pgfdeclarelayer{lines}
    \pgfsetlayers{background,lines,nodes}%s
    \begin{tikzpicture}
    \begin{pgfonlayer}{background}
''' % (tikzstyles, 
    r'''
    \begin{turn}{-90}''' if horizontal else '')

ordr = sorted([(lab, r['releases'] if 'releases' in r else [] + \
            [r['year']] if 'year' in r else [], r['label'] if 'label' in r else lab, 
            r['data_import'] if 'data_import' in r else [],
            'omnipath' if 'omnipath' in r and r['omnipath'] else 'others') \
        for lab, r in d.iteritems() \
        if 'year' in r or 'releases' in r], \
    key = lambda x: min(x[1]))

# the background grid and year labels
for i in years:
    tikz += r'''        \fill[anchor = south west, fill = %s, 
        inner sep = 0pt, outer sep = 0pt] '''\
        r'''(%f, %f) rectangle (%f, %f);
        \node[anchor = north west, rotate = 90, text width = %fcm, 
            fill = %s, inner sep = 0pt, outer sep = 0pt, align = center, 
            minimum size = %fcm] at (0.0, %f) {\small{\color{%s}%u}};
        \fill[fill = red] (%f, %f) circle (0.0pt);
''' % (
    rowbg, # background of row
    yearbarwidth + sepwidth, # left edge of row
    sum(lineheight[:i]), # top edge of row
    width, # right edge of the row
    sum(lineheight[:i + 1]) - sepwidth, # bottom edge of row
    lineheight[i] - sepwidth, # height of year label
    labelbg, # background of label
    yearbarwidth, # width of year label
    sum(lineheight[:i]), # top of year label
    labelfg, # text color of label
    firstyear + i, # year
    0.0, sum(lineheight[:i]) # red dot
    )

# new layer for nodes
tikz += r'''    \end{pgfonlayer}
    \begin{pgfonlayer}{nodes}
    '''

# horizontal distance between vertical columns
xdist = (width - yearbarwidth - sepwidth - xoffset) / float(len(ordr))
nodelabels = []

# drawing vertical dots, labels and connecting lines:
for i, r in enumerate(ordr):
    coox = xdist * i + yearbarwidth + sepwidth + xdist / 2.0 + xoffset
    ymax = max(r[1])
    ydots = [y for y in r[1] if y != ymax]
    ylaby = ymax - firstyear
    cooylab = sum(lineheight[:ylaby]) + lineheight[ylaby] / 2.0
    ydots = [sum(lineheight[:y - firstyear]) + lineheight[y - firstyear] / 2.0 for y in ydots]
    for j, cooy in enumerate(ydots):
        tikz += r'''        \node[circle, fill = %s, minimum size = %f, opacity = %f] 
            (%s) at (%f, %f) {};
        ''' % (
            dotcol, # fill color for dot
            dotsize, # size of the dot
            dotlineopacity, # opacity of dot
            '%s%u' % (r[0].lower(), j), # label
            coox, # x coordinate
            cooy # y coordinate
        )
    tikz += r'''        \node[%s] 
            (%s) at (%f, %f) 
            {\footnotesize %s};
    ''' % (
        r[4], # node style
        r[0].lower(), # node name
        coox, # node x coordinate
        cooylab, # node y coordinate
        r[0] # label text
    )
    nodelabels.append(r[0].lower())
    if len(r[1]) > 1:
        tikz += r'''        \draw[draw = %s, line width = %fpt, opacity = %f] (%s%s);
        ''' % (
            linecol, 
            linewidth,
            dotlineopacity, 
            '%s) -- (' % r[0].lower(), 
            ') -- ('.join( \
                ['%s%u' % (r[0].lower(), j) for j in xrange(len(ydots))])
        )

# new layer for crossing lines showing data transfers:
tikz += r'''\end{pgfonlayer}
    \begin{pgfonlayer}{lines}
    '''

# drawing data transfer lines:
for r in ordr:
    for s in r[3]:
        if r[0].lower() in nodelabels and s.lower() in nodelabels:
            tikz += r'''        \draw[-latex, draw = %s, line width = %fpt, opacity = %f] 
                (%s) -- (%s);
            ''' % (
                dataimportcol,
                linewidth,
                dotlineopacity,
                s.lower(),
                r[0].lower()
            )

# closing layer, tikzpicture and LaTeX document:
tikz += r'''    \end{pgfonlayer}
    \end{tikzpicture}%s
    \end{document}
''' % (r'''
    \end{turn}''' if horizontal else '')

# writing to file:
with open(tikzfname, 'w') as f:
    f.write(tikz)

# compiling with XeLaTeX:
subprocess.call([latex, tikzfname])
