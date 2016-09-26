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

# generic modules #

import sys
import os
import re
import cPickle as pickle
import copy
import igraph
import louvain
import cairo
import heapq
import operator
import textwrap
from collections import Counter, OrderedDict

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
from bioigraph import dataio
from bioigraph import go, gsea
from bioigraph.ig_drawing import DefaultGraphDrawerFFsupport


net = bioigraph.BioGraph(9606)
net.init_network(pfile = 'cache/default_plus_acsn.pickle')
net.genesymbol_labels()

psrc = dict((s, ([v['name'] for v in net.graph.vs if s in v['sources']])) for s in net.sources)

g = gsea.GSEA(mapper = net.mapper)
# hallmarks
g.load_collection('H')
# cancer modules
g.load_collection('CM')
# curated
g.load_collection('C2')
# oncogenic
g.load_collection('C6')

enr = gsea.GSEABinaryEnrichmentSet(basic_set = net.graph.vs['name'], gsea = g)

enr_db = {}
for s, uniprots in psrc.iteritems():
    enr.new_set(uniprots)
    enr_db[s] = enr.top_genesets(length = 100, min_set_size = 15, groups = ['H', 'C2'])

topcountr = Counter(flatList([x[:10] for x in enr_db.values()]))
topcountrall = Counter(flatList([x for x in enr_db.values()]))

mostcommon = OrderedDict(sorted(dict(topcountr).items(), 
    key = lambda i: i[1], reverse = True)[:20])
mostcomterm = set(mostcommon.keys())
mostspecial = dict([(db, [t for t in lst if topcountrall[t] == 1][:10]) \
    for db, lst in enr_db.iteritems()])

edges = []
for term in mostcommon.keys():
    for db, lst in enr_db.iteritems():
        if term in lst[:10]:
            edges.append((db, term))

for db, terms in mostspecial.iteritems():
    for term in terms[:3]:
        edges.append((db, term))

gg = igraph.Graph.TupleList(edges)

repl = {
    'signaling pathway': 'signaling',
    'Genes involved in ': '',
    'Toll-like receptor': 'TLR',
}
rerepl = re.compile('|'.join(repl.keys()))
regid = re.compile(r'[\[\(].*[\]\)]')

gg.vs['typ'] = ['gsetc' if v['name'] in mostcomterm \
    else 'gsets' if v['name'] not in net.sources \
    else 'db' for v in gg.vs]
gg.vs['color'] = ['#6EA945AA' if v['typ'] == 'db' \
    else '#007B7FAA' if v['typ'] == 'gsetc' \
        else '#FCCC06AA' for v in gg.vs]
gg.vs['shape'] = ['circle' if v['typ'] == 'db' else 'rectangle' for v in gg.vs]
gg.vs['label'] = [v['name'] if v['typ'] == 'db' else \
    regid.sub('', rerepl.sub(lambda x: repl[x.group()], v['name'])).split(';')[0] \
    for v in gg.vs]
gg.vs['label'] = [textwrap.wrap(v['label'], 16) \
    if v['typ'] != 'db' else v['label'] for v in gg.vs]

gg.vs['label_size'] = [14 if v['typ'] == 'db' else 7 for v in gg.vs]

width_unit = 5
height_unit = 9

gg.vs['size'] = [20 if v['typ'] == 'db' else 45 for v in gg.vs]
gg.vs['height'] = [len(v['label']) * height_unit \
    if v['typ'] != 'db' else 20 for v in gg.vs]
gg.vs['width'] = [max((len(x) for x in v['label'])) * width_unit \
    if v['typ'] != 'db' else 20 for v in gg.vs]

gg.vs['label'] = ['\n'.join(v['label']) \
    if type(v['label']) is list else v['label'] for v in gg.vs]

gg.vs['label_color'] = ['#FFFFFFFF' if v['typ'] == 'gsetc' else '#777777FF' for v in gg.vs]
gg.vs['label_family'] = ['Sentinel Book' \
    if v['typ'] == 'db' else 'HelveticaNeueLT Std Lt' for v in gg.vs]

gglo = gg.layout_fruchterman_reingold(repulserad = gg.vcount() ** 2.98, 
    maxiter = 1000, area = gg.vcount() ** 2.3)

sf = cairo.PDFSurface('gsea_enrichment_graph.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)

plot = igraph.plot(gg, vertex_label = gg.vs['label'],
    layout = gglo, 
    bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_frame_width = 0, 
    edge_color = '#007B7F55',
    edge_curved = False)

plot.redraw()
plot.save()

### GO ##
net.load_go()
enr = net.go_enrichment()
enr_db = {}
for s, uniprots in psrc.iteritems():
    enr.new_set(uniprots)
    enr_db[s] = enr.top_names(length = 100, min_set_size = 15)

topcountr = Counter(flatList([x[:10] for x in enr_db.values()]))
topcountrall = Counter(flatList([x for x in enr_db.values()]))

mostcommon = OrderedDict(sorted(dict(topcountr).items(), 
    key = lambda i: i[1], reverse = True)[:20])
mostcomterm = set(mostcommon.keys())
mostspecial = dict([(db, [t for t in lst if topcountrall[t] == 1][:10]) \
    for db, lst in enr_db.iteritems()])


topcountr = Counter(flatList([x[:10] for x in enr_db.values()]))
topcountrall = Counter(flatList([x for x in enr_db.values()]))

mostcommon = OrderedDict(sorted(dict(topcountr).items(), 
    key = lambda i: i[1], reverse = True)[:20])
mostcomterm = set(mostcommon.keys())
mostspecial = dict([(db, [t for t in lst if topcountrall[t] == 1][:10]) \
    for db, lst in enr_db.iteritems()])

edges = []
for term in mostcommon.keys():
    for db, lst in enr_db.iteritems():
        if term in lst[:10]:
            edges.append((db, term))

for db, terms in mostspecial.iteritems():
    for term in terms[:3]:
        edges.append((db, term))

gg = igraph.Graph.TupleList(edges)

repl = {
    'signaling pathway': 'signaling',
    'Genes involved in ': '',
    'toll-like receptor': 'TLR',
    'epidermal growth factor receptor': 'EGFR',
    'peptidyl-tyrosine': 'tyrosine',
    'signal transduction': 'signaling',
    'protein kinase A': 'PKA'
}
rerepl = re.compile('|'.join(repl.keys()))
regid = re.compile(r'[\[\(].*[\]\)]')

gg.vs['typ'] = ['gsetc' if v['name'] in mostcomterm \
    else 'gsets' if v['name'] not in net.sources \
    else 'db' for v in gg.vs]
gg.vs['color'] = ['#6EA945AA' if v['typ'] == 'db' \
    else '#007B7FAA' if v['typ'] == 'gsetc' \
        else '#FCCC06AA' for v in gg.vs]
gg.vs['shape'] = ['circle' if v['typ'] == 'db' else 'rectangle' for v in gg.vs]
gg.vs['label'] = [v['name'] if v['typ'] == 'db' else \
    regid.sub('', rerepl.sub(lambda x: repl[x.group()], v['name'])).split(';')[0] \
    for v in gg.vs]
gg.vs['label'] = [v['label'] if v['typ'] == 'db' else v['label'][0].upper() + v['label'][1:] for v in gg.vs]
gg.vs['label'] = [textwrap.wrap(v['label'], 16) \
    if v['typ'] != 'db' else v['label'] for v in gg.vs]

gg.vs['label_size'] = [14 if v['typ'] == 'db' else 7 for v in gg.vs]

width_unit = 5
height_unit = 9

gg.vs['size'] = [20 if v['typ'] == 'db' else 45 for v in gg.vs]
gg.vs['height'] = [len(v['label']) * height_unit \
    if v['typ'] != 'db' else 20 for v in gg.vs]
gg.vs['width'] = [max((len(x) for x in v['label'])) * width_unit \
    if v['typ'] != 'db' else 20 for v in gg.vs]

gg.vs['label'] = ['\n'.join(v['label']) \
    if type(v['label']) is list else v['label'] for v in gg.vs]

gg.vs['label_color'] = ['#FFFFFFFF' if v['typ'] == 'gsetc' else '#777777FF' for v in gg.vs]
gg.vs['label_family'] = ['Sentinel Book' \
    if v['typ'] == 'db' else 'HelveticaNeueLT Std Lt' for v in gg.vs]

gglo = gg.layout_fruchterman_reingold(repulserad = gg.vcount() ** 2.98, 
    maxiter = 1000, area = gg.vcount() ** 2.3)

sf = cairo.PDFSurface('go_enrichment_graph.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)

plot = igraph.plot(gg, vertex_label = gg.vs['label'],
    layout = gglo, 
    bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_frame_width = 0, 
    edge_color = '#007B7F55',
    edge_curved = False)

plot.redraw()
plot.save()
