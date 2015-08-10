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
import heapq
import operator

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
from matplotlib import gridspec

# from bioigraph #

import bioigraph
from bioigraph import chembl
from bioigraph.common import *
from bioigraph.data_formats import best, good, ugly, transcription
import _sensitivity as sens
from bioigraph import progress
from bioigraph import dataio
from bioigraph.ig_drawing import DefaultGraphDrawerFFsupport
from bioigraph import plot

omnipath = 'OmniPath'

net = bioigraph.BioGraph(9606)

#net.init_network(pfile = 'cache/plus_phospho.pickle')
#net.init_network()
#net.save_network(pfile = 'cache/default_network.pickle')

### phospho network + acsn + default :::
#net.load_resources(lst={'mimp': good['mimp']})
#net.load_resources(lst={'pnetworks': good['pnetworks']})
#net.load_resources(lst={'psite_noref': good['psite_noref']})
#net.save_network(pfile = 'cache/default_plus_acsn_phospho.pickle')

#net.get_directed(conv_edges = True, mutual = True)
#net.graph = net.dgraph
#net.save_network(pfile = 'cache/plus_phospho_directed.pickle')
#net.init_network(pfile = 'cache/default_network.pickle')
#net.load_resources(lst={'acsn': ugly['acsn']})
#net.save_network(pfile = 'cache/phospho_plus_acsn.pickle')
#net.save_network(pfile = 'cache/default_plus_acsn.pickle')
net.init_network(pfile = 'cache/default_plus_acsn.pickle')
net.load_resources(lst={'hsn': ugly['hsn']})

net.curation_tab()
net.curation_tab(fname = 'curation_stats_wang.tex')

net.genesymbol_labels()
net.set_tfs()
net.set_receptors()
net.load_corum()
sens.in_complex(net)
net.in_complex()
net.load_ptms()
net.read_list_file(bioigraph.data_formats.cgc)
net.read_list_file(bioigraph.data_formats.intogene_cancer)
# net.load_comppi()
sep = net.separate()

net.lists['rec'] = uniqList(flatList([net.mapper.map_name(rec, 'genesymbol', 'uniprot') \
    for rec in dataio.get_hpmr()]))
net.lists['tfs'] = uniqList(flatList([net.mapper.map_name(tf, 'ensg', 'uniprot') \
    for tf in dataio.get_tfcensus()['ensg']]))
# sep = dict([(s, net.get_network({'edge': {'sources': [s]}, 'node': {}})) for s in net.sources])

# source-vcount sens.barplot
d = zip(*[(s, len([v for v in net.graph.vs if s in v['sources']])) for s in net.sources] + \
    [(omnipath, net.graph.vcount())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'proteins-by-db.pdf', lab_size = 11, 
    ylab = 'Number of proteins', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y',
    y_break = (0.37, 0.15))

# source-ecount sens.barplot
d = zip(*[(s, len([e for e in net.graph.es if s in e['sources']])) for s in net.sources] + \
    [(omnipath, net.graph.ecount())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'interactions-by-db.pdf', lab_size = 11, 
    ylab = 'Number of interactions', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.19, 0.1))

# density sens.barplot
d = zip(*[(s, g.density()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.density())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'density-by-db.pdf', lab_size = 11, 
    ylab = 'Graph density', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.5, 0.1))

# transitivity sens.barplot
d = zip(*[(s, g.transitivity_undirected()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.transitivity_undirected())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'transitivity-by-db.pdf', lab_size = 11, 
    ylab = 'Graph global transitivity', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y')

# diameter sens.barplot
d = zip(*[(s, g.diameter()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.diameter())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'diameter-by-db.pdf', lab_size = 11, 
    ylab = 'Graph diameter', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y')

# receptors sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['rec']])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec']))])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptors-by-db.pdf', lab_size = 11, 
    ylab = 'Number of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', y_break = (0.6, 0.09))

# receptors prop sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['rec']])/float(g.vcount())*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec'])/float(net.graph.vcount())*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorprop-by-db.pdf', lab_size = 11, 
    ylab = 'Percentage of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', y_break = (0.39, 0.12))

# receptors prop sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['rec']])/float(len(net.lists['rec']))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec'])/float(len(net.lists['rec']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorcov-by-db.pdf', lab_size = 11, 
    ylab = 'Percentage of all\nhuman receptors covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y')

# transcription factors sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['tf']])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf']))])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfs-by-db.pdf', lab_size = 11, 
    ylab = 'Number of transcription factors', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y',
    y_break = (0.45, 0.2))

# transcription factor prop sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['tf']])/float(g.vcount())*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf'])/float(net.graph.vcount())*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfprop-by-db.pdf', lab_size = 11, 
    ylab = 'Percentage of transcription factors', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.53, 0.2))

d = zip(*[(s, len([v for v in g.vs if v['tf']])/float(len(net.lists['tfs']))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf'])/float(len(net.lists['tfs']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfcov-by-db.pdf', lab_size = 11, 
    ylab = 'Percentage of all human\ntranscription factors covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.5, 0.1))

# disease associations:
net.load_disgenet()
diss = dataio.get_disgenet()
dis = {}
for d in diss:
    ups = net.mapper.map_name(d['entrez'], 'entrez', 'uniprot')
    for u in ups:
        if u not in dis: dis[u] = []
        dis[u].append(d['disease'])

d = zip(*[(s, sum(len(x) for u, x in dis.iteritems() if u in g.vs['name'])/float(sum(len(x) for x in dis.values()))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(len(x) for x in net.graph.vs['dis'])/float(sum(len(x) for x in dis.values()))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'discov-by-db.pdf', lab_size = 11, 
    ylab = 'Percentage of disease-gene associations covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.53, 0.12))

# stacked sens.barplot example
#d = zip(*[(s, len([v for v in g.vs if v['tf']])/float(g.vcount())*100, 
    #len([v for v in g.vs if v['rec']])/float(g.vcount())*100,
    #len([v for v in g.vs if not v['tf'] and not v['rec']])/float(g.vcount())*100) \
    #for s, g in sep.iteritems()])
#sens.stacked_barplot(x = d[0], y = d[1:], names = ['Transcription factors', 'Receptors', 'Others'], 
    #data = None, fname = 'tfrec-by-db.pdf', lab_size = 11, 
    #ylab = 'Percentage of transcription factors', 
    #xlab = 'Pathway resources', order = 'y')

# complex/ptm/other stacked sens.barplot
d = zip(*[(s, 
    len([e for e in g.es if e['in_complex']]), 
    len([e for e in g.es if len(e['ptm']) > 0]),
    len([e for e in g.es if len(e['ptm']) == 0 and not e['in_complex']])) \
    for s, g in sep.iteritems()] + \
    [tuple([omnipath, 
        sum(net.graph.es['in_complex']), 
        len([e for e in net.graph.es if len(e['ptm']) > 0]),
        len([e for e in net.graph.es if len(e['ptm']) == 0 and not e['in_complex']])])]
    )
sens.stacked_barplot(x = d[0], y = d[1:], names = ['Within complex', 'Having PTM', 'None'], 
    data = None, fname = 'comp-ptm-by-db.pdf', lab_size = 11, 
    ylab = 'Interactions in complex / having known PTM / none of them', 
    xlab = 'Pathway resources', order = 'y')

# number of complexes sens.barplot
d = zip(*[(s, len(sens.complexes_in_network(g))) \
    for s, g in sep.iteritems()] + \
    [(omnipath, len(sens.complexes_in_network(net.graph)))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'complexes-by-db.pdf', lab_size = 11, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of complexes (of a total of %u CORUM complexes)'%\
        len(sens.complexes_in_network(net.graph)), 
    xlab = 'Pathway resources', order = 'y', y_break = (0.55, 0.05))

# number of ptms sens.barplot
d = zip(*[(s, sum([len(e['ptm']) for e in g.es])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum([len(e['ptm']) for e in net.graph.es]))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ptms-by-db.pdf', lab_size = 11, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of PTMs (of a total of %u PTMs)'%\
        sum([len(e['ptm']) for e in net.graph.es]), 
    xlab = 'Pathway resources', order = 'y', y_break = (0.6, 0.1))

# number of ptms sens.barplot
d = zip(*[(s, len([1 for e in g.es if len(e['ptm']) != 0])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, len([1 for e in net.graph.es if len(e['ptm']) != 0]))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'havingptm-by-db.pdf', lab_size = 11, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Interactions having PTM', 
    xlab = 'Pathway resources', order = 'y', y_break = (0.5, 0.15))

# cosmic cancer gene census coverage sens.barplot
d = zip(*[(s, len(set(net.lists['CancerGeneCensus']) & set(g.vs['name'])) / \
        float(len(net.lists['CancerGeneCensus']))*100) \
        for s, g in sep.iteritems()] + \
    [(omnipath, len(set(net.lists['CancerGeneCensus']) & set(net.graph.vs['name'])) / \
        float(len(net.lists['CancerGeneCensus']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ccgccov-by-db.pdf', lab_size = 11, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Percentage of\nCosmic Cancer Gene Census genes covered', 
    xlab = 'Pathway resources', order = 'y')

# intogene cancer driver coverage sens.barplot
d = zip(*[(s, len(set(net.lists['Intogene']) & set(g.vs['name'])) / \
        float(len(net.lists['Intogene']))*100) \
        for s, g in sep.iteritems()] + \
    [(omnipath, len(set(net.lists['Intogene']) & set(net.graph.vs['name'])) / \
        float(len(net.lists['Intogene']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov-by-db.pdf', lab_size = 11, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Percentage of\nIntogene cancer driver genes covered', 
    xlab = 'Pathway resources', order = 'y')

# dirs signes w/o omnipath
d = zip(*[(s, 
    sum([sum([s in e['dirs'].positive_sources[e['dirs'].straight],
        s in e['dirs'].positive_sources[e['dirs'].reverse]]) for e in g.es]),
    sum([sum([s in e['dirs'].negative_sources[e['dirs'].straight],
        s in e['dirs'].negative_sources[e['dirs'].reverse]]) for e in g.es]),
    sum([sum([s in ((e['dirs'].sources[e['dirs'].straight] - \
        e['dirs'].positive_sources[e['dirs'].straight]) - \
        e['dirs'].negative_sources[e['dirs'].straight]),
        s in ((e['dirs'].sources[e['dirs'].reverse] - \
        e['dirs'].positive_sources[e['dirs'].reverse]) - \
        e['dirs'].negative_sources[e['dirs'].reverse])]) for e in g.es]),
    sum([s in e['dirs'].sources['undirected'] for e in g.es])
        ) \
    for s, g in sep.iteritems()])
sens.stacked_barplot(x = d[0], y = d[1:], 
    names = ['positive', 'negative', 'unknown effect', 'unknown direction'],
    data = None, fname = 'dirs-signes-by-db-wo-op.pdf', lab_size = 11, 
    ylab = 'Interactions:\npositive/negative/directed unknown effect/undirected', 
    xlab = 'Pathway resources', order = 'y')

# number of complexes sens.barplot
d = zip(*[(s, 
    sum([sum([s in e['dirs'].positive_sources[e['dirs'].straight],
        s in e['dirs'].positive_sources[e['dirs'].reverse]]) for e in g.es]),
    sum([sum([s in e['dirs'].negative_sources[e['dirs'].straight],
        s in e['dirs'].negative_sources[e['dirs'].reverse]]) for e in g.es]),
    sum([sum([s in ((e['dirs'].sources[e['dirs'].straight] - \
        e['dirs'].positive_sources[e['dirs'].straight]) - \
        e['dirs'].negative_sources[e['dirs'].straight]),
        s in ((e['dirs'].sources[e['dirs'].reverse] - \
        e['dirs'].positive_sources[e['dirs'].reverse]) - \
        e['dirs'].negative_sources[e['dirs'].reverse])]) for e in g.es]),
    sum([s in e['dirs'].sources['undirected'] for e in g.es])
        ) \
    for s, g in sep.iteritems()] + \
        [tuple([omnipath, 
            sum([sum([e['dirs'].is_stimulation(e['dirs'].straight), 
                e['dirs'].is_stimulation(e['dirs'].reverse)]) for e in net.graph.es]),
            sum([sum([e['dirs'].is_inhibition(e['dirs'].straight), 
                e['dirs'].is_inhibition(e['dirs'].reverse)]) for e in net.graph.es]),
            sum([sum([e['dirs'].straight and not \
                e['dirs'].is_stimulation(e['dirs'].straight) and not \
                e['dirs'].is_inhibition(e['dirs'].straight),
                e['dirs'].reverse and not \
                e['dirs'].is_stimulation(e['dirs'].reverse) and not \
                e['dirs'].is_inhibition(e['dirs'].reverse)]) for e in net.graph.es]),
            sum([1 for e in net.graph.es if not e['dirs'].is_directed()])])])
sens.stacked_barplot(x = d[0], y = d[1:], 
    names = ['positive', 'negative', 'unknown effect', 'unknown direction'],
    data = None, fname = 'dirs-signes-by-db.pdf', lab_size = 11, 
    ylab = 'Interactions:\npositive/negative/directed unknown effect/undirected', 
    xlab = 'Pathway resources', order = 'y')

# comppi locations
comppi_locs = uniqList(flatList([x.keys() for x in net.graph.vs['comppi']]))
for v in net.graph.vs:
    for loc in comppi_locs:
        if loc not in v['comppi']:
            v['comppi'][loc] = 0.0

points = []
for v in net.graph.vs:
    for s in v['sources']:
        for l, c in v['comppi'].iteritems():
            points.append((s, l, c))

points = pd.DataFrame.from_records(points)
points.columns = ['database', 'localization', 'comppi_score']

sns.set_context('paper', rc = {'patch.linewidth': 0.0, 'axes.labelsize': 9, 
    'xtick.labelsize': 5, 'ytick.labelsize': 6, 'figure.figsize': [11.7, 8.3]})
grid = sns.FacetGrid(points, col="database", hue="database", col_wrap=8, size=2.5, 
    sharex = False, sharey = False)
grid.map(sns.sens.barplot, 'localization', 'comppi_score', estimator = np.mean, ci = None)
foo = [plt.setp(xax.xaxis.get_majorticklabels(), rotation=90) for xax in grid.axes]
for xax in grid.axes:
    start, end = xax.get_xlim()
    xax.set_title(xax.get_title().replace('database = ',''))

grid.set_ylabels('Cumulative ComPPI-score')
grid.set_xlabels('Localizations')
grid.fig.tight_layout()
grid.fig.savefig('proteins-locations.pdf')

# edge comppi score boxplot:
points = []
for e in net.graph.es:
    if e['comppi'] is not None:
        for s in e['sources']:
            points.append((s, e['comppi']))

points = pd.DataFrame.from_records(points)
points.columns = ['database', 'comppi_score']

fig, ax = plt.subplots()
sns.set(font = 'Helvetica Neue LT Std')
sns.boxplot(points['comppi_score'], groupby = points['database'], 
    color = sens.embl_colors, linewidth = 0.1, saturation = 0.66)

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(8)
    tick.label.set_color('black')

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(11)
    tick.label.set_color('black')

fig.savefig('comppi-edge-by-db.pdf')

# edge comppi sens.barplots:
points = []
for e in net.graph.es:
    for s in e['sources']:
        for l in e['loc']:
            points.append((s, l))

points = pd.DataFrame.from_records(points)
points.columns = ['database', 'localization']

def brp(data, **kwargs):
    globals()['sample'] = data
    ordr = list(reversed(data.value_counts().order().index))[:10]
    plt.bar(range(len(ordr)), list(reversed(data.value_counts().order()))[:10],
        **kwargs)
    plt.xticks(np.arange(0.35, len(ordr), 1.0), ordr)


sns.set_context('paper', rc = {'patch.linewidth': 0.0, 'axes.labelsize': 9, 
    'xtick.labelsize': 5, 'ytick.labelsize': 6, 'figure.figsize': [11.7, 8.3]})
grid = sns.FacetGrid(points, col="database", hue="database", col_wrap=8, size=2.5, 
    sharex = False, sharey = False)
grid.map(brp, 'localization')
foo = [plt.setp(xax.xaxis.get_majorticklabels(), rotation=90) for xax in grid.axes]
for xax in grid.axes:
    start, end = xax.get_xlim()
    xax.set_title(xax.get_title().replace('database = ',''))

grid.set_ylabels('Number of edges')
grid.set_xlabels('Localization')
grid.fig.tight_layout()
grid.fig.savefig('comppi-edge-by-db2.pdf')

# loc for vertices:
net.graph.vs['loc'] = [[] for _ in net.graph.vs]
for v in net.graph.vs:
    v['loc'] = list(set([i[0] for i in \
        heapq.nlargest(2 , v['comppi'].iteritems(), 
        operator.itemgetter(1))
    ]))

points = []
for v in net.graph.vs:
    for s in v['sources']:
        for l in v['loc']:
            points.append((s, l))

points = pd.DataFrame.from_records(points)
points.columns = ['database', 'localization']

def brp(data, **kwargs):
    globals()['sample'] = data
    ordr = list(reversed(data.value_counts().order().index))[:10]
    plt.bar(range(len(ordr)), list(reversed(data.value_counts().order()))[:10],
        **kwargs)
    plt.xticks(np.arange(0.35, len(ordr), 1.0), ordr)

sns.set_context('paper', rc = {'patch.linewidth': 0.0, 'axes.labelsize': 9, 
    'xtick.labelsize': 5, 'ytick.labelsize': 6, 'figure.figsize': [11.7, 8.3]})
grid = sns.FacetGrid(points, col="database", hue="database", col_wrap=8, size=2.5, 
    sharex = False, sharey = False)
grid.map(brp, 'localization')
foo = [plt.setp(xax.xaxis.get_majorticklabels(), rotation=90) for xax in grid.axes]
for xax in grid.axes:
    start, end = xax.get_xlim()
    xax.set_title(xax.get_title().replace('database = ',''))

grid.set_ylabels('Number of proteins')
grid.set_xlabels('Localization')
grid.fig.tight_layout()
grid.fig.savefig('comppi-node-by-db2.pdf')

# scatterplots:
topdata = {
'Proteins': [len([v for v in net.graph.vs if s in v['sources']]) for s in net.sources],
'Interactions': [len([e for e in net.graph.es if s in e['sources']]) for s in net.sources],
'Density': [sep[s].density() for s in net.sources],
'Transitivity': [sep[s].transitivity_undirected() for s in net.sources],
'Diameter': [sep[s].diameter() for s in net.sources]
}
topdf = pd.DataFrame(topdata, index = net.sources)
fig, ax = plt.subplots()
g = sns.pairplot(topdf)
#fig.tight_layout()
g.savefig('topology_pairplot.pdf')
plt.close(fig)

fig, ax = plt.subplots()
pc = plt.scatter(list(topdf['Density']), list(topdf['Diameter']), 
    s = [np.pi * (50.0 * x/float(max(topdf['Proteins'])))**2 \
        for x in list(topdf['Proteins'])], c = '#6ea945', alpha = 0.5, 
        edgecolors = 'none')
plt.xlabel('Density')
plt.ylabel('Diameter')
plt.axis('tight')
plt.xlim([-0.01, 0.07])
#ax.autoscale_view()
fig.savefig('dens-diam-vcount.pdf')
plt.close(fig)

fig, ax = plt.subplots()
pc = plt.scatter(list(topdf['Density']), list(topdf['Transitivity']), 
    s = [np.pi * (50.0 * x/float(max(topdf['Proteins'])))**2 \
        for x in list(topdf['Proteins'])], c = '#6ea945', alpha = 0.5, 
        edgecolors = 'none')
plt.xlabel('Density')
plt.ylabel('Transitivity')
#plt.xlim([-0.01, 0.07])
#ax.set_xscale('log')
#ax.set_yscale('log')
plt.axis('tight')
#ax.autoscale_view()
fig.savefig('dens-trans-vcount.pdf')
plt.close(fig)