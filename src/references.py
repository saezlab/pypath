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
import igraph

# stats and plotting modules #

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
import _sensitivity as sens
from bioigraph.common import *
from bioigraph import dataio
from bioigraph import plot

# constants
omnipath = 'OmniPath'
lab_size = (18, 21)
axis_lab_size = 36
htp_threshold = 50

net = bioigraph.BioGraph(9606)

#net.init_network()
net.init_network(pfile = 'cache/default_network_ltp.pickle')

net.htp_stats()

pubmeds = uniqList(flatList([[r.pmid for r in e['references']] for e in net.graph.es]))
pubmeds = set(pubmeds) - net.htp[htp_threshold]['htrefs']
##
for e in net.graph.es:
    for s, rs in e['refs_by_source'].iteritems():
        for r in rs:
            if not r.pmid.isdigit():
                print '%s\t%s' % (s, r.pmid)

notpmid = [i for i in pubmeds if not i.isdigit()]

#pmdata = dataio.get_pubmeds(pubmeds)

#pickle.dump(pmdata, open('cache/pubmed2.pickle', 'wb'))
pmdata = pickle.load(open('cache/pubmed2.pickle', 'rb'))

points = []
for e in net.graph.es:
    for s, rs in e['refs_by_source'].iteritems():
        for r in rs:
            if r.pmid not in net.htp[htp_threshold]['htrefs'] and r.pmid in pmdata and 'pubdate' in pmdata[r.pmid]:
                points.append((s, r.pmid, int(pmdata[r.pmid]['pubdate'][:4]), 
                    pmdata[r.pmid]['source']))

points = uniqList(points)

points = pd.DataFrame.from_records(points)
points.columns = ['database', 'pmid', 'year', 'journal']

# ## References by Year
ordr = sorted(uniqList(list(points.year)))
fig, ax = plt.subplots()
sns.set(font = 'Helvetica Neue LT Std')
sns.set_context('poster', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0})
ax = sns.countplot(data = points, x = 'year', color = '#007b7f', order = ordr)
sns.set_context('poster', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0})

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(9)

ax.set_ylabel('Number of PubMed IDs')
ax.set_xlabel('Years')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)
for t in ax.xaxis.get_major_ticks():
    t.label.set_fontsize(lab_size[0])

for t in ax.yaxis.get_major_ticks():
    t.label.set_fontsize(lab_size[1])

ax.xaxis.label.set_size(axis_lab_size)
ax.yaxis.label.set_size(axis_lab_size)

fig.tight_layout()
fig.savefig('references-by-year.pdf')

# ## References by Database
ordr = reversed(points.database.value_counts().order().index)
ordr = [omnipath] + list(ordr)

bp = plot.Barplot(x = ordr, y = [len(pubmeds)] + list(points.database.value_counts()), 
    data = None, fname = 'references-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, desc = False, 
    ylab = 'Number of PubMed IDs', color = ['#6EA945'] + ['#007B7F'] * len(points.database.unique()),
    xlab = 'Resources', order = ordr, y_break = (0.3, 0.07))

# ## References by Journal
for i in [50, 100]:
    ordr = list(reversed(points.journal.value_counts().order().index))[:i]
    fig, ax = plt.subplots()
    sns.set(font = 'Helvetica Neue LT Std')
    sns.set_context('poster', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
    ax.bar(range(len(ordr)), 
        np.log10(list(reversed(points.journal.value_counts().order()))[:i]), 
        label = ordr, color = '#007b7f')
    plt.xticks(np.arange(0.35, len(ordr), 1.0), ordr)
    plt.yticks(np.arange(0,5), [str(int(10.0**yt)) for yt in np.arange(0, 5)])
    sns.set_context('poster', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8 * 100 / i)
    
    for t in ax.yaxis.get_major_ticks():
        t.label.set_fontsize(lab_size[1])
    
    ax.xaxis.label.set_size(axis_lab_size)
    ax.yaxis.label.set_size(axis_lab_size * 0.66)
    
    ax.set_ylabel('Number of PubMed IDs (log)')
    ax.set_xlabel('Journals')
    nul = plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
    fig.tight_layout()
    fig.savefig('references-by-journal-%u.pdf' % i)

# ## References by Database & Journal
def brp(data, **kwargs):
    globals()['sample'] = data
    ordr = list(reversed(data.value_counts().order().index))[:10]
    plt.bar(range(len(ordr)), list(reversed(data.value_counts().order()))[:10],
        **kwargs)
    plt.xticks(np.arange(0.35, len(ordr), 1.0), ordr)

sns.set_context('paper', rc = {'patch.linewidth': 0.0, 'axes.labelsize': 18, 'axes.titlesize': 27,
    'xtick.labelsize': 13, 'ytick.labelsize': 15, 'figure.figsize': [11.7, 17.3]})
grid = sns.FacetGrid(points, col="database", hue="database", col_wrap=8, size=5.5, 
    sharex = False, sharey = False, 
    col_order = [j[1] for j in sorted([(i.upper(), i) for i in points.database.unique()], key = lambda x: x[0])])
grid.map(brp, 'journal')
foo = [plt.setp(xax.xaxis.get_majorticklabels(), rotation=90) for xax in grid.axes]
for xax in grid.axes:
    start, end = xax.get_xlim()
    nul = xax.set_title(xax.get_title().replace('database = ',''))

grid.set_ylabels('PubMed IDs')
grid.set_xlabels('Journals')
grid.fig.tight_layout()
grid.fig.savefig('references-by-db-journal.pdf')

# ## References by Database & Year
def brp(data, **kwargs):
    years = [str(y) for y in xrange(min(data), max(data) + 1, 1)]
    globals()['sample'] = data
    data = [str(i) for i in data]
    plt.bar(range(len(years)), [list(data).count(y) if y in data else 0 for y in years],
        **kwargs)
    plt.xticks(np.arange(0.35, len(years), 1.0), years)

sns.set_context('paper', rc = {'patch.linewidth': 0.0, 'axes.labelsize': 18, 'axes.titlesize': 27,
    'xtick.labelsize': 13, 'ytick.labelsize': 15, 'figure.figsize': [11.7, 17.3]})
grid = sns.FacetGrid(points, col="database", hue="database", col_wrap=5, size=3.5, 
    sharex = False, sharey = False, 
    col_order = [j[1] for j in sorted([(i.upper(), i) for i in points.database.unique()], key = lambda x: x[0])])
grid.map(brp, 'year')
foo = [plt.setp(xax.xaxis.get_majorticklabels(), rotation=90) for xax in grid.axes]
for xax in grid.axes:
    start, end = xax.get_xlim()
    newpos = np.arange(start, end, 3)
    vislabels = range(0, len(newpos)*3, 3)
    foo = xax.set_xticklabels([t._text for i, t in \
        enumerate(list(xax.get_xticklabels())) \
        if i in vislabels])
    foo = xax.xaxis.set_ticks([x + 0.5 for x in newpos])
    xax.set_title(xax.get_title().replace('database = ',''))

grid.set_ylabels('PubMed IDs')
grid.set_xlabels('Years')
grid.fig.tight_layout()
grid.fig.savefig('references-by-db-year-new.pdf')

# ## References by Database & Year # boxplot
op_refs = uniqList(flatList([[rr.pmid for rr in e['references']] for e in net.graph.es]))
medpubyr = [(s, np.median(points[(points.database == s)]['year'])) \
    for s in net.sources if s != 'ACSN'] + [(omnipath, np.median([int(pmdata[r]['pubdate'][:4]) \
        for r in op_refs if r in pmdata and 'pubdate' in pmdata[r]]))]
medpubyr = sorted(medpubyr, key = lambda x: x[1])
colors = []

oppoints = []
for r in op_refs:
    if r not in net.htp[htp_threshold]['htrefs'] and r in pmdata and 'pubdate' in pmdata[r]:
                oppoints.append((omnipath, r, int(pmdata[r]['pubdate'][:4]), 
                    pmdata[r]['source']))

oppoints = uniqList(oppoints)
oppoints = pd.DataFrame.from_records(oppoints)
oppoints.columns = ['database', 'pmid', 'year', 'journal']

allpoints = points.append(oppoints, ignore_index = True)

sns.set_context('poster', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0, 'figure.figsize': [12.8, 8.8]})
fig, ax = plt.subplots()
ax = sns.boxplot('database', 'year', data = allpoints, 
    palette = ['#6EA945' if i[0] == omnipath else '#007B7F' for i in  medpubyr], 
    linewidth = 0.1, saturation = 0.66, order = [i[0] for i in medpubyr])
ax.set_xlabel('Resources', weight = 'light', fontsize = 12, 
        variant = 'normal', color = '#000000', stretch = 'normal')
ax.set_ylabel('Year', weight = 'light', fontsize = 12, 
    variant = 'normal', color = '#000000', stretch = 'normal')
plt.setp(ax.xaxis.get_majorticklabels(), rotation = 90)
for t in ax.xaxis.get_major_ticks():
    t.label.set_fontsize(lab_size[0])

for t in ax.yaxis.get_major_ticks():
    t.label.set_fontsize(lab_size[1])

ax.xaxis.label.set_size(axis_lab_size)
ax.yaxis.label.set_size(axis_lab_size)

fig.tight_layout()
fig.savefig('pubyear-boxplot.pdf')
plt.close(fig)

#points.groupby(['database','year']).count()