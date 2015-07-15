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

net = bioigraph.BioGraph(9606)

#net.init_network()
net.init_network(pfile = 'cache/default_network.pickle')

pubmeds = uniqList(flatList([[r.pmid for r in e['references']] for e in net.graph.es]))

pmdata = dataio.get_pubmeds(pubmeds)

pickle.dump(pmdata, open('cache/pubmed2.pickle', 'wb'))
points = []
for e in net.graph.es:
    for s, rs in e['refs_by_source'].iteritems():
        for r in rs:
            if r.pmid in pmdata and 'pubdate' in pmdata[r.pmid]:
                points.append((s, r.pmid, int(pmdata[r.pmid]['pubdate'][:4]), 
                    pmdata[r.pmid]['source']))

points = uniqList(points)

points = pd.DataFrame.from_records(points)
points.columns = ['database', 'pmid', 'year', 'journal']



# fig = points['year'].hist(by = points['database'])

fig, ax = plt.subplots()
sns.set(font = 'Helvetica Neue LT Std')
sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0})
# sns.axes_style({'lines.linewidth': 0.0})
ax = sns.barplot(data = points, x = 'year', color = '#007b7f')
sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0})

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(9)

ax.set_ylabel('Number of PubMed IDs')
ax.set_xlabel('Years')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)
fig.tight_layout()
fig.savefig('references-by-year.pdf')

# ## #

ordr = points.database.value_counts().order().index
fig, ax = plt.subplots()
sns.set(font = 'Helvetica Neue LT Std')
sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0})
# sns.axes_style({'lines.linewidth': 0.0})
ax = sns.barplot(data = points, x = 'database', color = '#007b7f', x_order = ordr)
sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0})

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(9)

ax.set_ylabel('Number of PubMed IDs')
ax.set_xlabel('Resources')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)
fig.tight_layout()
fig.savefig('references-by-db.pdf')

# ## #

ordr = list(reversed(points.journal.value_counts().order().index))[:100]
fig, ax = plt.subplots()
sns.set(font = 'Helvetica Neue LT Std')
sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0})
# sns.axes_style({'lines.linewidth': 0.0})
ax.bar(range(len(ordr)), 
    np.log10(list(reversed(points.journal.value_counts().order()))[:100]), 
    label = ordr, color = '#007b7f')
#sns.barplot(x = ordr,
    #y = list(reversed(points.journal.value_counts().order()))[:50], 
    #color = '#007b7f', label = ordr)
plt.xticks(np.arange(0.35, len(ordr), 1.0), ordr)
plt.yticks(np.arange(0,5), [str(int(10.0**yt)) for yt in np.arange(0, 5)])
sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0})

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(6)

ax.set_ylabel('Number of PubMed IDs (log)')
ax.set_xlabel('Journals')
plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
fig.tight_layout()
fig.savefig('references-by-journal.pdf')

# ## #



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
grid.map(brp, 'journal')
foo = [plt.setp(xax.xaxis.get_majorticklabels(), rotation=90) for xax in grid.axes]
for xax in grid.axes:
    start, end = xax.get_xlim()
    xax.set_title(xax.get_title().replace('database = ',''))

grid.set_ylabels('PubMed IDs')
grid.set_xlabels('Journals')
grid.fig.tight_layout()
grid.fig.savefig('references-by-db-journal.pdf')

# ## #


sns.set_context('paper', rc = {'patch.linewidth': 0.0, 'axes.labelsize': 9, 
    'xtick.labelsize': 4})
grid = sns.FacetGrid(points, col="database", hue="database", col_wrap=5, size=1.5, 
    sharex = False, sharey = False)
grid.map(sns.barplot, 'year')
foo = [plt.setp(xax.xaxis.get_majorticklabels(), rotation=90) for xax in grid.axes]
for xax in grid.axes:
    start, end = xax.get_xlim()
    newpos = np.arange(start, end, 3)
    vislabels = range(0, len(newpos)*3, 3)
    print vislabels
    print [t._text for i, t in \
        enumerate(list(xax.get_xticklabels())) \
        if i in vislabels]
    foo = xax.set_xticklabels([t._text for i, t in \
        enumerate(list(xax.get_xticklabels())) \
        if i in vislabels])
    foo = xax.xaxis.set_ticks([x + 0.5 for x in newpos])
    xax.set_title(xax.get_title().replace('database = ',''))

grid.set_ylabels('PubMed IDs')
grid.set_xlabels('Years')
grid.fig.tight_layout()
grid.fig.savefig('references-by-db-year.pdf')

#points.groupby(['database','year']).count()