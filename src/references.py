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
import matplotlib.gridspec as gridspec
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

net.read_list_file(bioigraph.data_formats.cgc)
net.read_list_file(bioigraph.data_formats.intogene_cancer)

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
    
    if i == 50:
        ax.set_ylim(1.8, 4)
    
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
    'xtick.labelsize': 13, 'ytick.labelsize': 15, 'figure.figsize': [9.7, 17.3]})
grid = sns.FacetGrid(points, col="database", hue="database", col_wrap=5, size=2.5, 
    sharex = False, sharey = False, aspect = 1.33, 
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
grid.fig.savefig('references-by-db-year.pdf')

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

sns.set(font = 'Helvetica Neue LT Std')
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

# ## References by Database & Year & number of refs # boxplot & barplot
# this works if the boxplot code above has been run:
boxplot_ordr = [l._text for l in ax.get_xticklabels()]

fig = plt.figure()
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 8], width_ratios = [1])
#, wspace = 0.05, hspace = 0.05)
ax0 = plt.subplot(gs[0])
# here comes the top barplot
refc_by_db = points.database.value_counts()
refc = [refc_by_db[s] if s != omnipath else len(pubmeds) for s in boxplot_ordr]
barplot = ax0.bar(range(len(refc)), refc, 
    color = ['#6EA945' if s == omnipath else '#007B7F' for s in boxplot_ordr])
tick_loc = np.arange(len(boxplot_ordr)) + 0.3
plt.xticks(tick_loc)
ax0.set_xticklabels(boxplot_ordr, rotation = 90, fontsize = lab_size[0] * 0.66)
ax0.set_yscale('log')
ax0.set_xlim([-0.2, 27.8])
ax0.set_xlabel('')
ax0.set_ylabel('Number of\nPubMed IDs')

# here the boxplot
ax1 = plt.subplot(gs[1])
sns.set(font = 'Helvetica Neue LT Std')
sns.set_context('poster', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0, 'figure.figsize': [12.8, 8.8]})
ax1 = sns.boxplot('database', 'year', data = allpoints, 
    palette = ['#6EA945' if i[0] == omnipath else '#007B7F' for i in  medpubyr], 
    linewidth = 0.1, saturation = 0.66, order = [i[0] for i in medpubyr], 
    ax = ax1)
ax1.set_xlabel('Resources', weight = 'light', fontsize = 12, 
        variant = 'normal', color = '#000000', stretch = 'normal')
ax1.set_ylabel('Year', weight = 'light', fontsize = 12, 
    variant = 'normal', color = '#000000', stretch = 'normal')
plt.setp(ax1.xaxis.get_majorticklabels(), rotation = 90)
for t in ax1.xaxis.get_major_ticks():
    t.label.set_fontsize(lab_size[0])

for t in ax1.yaxis.get_major_ticks():
    t.label.set_fontsize(lab_size[1])

ax1.set_xticklabels([''])

ax1.xaxis.label.set_size(axis_lab_size)
ax1.yaxis.label.set_size(axis_lab_size)

fig.tight_layout()
plt.savefig('refs-year-db.pdf')
plt.close()

# ## References by Database & Year & number of refs # boxplot & barplot
# this works if the boxplot code above has been run:

fig = plt.figure()
gs = gridspec.GridSpec(2, 2, height_ratios=[2, 8], width_ratios = [2,8])
#, wspace = 0.05, hspace = 0.05)
# here comes the top barplot
ax0 = plt.subplot(gs[1])
refc_by_db = points.database.value_counts()
refc = [refc_by_db[s] if s != omnipath else len(pubmeds) for s in boxplot_ordr]
barplot = ax0.bar(range(len(refc)), refc, 
    color = ['#6EA945' if s == omnipath else '#007B7F' for s in boxplot_ordr])
tick_loc = np.arange(len(boxplot_ordr)) + 0.3
plt.xticks(tick_loc)
ax0.set_xticklabels(boxplot_ordr, rotation = 90, fontsize = lab_size[0] * 0.66)
ax0.set_yscale('log')
ax0.set_xlim([-0.2, 27.8])
ax0.set_xlabel('')
ax0.set_ylabel('Number of\nPubMed IDs')

# the refs by year barplot
ax2 = plt.subplot(gs[2])
refc_by_y = points.year.value_counts()
refc = [refc_by_y[y] if y in refc_by_y else 0 \
    for y in xrange(min(refc_by_y.index) - 1,max(refc_by_y.index) + 1)]
barplot = ax2.barh(np.arange(len(refc)) - 0.5, refc,
    color = ['#007B7F'] * len(refc))
ax2.set_xscale('log')
ax2.set_ylim([1.0, len(refc) - 0.2])
ax2.set_xlim([10000, 0])
ax2.set_yticklabels([''])
ax2.set_xlabel('Number of\nPubMed IDs')

# here the boxplot
ax1 = plt.subplot(gs[3])
sns.set(font = 'Helvetica Neue LT Std')
sns.set_context('poster', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0, 'figure.figsize': [12.8, 8.8]})
ax1 = sns.boxplot('database', 'year', data = allpoints, 
    palette = ['#6EA945' if i[0] == omnipath else '#007B7F' for i in  medpubyr], 
    linewidth = 0.1, saturation = 0.66, order = [i[0] for i in medpubyr], 
    ax = ax1)
ax1.set_xlabel('Resources', weight = 'light', fontsize = 12, 
        variant = 'normal', color = '#000000', stretch = 'normal')
ax1.set_ylabel('Year', weight = 'light', fontsize = 12, 
    variant = 'normal', color = '#000000', stretch = 'normal')
#tick_loc = np.arange(min(refc_by_y.index) - 1, max(refc_by_y.index) + 1)
#ax1.set_yticks(tick_loc)
#ax1.set_ylim([min(tick_loc) - 1, max(tick_loc) + 1])
#ax1.set_yticklabels(['%u'%y for y in tick_loc], fontsize = lab_size[0]*0.3)
plt.setp(ax1.xaxis.get_majorticklabels(), rotation = 90)
for t in ax1.xaxis.get_major_ticks():
    t.label.set_fontsize(lab_size[0])

#for t in ax1.yaxis.get_major_ticks():
#    #t.label.set_fontsize(lab_size[1])

ax1.set_xticklabels([''])

ax1.xaxis.label.set_size(axis_lab_size*0.66)
ax1.yaxis.label.set_size(axis_lab_size*0.66)

fig.tight_layout()
plt.savefig('refs-year-db-y.pdf')
plt.close()


# ## # ## # ## # ## # ##
# Rolland 2014 like visualization
a = net.graph.get_adjacency()
a = list(a)
ordr = [j[0] for j in sorted([(i, sum(l)) for i, l in enumerate(a)], key = lambda x: x[1], reverse = True)]
refs_per_protein = [len(net.graph.vs[i]['references']) for i in ordr]
an = np.array(a)
# sorting by cols and rows
an = an[an.sum(axis = 1).argsort()[::-1],:][:,an.sum(axis = 0).argsort()[::-1]]
bincount = 50
binsize = int(math.ceil(len(an) / float(bincount)))
abins = []
refs_pprot_bins = []
for i in xrange(0, bincount * binsize, binsize):
    thisRow = []
    for j in xrange(0, bincount * binsize, binsize):
        thisRow.append(sum(sum(an[i:i + binsize, j:j + binsize])))
    abins.append(thisRow)
    refs_pprot_bins.append(sum(refs_per_protein[i:i + binsize]))

anbins = np.array(abins)
norm_anbins = np.tril(anbins, k = -1)
norm_anbins = np.log10(norm_anbins)
norm_anbins[norm_anbins == -np.inf] = 0.0
norm_anbins = norm_anbins / norm_anbins.max()
# rgb(0, 123, 127)
# (0.0, 0.4823529411764706, 0.4980392156862745)
cdict = {'red':   [(0.0,  1.0, 1.0),
                   (1.0,  0.0, 0.0)],
         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.4823529411764706, 0.4823529411764706)],
         'blue':  [(0.0,  1.0, 1.0),
                   (1.0,  0.4980392156862745, 0.4980392156862745)]}
cmap = mpl.colors.LinearSegmentedColormap('emblpetrol', cdict)
sns.set(font = 'Helvetica Neue LT Std')
sns.set_context('poster', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0, 'figure.figsize': [8.8, 8.8]})
fig, ax = plt.subplots()
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
heatmap = ax.pcolor(norm_anbins.T[::-1], cmap = cmap, color = '#FFFFFF')
fig.savefig('refs-heatmap.pdf')
plt.close(fig)

# Rolland 2014 like visualization
a = net.graph.get_adjacency()
a = list(a)
ordr = [j[0] for j in sorted([(i, sum(l)) for i, l in enumerate(a)], 
    key = lambda x: x[1], reverse = True)]
refs_per_protein = [len(net.graph.vs[i]['references']) for i in ordr]
an = np.array(a)
# sorting by cols and rows
an = an[an.sum(axis = 1).argsort()[::-1],:][:,an.sum(axis = 0).argsort()[::-1]]
bincount = 50
binsize = int(math.ceil(len(an) / float(bincount)))
abins = []
prot_bins = []
refs_pprot_bins = []
for i in xrange(0, bincount * binsize, binsize):
    thisRow = []
    for j in xrange(0, bincount * binsize, binsize):
        thisRow.append(sum(sum(an[i:i + binsize, j:j + binsize])))
    abins.append(thisRow)
    refs_pprot_bins.append(sum(refs_per_protein[i:i + binsize]))
    prot_bins.append(set(ordr[i:i + binsize]))

refs_int_bins = np.array([len(uniqList(flatList([[r.pmid for r in e['references']] \
    for e in net.graph.es if e.source in p or e.target in p]))) for p in prot_bins])
int_bins = np.array([len([e \
    for e in net.graph.es if e.source in p or e.target in p]) for p in prot_bins])
refs_per_int = refs_int_bins / int_bins.astype(float)

itg_bins = [len(set([net.graph.vs[p]['name'] for p in ps]) & set(net.lists['Intogene'])) \
    for ps in prot_bins]

cgc_bins = [len(set([net.graph.vs[p]['name'] for p in ps]) & \
        set(net.lists['CancerGeneCensus'])) \
    for ps in prot_bins]

anbins = np.array(abins)
norm_anbins = np.tril(anbins, k = -1) + np.diagflat(np.diagonal(anbins))
norm_anbins = np.log10(norm_anbins)
norm_anbins[norm_anbins == -np.inf] = 0.0
norm_anbins = norm_anbins / norm_anbins.max()
# rgb(0, 123, 127)
# (0.0, 0.4823529411764706, 0.4980392156862745)
cdict = {'red':   [(0.0,  1.0, 1.0),
                   (1.0,  0.0, 0.0)],
         'green': [(0.0,  1.0, 1.0),
                   (1.0,  0.4823529411764706, 0.4823529411764706)],
         'blue':  [(0.0,  1.0, 1.0),
                   (1.0,  0.4980392156862745, 0.4980392156862745)]}
cmap = mpl.colors.LinearSegmentedColormap('emblpetrol', cdict)
# sns.set_style('white')
sns.set(font = 'Helvetica Neue LT Std')
sns.set_context('poster', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
    'grid.linewidth': 1.0, 'figure.figsize': [10, 16]})
sns.set_style({'axes.facecolor': '#FFFFFF', 'figure.facecolor': '#FFFFFF'})

fig = plt.figure()
gs = gridspec.GridSpec(6, 2, height_ratios=[2, 8, 2, 2, 1, 1], width_ratios = [4, 1])
#, wspace = 0.05, hspace = 0.05)
ax0 = plt.subplot(gs[0])
barplot = ax0.bar(range(len(refs_int_bins)), refs_int_bins, color = '#6EA945')
plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'off', labelbottom = 'off')
ax0.set_ylabel('Number of\npublications', fontsize = axis_lab_size * 0.33)
ax1 = plt.subplot(gs[2])
heatmap = ax1.pcolor(norm_anbins.T[::-1], cmap = cmap, color = '#FFFFFF')
plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'off', labelbottom = 'off')
plt.tick_params(axis = 'y', which = 'both', left = 'off', right = 'off', labelleft = 'off')
ax2 = plt.subplot(gs[3])
barplot = ax2.barh(range(len(refs_int_bins)), list(reversed(refs_int_bins)), color = '#6EA945')
plt.tick_params(axis = 'y', which = 'both', left = 'off', right = 'off', labelleft = 'off')
plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'on', 
    labelbottom = 'off', labeltop = 'on')
ax2.set_xlabel('Number of\npublications', fontsize = axis_lab_size * 0.33)
ax2.xaxis.set_label_position('top') 
plt.setp(ax2.xaxis.get_majorticklabels(), rotation = -90)
ax4 = plt.subplot(gs[4])
barplot = ax4.bar(range(len(int_bins)), int_bins, color = '#6EA945')
plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'off', labelbottom = 'off')
ax4.set_ylabel('Number of\ninteractions', fontsize = axis_lab_size * 0.33)
ax6 = plt.subplot(gs[6])
barplot = ax6.bar(range(len(refs_per_int)), refs_per_int, color = '#6EA945')
plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'off', labelbottom = 'off')
ax6.set_ylabel('Reference per\ninteraction', fontsize = axis_lab_size * 0.33)
ax8 = plt.subplot(gs[8])
sctplot = ax8.scatter(np.array(range(len(itg_bins))) + 0.5, [0] * len(itg_bins), 
    s = [r*4.7 for r in itg_bins],
    color = '#FCCC06', alpha = 0.5, clip_on = False)
plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'off', labelbottom = 'off')
plt.tick_params(axis = 'y', which = 'both', left = 'off', right = 'off', labelleft = 'off')
ax8.set_xlim([0.0, 50.0])
ax8.set_ylabel('Intogene\ncancer drivers', fontsize = axis_lab_size * 0.33)
ax10 = plt.subplot(gs[10])
sctplot = ax10.scatter(np.array(range(len(cgc_bins))) + 0.5, [0] * len(cgc_bins), 
    s = [r*4.7 for r in cgc_bins],
    color = '#FCCC06', alpha = 0.5, clip_on = False)
plt.tick_params(axis = 'x', which = 'both', bottom = 'off', top = 'off', labelbottom = 'off')
plt.tick_params(axis = 'y', which = 'both', left = 'off', right = 'off', labelleft = 'off')
ax10.set_xlim([0.0, 50.0])
ax10.set_ylabel('CGC\ncancer drivers', fontsize = axis_lab_size * 0.33)
plt.tight_layout(fig)
plt.savefig('refs-heatmap-2.pdf')
plt.close()

