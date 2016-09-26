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

import math

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
from bioigraph.common import *

# literature curated network
net = bioigraph.BioGraph(9606)
net.init_network(pfile = 'cache/default_network_wo-intact_ltp-only.pickle')
net.read_list_file(bioigraph.data_formats.cgc)
net.read_list_file(bioigraph.data_formats.intogene_cancer)
# high confidence interactome from Vidal Lab, 2nd release, 2014
hi2 = bioigraph.BioGraph(9606)
hi2.init_network(lst = {'hi2': bioigraph.data_formats.ugly['hi2']})

## Rolland 2014 like visualization
# getting adjacency matrices
a = net.graph.get_adjacency()
b = hi2.graph.get_adjacency()
a = list(a)
b = list(b)
# ordering literature curated network's vertex ids by degree
ordr = [j[0] for j in sorted([(i, sum(l)) for i, l in enumerate(a)], 
    key = lambda x: x[1], reverse = True)]

# pubmed ids for each protein
refs_per_protein = [len(net.graph.vs[i]['references']) for i in ordr]
an = np.array(a)
bn = np.array(b)
# sorting by cols and rows
an = an[an.sum(axis = 1).argsort()[::-1],:][:,an.sum(axis = 0).argsort()[::-1]]
bincount = 50
binsize = int(math.ceil(len(an) / float(bincount)))
abins = []
bbins = []
prot_bins = []
hi2_prot_bins = []
refs_pprot_bins = []
for i in xrange(0, bincount * binsize, binsize):
    thisRow = []
    thisRowHi2 = []
    for j in xrange(0, bincount * binsize, binsize):
        thisRow.append(sum(sum(an[i:i + binsize, j:j + binsize])))
    abins.append(thisRow)
    refs_pprot_bins.append(sum(refs_per_protein[i:i + binsize]))
    prot_bins.append(set(ordr[i:i + binsize]))
    # vertex ids of proteins in hi2 for each bin
    hi2_prot_bins.append(set([hi2.graph.vs.find(name = net.graph.vs[p]['name']).index \
        for p in ordr[i:i + binsize] \
            if net.graph.vs[p]['name'] in hi2.graph.vs['name']]))

for i in xrange(len(hi2_prot_bins)):
    thisRow = []
    ivids = np.array(sorted(list(hi2_prot_bins[i])))
    for j in xrange(len(hi2_prot_bins)):
        jvids = np.array(sorted(list(hi2_prot_bins[j])))
        thisRow.append(np.sum(bn[ivids[:, None], jvids]) if len(ivids) > 0 and len(jvids) > 0 else 0)
    bbins.append(thisRow)

refs_int_bins = np.array([len(uniqList(flatList([[r.pmid for r in e['references']] \
    for e in net.graph.es if e.source in p or e.target in p]))) for p in prot_bins])
int_bins = np.array([len([e \
    for e in net.graph.es if e.source in p or e.target in p]) for p in prot_bins])
hi2_int_bins = np.array([len([e \
    for e in hi2.graph.es if e.source in p or e.target in p]) for p in hi2_prot_bins])
refs_per_int = refs_int_bins / int_bins.astype(float)

itg_bins = [len(set([net.graph.vs[p]['name'] for p in ps]) & set(net.lists['Intogene'])) \
    for ps in prot_bins]

cgc_bins = [len(set([net.graph.vs[p]['name'] for p in ps]) & \
        set(net.lists['CancerGeneCensus'])) \
    for ps in prot_bins]

anbins = np.array(abins)
bnbins = np.array(bbins)

norm_anbins = np.tril(anbins, k = -1) + np.diagflat(np.diagonal(anbins))
norm_anbins = np.log10(norm_anbins)
norm_anbins[norm_anbins == -np.inf] = 0.0
norm_anbins = norm_anbins / norm_anbins.max()

norm_bnbins = np.tril(bnbins, k = -1) + np.diagflat(np.diagonal(bnbins))
norm_bnbins = np.log10(norm_bnbins)
norm_bnbins[norm_bnbins == -np.inf] = 0.0
norm_bnbins = norm_bnbins / norm_bnbins.max()
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
plt.savefig('refs-heatmap-3.pdf')
plt.close()
