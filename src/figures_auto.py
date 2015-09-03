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
    # This can be run non interactively
    # to regenerate the article figures.
#

# generic modules #

from collections import Counter

# stats and plotting modules #

import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import _sensitivity as sens

# from bioigraph #

import bioigraph
from bioigraph.common import *
from bioigraph.data_formats import best, good, ugly, transcription
from bioigraph import dataio
from bioigraph import plot

# parameters
omnipath = 'OmniPath'
lab_size = (18, 21)
axis_lab_size = 36

net = bioigraph.BioGraph(9606)

#net.init_network(pfile = 'cache/default_plus_acsn_wo-intact.pickle')
net.init_network(pfile = 'cache/default_network_wo-intact_ltp-only.pickle')

net.curation_tab(latex_hdr = False, fname = 'curation_stats_stripped.tex')

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

cdrivers_notcov = uniqList(flatList([net.mapper.map_name(n, 'uniprot', 'genesymbol') \
    for n in (set(net.lists['CancerGeneCensus']) | set(net.lists['Intogene'])) - \
        set(net.graph.vs['name'])]))

cdrivers_notcov = [
    (n,
    ';'.join(net.mapper.map_name(n, 'uniprot', 'genesymbol')),
    ';'.join(net.mapper.map_name(n, 'uniprot', 'protein-name'))) \
    for n in list((set(net.lists['CancerGeneCensus']) | set(net.lists['Intogene'])) - \
        set(net.graph.vs['name']))]

with open('cdrivers_notcov', 'w') as f:
    f.write('\n'.join(['\t'.join(i) for i in cdrivers_notcov]))

# source-vcount sens.barplot
d = zip(*[(s, len([v for v in net.graph.vs if s in v['sources']])) for s in net.sources] + \
    [(omnipath, net.graph.vcount())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'proteins-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Number of proteins', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y',
    y_break = (0.37, 0.15))

vcount_ordr = list(bp.ordr)

# source-ecount sens.barplot

d = zip(*[(s, len([e for e in net.graph.es if s in e['sources']])) for s in net.sources] + \
    [(omnipath, net.graph.ecount())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'interactions-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Number of interactions', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.19, 0.1))

# density sens.barplot
d = zip(*[(s, g.density()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.density())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'density-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Graph density', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.5, 0.1))

# density with same order as protein count:
d = zip(*[(s, g.density()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.density())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'density-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Graph density', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.5, 0.1))

# ecount - cce scatterplot:
grp = {
    'CancerCellMap': 'p',
    'InnateDB': 'i',
    'SPIKE': 'p',
    'LMPID': 'm',
    'DIP': 'i',
    'HPRD': 'm',
    'PDZBase': 'p',
    'dbPTM': 'm',
    'MatrixDB': 'i',
    'DOMINO': 'm',
    'Signor': 'p',
    'Macrophage': 'p',
    'NetPath': 'r',
    'ELM': 'm',
    'SignaLink2': 'p',
    'NRF2ome': 'p',
    'DEPOD': 'm',
    'phosphoELM': 'm',
    'MPPI': 'i',
    'Guide2Pharmacology': 'p',
    'TRIP': 'p',
    'AlzPathway': 'r',
    'PhosphoSite': 'm',
    'CA1': 'p',
    'NCI-PID': 'r',
    'DeathDomain': 'p',
    'ARN': 'p'
}

cl = {
    'p': '#6ea945',
    'm': '#007B7F',
    'i': '#FCCC06',
    'r': '#646567'
}

labs = {
    'p': 'Pathway',
    'm': 'PTM',
    'r': 'Reaction',
    'i': 'Interaction'
}

cs = net.curation_stats()
fig, ax = plt.subplots()
font_family = 'Helvetica Neue LT Std'
sns.set(font = font_family)
ax.set_yscale('log')
ax.set_xscale('log')
p = plt.scatter([sep[s].ecount() for s in sorted(sep.keys())], 
    [cs[s]['corrected_curation_effort'] for s in sorted(sep.keys())], 
    c = [cl[grp[s]] for s in sorted(sep.keys())], edgecolors = 'none')
for g, c in cl.iteritems():
    xxx = np.log10(np.array([sep[s].ecount() for s in sorted(sep.keys()) if grp[s] == g]))
    m, b = np.polyfit(xxx, 
        np.log10(np.array([cs[s]['corrected_curation_effort'] for s in sorted(sep.keys()) if grp[s] == g])), 1)
    plt.plot([10**xx for xx in xxx], [10**y for y in m * xxx + b], '-', color = c, label = labs[g])

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.xlabel('Number of interactions')
plt.ylabel('Corrected curation effort')
ax.xaxis.label.set_size(axis_lab_size*0.66)
ax.yaxis.label.set_size(axis_lab_size*0.66)
fig.tight_layout()
fig.savefig('ecount-cce.pdf')
plt.close(fig)

# vcount - ecount scatterplot:
topdata = {
'Proteins': [len([v for v in net.graph.vs if s in v['sources']]) for s in net.sources],
'Interactions': [len([e for e in net.graph.es if s in e['sources']]) for s in net.sources],
'Density': [sep[s].density() for s in net.sources],
'Transitivity': [sep[s].transitivity_undirected() for s in net.sources],
'Diameter': [sep[s].diameter() for s in net.sources],
'Database': net.sources
}
topdf = pd.DataFrame(topdata, index = net.sources)
topdf.sort(columns = ['Proteins'], inplace = True)
def rotate_labels(angles = (0, -90, -135, -180, -270, -315)):
    i = 0
    while True:
        yield angles[i % len(angles)]
        i += 1

def move_labels(dist = (0, 10, 20, 30, 40, 50)):
    i = 0
    while True:
        yield dist[i % len(dist)] if i < 20 else 10
        i += 1

dists = move_labels()
fig, ax = plt.subplots()
font_family = 'Helvetica Neue LT Std'
sns.set(font = font_family)
rads = [30.0 * xi/float(max(topdf['Interactions'])) + 5 \
    for xi in list(topdf['Proteins'])]
ax.set_yscale('log')
ax.set_xscale('log')
pc = plt.scatter(list(topdf['Proteins']), list(topdf['Interactions']), \
    s = [np.pi * r**2 for r in rads], c = '#6ea945', alpha = 0.5, \
        edgecolors = 'none')

ax.set_xlim([-20.0, 3500.0])
ax.set_ylim([-20.0, 11000.0])

tick_loc = [10, 20, 50 ,100, 200, 500, 1000, 2000, 5000, 10000, 20000]
plt.xticks(tick_loc)
plt.yticks(tick_loc)

for label, x, y, o in \
    zip(topdf['Database'], topdf['Proteins'], topdf['Interactions'], rads):
    dst = dists.next()
    coo = (-7 - dst / float(3), 21 + dst)
    plt.annotate(
        label, 
        xy = (x, y), xytext = coo,
        xycoords = 'data',
        textcoords = 'offset points', ha = 'center', va = 'bottom', color = '#007B7F',
        arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc,rad=.0',
            color = '#007B7F', edgecolor = '#007B7F', alpha = 1.0, 
            visible = True, linewidth = 0.2), 
            )

plt.xlabel('Number of proteins')
plt.ylabel('Number of interacting pairs')

ax.set_xlim([-20.0, 3500.0])
ax.set_ylim([-20.0, 11000.0])

tlabs = []
for i, t in enumerate(list(ax.xaxis.get_major_locator().locs)):
    tlabs.append(str(t))

ax.set_xticklabels(tlabs)

tlabs = []
for i, t in enumerate(list(ax.yaxis.get_major_locator().locs)):
    tlabs.append(str(t))

ax.set_yticklabels(tlabs)

for t in ax.xaxis.get_major_ticks():
    t.label.set_fontsize(lab_size[0]*0.66)

for t in ax.yaxis.get_major_ticks():
    t.label.set_fontsize(lab_size[1]*0.66)

ax.xaxis.label.set_size(axis_lab_size*0.66)
ax.yaxis.label.set_size(axis_lab_size*0.66)
fig.tight_layout()
#ax.autoscale_view()
fig.savefig('vcount-ecount-log.pdf')
plt.close(fig)


# transitivity sens.barplot
d = zip(*[(s, g.transitivity_undirected()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.transitivity_undirected())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'transitivity-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Graph global transitivity', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'transitivity-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Graph global transitivity', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

# diameter sens.barplot
d = zip(*[(s, g.diameter()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.diameter())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'diameter-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Graph diameter', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'diameter-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Graph diameter', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

# receptors sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['rec']])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec']))])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptors-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Number of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', y_break = (0.6, 0.09))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptors-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Number of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, y_break = (0.6, 0.09))

# receptors prop sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['rec']])/float(g.vcount())*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec'])/float(net.graph.vcount())*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorprop-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', y_break = (0.39, 0.12))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorprop-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, y_break = (0.39, 0.12))

# receptors prop sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['rec']])/float(len(net.lists['rec']))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec'])/float(len(net.lists['rec']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorcov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of receptors covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorcov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of receptors covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

# transcription factors sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['tf']])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf']))])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfs-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Number of transcription factors', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y',
    y_break = (0.45, 0.2))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfs-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Number of transcription factors', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.45, 0.2))

# transcription factor prop sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['tf']])/float(g.vcount())*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf'])/float(net.graph.vcount())*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfprop-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of transcription factors', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.53, 0.2))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfprop-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of transcription factors', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False,
    y_break = (0.53, 0.2))

d = zip(*[(s, len([v for v in g.vs if v['tf']])/float(len(net.lists['tfs']))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf'])/float(len(net.lists['tfs']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfcov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of TFs covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.5, 0.1))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfcov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of TFs covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
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
    data = None, fname = 'discov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of disease-gene' + '\nassociations covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.53, 0.12))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'discov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = r'% of disease-gene' + '\nassociations covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False,
    y_break = (0.53, 0.12))


# number of complexes sens.barplot
d = zip(*[(s, len(sens.complexes_in_network(g))) \
    for s, g in sep.iteritems()] + \
    [(omnipath, len(sens.complexes_in_network(net.graph)))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'complexes-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of complexes',
    #\n(out of %u)'%\
        #len(sens.complexes_in_network(net.graph)), 
    xlab = 'Pathway resources', order = 'y', y_break = (0.55, 0.05))
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'complexes-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of complexes',
    #\n(out of %u)'%\
        #len(sens.complexes_in_network(net.graph)), 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.55, 0.05))


# number of ptms sens.barplot
d = zip(*[(s, sum([len(e['ptm']) for e in g.es])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum([len(e['ptm']) for e in net.graph.es]))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ptms-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of PTMs',
    #\n(out of %u PTMs)'%\
        #sum([len(e['ptm']) for e in net.graph.es]), 
    xlab = 'Pathway resources', order = 'y', y_break = (0.6, 0.1))
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ptms-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of PTMs',
    #\n(out of %u PTMs)'%\
        #sum([len(e['ptm']) for e in net.graph.es]), 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.6, 0.1))


# number of ptms sens.barplot
d = zip(*[(s, len([1 for e in g.es if len(e['ptm']) != 0])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, len([1 for e in net.graph.es if len(e['ptm']) != 0]))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'havingptm-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Interactions having PTM', 
    xlab = 'Pathway resources', order = 'y', y_break = (0.5, 0.15))
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'havingptm-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Interactions having PTM', 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.5, 0.15))


# TODO difference between COSMIC and Intogene

# cosmic cancer gene census coverage sens.barplot
d = zip(*[(s, len(set(net.lists['CancerGeneCensus']) & set(g.vs['name'])) / \
        float(len(net.lists['CancerGeneCensus']))*100) \
        for s, g in sep.iteritems()] + \
    [(omnipath, len(set(net.lists['CancerGeneCensus']) & set(net.graph.vs['name'])) / \
        float(len(net.lists['CancerGeneCensus']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ccgccov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of CGC genes', 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ccgccov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of CGC genes', 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

# intogene cancer driver coverage sens.barplot
d = zip(*[(s, len(set(net.lists['Intogene']) & set(g.vs['name'])) / \
        float(len(net.lists['Intogene']))*100) \
        for s, g in sep.iteritems()] + \
    [(omnipath, len(set(net.lists['Intogene']) & set(net.graph.vs['name'])) / \
        float(len(net.lists['Intogene']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of Intogene genes', 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of Intogene genes', 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

# dirs signes w/o omnipath
vcount_ordr.remove(omnipath)
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
stacked_barplot(x = d[0], y = d[1:], 
    names = ['positive', 'negative', 'unknown effect', 'unknown direction'],
    data = None, fname = 'dirs-signes-by-db-wo-op-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Interactions', 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

stacked_barplot(x = d[0], y = d[1:], 
    names = ['positive', 'negative', 'unknown effect', 'unknown direction'],
    data = None, fname = 'dirs-signes-by-db-wo-op.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    ylab = 'Interactions', 
    xlab = 'Pathway resources', order = 'y')

# ## ##
net.init_network(pfile = 'cache/default_network.pickle')

# ref ecount barplot
refc = Counter(flatList((r.pmid for r in e['references']) for e in net.graph.es))

# percentage of high throughput interactions
htdata = {}
for htlim in reversed(xrange(5, 201)):
    htrefs = set([i[0] for i in refc.most_common() if i[1] > htlim])
    htedgs = [e.index for e in net.graph.es if \
        len(set([r.pmid for r in e['references']]) - htrefs) == 0]
    htsrcs = uniqList(flatList([net.graph.es[e]['sources'] for e in htedgs]))
    htdata[htlim] = {
        'rnum': len(htrefs),
        'enum': len(htedgs),
        'snum': len(htsrcs),
        'htrefs': htrefs
    }

for htlim in reversed(xrange(5, 201)):
    htedgs = [e.index for e in net.graph.es if \
        len(set([r.pmid for r in e['references']]) - htdata[htlim]['htrefs']) == 0]
    net.graph.delete_edges(htedgs)
    zerodeg = [v.index for v in net.graph.vs if v.degree() == 0]
    net.graph.delete_vertices(zerodeg)
    ltvcnt = net.graph.vcount()
    ltecnt = net.graph.ecount()
    htdata[htlim]['lenum'] = ltecnt
    htdata[htlim]['lvnum'] = ltvcnt

srcs10 = uniqList(flatList([net.graph.es[e]['sources'] for e in e10]))
srcs50 = uniqList(flatList([net.graph.es[e]['sources'] for e in e50]))
srcs100 = uniqList(flatList([net.graph.es[e]['sources'] for e in e100]))

fig, axs = plt.subplots(5, figsize = (10, 20), sharex=True)
axs[0].plot(sorted(htdata.keys()), [htdata[h]['rnum'] for h in sorted(htdata.keys())], 
    '-', color = '#007B7F')
axs[0].set_ylabel('Number of references', fontsize = axis_lab_size * 0.45)
axs[1].plot(sorted(htdata.keys()), [htdata[h]['enum'] for h in sorted(htdata.keys())], 
    '-', color = '#6EA945')
axs[1].set_ylabel('Number of edges', fontsize = axis_lab_size * 0.45)
axs[2].plot(sorted(htdata.keys()), [htdata[h]['snum'] for h in sorted(htdata.keys())], 
    '-', color = '#DA0025')
axs[2].set_ylabel('Number of resources', fontsize = axis_lab_size * 0.45)
axs[3].plot(sorted(htdata.keys()), [htdata[h]['lenum'] for h in sorted(htdata.keys())], 
    '-', color = '#996A44')
axs[3].set_ylabel('LT network edge count', fontsize = axis_lab_size * 0.45)
axs[4].plot(sorted(htdata.keys()), [htdata[h]['lvnum'] for h in sorted(htdata.keys())], 
    '-', color = '#FCCC06')
axs[4].set_xlabel('HT limit', fontsize = axis_lab_size * 0.66)
axs[4].set_ylabel('LT network node count', fontsize = axis_lab_size * 0.45)
fig.savefig('ht-limits.pdf')
plt.close(fig)



d = zip(*[(s, len(set(net.lists['Intogene']) & set(g.vs['name'])) / \
        float(len(net.lists['Intogene']))*100) \
        for s, g in sep.iteritems()] + \
    [(omnipath, len(set(net.lists['Intogene']) & set(net.graph.vs['name'])) / \
        float(len(net.lists['Intogene']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of Intogene genes', 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of Intogene genes', 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)