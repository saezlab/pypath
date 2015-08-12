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

# stats and plotting modules #

import numpy as np
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

net.init_network(pfile = 'cache/default_plus_acsn.pickle')

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
rads = [30.0 * xi/float(max(topdf['Interactions'])) + 5 \
    for xi in list(topdf['Proteins'])]
pc = plt.scatter(list(topdf['Proteins']), list(topdf['Interactions']), \
    s = [np.pi * r**2 for r in rads], c = '#6ea945', alpha = 0.5, \
        edgecolors = 'none')

for label, x, y, o in \
    zip(topdf['Database'], topdf['Proteins'], topdf['Interactions'], rads):
    ang = angles.next()
    print 'angle of %s: %u' % (label, ang)
    coo = rotate((-1 * max(o, 30), -1 * max(o, 30)), ang)
    dst = dists.next()
    coo = (-7 - dst / float(3), 21 + dst)
    print coo
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
plt.axis('tight')
# plt.xlim([-0.01, 0.07])
#ax.autoscale_view()
fig.savefig('vcount-ecount.pdf')
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
