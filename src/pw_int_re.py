#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  Copyright (c) 2014-2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  This code is not for public use.
#  Please do not redistribute.
#  For permission please contact me.
#
#  Website: http://www.ebi.ac.uk/~denes
#

import pypath
from pypath import data_formats
from pypath import plot
from collections import Counter

import matplotlib as mpl
import seaborn as sns

import scipy.stats as stats

pw = pypath.PyPath()
pw.init_network(data_formats.pathway)

it = pypath.PyPath()
it.init_network(data_formats.interaction_htp)

def significance(p):
    return 'ns' if p > 0.05 \
        else '*' if p > 0.01 \
        else '**' if p > 0.001 \
        else '***' if p > 0.0001 \
        else '****'

labels = [
    'Activitiy flow\npathways',
    'Interaction\nnetworks'
]

degree = [
    pw.graph.vs.degree(),
    it.graph.vs.degree()
]

nbtw = [
    np.array(pw.graph.vs.betweenness()),
    np.array(it.graph.vs.betweenness())
]

nbtw = [
    (nbtw[0] - min(nbtw[0])) / (max(nbtw[0]) - min(nbtw[0])),
    (nbtw[1] - min(nbtw[1])) / (max(nbtw[1]) - min(nbtw[1]))
]

ebtw = [
    np.array(pw.graph.es.edge_betweenness()),
    np.array(it.graph.es.edge_betweenness())
]

trans = [
    np.array(pw.graph.transitivity_local_undirected(mode = "zero")),
    np.array(it.graph.transitivity_local_undirected(mode = "zero"))
]

refc = [
    Counter(
        common.flatList(
            map(
                lambda e:
                    map(
                        lambda r:
                            r.pmid,
                        e['references']
                    ),
                pw.graph.es
            )
        )
    ).values(),
    Counter(
        common.flatList(
            map(
                lambda e:
                    map(
                        lambda r:
                            r.pmid,
                        e['references']
                    ),
                it.graph.es
            )
        )
    ).values()
]

refc2 = [
    map(lambda e: len(e['references']), pw.graph.es),
    map(lambda e: len(e['references']), it.graph.es)
]

fig = mpl.figure.Figure(figsize = (9.0, 6.0))

mw = stats.mannwhitneyu(degree[0], degree[1])
ax = fig.add_subplot(231)
vi = ax.boxplot(degree)
ax.set_yscale('log')
ax.set_xticks([1.0, 2.0])
xtlabs = ax.xaxis.set_ticklabels(labels)
ax.set_title('Degrees')
ylim = ax.get_ylim()
ax.text(1.5, ylim[1] * 0.7, '%s, p = %g' % (significance(mw.pvalue), mw.pvalue),
        horizontalalignment = 'center', fontsize = 10)
ax.xaxis.grid(False)

mw = stats.mannwhitneyu(nbtw[0], nbtw[1])
ax = fig.add_subplot(232)
vi = ax.boxplot(nbtw)
ax.set_yscale('log')
ax.set_xticks([1.0, 2.0])
xtlabs = ax.xaxis.set_ticklabels(labels)
ax.set_title('Node betweenness')
ylim = ax.get_ylim()
ax.text(1.5, ylim[1] * 0.6, '%s, p = %g' % (significance(mw.pvalue), mw.pvalue),
        horizontalalignment = 'center', fontsize = 10)
ax.xaxis.grid(False)

mw = stats.mannwhitneyu(ebtw[0], ebtw[1])
ax = fig.add_subplot(233)
vi = ax.boxplot(ebtw)
ax.set_yscale('log')
ax.set_xticks([1.0, 2.0])
xtlabs = ax.xaxis.set_ticklabels(labels)
ax.set_title('Edge betweenness')
ylim = ax.get_ylim()
ax.text(1.5, ylim[1] * 0.5, '%s, p = %g' % (significance(mw.pvalue), mw.pvalue),
        horizontalalignment = 'center', fontsize = 10)
ax.xaxis.grid(False)

mw = stats.mannwhitneyu(trans[0], trans[1])
ax = fig.add_subplot(234)
vi = ax.boxplot(trans)
ax.set_yscale('log')
ax.set_xticks([1.0, 2.0])
xtlabs = ax.xaxis.set_ticklabels(labels)
ax.set_title('Local transitivity')
ylim = ax.get_ylim()
ax.text(1.5, ylim[1] * 0.7, '%s, p = %g' % (significance(mw.pvalue), mw.pvalue),
        horizontalalignment = 'center', fontsize = 10)
ax.xaxis.grid(False)

mw = stats.mannwhitneyu(refc[0], refc[1])
ax = fig.add_subplot(235)
vi = ax.boxplot(refc)
ax.set_yscale('log')
ax.set_xticks([1.0, 2.0])
xtlabs = ax.xaxis.set_ticklabels(labels)
ax.set_title('References: edge count')
ylim = ax.get_ylim()
ax.text(1.5, ylim[1] * 0.5, '%s, p = %g' % (significance(mw.pvalue), mw.pvalue),
        horizontalalignment = 'center', fontsize = 10)
ax.xaxis.grid(False)

mw = stats.mannwhitneyu(refc2[0], refc2[1])
ax = fig.add_subplot(236)
vi = ax.boxplot(refc2)
ax.set_yscale('log')
ax.set_xticks([1.0, 2.0])
xtlabs = ax.xaxis.set_ticklabels(labels)
ax.set_title('Edges: ref count')
ylim = ax.get_ylim()
ax.text(1.5, ylim[1] * 0.6, '%s, p = %g' % (significance(mw.pvalue), mw.pvalue),
        horizontalalignment = 'center', fontsize = 10)
ax.xaxis.grid(False)

fig.tight_layout()
fig.savefig('boxplots.pdf')

reload(plot)

hist = plot.Histogram(data, labels,
                      fname = 'degree_hist.pdf',
                      xlab = 'Degree',
                      ylab = 'Nodes',
                      tone = 0,
                      legend_size = 7,
                      x_log = True,
                      y_log = True,
                      range = (0, 200),
                      normed = True,
                      finish = True)

hist.finish()





reload(plot)

hist = plot.Histogram(data, labels,
                      fname = 'betweenness_hist.pdf',
                      xlab = 'Betweenness',
                      ylab = 'Nodes',
                      tone = 0,
                      legend_size = 7,
                      kde_base = .05,
                      x_log = True,
                      y_log = False,
                      normed = True)

hist.finish()




reload(plot)

hist = plot.Histogram(data, labels,
                      fname = 'ebetweenness_hist.pdf',
                      xlab = 'Betweenness',
                      ylab = 'Edges',
                      tone = 0,
                      legend_size = 7,
                      kde_base = .05,
                      #xlim = (1, 10),
                      nbins = 10,
                      x_log = True,
                      y_log = False,
                      normed = True)

hist.finish()






reload(plot)

hist = plot.Histogram(data, labels,
                      fname = 'transitivity_hist.pdf',
                      xlab = 'Local transitivity',
                      ylab = 'Nodes',
                      tone = 0,
                      legend_size = 7,
                      kde_base = .05,
                      #xlim = (1, 10),
                      nbins = 10,
                      x_log = True,
                      y_log = False,
                      normed = True)

hist.finish()


reload(plot)

hist = plot.Histogram(data, labels,
                      fname = 'htp_hist.pdf',
                      xlab = 'Interactions',
                      ylab = 'References',
                      tone = 0,
                      nbins = 100,
                      legend_size = 7,
                      kde_base = .05,
                      #xlim = (1, 10),
                      x_log = True,
                      y_log = True,
                      normed = True)

hist.finish()


# ### ### ### HTP analysis

refcc = [
    Counter(
        common.flatList(
            map(
                lambda e:
                    map(
                        lambda r:
                            r.pmid,
                        e['references']
                    ),
                pw.graph.es
            )
        )
    ),
    Counter(
        common.flatList(
            map(
                lambda e:
                    map(
                        lambda r:
                            r.pmid,
                        e['references']
                    ),
                it.graph.es
            )
        )
    )
]

htp_refs_pw = \
set(
    map(
        lambda (pubmed, n):
            pubmed,
        filter(
            lambda (pubmed, n):
                n > 25,
            refcc[0].iteritems()
        )
    )
)

htp_refs_it = \
set(
    map(
        lambda (pubmed, n):
            pubmed,
        filter(
            lambda (pubmed, n):
                n > 25,
            refcc[1].iteritems()
        )
    )
)

htp_edges_it = \
set(
    map(
        lambda e:
            e.index,
        filter(
            lambda e:
                len(set(
                    map(
                        lambda r:
                            r.pmid,
                        e['references']
                    )
                ) - htp_refs_it) == 0,
            it.graph.es
        )
    )
)

len(htp_edges_it) / float(it.graph.ecount())

htp_edges_pw = \
set(
    map(
        lambda e:
            e.index,
        filter(
            lambda e:
                len(set(
                    map(
                        lambda r:
                            r.pmid,
                        e['references']
                    )
                ) - htp_refs_pw) == 0,
            pw.graph.es
        )
    )
)

len(htp_edges_pw) / float(pw.graph.ecount())

htp_sources_it = {}

for e in it.graph.es:
    for s, rs in e['refs_by_source'].iteritems():
        if not len(set([r.pmid for r in rs]) - htp_refs_it):
            if s not in htp_sources_it:
                htp_sources_it[s] = 0
            htp_sources_it[s] += 1

htp_sources_pw = {}

for e in pw.graph.es:
    for s, rs in e['refs_by_source'].iteritems():
        if not len(set([r.pmid for r in rs]) - htp_refs_pw):
            if s not in htp_sources_it:
                htp_sources_pw[s] = 0
            htp_sources_pw[s] += 1

it.graph.delete_edges(htp_edges_it)
list(np.where(np.array(it.graph.vs.degree()) == 0)[0])

len(list(np.where(np.array(it.graph.vs.degree()) == 0)[0])) / \
    float(it.graph.vcount())

pw.graph.delete_edges(htp_edges_pw)
list(np.where(np.array(pw.graph.vs.degree()) == 0)[0])

len(list(np.where(np.array(pw.graph.vs.degree()) == 0)[0])) / \
    float(pw.graph.vcount())
