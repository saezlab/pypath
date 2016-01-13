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
import locale
locale.setlocale('en_GB')

# stats and plotting modules #

import math
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import _sensitivity as sens
from scipy import stats

# from pypath #

import pypath
from pypath.common import *
from pypath.data_formats import best, good, ugly, transcription
from pypath import dataio
from pypath import plot

# parameters
omnipath = 'OmniPath'
lab_size = (18, 21)
axis_lab_size = 36

net = pypath.Pypath(9606)

#net.init_network(exclude = ['intact'])
#net.third_source_directions()
#net.remove_htp()
#net.save_network('default_network_20151230.pickle')
net.init_network('default_network_20151230.pickle')
#net.init_network(pfile = 'cache/default_plus_acsn_wo-intact.pickle')
#net.init_network(pfile = 'cache/default_network_wo-intact_ltp-only.pickle')

net.curation_tab(latex_hdr = False, fname = 'curation_stats_stripped.tex')

net.genesymbol_labels()
net.set_tfs()
net.set_receptors()
net.set_chembl_mysql('chembl_ebi')
net.set_drugtargets()
net.set_kinases()
net.set_druggability()
net.load_disgenet()
net.load_corum()
sens.in_complex(net)
net.in_complex()
net.load_ptms()
net.read_list_file(pypath.data_formats.cgc)
net.read_list_file(pypath.data_formats.intogene_cancer)
# net.load_comppi()
sep = net.separate()

net.lists['rec'] = uniqList(flatList([net.mapper.map_name(rec, 'genesymbol', 'uniprot') \
    for rec in dataio.get_hpmr()]))
net.lists['dgb'] = uniqList(flatList([net.mapper.map_name(dgb, 'genesymbol', 'uniprot') \
    for dgb in dataio.get_dgidb()]))
net.lists['kin'] = uniqList(flatList([net.mapper.map_name(kin, 'genesymbol', 'uniprot') \
    for kin in dataio.get_kinases()]))
net.lists['tfs'] = uniqList(flatList([net.mapper.map_name(tf, 'ensg', 'uniprot') \
    for tf in dataio.get_tfcensus()['ensg']]))
net.lists['dis'] = uniqList(flatList([net.mapper.map_name(dis['genesymbol'], 'genesymbol', 'uniprot') \
    for dis in dataio.get_disgenet()]))

proteome = dataio.all_uniprots(swissprot = 'yes')

contDisg = np.array([[len(proteome), net.graph.vcount()], [len(net.lists['dis']), len([1 for v in net.graph.vs if len(v['dis']) > 0])]])
stats.fisher_exact(contDisg)

contCanc = np.array([[len(proteome), net.graph.vcount()], [len(uniqList(net.lists['CancerGeneCensus'] + net.lists['Intogene'])), 
    len(set(net.lists['CancerGeneCensus'] + net.lists['Intogene']) & set(net.graph.vs['name']))]])
stats.fisher_exact(contCanc)

contDgb = np.array([[len(proteome), net.graph.vcount()], [len(net.lists['dgb']), len([1 for v in net.graph.vs if v['dgb']])]])
stats.fisher_exact(contDgb)

contRec = np.array([[len(proteome), net.graph.vcount()], [len(net.lists['rec']), len([1 for v in net.graph.vs if v['rec']])]])
stats.fisher_exact(contRec)

contTf = np.array([[len(proteome), net.graph.vcount()], [len(net.lists['tfs']), len([1 for v in net.graph.vs if v['tf']])]])
stats.fisher_exact(contTf)

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
fp = mpl.font_manager.FontProperties(family = 'Helvetica Neue LT Std', style = 'normal',
    variant = 'normal', weight = 'normal', stretch = 'normal')

d = zip(*[(s, len([v for v in net.graph.vs if s in v['sources']])) for s in net.sources] + \
    [(omnipath, net.graph.vcount())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'proteins-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Number of proteins', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y',
    y_break = (0.37, 0.15))

vcount_ordr = list(bp.ordr)

# source-ecount sens.barplot

d = zip(*[(s, len([e for e in net.graph.es if s in e['sources']])) for s in net.sources] + \
    [(omnipath, net.graph.ecount())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'interactions-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Number of interactions', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = vcount_ordr, fin = False, desc = False, 
    y_break = (0.23, 0.07))

bp.ax.yaxis.labelpad = 15
bp.finish()

# density sens.barplot
d = zip(*[(s, g.density()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.density())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'density-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Graph density', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.5, 0.1))

# density with same order as protein count:
d = zip(*[(s, g.density()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.density())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'density-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
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

cl2 = {
    'p': '#DFEAD2',
    'm': '#CADADB',
    'i': '#FFF1BC',
    'r': '#CFD0D1'
}

labs = {
    'p': 'Pathway',
    'm': 'PTM',
    'r': 'Reaction',
    'i': 'Interaction'
}

# thanks for:
# http://stackoverflow.com/a/28336695/854988

cs = net.curation_stats()
fig, ax = plt.subplots()
font_family = 'Helvetica Neue LT Std'
sns.set(font = font_family)
ax.set_yscale('log')
ax.set_xscale('log')
x = np.array([sep[s].ecount() for s in sorted(sep.keys())])
y = np.array([cs[s]['corrected_curation_effort'] for s in sorted(sep.keys())])
cols = [cl[grp[s]] for s in sorted(sep.keys())]
p = plt.scatter(x, y, c = cols, edgecolors = 'none')
for g, c in cl.iteritems():
    xxx = np.log10(np.array([sep[s].ecount() for s in sorted(sep.keys()) if grp[s] == g]))
    yyy = np.log10(np.array([cs[s]['corrected_curation_effort'] for s in sorted(sep.keys()) if grp[s] == g]))
    (m, b), V = np.polyfit(xxx, yyy, 1, cov = True)
    n = xxx.size
    y_fit = np.polyval((m, b), xxx)
    df = n - 2
    t = stats.t.ppf(0.95, df)
    resid = yyy - y_fit
    chi2 = np.sum((resid / y_fit)**2)
    chi2_red = chi2 / df
    s_err = np.sqrt(np.sum(resid**2) / df)
    plt.plot([10**xx for xx in xxx], [10**y for y in m * xxx + b], '-', color = c, label = labs[g])
    x2 = np.linspace(np.min(xxx), np.max(xxx), 100)
    y2 = np.linspace(np.min(y_fit), np.max(y_fit), 100)
    # confidence interval
    ci = t * s_err * np.sqrt(1 / n + (x2 - np.mean(xxx))**2 / np.sum((xxx - np.mean(xxx))**2))
    # prediction interval
    pi = t * s_err * np.sqrt(1 + 1 / n + (x2 - np.mean(xxx))**2 / np.sum((xxx - np.mean(xxx))**2))
    ax.fill_between([10**xx for xx in x2], [10**yy for yy in (y2 + ci)], [10**yy for yy in (y2 - ci)], color = c, edgecolor = '', alpha = 0.2)
    print pi
    print [10**yy for yy in (y2 - pi)]
    print [10**yy for yy in (y2 + pi)]
    ax.plot([10**xx for xx in x2], [10**yy for yy in (y2 - pi)], '--', color = c, linewidth = 0.5)
    ax.plot([10**xx for xx in x2], [10**yy for yy in (y2 + pi)], '--', color = c, linewidth = 0.5)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.xlabel('Number of interactions')
plt.ylabel('Corrected curation effort')
ax.xaxis.label.set_size(axis_lab_size*0.66)
ax.yaxis.label.set_size(axis_lab_size*0.66)
fig.tight_layout()
fig.savefig('ecount-cce.pdf')
plt.close(fig)

# ## ## ## ##

# vcount - ecount scatterplot:
topdata = {
'Proteins': [len([v for v in net.graph.vs if s in v['sources']]) for s in net.sources],
'Interactions': [len([e for e in net.graph.es if s in e['sources']]) for s in net.sources],
'Density': [sep[s].density() for s in net.sources],
'Redensity': [1.0 / sep[s].density() for s in net.sources],
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

def move_labels(dist = (0, 10, 20, 30, 40, 50, 60, 70)):
    i = 0
    while True:
        yield dist[i % len(dist)] if i < 20 else 20
        i += 1

def overlap(bbox1, bbox2):
    return (bbox1._points[0][0] > bbox2._points[0][0] and \
        bbox1._points[0][0] < bbox2._points[1][0] or \
        bbox1._points[1][0] > bbox2._points[0][0] and \
        bbox1._points[1][0] < bbox2._points[1][0] or \
        bbox2._points[0][0] > bbox1._points[0][0] and \
        bbox2._points[0][0] < bbox1._points[1][0] or \
        bbox2._points[1][0] > bbox1._points[0][0] and \
        bbox2._points[1][0] < bbox1._points[1][0]) and \
        (bbox1._points[0][1] > bbox2._points[0][1] and \
        bbox1._points[0][1] < bbox2._points[1][1] or \
        bbox1._points[1][1] > bbox2._points[0][1] and \
        bbox1._points[1][1] < bbox2._points[1][1] or \
        bbox2._points[0][1] > bbox1._points[0][1] and \
        bbox2._points[0][1] < bbox1._points[1][1] or \
        bbox2._points[1][1] > bbox1._points[0][1] and \
        bbox2._points[1][1] < bbox1._points[1][1])

def get_moves(bbox1, bbox2):
    xmove = 0
    ymove = 0
    if bbox1._points[0][0] > bbox2._points[0][0] and \
        bbox1._points[0][0] < bbox2._points[1][0] or \
        bbox1._points[1][0] > bbox2._points[0][0] and \
        bbox1._points[1][0] < bbox2._points[1][0] or \
        bbox2._points[0][0] > bbox1._points[0][0] and \
        bbox2._points[0][0] < bbox1._points[1][0] or \
        bbox2._points[1][0] > bbox1._points[0][0] and \
        bbox2._points[1][0] < bbox1._points[1][0]:
        if (bbox1._points[0][0] + bbox1._points[1][0]) / 2.0 < \
            (bbox2._points[0][0] + bbox2._points[1][0]) / 2.0:
            xmove = bbox1._points[1][0] - bbox2._points[0][0]
        else:
            xmove = bbox1._points[0][0] - bbox2._points[1][0]
    if bbox1._points[0][1] > bbox2._points[0][1] and \
        bbox1._points[0][1] < bbox2._points[1][1] or \
        bbox1._points[1][1] > bbox2._points[0][1] and \
        bbox1._points[1][1] < bbox2._points[1][1] or \
        bbox2._points[0][1] > bbox1._points[0][1] and \
        bbox2._points[0][1] < bbox1._points[1][1] or \
        bbox2._points[1][1] > bbox1._points[0][1] and \
        bbox2._points[1][1] < bbox1._points[1][1]:
        if (bbox1._points[0][1] + bbox1._points[1][1]) / 2.0 < \
            (bbox2._points[0][1] + bbox2._points[1][1]) / 2.0:
            ymove = bbox1._points[1][1] - bbox2._points[0][1]
        else:
            ymove = bbox1._points[0][1] - bbox2._points[1][1]
    return (xmove, ymove)

dists = move_labels()
fig, ax = plt.subplots()
font_family = 'Helvetica Neue LT Std'
sns.set(font = font_family)
rads = [30.0 * xi/float(max(topdf['Interactions'])) + 5 \
    for xi in list(topdf['Proteins'])]
rads = 20.0 * (topdf['Redensity'] / max(topdf['Redensity']))**0.66
ax.set_yscale('log')
ax.set_xscale('log')
# the points:
x = np.array(topdf['Proteins'])
y = np.array(topdf['Interactions'])
pc = plt.scatter(x, y, \
    s = [np.pi * r**2 for r in rads], c = '#6ea945', alpha = 0.5, \
        edgecolors = 'none')
x = np.log10(x)
y = np.log10(y)
# (log)linear fit with confidence and prediction interval:
(m, b), V = np.polyfit(x, y, 1, cov = True)
n = x.size
y_fit = np.polyval((m, b), x)
df = n - 2
t = stats.t.ppf(0.95, df)
resid = y - y_fit
chi2 = np.sum((resid / y_fit)**2)
chi2_red = chi2 / df
s_err = np.sqrt(np.sum(resid**2) / df)
x2 = np.linspace(np.min(x), np.max(x), 100)
y2 = np.linspace(np.min(y_fit), np.max(y_fit), 100)
# confidence interval
ci = t * s_err * np.sqrt(1 / n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))
# prediction interval
pi = t * s_err * np.sqrt(1 + 1 / n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))
# regression line
plt.plot([10**xx for xx in x], [10**yy for yy in y_fit], 
    '-', color = '#B6B7B9', alpha = 0.5, rc = {'font.weight': 'light'})
# confidence interval
ax.fill_between([10**xx for xx in x2], [10**yy for yy in (y2 + ci)], 
    [10**yy for yy in (y2 - ci)], 
    color = '#B6B7B9', edgecolor = '', alpha = 0.2)
# prediction intreval
ax.plot([10**xx for xx in x2], [10**yy for yy in (y2 - pi)], 
    '--', color = '#B6B7B9', linewidth = 0.5)
ax.plot([10**xx for xx in x2], [10**yy for yy in (y2 + pi)], 
    '--', color = '#B6B7B9', linewidth = 0.5)

ax.set_xlim([-20.0, 5000.0])
ax.set_ylim([-20.0, 50000.0])

tick_loc = [10, 20, 50 ,100, 200, 500, 1000, 2000, 5000, 10000, 20000]
plt.xticks(tick_loc)
plt.yticks(tick_loc)

annots = []
for label, xx, yy, yf, o in \
    zip(topdf['Database'], topdf['Proteins'], topdf['Interactions'], y_fit, rads):
    dst = dists.next()
    d = 1.0 if yy > 10**yf else -1.0
    coo = ((-7 - dst) * d / 3.0, (21 + dst) * d)
    annots.append(plt.annotate(
        label, 
        xy = (xx, yy), xytext = coo,
        xycoords = 'data',
        textcoords = 'offset points', ha = 'center', va = 'bottom', color = '#007B7F',
        arrowprops = dict(arrowstyle = '-', connectionstyle = 'arc,rad=.0',
            color = '#007B7F', edgecolor = '#007B7F', alpha = 1.0, 
            visible = True, linewidth = 0.2), 
            ))

# legend:
lhandles = []
llabels = []
sleg = [0.001, 0.005, 0.010, 0.02, 0.05]
for s in sleg:
    lhandles.append(mpl.legend.Line2D(range(1), range(1), 
        color = 'none', marker = 'o', 
        markersize = 20.0 * (1 / s / max(topdf['Redensity']))**0.66, 
        markerfacecolor = '#6ea945', markeredgecolor = 'none', alpha = .5, 
        label = '%.3f' % s))
    llabels.append('%.3f' % s)

#xleg = [10**(max(x) * 0.95)] * len(sleg)
#ymin = min(y)
#yrng = max(y) - ymin
#ymin = ymin + yrng * 0.02
#sleg = 20.0 * (sleg/max(topdf['Redensity']))**0.66
#yleg = (np.array([np.log10(sum(sleg[:yy + 1])) for yy in xrange(len(sleg))])**13.6) + 80
#ax.scatter(xleg, yleg, marker = 'o', s = [np.pi * s**2 for s in sleg], 
    #c = '#6ea945', alpha = 0.5, edgecolors = 'none')

leg = plt.legend(lhandles, llabels, loc = 4, title = 'Density', labelspacing = .9,
    borderaxespad = .9)

leg.title().set_fontproperties(weight = 'light')

plt.xlabel('Number of proteins')
plt.ylabel('Number of interacting pairs')

ax.set_xlim([-20.0, 5000.0])
ax.set_ylim([-20.0, 20000.0])

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
fig.savefig('vcount-ecount-log-3.pdf')

steps = [0] * len(annots)
for i, a2 in enumerate(annots):
    overlaps = False
    for z in xrange(100):
        for a1 in annots[:i]:
            if overlap(a1.get_window_extent(), a2.get_window_extent()):
                print '\tOverlapping labels: %s and %s' % (a1._text, a2._text)
                mv = get_moves(a1.get_window_extent(), a2.get_window_extent())
                if steps[i] % 2 == 0:
                    a2.xyann = (a2.xyann[0] + mv[0] * 1.1 * (z / 2 + 1), a2.xyann[1])
                else:
                    a2.xyann = (a2.xyann[0], a2.xyann[1] + mv[1] * 1.1* (z / 2 + 1))
                steps[i] += 1
            else:
                print '\tOK, these do not overlap: %s and %s' % (a1._text, a2._text)
        if not overlaps:
            break

fig.savefig('vcount-ecount-log-2.pdf')

plt.close(fig)

# transitivity sens.barplot
d = zip(*[(s, g.transitivity_undirected()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.transitivity_undirected())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'transitivity-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Graph global transitivity', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'transitivity-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Graph global transitivity', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

# diameter sens.barplot
d = zip(*[(s, g.diameter()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.diameter())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'diameter-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Graph diameter', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'diameter-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Graph diameter', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

# receptors sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['rec']])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec']))])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptors-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Number of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', y_break = (0.6, 0.09))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptors-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Number of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, y_break = (0.6, 0.09))

# receptors prop sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['rec']])/float(g.vcount())*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec'])/float(net.graph.vcount())*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorprop-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', y_break = (0.39, 0.12))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorprop-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, y_break = (0.39, 0.12))

# receptors coverage sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['rec']])/float(len(net.lists['rec']))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec'])/float(len(net.lists['rec']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorcov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of receptors covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorcov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of receptors covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

# transcription factors sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['tf']])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf']))])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfs-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Number of TFs', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y',
    y_break = (0.45, 0.2))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfs-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Number of TFs', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.45, 0.2))

# transcription factor prop sens.barplot
d = zip(*[(s, len([v for v in g.vs if v['tf']])/float(g.vcount())*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf'])/float(net.graph.vcount())*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfprop-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of transcription factors', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfprop-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of transcription factors', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

d = zip(*[(s, len([v for v in g.vs if v['tf']])/float(len(net.lists['tfs']))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf'])/float(len(net.lists['tfs']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfcov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of TFs covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.5, 0.1))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfcov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of TFs covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.5, 0.1))

# receptors & TFs coverage barplot together
d = zip(*[(omnipath, 
        sum(net.graph.vs['rec'])/float(len(net.lists['rec']))*100,
        sum(net.graph.vs['tf'])/float(len(net.lists['tfs']))*100,
        sum(net.graph.vs['kin'])/float(len(net.lists['kin']))*100,
        sum(net.graph.vs['dgb'])/float(len(net.lists['dgb']))*100
    )] + \
    [(s, 
        len([v for v in sep[s].vs if v['rec']])/float(len(net.lists['rec']))*100,
        len([v for v in sep[s].vs if v['tf']])/float(len(net.lists['tfs']))*100,
        len([v for v in sep[s].vs if v['kin']])/float(len(net.lists['kin']))*100,
        len([v for v in sep[s].vs if v['dgb']])/float(len(net.lists['dgb']))*100
    )\
    for s in vcount_ordr if s != omnipath])

d = dict(zip(['Database', 'Receptors', 'TFs', 'Kinases', 'Druggables'], 
    [list(dd) for dd in d]))
d = pd.DataFrame(d, index = d['Database'])

fig, ax = plt.subplots()
font_family = 'Helvetica Neue LT Std'
sns.set(font = font_family)
x = np.arange(len(d['Database']))
w = 0.2
c = {
    'rec': '#7AA0A1',
    'tf': '#C6909C',
    'kin': '#C5B26E',
    'dgb': '#9D8BB7'
}
rectsa = ax.bar(x, d['Receptors'], w, color = c['rec'], edgecolor = 'none')
rectsn = ax.bar(x + w, d['TFs'], w, color = c['tf'], edgecolor = 'none')
rectsn = ax.bar(x + 2 * w, d['Kinases'], w, color = c['kin'], edgecolor = 'none')
rectsn = ax.bar(x + 3 * w, d['Druggables'], w, color = c['dgb'], edgecolor = 'none')

lhandles = [mpl.patches.Patch(color=c['rec'], label='Receptors (all: %s)' %\
        locale.format('%d', len(net.lists['rec']), grouping = True)), 
    mpl.patches.Patch(color=c['tf'], label='Transcription factors (all: %s)' %\
        locale.format('%d', len(net.lists['tfs']), grouping = True)),
    mpl.patches.Patch(color=c['kin'], label='Kinases (all: %s)' %\
        locale.format('%d', len(net.lists['kin']), grouping = True)),
    mpl.patches.Patch(color=c['dgb'], label='Druggable proteins (all: %s)' %\
        locale.format('%d', len(net.lists['dgb']), grouping = True))]
leg = ax.legend(handles = lhandles)
ax.set_xlim(-0.5, 9.5)
ax.set_ylabel(r'% of proteins covered', fontsize = axis_lab_size * 0.56)
ax.set_xlabel('Pathway resources', fontsize = axis_lab_size * 0.56)
ax.set_xticks(x + 2 * w)
ax.set_xticklabels(d['Database'], rotation = 90, fontsize = lab_size[0] * 0.76)

plt.tight_layout()

fig.savefig('interesting-proteins-cov-by-db.pdf')

plt.close()


#d = d.sort(columns = 'Receptors', ascending = False)

# disease associations:
net.load_disgenet()
diss = dataio.get_disgenet()
dis = {}
for d in diss:
    ups = net.mapper.map_name(d['entrez'], 'entrez', 'uniprot')
    for u in ups:
        if u not in dis: dis[u] = []
        dis[u].append(d['disease'])

d = zip(*[(s, sum(len(x) for u, x in dis.iteritems() \
    if u in g.vs['name'])/float(sum(len(x) \
    for x in dis.values()))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(len(x) for x in net.graph.vs['dis']) / \
        float(sum(len(x) for x in dis.values()))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'discov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of disease-gene' + '\nassociations covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = 'y', 
    y_break = (0.65, 0.12))
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'discov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = r'% of disease-gene' + '\nassociations covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False,
    y_break = (0.65, 0.12))

# number of complexes sens.barplot
d = zip(*[(s, len(sens.complexes_in_network(g))) \
    for s, g in sep.iteritems()] + \
    [(omnipath, len(sens.complexes_in_network(net.graph)))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'complexes-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of complexes',
    #\n(out of %u)'%\
        #len(sens.complexes_in_network(net.graph)), 
    xlab = 'Pathway resources', order = 'y', y_break = (0.55, 0.15))
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'complexes-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of complexes',
    #\n(out of %u)'%\
        #len(sens.complexes_in_network(net.graph)), 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.55, 0.15))

# number of ptms sens.barplot
d = zip(*[(s, sum([len(e['ptm']) for e in g.es])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum([len(e['ptm']) for e in net.graph.es]))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ptms-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of PTMs',
    #\n(out of %u PTMs)'%\
        #sum([len(e['ptm']) for e in net.graph.es]), 
    xlab = 'Pathway resources', order = 'y', y_break = (0.55, 0.15))
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ptms-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of PTMs',
    #\n(out of %u PTMs)'%\
        #sum([len(e['ptm']) for e in net.graph.es]), 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.55, 0.15))

# number of ptms sens.barplot
d = zip(*[(s, len([1 for e in g.es if len(e['ptm']) != 0])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, len([1 for e in net.graph.es if len(e['ptm']) != 0]))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'havingptm-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Interactions having PTM', 
    xlab = 'Pathway resources', order = 'y', y_break = (0.49, 0.1))
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'havingptm-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Interactions having PTM', 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False, 
    y_break = (0.49, 0.1))

# TODO difference between COSMIC and Intogene

# cosmic cancer gene census coverage sens.barplot
d = zip(*[(s, len(set(net.lists['CancerGeneCensus']) & set(g.vs['name'])) / \
        float(len(net.lists['CancerGeneCensus']))*100) \
        for s, g in sep.iteritems()] + \
    [(omnipath, len(set(net.lists['CancerGeneCensus']) & set(net.graph.vs['name'])) / \
        float(len(net.lists['CancerGeneCensus']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ccgccov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of CGC genes', 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ccgccov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
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
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of Intogene genes', 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
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
sens.stacked_barplot(x = d[0], y = d[1:], 
    names = ['Positive', 'Negative', 'Unknown effect', 'Unknown direction'],
    data = None, fname = 'dirs-signes-by-db-wo-op-o-st.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Interactions', 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

sens.stacked_barplot(x = d[0], y = d[1:], 
    names = ['Positive', 'Negative', 'Unknown effect', 'Unknown direction'],
    data = None, fname = 'dirs-signes-by-db-wo-op-st.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Interactions', 
    xlab = 'Pathway resources', order = 'y')

# dirs & signs, grouped
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
d = dict(zip(['Database', 'Positive', 'Negative', 'Directed', 'Undirected'], 
    [list(dd) for dd in d]))
d = pd.DataFrame(d, index = vcount_ordr)
fig, ax = plt.subplots()
font_family = 'Helvetica Neue LT Std'
sns.set(font = font_family)
x = np.arange(len(d['Database']))
w = 0.2
c = {
    'positive': '#7AA0A1',
    'negative': '#C6909C',
    'directed': '#92C1D6',
    'undirected': '#C5B26E'
}
rectsp = ax.bar(x, d['Positive'], w, color = c['positive'], edgecolor = 'none')
rectsn = ax.bar(x + w, d['Negative'], w, color = c['negative'], edgecolor = 'none')
rectsd = ax.bar(x + 2 * w, d['Directed'], w, color = c['directed'], edgecolor = 'none')
rectsd = ax.bar(x + 3 * w, d['Undirected'], w, color = c['undirected'], edgecolor = 'none')
lhandles = [mpl.patches.Patch(color=c['positive'], label='Stimulation'), 
    mpl.patches.Patch(color=c['negative'], label='Inhibition'), 
    mpl.patches.Patch(color=c['directed'], label='Unknown effect'), 
    mpl.patches.Patch(color=c['undirected'], label='Undirected')]
ax.legend(handles = lhandles)
ax.set_xlim(-1, 27)
ax.set_ylabel('Interactions', fontsize = axis_lab_size * 0.66)
ax.set_xlabel('Pathway resources', fontsize = axis_lab_size * 0.66)
ax.set_xticks(x + 2 * w)
ax.set_xticklabels(d['Database'], rotation = 90, fontsize = lab_size[0] * 0.66)

plt.tight_layout()

fig.savefig('dirs-signes-by-db-wo-op-o.pdf')

plt.close()

stacked_barplot(x = d[0], y = d[1:], 
    names = ['positive', 'negative', 'unknown effect', 'unknown direction'],
    data = None, fname = 'dirs-signes-by-db-wo-op-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    ylab = 'Interactions', 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

stacked_barplot(x = d[0], y = d[1:], 
    names = ['positive', 'negative', 'unknown effect', 'unknown direction'],
    data = None, fname = 'dirs-signes-by-db-wo-op.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
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

# ## #

d = zip(*[(s, len(set(net.lists['Intogene']) & set(g.vs['name'])) / \
        float(len(net.lists['Intogene']))*100) \
        for s, g in sep.iteritems()] + \
    [(omnipath, len(set(net.lists['Intogene']) & set(net.graph.vs['name'])) / \
        float(len(net.lists['Intogene']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of IntOGen genes', 
    xlab = 'Pathway resources', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = 'light',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of IntOGen genes', 
    xlab = 'Pathway resources', order = vcount_ordr, desc = False)

# ## #
# PTMs in PTM resources and in the network
# ## #

ptms = {
    'Signor': dataio.load_signor_ptms(),
    'MIMP': dataio.get_mimp(),
    'PhosphoNetworks': dataio.get_phosphonetworks(),
    'phosphoELM': dataio.get_phosphoelm(),
    'dbPTM': dataio.get_dbptm(),
    'PhosphoSite': dataio.get_psite_phos(),
    'HPRD': dataio.get_hprd_ptms(),
    'Li2012': dataio.li2012_phospho()
}
ptms = {
    #'DEPOD': net.load_depod_dmi(return_raw = True),
    'Signor': net.load_signor_ptms(return_raw = True),
    'Li2012': net.load_li2012_ptms(return_raw = True),
    'HPRD': net.load_hprd_ptms(return_raw = True),
    'MIMP': net.load_mimp_dmi(return_raw = True),
    'PhosphoNetworks': net.load_pnetworks_dmi(return_raw = True),
    'PhosphoELM': net.load_phosphoelm(return_raw = True),
    'dbPTM': net.load_dbptm(return_raw = True),
    'PhosphoSite': net.load_psite_phos(return_raw = True)
}

d = zip(*[(s, len(uniqList(m)), 
            len([p for e in net.graph.es for p in e['ptm'] if s.lower() in [ps.lower() for ps in p.sources]])
        ) for s, m in ptms.items()] + \
        [('All', len(uniqList(flatList(ptms.values()))), len(uniqList(flatList(net.graph.es['ptm']))))])

d = dict(zip(['Database', 'All', 'In_network'], 
    [list(dd) for dd in d]))
d = pd.DataFrame(d, index = d['Database'])

d = d.sort(columns = 'All', ascending = False)

fig, ax = plt.subplots()
font_family = 'Helvetica Neue LT Std'
sns.set(font = font_family)
x = np.arange(len(d['Database']))
w = 0.44
c = {
    'all': '#7AA0A1',
    'inn': '#C6909C'
}
rectsa = ax.bar(x, d['All'], w, color = c['all'], edgecolor = 'none')
rectsn = ax.bar(x + w, d['In_network'], w, color = c['inn'], edgecolor = 'none')

lhandles = [mpl.patches.Patch(color=c['all'], label='All PTMs'), 
    mpl.patches.Patch(color=c['inn'], label='Mapped on OmniPath')]
ax.legend(handles = lhandles)
ax.set_xlim(-0.5, 9.5)
ax.set_ylabel('PTMs', fontsize = axis_lab_size * 0.66)
ax.set_xlabel('PTM resources', fontsize = axis_lab_size * 0.66)
ax.set_xticks(x + w)
ax.set_xticklabels(d['Database'], rotation = 90, fontsize = lab_size[0] * 0.66)

plt.tight_layout()

fig.savefig('ptms-by-ptmdb.pdf')

plt.close()
