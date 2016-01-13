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

from collections import Counter, OrderedDict
from itertools import groupby
import pypath
import copy
from itertools import chain
from pypath.data_formats import best, good, ugly, transcription
from pypath import dataio
from pypath.common import uniqList, flatList

import _sensitivity as sens

import numpy as np
import numpy
import numpy.ma as ma
from numpy.random import randn
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as hc
import hcluster as hc2
import matplotlib.patches as mpatches
from matplotlib.mlab import PCA

def prdb_tissue_expr(tissue, prdb, occurrence = 1, 
    group_function = lambda x: sum(x) / float(len(x))):
    nsamples = len(prdb.samples[tissue])
    occurrence = min(nsamples, occurrence) \
        if type(occurrence) is int \
        else nsamples * occurrence
    proteins_present = set([uniprot for uniprot, cnt in \
        Counter([uniprot for uniprots in \
            [expr.keys() for expr in \
                [prdb.expression[sample] for sample in prdb.samples[tissue]] \
            ] for uniprot in uniprots]).iteritems() \
        if cnt >= occurrence])
    expressions = dict([(uniprot, group_function(
        [prdb.expression[sample][uniprot] \
            for sample in prdb.samples[tissue] \
            if uniprot in prdb.expression[sample]])) \
        for uniprot in proteins_present])
    return expressions

def read_proteinatlas(mapper, fname = '/home/denes/Dokumentumok/pw/data/normal_tissue.csv'):
    result = {}
    with open(fname, 'r') as f:
        l = f.readline()
        for l in f:
            l = l.strip().replace('"', '').split(',')
            uniprots = mapper.map_name(l[0], 'ensg', 'uniprot')
            tissue = (l[1], l[2])
            if tissue not in result:
                result[tissue] = {}
            for u in uniprots:
                result[tissue][u] = l[-3:]
    return result

def clusterd(data, xlab = 'X label', ylab = 'Y label', method = 'ward', 
    fname = 'hc-test.pdf', fontface = 'sans-serif', legend = None,
    textcol = 'black', labelcols = None, orientation = 'right', **kwargs):
    mpl.rcParams['lines.linewidth'] = 0.5
    # fig, ax = plt.subplots()
    link = hc2.linkage(data, method = method)
    fig = plt.figure(figsize=(11.3, data.shape[0] * 0.17))
    ax = fig.gca()
    dend = hc2.dendrogram(link, orientation = orientation, **kwargs)
    ax.set_xlabel(xlab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    ax.set_ylabel(ylab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    # ax.yaxis.tick_right()
    ax.grid(b = None, visible = False, axis = 'y')
    if labelcols is not None:
        lbls = ax.get_ymajorticklabels() if orientation in ['left', 'right'] else \
            ax.get_xmajorticklabels()
        for lbl in lbls:
            lbl_text = lbl.get_text()
            if lbl_text in labelcols:
                lbl.set_color(labelcols[lbl_text])
        if type(legend) is dict:
            legend_patches = [mpatches.Patch(color=col, label=lab) \
                for lab, col in legend.iteritems()]
            ax.legend(handles = legend_patches)
    fig.tight_layout()
    fig.savefig(fname)

def rgb2hex(rgb):
    return '#%02x%02x%02x' % rgb

def hex2rgb(self, rgbhex):
    rgbhex = rgbhex.lstrip('#')
    lv = len(rgbhex)
    return tuple(int(rgbhex[i:i + 2], 16) for i in range(0, lv, 2))

def rgb1(self, rgb256):
    return rgb256 if not any([i > 1 for i in rgb256]) \
        else tuple([x / float(255) for x in rgb256])

def rgb256(rgb1):
    return rgb1 if any([i > 1.0 for i in rgb1]) \
        else tuple([x * 255.0 for x in rgb1])

mysql_gelati = (None,'mapping_gelati')
mysql_chembl = (None,'chembl_ebi')

net = pypath.Pypath(9606, mysql=mysql_gelati, name="demo")

#net.load_resources({'netpath': best['slk']})
#net.load_resources(lst={'ccmap': best['ccmap2']})
net.init_network(pfile = 'cache/default_plus_acsn_phospho.pickle')

# net.edge_loc()

ltps = OrderedDict()
with open('ltplist.csv', 'r') as f:
    for l in f:
        l = l.strip().split('\t')
        ltps[l[2]] = l

len(set(ltps.keys()) & set(net.graph.vs['name']))
len(set(ltps.keys()) - set(net.graph.vs['name']))

with open('ltpcols.csv', 'r') as f:
    ltpfcols = dict([tuple(l.split('\t')) for l in f.read().split('\n') if len(l) > 0])

ltpcols = dict([(l[1], ltpfcols[l[0]]) for l in ltps.values()])

with open('ltpcols2.csv', 'w') as f:
    f.write('\n'.join(['"%s"\t"%s"'%l for l in ltpcols.items()]))

net.graph.vs['ltp'] = [v['name'] in ltps for v in net.graph.vs]
net.genesymbol_labels()


for ltp, syn in ltps.iteritems():
    if ltp in net.graph.vs['name']:
        ltpnet = net.neighbourhood_network(net.graph.vs.find(name=ltp).index, second = True)
        trace = net.small_plot(ltpnet, name = 'ltp-%s'%syn[1])

net.get_proteomicsdb('denes', '33tengerimalaC')
prdb = net.proteomicsdb

# LTP neigborhood networks...
for tissue_id in net.proteomicsdb.tissues_loaded:
    tissue_name = [t['TISSUE_NAME'] for t in net.proteomicsdb.tissues \
        if t['TISSUE_ID'] == tissue_id][0]
    sys.stdout.write('\t:: Generating tissue specific network: %s\n' % tissue_name)
    sys.stdout.flush()
    tissuenet = net.prdb_tissue_network(tissue_id)
    sys.stdout.write('\t:: Signaling network of %s consists of '\
        '%u proteins and %u edges\n' % \
        (tissue_name, tissuenet.vcount(), tissuenet.ecount()))
    sys.stdout.flush()
    for ltp, syn in ltps.iteritems():
        if ltp in tissuenet.vs['name']:
            sys.stdout.write('\t\t:: Next protein: %s\n' % syn[1].split()[0])
            sys.stdout.flush()
            ltpnet = net.neighbourhood_network(tissuenet.vs.find(name=ltp).index, 
                second = True)
            ltpname = net.mapper.map_name(ltp, 'uniprot', 'genesymbol')
            ltpname = ltp if len(ltpname) == 0 else ltpname[0]
            trace = net.small_plot(ltpnet, name = 'ltp-%s-%s' % \
                (syn[1].split()[0], tissue_id), 
                title_text = 'Protein %s in %s' % (ltpname, tissue_name))

# building ltp expression by tissue dict
# and tissue id -- name conversion dict
cellines = [t['TISSUE_NAME'] for t in prdb.tissues if t['TISSUE_CATEGORY'] == 'cell line']
pd = {}
pdcl = []
ltpexp = {}
tissnm = {}
for tissue_id in net.proteomicsdb.tissues_loaded:
    tissue_name = [t['TISSUE_NAME'] for t in net.proteomicsdb.tissues \
        if t['TISSUE_ID'] == tissue_id][0]
    expr = prdb_tissue_expr(tissue_id, net.proteomicsdb)
    if tissue_name in cellines and tissue_name != 'Mixed':
        pdcl.append(tissue_id)
        pd[tissue_id] = dict([(u, e) \
            for u, e in expr.iteritems() if u in ltps])
    tissnm[tissue_id] = tissue_name
    ltpexp[tissue_id] = [u for u, e in expr.iteritems() \
        if u in ltps and e > 0.0]

tissid = dict([(v, k) for k,v in tissnm.iteritems()])

# Human Protein Atlas

pa = read_proteinatlas(net.mapper)

pavals = {
    'High': 3, 
    'Not detected': 0, 
    'Medium': 2, 
    'Low': 1
}

# gene ontology
gos = dataio.get_go_quick(slim = True)
gos2 = dataio.get_go_quick(
    slim = 'http://www.geneontology.org/ontology/subsets/goslim_pir.obo')


ltpfam = list(set([l[0] for l in ltps.values()]))
ltpnm = dict([(l[1], u) for u, l in ltps.iteritems()])
families = dict([(fam, []) for fam in ltpfam])
for l in ltps.values():
    families[l[0]].append(l[1])

which_family = dict([(l[1], l[0]) for l in ltps.values()])

mterms = ['metabol', 'transport', 'biosynth', 'catabol']
sterms = ['signal', 'communicat', 'response', 'immune', 'regulat', 
    'cell cycl', 'cell diff']

### exporting ProteomicsDB expression table:
ltp_prdb_cellines = 'ltp_prdb_celllines.tab'
ltp_prdb_hdr = ['GeneSymbol'] + ['"%s"'%tissnm[t] \
    for t in sorted(pd.keys()) if t in cellineswltp]
with open(ltp_prdb_cellines, 'w') as f:
    f.write('\t'.join(ltp_prdb_hdr) + '\n')
    for u, ltp in ltps.iteritems():
        prdbvals = ['0.0' if u not in pd[t] else \
                str(pd[t][u]) for t in sorted(pd.keys()) if t in cellineswltp]
        if u in ltps_in_pd:
            f.write('\t'.join(['"%s"'%ltp[1]] + prdbvals) + '\n')

### exporting expression table
ltp_prdb_file = 'ltp_proteomicsdb_v2.tab'
ltp_hpa_normal = 'ltp_hpa_normal.tab'
ltp_prdb_hdr = ['UniProt', 'GeneSymbol', 'Family', 'GOslim BP', 'GOslim BP PIR',
    'Signaling', 'Metabolism', 'In pathways', 'Samples of occurrence in ProteomicsDB', 
    'Samples of occurrence in HPA'] + sorted(tissid.keys()) + \
    ['%s :: %s'%t for t in sorted(pa.keys())]
ltp_hpa_hdr = ['GeneSymbol'] + ['"%s :: %s"' % \
    tuple([t[0].capitalize(), t[1].replace('endometrial ', '').\
        replace('/', ' ').replace('-', ' ')]) \
    for t in sorted(pa.keys()) if t in tis_mes_hc]
with open(ltp_prdb_file, 'w') as f:
    with open(ltp_hpa_normal, 'w') as h:
        # header:
        f.write('\t'.join(ltp_prdb_hdr) + '\n')
        h.write('\t'.join(ltp_hpa_hdr) + '\n')
        # data rows per ltp:
        for u, ltp in ltps.iteritems():
            go = '' if u not in gos['terms']['P'] \
                else ';'.join([gos['names'][g] \
                    for g in gos['terms']['P'][u]])
            go2 = '' if u not in gos2['terms']['P'] \
                else ';'.join([gos2['names'][g] \
                    for g in gos2['terms']['P'][u]])
            s = 1 if len([i for i in sterms if i in go or i in go2]) > 0 else 0
            m = 1 if len([i for i in mterms if i in go or i in go2]) > 0 else 0
            p = 1 if u in net.graph.vs['name'] else 0
            prdbvals = [1 if u in ltpexp[tissid[t]] else 0 for t in sorted(tissid.keys())]
            hpavals_hc = ['na' if u not in pa[t] else 'uc' if pa[t][u][2] == 'Uncertain' else \
                str(pavals[pa[t][u][0]]) for t in sorted(pa.keys()) if t in tis_mes_hc]
            hpavals = ['na' if u not in pa[t] else 'uc' if pa[t][u][2] == 'Uncertain' else \
                str(pavals[pa[t][u][0]]) for t in sorted(pa.keys())]
            if ltp[2] in mes_hc:
                h.write('\t'.join(['"%s"'%ltp[1]] + hpavals) + '\n')
            f.write('\t'.join([u, ltp[1], ltp[0], go, go2, 
                '%u'%s, '%u'%m, '%u'%p, '%u'%sum(prdbvals), 
                # expressed = supportive & detected
                '%s'%len([1 for t in pa.values() if u in t and t[u][0] != 'Not detected' \
                    and t[u][2] == 'Supportive'])] + \
                ['%s'%v for v in prdbvals] + ['%s'%v for v in hpavals]) + '\n')
        # bottom sum row:
        f.write('\t'.join(['All LTPs expressed', '', '', '', '', '', '', '', 
            '%u'%len(set([it for sl in ltpexp.values() for it in sl])),
            '%u'%len([u for u in ltps.keys() \
                if len([1 for t in pa.values() \
                    if u in t and t[u][0] != 'Not detected' \
                    and t[u][2] == 'Supportive']) != 0])] + \
            ['%u'%len(ltpexp[tissid[t]]) for t in sorted(tissid.keys())] + \
            ['%u'%len([1 for u, e in pa[t].iteritems() \
                if u in ltps and e[0] != 'Not detected' and e[2] == 'Supportive']) \
                for t in sorted(pa.keys())] \
            ) + '\n')

### building numpy arrays:
prdb_ltp_lst = []
hpa_ltp_lst = []
labels = []
# data rows per ltp:
for u, ltp in ltps.iteritems():
    prdb_ltp_lst.append([1 if u in ltpexp[tissid[t]] else 0 \
        for t in sorted(tissid.keys())])
    hpa_ltp_lst.append([0 if u not in pa[t] else pavals[pa[t][u][0]] \
        for t in sorted(pa.keys())])
    labels.append(ltp[1])

prdb_ltp = np.array(prdb_ltp_lst)
hpa_ltp = np.array(hpa_ltp_lst)

pal = sns.color_palette('husl', 10)
groups = list(set([l[0] for l in ltps.values()]))
colors = dict([(l[1], pal[groups.index(l[0])]) for l in ltps.values()])
legcol = dict([(g, pal[i]) for i, g in enumerate(groups)])

### Human Protein Atlas numbers:
def exp_hc(uniprot, padata):
    nod = 'Not detected'
    sup = 'Supportive'
    return uniprot in padata and padata[uniprot][0] != nod and padata[uniprot][2] == sup

def measured(uniprot, padata):
    sup = 'Supportive'
    return uniprot in padata and padata[uniprot][2] == sup

nod = 'Not detected'
sup = 'Supportive'
unc = 'Uncertain'
in_pa = uniqList(flatList([p.keys() for p in pa.values()]))
# LTPs not in HPA:
len(set(ltps.keys()) - set(in_pa)) # 22
ltps_in_pa = list(set(ltps.keys()) & set(in_pa))
# LTPs with only uncertain values:
[len(set([t[p][2] for t in pa.values() if p in t]) \
    - set([unc])) for p in ltps_in_pa].count(0) # 57
# LTPs with at least one high confidence expression data:
[len(set([t[p][2] for t in pa.values() if p in t]) \
    - set([unc])) for p in ltps_in_pa].count(1) # 45
mes_hc = [p for p in ltps_in_pa if len(set([t[p][2] for t in pa.values() if p in t]) \
    - set([unc])) > 0]
# tissues with any LTPs measured with high confidence:
len([t for t, d in pa.iteritems() if \
    len([p for p in ltps.keys() if p in d and d[p][2] == sup]) > 0]) # 82
tis_mes_hc = [t for t, d in pa.iteritems() if \
    len([p for p in ltps.keys() if p in d and d[p][2] == sup]) > 0]
# tissues with any LTPs expressed with high confidence:
len([t for t, d in pa.iteritems() if \
    len([p for p in ltps.keys() if \
        p in d and d[p][2] == sup and d[p][0] != nod]) > 0]) # 82

# number of LTPs expressed with high confidence by tissue:
nltpbytis = zip(*[(('%s, %s'%t).capitalize(), len([p for p in ltps_in_pa \
        if exp_hc(p, d)])) \
    for t, d in pa.iteritems()])
barplot(x = nltpbytis[0], y = nltpbytis[1], 
    data = None, fname = 'numof-ltps-by-tissue.pdf', lab_size = 7, 
    ylab = 'Number of LTPs expressed', 
    xlab = 'Tissues in HPA', order = 'y')

### normalized by measured
# percentage of LTPs expressed with high confidence by tissue:
pltpmesbytis = zip(*[(('%s, %s'%t).capitalize(), len([p for p in ltps_in_pa \
        if exp_hc(p, d)]) / \
        float(len([p for p in ltps_in_pa if measured(p, d)])) * 100) \
    for t, d in pa.iteritems()])
barplot(x = pltpmesbytis[0], y = pltpmesbytis[1], 
    data = None, fname = 'pctof-ltps-measured-by-tissue.pdf', lab_size = 7, 
    ylab = 'Percentage of measured LTPs found expressed', 
    xlab = 'Tissues in HPA', order = 'y')

# number of LTPs expressed with high confidence by organ:
organs = uniqList([t[0] for t in pa.keys()])
nltpbyorg = zip(*[(o.capitalize(), \
    len(uniqList(flatList(\
        [[p for p in ltps_in_pa if exp_hc(p, d)] for t, d in pa.iteritems() if t[0] == o]))))\
    for o in organs])
barplot(x = nltpbyorg[0], y = nltpbyorg[1], 
    data = None, fname = 'numof-ltps-by-organ.pdf', lab_size = 11, 
    ylab = 'Number of LTPs expressed', 
    xlab = 'Organs in HPA', order = 'y')

# number of LTPs expressed with high confidence by category:
nltpbycat = zip(*[(grp.capitalize(), \
    len(uniqList(flatList(\
        [[p for p in ltps_in_pa if exp_hc(p, d)] for t, d in pa.iteritems() if t in tiss]))))\
    for grp, tiss in tiss_grp.iteritems()])
sens.barplot(x = nltpbycat[0], y = nltpbycat[1], 
    data = None, fname = 'numof-ltps-by-category.pdf', lab_size = 11, 
    ylab = 'Number of LTPs expressed', 
    xlab = 'Tissue categories', order = 'y')

### measured
# number of LTPs expressed with high confidence by tissue:
nltpmesbytis = zip(*[(('%s, %s'%t).capitalize(), len([p for p in ltps_in_pa \
        if exp_hc(p, d)]), \
        len([p for p in ltps_in_pa if measured(p, d)]) - \
        len([p for p in ltps_in_pa if exp_hc(p, d)])) \
    for t, d in pa.iteritems()])
stacked_barplot(x = nltpmesbytis[0], y = nltpmesbytis[1:], names = ['Expressed', 'Measured'],
    data = None, fname = 'numof-ltps-measured-by-tissue.pdf', lab_size = 7, 
    ylab = 'Number of LTPs measured / found expressed', 
    xlab = 'Tissues in HPA', order = 0)
# number of LTPs expressed with high confidence by organ:
organs = uniqList([t[0] for t in pa.keys()])
nltpmesbyorg = zip(*[(o.capitalize(), \
    len(uniqList(flatList(\
        [[p for p in ltps_in_pa if exp_hc(p, d)] for t, d in pa.iteritems() if t[0] == o]))),
    len(uniqList(flatList(\
        [[p for p in ltps_in_pa if measured(p, d)] \
            for t, d in pa.iteritems() if t[0] == o]))) -\
    len(uniqList(flatList(\
        [[p for p in ltps_in_pa if exp_hc(p, d)] for t, d in pa.iteritems() if t[0] == o]))))\
    for o in organs])
stacked_barplot(x = nltpmesbyorg[0], y = nltpmesbyorg[1:], names = ['Expressed', 'Measured'],
    data = None, fname = 'numof-ltps-measured-by-organ.pdf', lab_size = 7, 
    ylab = 'Number of LTPs measured / found expressed', 
    xlab = 'Organs in HPA', order = 0)

# normalized (percentage) per tissue:
pltpmesbyorg = zip(*[(o.capitalize(), \
    len(uniqList(flatList(\
        [[p for p in ltps_in_pa if exp_hc(p, d)] for t, d in pa.iteritems() if t[0] == o]))) / \
    float(len(uniqList(flatList(\
        [[p for p in ltps_in_pa if measured(p, d)] \
            for t, d in pa.iteritems() if t[0] == o])))) * 100.0)\
    for o in organs])
barplot(x = pltpmesbyorg[0], y = pltpmesbyorg[1],
    data = None, fname = 'pctof-ltps-measured-by-organ.pdf', lab_size = 11, 
    ylab = 'Percentage of measured LTPs found expressed', 
    xlab = 'Organs in HPA', order = 'y')

# normalized by category:
pltpbycat = zip(*[(grp.capitalize(), \
    len(uniqList(flatList(\
        [[p for p in ltps_in_pa if exp_hc(p, d)] \
        for t, d in pa.iteritems() if t in tiss]))) / \
    float(len(uniqList(flatList(\
        [[p for p in ltps_in_pa if measured(p, d)] \
        for t, d in pa.iteritems() if t in tiss])))) * 100.0 \
    )\
    for grp, tiss in tiss_grp.iteritems()])
sens.barplot(x = pltpbycat[0], y = pltpbycat[1], 
    data = None, fname = 'pctof-ltps-by-category-b.pdf', lab_size = 11, 
    ylab = 'Percentage of measured LTPs\nfound expressed', 
    xlab = 'Tissue categories', order = 'y')

def mostCommon(lst):
    cnt = Counter(lst)
    mxo = max(cnt.values())
    return [v for v, c in cnt.iteritems() if c == mxo]

with open('ltp-hpa-cat.tab', 'w') as f:
    f.write('\t'.join(['"GeneSymbol"'] + ['"%s"'%t.capitalize() \
        for t in sorted(tiss_grp.keys())]) + '\n')
    for ltp in ltps.values():
        if ltp[2] in mes_hc:
            grp_exp = [max(mostCommon(
                [pavals[pa[tis][ltp[2]][0]] if measured(ltp[2], pa[tis]) \
                    else -1 if ltp[2] in pa[tis] and pa[tis][ltp[2]][2] == unc \
                    else -2 \
                    for tis in tiss])) \
            for grp, tiss in tiss_grp.iteritems()]
            if any([i > 0 for i in grp_exp]):
                f.write('\t'.join(['"%s"'%ltp[1].strip()] + \
                    [str(i) for i in grp_exp]) + '\n')

### pleiotropy
# pleiotropy barplot, HPA

ltpplebytis = sorted([(ltp[1], which_family[ltp[1]], 
    len(uniqList([t for t, d in pa.iteritems() if exp_hc(p, d)]))/\
    float(len(tis_mes_hc)) * 100.0)\
    for p, ltp in ltps.iteritems() if p in mes_hc], 
    key = lambda x: (x[1], x[2]), reverse = True)

fig, ax = plt.subplots()
bar_width = 0.9
lab_angle = 90
xvals = np.arange(len(ltpplebytis))
sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
xticklabs = plt.bar(xvals, zip(*ltpplebytis)[2], bar_width, 
    color = [ltpcols[ltp] for ltp in zip(*ltpplebytis)[0]])
plt.xticks(xvals + bar_width, zip(*ltpplebytis)[0])
plt.setp(ax.xaxis.get_majorticklabels(), rotation = lab_angle)
plt.xlabel('LTPs by families')
plt.ylabel('Percentage of tissues in HPA')
plt.title('Pleiotropy of LTPs\nmeasured by high confidence in HPA normal tissues')
ax.grid(axis = 'x', alpha = 0.0)
legend_patches = [mpatches.Patch(color=col, label=lab) \
    for lab, col in ltpfcols.iteritems()]
ax.legend(handles = legend_patches, loc = (0.6, 0.2))
plt.axis('tight')
fig.tight_layout()
fig.savefig('ltp-hpa-pleiotropy.pdf')
plt.close(fig)

# pleiotropy in ProteomicsDB cell lines:
ltpplebytis = sorted([(ltp[1], which_family[ltp[1]], 
    len(uniqList([t for t, d in pd.iteritems() if p in d and d[p] > 0.0]))/\
    float(len(cellineswltp)) * 100.0)\
    for p, ltp in ltps.iteritems() if p in ltps_in_pd], 
    key = lambda x: (x[1], x[2]), reverse = True)

fig, ax = plt.subplots(figsize = (12,6))
bar_width = 0.9
lab_angle = 90
lab_size = 5
xvals = np.arange(len(ltpplebytis))
sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
xticklabs = plt.bar(xvals, zip(*ltpplebytis)[2], bar_width, 
    color = [ltpcols[ltp] for ltp in zip(*ltpplebytis)[0]])
plt.xticks(xvals + bar_width, zip(*ltpplebytis)[0])
plt.setp(ax.xaxis.get_majorticklabels(), rotation = lab_angle)
plt.xlabel('LTPs by families')
plt.ylabel('Percentage of\ncell lines in ProteomicsDB')
plt.title('Pleiotropy of LTPs\nmeasured in ProteomicsDB cell lines')
ax.grid(axis = 'x', alpha = 0.0)
legend_patches = [mpatches.Patch(color=col, label=lab) \
    for lab, col in ltpfcols.iteritems()]

for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(lab_size)

ax.legend(handles = legend_patches, loc = (0.5, 0.2), prop = {'size': 10})
plt.axis('tight')
fig.tight_layout()
fig.savefig('ltp-proteomicsdb-pleiotropy.pdf')
plt.close(fig)

# pie charts:
organ_tis = dict([(org if not org.split()[-1].isdigit() else ' '.join(org.split()[:-1]), []) \
    for org in organs])

for t in pa.keys():
    org = t[0] if not t[0].split()[-1].isdigit() else ' '.join(t[0].split()[:-1])
    if org in organ_tis:
        organ_tis[org].append(t)

org_prop = dict([(org, 
    OrderedDict(dict([(ltpf, sum([int(bool(sum([int(exp_hc(ltpnm[ltp], pa[tis])) \
                for tis in tiss]))) \
            for ltp in families[ltpf]])) \
        for ltpf in sorted(families.keys())]))) \
    for org, tiss in organ_tis.iteritems()])

for org, prop in org_prop.iteritems():
    fig = plt.figure(figsize = (6, 6))
    sns.set_context('paper')
    sns.set(font = 'Helvetica Neue LT Std')
    patches, texts = plt.pie(prop.values(), 
        labels = [ltpf if num > 0 else '' for ltpf, num in prop.iteritems()],
        colors = [ltpfcols[ltpf] for ltpf in prop.keys()], 
        wedgeprops = { 'linewidth' : 0.0 })
    nul = [txt.set_fontsize(14) for txt in texts]
    titl = plt.title(org.capitalize())
    titl.set_fontsize(21)
    fig.tight_layout()
    plt.savefig('pie-%s.pdf'%org.replace(' ', '_'))
    plt.close(fig)
    # ## #
    fig = plt.figure(figsize = (6, 6))
    sns.set_context('paper')
    sns.set(font = 'Helvetica Neue LT Std')
    patches, texts = plt.pie(prop.values(), 
        labels = None,
        colors = [ltpfcols[ltpf] for ltpf in prop.keys()], 
        wedgeprops = { 'linewidth' : 0.0 })
    nul = [txt.set_fontsize(14) for txt in texts]
    fig.tight_layout()
    plt.savefig('pie-%s-%u.svg' % (org.replace(' ', '_'), sum(prop.values())))
    plt.close(fig)


#### #######################################
#### #######################################

### ProteomicsDB numbers:

# LTPs in ProteomicsDB cell lines:
ltps_in_pd = uniqList(flatList([p.keys() for p in pd.values()])) # 103
# number of cell lines:
len(uniqList(cellines)) # 141
# LTPs not in ProteomicsDB:
len(set(ltps.keys()) - set(ltps_in_pd)) # 22
# cell lines with LTP detected:
cellineswltp = [t for t, d in pd.iteritems() if len(d) > 0]
len(cellineswltp) # 98

# number of LTPs expressed with high confidence by tissue:
nltpbycel = zip(*[(tissnm[t], len(d)) \
    for t, d in pd.iteritems()])
sens.barplot(x = nltpbycel[0], y = nltpbycel[1], 
    data = None, fname = 'numof-ltps-by-cel.pdf', lab_size = 5, 
    ylab = 'Number of LTPs expressed', 
    xlab = 'Cell lines in ProteomicsDB', order = 'y')

nltpbycel = zip(*[(tissnm[t], len(d)) \
    for t, d in pd.iteritems() if len(d) > 0])
sens.barplot(x = nltpbycel[0], y = nltpbycel[1], 
    data = None, fname = 'numof-ltps-by-cel-nonzero.pdf', lab_size = 5, 
    ylab = 'Number of LTPs expressed', 
    xlab = 'Cell lines in ProteomicsDB', order = 'y')


### pleiotropy
# pleiotropy barplot, HPA
ltpplebytis = zip(*[(ltp[1], len(uniqList([t for t, d in pd.iteritems() if d > 0.0]))/\
    float(len(tis_mes_hc)))\
    for p, ltp in ltps.iteritems() if p in mes_hc])

barplot(ltpplebytis[0], ltpplebytis[1], 
    data = None, xlab = 'Lipid Transfer Proteins', 
    ylab = 'Pleiotropy in HPA tissues', 
    fname = 'LTPs_HPA_tissues_pleiotropy.pdf', lab_size = 9, 
    order = 'y', legend = legcol)

fig, ax = barplot(ltpplebytis[0], ltpplebytis[1], 
    hue = [which_family[l] for l in ltpplebytis[0]], 
    data = None, xlab = 'Lipid Transfer Proteins', 
    ylab = 'Pleiotropy in HPA tissues', 
    palette = dict([(which_family[lab], colors[lab]) for lab in ltpplebytis[0]]), 
    fname = 'LTPs_HPA_tissues_pleiotropy2.pdf', lab_size = 9, 
    order = None, legend = legcol, fin = False)


####
####

### clustering HPA by tissue:
data = ma.masked_invalid(
data = np.array([[pavals[d[p][0]] \
        if measured(p, d) else 0.0 for t, d in pa.iteritems()] \
    for p in mes_hc])
    )
clusterd(data, 
    labels = [ltps[p][1] for p in mes_hc], orientation = 'left', labelcols = colors, 
    legend = legcol, ylab = 'Lipid transfer proteins', method = 'ward', 
    xlab = 'Similarity of LTP expression patterns HPA tissues',
    fname = 'LTPs_HPA_ward_hc.pdf')

### old clusterings:
clusterd(prdb_ltp, labels = labels, orientation = 'left', labelcols = colors, 
    legend = legcol, ylab = 'Lipid transfer proteins', method = 'ward', 
    xlab = 'Similarity of cell/tissue expression patterns in ProteomicsDB')

def prdb_filter(prdb, ltpexp, ltps, categories = ['tissue'], 
        hpa = None, hpa_levels = None):
    hpa_levels = hpa_levels if hpa_levels is not None else ['Low', 'Medium', 'High']
    tissue_cat = dict([(t['TISSUE_ID'], t['TISSUE_CATEGORY']) for t in prdb.tissues])
    prdb_ltp_lst = []
    labels = []
    # data rows per ltp:
    for u, ltp in ltps.iteritems():
        this_ltp = [1 if u in ltpexp[tissid[t]] else 0 \
            for t, tid in sorted(tissid.iteritems()) \
            if tissue_cat[tid] in categories] if ltpexp is not None else []
        if hpa is not None:
            this_ltp += [0 if u not in hpa[t] or hpa[t][u][0] in hpa_levels else 1 \
                for t in sorted(hpa.keys())]
        if sum(this_ltp) > 1:
            prdb_ltp_lst.append(this_ltp)
            labels.append(ltp[1])
    cnonzero = []
    for i in xrange(len(prdb_ltp_lst[0])):
        if sum([t[i] for t in prdb_ltp_lst]) != 0:
            cnonzero.append(i)
    prdb_ltp_lst = [[t[i] for i in cnonzero] for t in prdb_ltp_lst]
    return numpy.array(prdb_ltp_lst), labels

def prdb_by_tissue(prdb, ltpexp, ltps, categories = ['tissue']):
    tissue_cat = dict([(t['TISSUE_ID'], t['TISSUE_CATEGORY']) for t in prdb.tissues])
    tissue_nam = dict([(t['TISSUE_ID'], t['TISSUE_NAME']) for t in prdb.tissues])
    prdb_tis_lst = []
    labels = []
    for tis in tissnm.keys():
        this_tis = [1 if u in ltpexp[tis] else 0 for u in ltps.keys()]
        if sum(this_tis) > 1 and tissue_cat[tis] in categories:
            prdb_tis_lst.append(this_tis)
            labels.append(tissue_nam[tis])
    return numpy.array(prdb_tis_lst), labels

def hpa_by_tissue(hpa, pavals, ltps):
    hpa_tis_lst = []
    labels = []
    for tis, exp in hpa.iteritems():
        this_tis = [0 if u not in exp else pavals[exp[u][0]] for u in ltps.keys()]
        if sum(this_tis) > 1:
            hpa_tis_lst.append(this_tis)
            labels.append(' :: '.join(tis))
    return numpy.array(hpa_tis_lst), labels
    #return hpa_tis_lst, labels

data, labels = prdb_filter(net.proteomicsdb, ltpexp, ltps, hpa = pa)
clusterd(data, labels = labels, orientation = 'left', labelcols = colors, 
    legend = legcol, ylab = 'Lipid transfer proteins', method = 'ward', 
    xlab = 'Similarity of tissue expression patterns in ProteomicsDB & HPA tissues',
    fname = 'LTPs_ProteomicsDB_tissues_HPA_ward.pdf')

data, labels = prdb_filter(net.proteomicsdb, ltpexp, ltps)
clusterd(data, labels = labels, orientation = 'left', labelcols = colors, 
    legend = legcol, ylab = 'Lipid transfer proteins', method = 'ward', 
    xlab = 'Similarity of tissue expression patterns in ProteomicsDB tissues',
    fname = 'LTPs_ProteomicsDB_tissues_ward.pdf')

# pleiotropy barplot, ProteomicsDB
prdb_ple = np.array([x/float(data.shape[1]) for x in np.sum(data, axis = 1)])
barplot(labels, prdb_ple, data = None, xlab = 'Lipid Transfer Proteins', 
    ylab = 'Pleiotropy in ProteomicsDB tissues', color = [colors[lab] for lab in labels], 
    fname = 'LTPs_ProteomicsDB_tissues_pleiotropy.pdf', lab_size = 7, 
    order = 'y', legend = legcol)

# PCA
prdb_pca = PCA(data)
scatterplot([prdb_pca.Y.take(0, axis = 1)], [prdb_pca.Y.take(1, axis = 1)], 
    'PCA 1', 'PCA 2', 'LTPs_ProteomicsDB_tissues_PCA.pdf', [[colors[lab] for lab in labels]], 
    fontface = 'Helvetica Neue LT Std', legend = legcol, 
    textcol = 'black', xlog = False, ylog = False, ylim = None)

# pleiotropy barplot, HPA
data, labels = prdb_filter(net.proteomicsdb, None, ltps, hpa = pa)
hpa_ple = np.array([x/float(data.shape[1]) for x in np.sum(data, axis = 1)])
barplot(labels, hpa_ple, data = None, xlab = 'Lipid Transfer Proteins', 
    ylab = 'Pleiotropy in HPA tissues', color = [colors[lab] for lab in labels], 
    fname = 'LTPs_HPA_tissues_pleiotropy.pdf', lab_size = 7, 
    order = 'y', legend = legcol)

# PCA, HPA
hpa_pca = PCA(data)
scatterplot([hpa_pca.Y.take(0, axis = 1)], [hpa_pca.Y.take(1, axis = 1)], 
    'PCA 1', 'PCA 2', 'LTPs_HPA_tissues_PCA.pdf', [[colors[lab] for lab in labels]], 
    fontface = 'Helvetica Neue LT Std', legend = legcol, 
    textcol = 'black', xlog = False, ylog = False, ylim = None)

hpa_ltp_lst = []
labels = []
# data rows per ltp:
for u, ltp in ltps.iteritems():
    this_ltp = [0 if u not in pa[t] else pavals[pa[t][u][0]] \
        for t in sorted(pa.keys())]
    if sum(this_ltp) > 1:
        hpa_ltp_lst.append(this_ltp)
        labels.append(ltp[1])

hpa_ltp = numpy.array(hpa_ltp_lst)

clusterd(hpa_ltp, labels = labels, orientation = 'left', labelcols = colors, 
    legend = legcol, ylab = 'Lipid transfer proteins', method = 'ward', 
    xlab = 'Similarity of tissue expression patterns in HPA tissues',
    fname = 'LTPs_HPA_ward.pdf')

prdb_ltp_lst = []
labels = []
# data rows per ltp:
for u, ltp in ltps.iteritems():
    this_ltp = [1 if u in ltpexp[tissid[t]] else 0 \
        for t in sorted(tissid.keys())]
    if sum(this_ltp) > 1:
        prdb_ltp_lst.append(this_ltp)
        labels.append(ltp[1])

prdb_ltp = numpy.array(prdb_ltp_lst)

clusterd(prdb_ltp, labels = labels, orientation = 'left', labelcols = colors, 
    legend = legcol, ylab = 'Lipid transfer proteins', method = 'ward', 
    xlab = 'Similarity of tissue expression patterns in ProteomicsDB tissues and cell lines',
    fname = 'LTPs_ProteomicsDB_all_ward.pdf')

### by tissue:
data, labels = prdb_by_tissue(net.proteomicsdb, ltpexp, ltps)
clusterd(data, labels = labels, orientation = 'left', ylab = 'Tissues', method = 'ward', 
    xlab = 'Similarity of LTP expression patterns in ProteomicsDB tissues',
    fname = 'tissues_ProteomicsDB_tissues_ward.pdf')

data, labels = hpa_by_tissue(pa, pavals, ltps)
clusterd(data, labels = labels, orientation = 'left', ylab = 'Tissues', method = 'ward', 
    xlab = 'Similarity of LTP expression patterns in HPA tissues',
    fname = 'tissues_HPA_tissues_ward.pdf')

# grouping of cell types / tissues:
tiss_grp = OrderedDict([
('connective tissue', [
    ('skin 1', 'fibroblasts'),
    ('soft tissue 2', 'fibroblasts'),
    ('soft tissue 1', 'fibroblasts'),
    ('soft tissue 2', 'chondrocytes'),
    ('soft tissue 1', 'chondrocytes'),
    ('breast', 'myoepithelial cells'),
    ('ovary', 'ovarian stroma cells'),
    ('endometrium 1', 'cells in endometrial stroma'),
    ('endometrium 2', 'cells in endometrial stroma')]),
('glandular cells', [
    ('epididymis', 'glandular cells'),
    ('prostate', 'glandular cells'),
    ('endometrium 1', 'glandular cells'),
    ('endometrium 2', 'glandular cells'),
    ('fallopian tube', 'glandular cells'),
    ('breast', 'glandular cells'),
    ('seminal vesicle', 'glandular cells')]),
('lymphatic and immune tissue', [
    ('tonsil', 'germinal center cells'),
    ('tonsil', 'non-germinal center cells'),
    ('spleen', 'cells in red pulp'),
    ('spleen', 'cells in white pulp'),
    ('tonsil', 'squamous epithelial cells'),
    ('lymph node', 'germinal center cells'),
    ('lymph node', 'non-germinal center cells'),
    ('appendix', 'lymphoid tissue'),
    ('bone marrow', 'hematopoietic cells'),
    ('lung', 'macrophages'),
    ('skin 1', 'Langerhans')]),
('gastrointestinal glandules', [
    ('duodenum', 'glandular cells'),
    ('colon', 'glandular cells'),
    ('stomach 2', 'glandular cells'),
    ('rectum', 'glandular cells'),
    ('appendix', 'glandular cells'),
    ('small intestine', 'glandular cells'),
    ('stomach 1', 'glandular cells'),
    ('salivary gland', 'glandular cells'),
    ('pancreas', 'exocrine glandular cells'),
    ('gallbladder', 'glandular cells')]),
('nervous tissue', [
    ('cerebral cortex', 'neuropil'),
    ('lateral ventricle', 'neuronal cells'),
    ('cerebellum', 'Purkinje cells'),
    ('cerebellum', 'cells in granular layer'),
    ('cerebellum', 'cells in molecular layer'),
    ('cerebral cortex', 'neuronal cells'),
    ('hippocampus', 'neuronal cells'),
    ('soft tissue 2', 'peripheral nerve'),
    ('soft tissue 1', 'peripheral nerve'),
    ('colon', 'peripheral nerve/ganglion')]),
('endocrine cells', [
    ('parathyroid gland', 'glandular cells'),
    ('thyroid gland', 'glandular cells'),
    ('adrenal gland', 'glandular cells'),
    ('pancreas', 'islets of Langerhans'),
    ('testis', 'Leydig cells')]),
('mucosal epithelia', [
    ('esophagus', 'squamous epithelial cells'),
    ('vagina', 'squamous epithelial cells'),
    ('cervix', ' uterine'),
    ('oral mucosa', 'squamous epithelial cells')]),
('lung', [
    ('lung', 'pneumocytes'),
    ('bronchus', 'respiratory epithelial cells'),
    ('nasopharynx', 'respiratory epithelial cells')]),
('internal epithelia', [
    ('liver', 'bile duct cells'),
    ('colon', 'endothelial cells'),
    ('urinary bladder', 'urothelial cells'),
    ('cerebral cortex', 'endothelial cells')]),
('testis', [
    ('testis', 'cells in seminiferous ducts')]),
('ovary', [
    ('ovary', 'follicle cells')]),
('epidermis', [
    ('skin 2', 'epidermal cells'),
    ('skin 1', 'keratinocytes'),
    ('skin 1', 'melanocytes')]),
('adipose tissue', [
    ('soft tissue 1', 'adipocytes'),
    ('soft tissue 2', 'adipocytes'),
    ('breast', 'adipocytes')]), 
('liver', [
    ('liver', 'hepatocytes')]),
('glial tissue', [
    ('cerebral cortex', 'glial cells'),
    ('hippocampus', 'glial cells'),
    ('lateral ventricle', 'glial cells')]),
('myocytes', [
    ('skeletal muscle', 'myocytes'),
    ('heart muscle', 'myocytes'),
    ('smooth muscle', 'smooth muscle cells')]),
('placenta', [
    ('placenta', 'trophoblastic cells'),
    ('placenta', 'decidual cells')]),
('kidney', [
    ('kidney', 'cells in tubules'),
    ('kidney', 'cells in glomeruli')])])

which_grp = dict([(t, g) for g, ts in tiss_grp.iteritems() for t in ts])

pal24 = []
with open('/home/denes/24maxcontrast', 'r') as f:
    for l in f:
        l = [float(x.strip()) for x in l.split(',')[:3]]
        pal24.append(rgb2hex(tuple(rgb256(l))))

tisscols = []
tissgcols = {}
for i, grp in enumerate(tiss_grp.keys()):
    colnum = i
    tissgcols[grp] = pal24[colnum]
    for tis in tiss_grp[grp]:
        tisscols.append(['"%s"'%x for x in ['%s%s, %s' % \
            tuple([tis[0][0].capitalize(), tis[0][1:], tis[1]]), 
            grp[0].capitalize() + grp[1:], pal24[colnum]]])

with open('tissuecolors', 'w') as f:
    for l in tisscols:
        f.write('\t'.join(l).replace('endometrial ', '')\
            .replace('/', ' ').replace('-', ' ') + '\n')

with open('tissuegcolors', 'w') as f:
    for g in sorted(tissgcols.keys()):
        f.write('\t'.join(['"%s"'%x \
            for x in [''.join([g[0].capitalize(), g[1:]]), tissgcols[g]]]) + '\n')

# tables:
pa_table = []
for ltp, data in ltps.iteritems():
    if ltp in mes_hc:
        for tis, d in pa.iteritems():
            if measured(ltp, d):
                pa_table.append([data[1], data[0], ltp,  '%s%s, %s' % \
                    (tis[0][0].capitalize(), tis[0][1:], tis[1]), d[ltp][0]])

pa_table = sorted(pa_table, key = lambda x: (x[0], x[3]))
with open('ltp-hpa-table.tsv', 'w') as f:
    f.write('\n'.join('\t'.join(c for c in r) for r in pa_table))


pd_table = []
for tis, data in pd.iteritems():
    for u, val in data.iteritems():
        pd_table.append([ltps[u][1], ltps[u][0], u, 
            '%s%s' % (tissnm[tis][0].capitalize(), tissnm[tis][1:]), str(val)])

pd_table = sorted(pd_table, key = lambda x: (x[0], x[1]))
with open('ltp-prdb-table.tsv', 'w') as f:
    f.write('\n'.join('\t'.join(c for c in r) for r in pd_table))


