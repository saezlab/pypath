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
    # Some functions to use for processing the ANOVA results
    # from GDSC drug sensitivity screening.
#

import re
import random
import sys

import math
import numpy as np
from numpy.random import randn
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as hc
import hcluster as hc2
import matplotlib.gridspec as gridspec

from bioigraph.progress import Progress
from bioigraph.common import uniqList

def gdsc_chembls(infile = '/home/denes/Dokumentumok/gdsc/drugs_all_unif_manual_mapped2.tsv'):
    result = {}
    with open(infile, 'r') as f:
        for l in f:
            l = l.split('\t')
            result[l[0]] = select_chembls(l)
    return result

def select_chembls(l):
    chembls = l[24].replace('"', '').split(';')
    if len(chembls) > 9:
        chembls = l[25].split(';') + l[26].split(';')
    chembls = uniqList(chembls)
    if '' in chembls:
        chembls.remove('')
    return chembls

def read_sensitivity(infile = '/home/denes/Dokumentumok/gdsc/sensitivity/ANOVA_results.tsv'):
    regs = re.compile(r'[A-Z0-9]*')
    gsec = ['gamma-secretase', 'g-secretase']
    gsec_sub = ['PSEN1', 'NCSTN', 'APH1B', 'APH1A', 'PSENEN']
    chembls = gdsc_chembls()
    result = []
    with open(infile, 'r') as f:
        r = f.readline()
        for r in f:
            r = [i.strip() for i in r.split('\t')]
            result.append({
                'id': r[0],
                'drugid': r[2],
                'mutatedStr': r[1],
                'mutated': [i for i in regs.findall(r[1]) \
                    if len(i) > 0 and not i .startswith('PANCAN')],
                'targetStr': r[4],
                'target': [i for i in regs.findall(r[4]) if len(i) > 0] + \
                    [] if not sum([gs in r[4] for gs in gsec]) else gsec_sub,
                'drug': r[3].strip(),
                'pMaxC': float(r[7]),
                'pIC50pos': float(r[9]),
                'pIC50neg': float(r[10]),
                'sdIC50pos': float(r[12]),
                'sdIC50neg': float(r[13]),
                'effectpIC50': float(r[14]),
                'deltapIC50': float(r[11]),
                'fdr': float(r[21]),
                'pvalue': float(r[17]),
                'chembls': chembls[r[2]]
            })
    return result

def apply_sensitivity(graph, mapper, nominal_tgts, condition = lambda x: x['fdr'] <= 20.0):
    pairs = {}
    sens = read_sensitivity()
    in_netw = set(graph.vs['name'])
    drna = ['DNA', 'RNA']
    graph.vs['gdsc_mut'] = [{} for _ in graph.vs]
    graph.vs['gdsc_tgt'] = [{} for _ in graph.vs]
    prg = Progress(len(sens), 'Applying sensitivity data on network', 1)
    for s in sens:
        prg.step()
        if condition(s):
            target_up = set(nominal_tgts[s['drugid']]) & in_netw
            mutated_up = set([it for sl in \
                [mapper.map_name(g, 'genesymbol', 'uniprot') for g in s['mutated']] \
                for it in sl]) & in_netw
            if len(mutated_up) > 0 and len(target_up) > 0:
                for tgt in target_up:
                    v = name2vid(graph, tgt)
                    graph.vs[v]['gdsc_tgt'][s['id']] = (s['deltapIC50'], s['drug'], 
                        s['pvalue'], mutated_up, s['chembls'], s['drugid'])
                for mut in mutated_up:
                    v = name2vid(graph, mut)
                    graph.vs[v]['gdsc_mut'][s['id']] = (s['deltapIC50'], s['drug'], 
                        s['pvalue'], target_up, s['chembls'], s['drugid'])
                for mut in mutated_up:
                    for tgt in target_up:
                        if s['id'] not in pairs:
                            pairs[s['id']] = []
                        pairs[s['id']].append((mut, tgt, s['drug'], 
                            s['deltapIC50'], s['chembls']))
    prg.terminate()
    return pairs

def name2vid(graph, name):
    return name if type(name) is int else graph.vs.find(name=name).index

def find_all_paths2(graph, start, end, mode = 'OUT', maxlen = 2, 
        psize = 100, silent = False):
    def one_step(paths, adjlist):
        # extends all paths by one step using all neighbors in adjacency list
        return [i for ii in [[s + [a] for a in adjlist[s[-1]] if a not in s] \
            for s in paths] for i in ii]
    def parts(paths, targets, adjlist, maxlen = 2, psize = 100, depth = 1):
        complete_paths = [p for p in paths if p[-1] in targets]
        if len(paths) > 0 and len(paths[0]) <= maxlen:
            for i in xrange(0, len(paths), psize):
                new_paths = one_step(paths[i:i+psize], adjlist)
                complete_paths += parts(new_paths, targets, adjlist, maxlen, psize, 
                    depth = depth + 1)
                if not silent:
                    sys.stdout.write("\r"+" "*90)
                    sys.stdout.write('\r\tDepth: %u :: Paths found: %u' % \
                        (depth, len(complete_paths)))
                    sys.stdout.flush()
        return complete_paths
    all_paths = []
    start = [name2vid(graph, n) for n in (start if type(start) is list else [start])]
    end = [name2vid(graph, n) for n in (end if type(end) is list else [end])]
    send = set(end)
    if not silent:
        sys.stdout.write('\n')
    adjlist = [set(graph.neighbors(node, mode = mode)) \
        for node in xrange(graph.vcount())]
    paths = [[s] for s in start]
    all_paths = parts(paths, end, adjlist, maxlen, psize)
    if not silent:
        sys.stdout.write('\n')
    return all_paths

def shortest_paths(graph, pairs, neg = False, rev = False, maxflow = False):
    pairs_list = set([(it[0], it[1]) for sl in pairs.values() for it in sl])
    all_targets = set([i[1] for i in pairs_list])
    all_mutated = set([i[0] for i in pairs_list])
    spaths = {}
    prg = Progress(len(pairs), 'Looking up paths', 1)
    for ps in pairs.values():
        prg.step()
        for p in ps:
            if (neg and p[3] < 0.0) or (not neg and p[3] > 0.0):
                if maxflow:
                    spaths[p] = graph.maxflow_value(
                        name2vid(graph, p[0] if rev else p[1]), 
                        name2vid(graph, p[1] if rev else p[0]))
                else:
                    spaths[p] = graph.shortest_paths(
                        [name2vid(graph, p[0] if rev else p[1])], 
                        [name2vid(graph, p[1] if rev else p[0])])
    prg.terminate()
    rspaths = {}
    prg = Progress(len(all_targets) * len(all_mutated), 'Looking up paths', 1)
    for tgt in (all_mutated if rev else all_targets):
        for mut in (all_targets if rev else all_mutated):
            prg.step()
            p = (mut, tgt)
            if p not in pairs_list:
                if maxflow:
                    rspaths[p] = graph.maxflow_value(
                        name2vid(graph, mut if rev else tgt), 
                        name2vid(graph, tgt if rev else mut))
                else:
                    rspaths[p] = graph.shortest_paths(
                        [name2vid(graph, mut if rev else tgt)], 
                        [name2vid(graph, tgt if rev else mut)])
    prg.terminate()
    return spaths, rspaths

def shortest_path_lengths(spaths, rspaths):
    splens = [i[0][0] for i in spaths.values() if i[0][0] != float('inf')]
    rsplens = [i[0][0] for i in rspaths.values() if i[0][0] != float('inf')]
    return [(i, splens.count(i), rsplens.count(i) * len(splens) / len(rsplens)) \
        for i in xrange(max(splens + rsplens))]

def scatterplot(x, y, xlab, ylab, fname, colors, 
    fontface = 'Helvetica Neue LT Std', legend = None, 
    textcol = 'black', xlog = False, ylog = False, ylim = None, **kwargs):
    fig, ax = plt.subplots()
    sns.set(font = fontface)
    # ax = sns.regplot(x, y, fit_reg = False, 
    #    color = embl_colors, saturation = 0.66)
    for i, xi in enumerate(x):
        xx = [0.0 if n <= 0.0 else math.log10(n) for n in xi] if xlog else xi
        yy = [0.0 if n <= 0.0 else math.log10(n) for n in y[i]] if ylog else y[i]
        ax.scatter(xx, yy, color = colors[i], **kwargs)
    ax.set_xlabel(xlab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    ax.set_ylabel(ylab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    if ylim is not None:
        plt.ylim(ylim)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(9)
        tick.label.set_color(textcol)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(9)
        tick.label.set_color(textcol)
    if type(legend) is dict:
        legend_patches = [mpatches.Patch(color=col, label=lab) \
            for lab, col in legend.iteritems()]
        ax.legend(handles = legend_patches)
        xlm = ax.get_xlim()
        ax.set_xlim((xlm[0], xlm[1]*1.2))
    fig.tight_layout()
    fig.savefig(fname)

def embl_palette(textf = '/home/denes/embl_colors'):
    cols = []
    with open(textf, 'r') as f:
        series = []
        for i, l in enumerate(f):
            l = [x.strip() for x in l.split(',')]
            series.append(rgb2hex(tuple([256 * float(x) for x in l[0:3]])))
            if len(series) == 7:
                cols.append(series)
                series = []
    return cols

def boxplot(data, labels, xlab, ylab, fname, fontfamily = 'Helvetica Neue LT Std',
    textcol = 'black', violin = False):
    fig, ax = plt.subplots()
    sns.set(font = fontfamily)
    if violin:
        ax = sns.violinplot(data, names = labels, 
            color = embl_colors, linewidth = 0.1, saturation = 0.66)
    else:
        ax = sns.boxplot(data, names = labels, 
            color = embl_colors, linewidth = 0.1, saturation = 0.66)
    ax.set_xlabel(xlab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    ax.set_ylabel(ylab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
        tick.label.set_color(textcol)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(11)
        tick.label.set_color(textcol)
    fig.savefig(fname)

def stacked_barplot(x, y, data, fname, names, font_family = 'Helvetica Neue LT Std', 
    xlab = '', ylab = '', lab_angle = 90, lab_size = (18, 21), axis_lab_size = 36, 
    legend = True, 
    colors = ['#7AA0A1', '#C6909C', '#92C1D6', '#C5B26E', '#da0025'], 
    order = False, desc = True):
    if type(x) is list or type(x) is tuple:
        x = np.array(x)
    for i, yi in enumerate(y):
        if type(yi) is list or type(yi) is tuple:
            y[i] = np.array(yi)
    total = np.array([sum([y[j][i] for j in xrange(len(y))]) for i in xrange(len(y[0]))])
    if order == 'x':
        ordr = np.array([x[i] for i in x.argsort()])
    elif order == 'y':
        ordr = np.array([x[i] for i in total.argsort()])
    elif type(order) is int:
        ordr = np.array([x[i] for i in y[order].argsort()])
    elif len(set(order) & set(x)) == len(x):
        ordr = order
    else:
        ordr = x
    if desc:
        ordr = ordr[::-1]
    fig, ax = plt.subplots()
    sns.set(font = font_family)
    sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0, 'axes.labelsize': axis_lab_size})
    for j in xrange(len(y), 0, -1):
        this_level = np.array([sum([y[jj][i] for jj in xrange(j)]) \
            for i in xrange(len(y[0]))])
        ax = sns.barplot(x, y = this_level, data = data, 
            color = colors[j-1], order = ordr)
    #ax = sns.barplot(x, y = y, data = data, color = color, x_order = ordr)
    #plt.bar(range(len(ordr)), [y[i] for i in ordr], align = 'center')
    #plt.xticks(list(ordr), [x[i] for i in ordr])
    sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0, 'axes.labelsize': axis_lab_size})
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(lab_size[0])
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(lab_size[1])
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = lab_angle)
    if legend:
        lhandles = [mpl.patches.Patch(color = colors[i], label = names[i]) for i in xrange(len(y))]
    ax.legend(handles = lhandles)
    fig.tight_layout()
    fig.savefig(fname)
    plt.close(fig)

def barplot(x, y, data, fname, font_family = 'Helvetica Neue LT Std', 
    xlab = '', ylab = '', lab_angle = 90, lab_size = 9, color = '#007b7f', 
    order = False, desc = True, legend = None, fin = True, 
    y_break = None, **kwargs):
    '''
    y_break : tuple
    If not None, the y-axis will have a break. 2 floats in the tuple, < 1.0, 
    mean the lower and upper proportion of the plot shown. The part between
    them will be hidden. E.g. y_break = (0.3, 0.1) shows the lower 30% and 
    upper 10%, but 60% in the middle will be cut out.
    '''
    ax2 = None
    if type(x) is list or type(x) is tuple:
        x = np.array(x)
    if type(y) is list or type(y) is tuple:
        y = np.array(y)
    if order == 'x':
        ordr = np.array([x[i] for i in x.argsort()])
    elif order == 'y':
        ordr = np.array([x[i] for i in y.argsort()])
    else:
        ordr = x
    if desc:
        ordr = ordr[::-1]
    if type(color) is list and len(color) == len(ordr):
        xl = list(x)
        color = [color[xl.index(xi)] for xi in ordr]
    if y_break:
        gs = gridspec.GridSpec(2, 1, height_ratios = [y_break[1] / sum(y_break), y_break[0] / sum(y_break)])
        fig = plt.figure()
        ax2 = fig.add_subplot(gs[0])
        ax = fig.add_subplot(gs[1])
    else:
        fig, ax = plt.subplots()
    sns.set(font = font_family)
    sns.set_context('poster', rc = {'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
    ax = sns.barplot(x, y = y, data = data, color = color, order = ordr, ax = ax, **kwargs)
    if y_break:
        ax2 = sns.barplot(x, y = y, data = data, color = color, order = ordr, ax = ax2, **kwargs)
        ax2.yaxis.set_major_locator(MaxNLocator(nbins = int(9/sum(y_break) + 1), steps = [1, 2, 5, 10]))
        originalYticks = ax2.get_yticks()
        ymin, ymax = ax.get_ylim()
        ymax = min(ytick for ytick in ax.get_yticks() if ytick > max(y))
        ax.set_ylim((ymin, ymax * y_break[0]))
        ax2.set_ylim((ymax - ymax * y_break[1], ymax))
        axmin, axmax = ax.get_ylim()
        ax2min, ax2max = ax2.get_ylim()
        plt.subplots_adjust(hspace = 0.08)
        ax2.spines['bottom'].set_visible(False)
        plt.setp(ax2.xaxis.get_majorticklabels(), visible = False)
        ax.spines['top'].set_visible(False)
        ax.set_yticks([yt for yt in originalYticks if yt >= axmin and yt <= axmax])
        ax2.set_yticks([yt for yt in originalYticks if yt >= ax2min and yt <= ax2max])
    sns.set_context('poster', rc = {'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(lab_size)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = lab_angle)
    if type(legend) is dict:
        legend_patches = [mpatches.Patch(color = col, label = lab) \
            for lab, col in legend.iteritems()]
        ax.legend(handles = legend_patches)
    if fin:
        finish(fig, fname)
    else:
        return fig, ax, ax2

def finish(fig, fname):
    fig.tight_layout()
    fig.savefig(fname)
    plt.close(fig)

def in_complex(self, csources = ['corum']):
    self.graph.es['in_complex'] = \
    [sum([len(set(self.graph.vs[e.source]['complexes'][cs].keys()) & \
    set(self.graph.vs[e.target]['complexes'][cs].keys())) for cs in csources]) > 0 \
    for e in self.graph.es]

def complexes_in_network(g, csource = 'corum'):
    cdict = {}
    allv = set(g.vs['name'])
    for v in g.vs:
        for c, cdata in v['complexes'][csource].iteritems():
            if c not in cdict:
                cdict[c] = set(cdata['all_members'])
    return [c for c, memb in cdict.iteritems() if len(memb - allv) == 0]

def rgb2hex(rgb):
    return '#%02x%02x%02x' % rgb

def hex2rgb(self, rgbhex):
    rgbhex = rgbhex.lstrip('#')
    lv = len(rgbhex)
    return tuple(int(rgbhex[i:i + 2], 16) for i in range(0, lv, 2))

def rgb1(self, rgb256):
    return rgb256 if not any([i > 1 for i in rgb256]) \
        else tuple([x / float(255) for x in rgb256])

embl_pal = embl_palette()
embl_colors = sns.blend_palette([x[0] for x in embl_palette()][:-1], 14)