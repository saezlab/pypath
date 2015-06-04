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
    chembls = l[24].split(';')
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
    fontface = 'sans-serif',
    textcol = 'black', xlog = False, ylog = False, ylim = None, **kwargs):
    fig, ax = plt.subplots()
    sns.set(font = 'Helvetica Neue LT Std')
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

def barplot(x, y, data, fname, font_family = 'Helvetica Neue LT Std', 
    xlab = '', ylab = '', lab_angle = 90, lab_size = 9, color = '#007b7f'):
    if type(x) is list:
        x = np.array(x)
    if type(y) is list:
        y = np.array(y)
    fig, ax = plt.subplots()
    sns.set(font = font_family)
    sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
    ax = sns.barplot(x, y = y, data = data, color = color)
    sns.set_context('talk', rc={'lines.linewidth': 1.0, 'patch.linewidth': 0.0,
        'grid.linewidth': 1.0})
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(lab_size)
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation = lab_angle)
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

embl_pal = embl_palette()
embl_colors = sns.blend_palette([x[0] for x in embl_palette()][:-1], 14)