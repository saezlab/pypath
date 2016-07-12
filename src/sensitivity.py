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

from future.utils import iteritems

import re
import random
import sys

import pypath
import copy
from itertools import chain
from pypath.data_formats import best, good, ugly, transcription
from pypath import dataio
from pypath.progress import Progress

def gdsc_chembls(infile = '/home/denes/Dokumentumok/gdsc/drugs_chembl_pubchem.list'):
    result = {}
    with open(infile, 'r') as f:
        for l in f:
            l = l.split('\t')
            result[l[0]] = l[2].strip().split(';') if len(l[2]) > 0 else []
    return result

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

def apply_sensitivity(graph, mapper, condition = lambda x: x['fdr'] <= 20.0):
    pairs = {}
    sens = read_sensitivity()
    in_netw = set(graph.vs['name'])
    drna = ['DNA', 'RNA']
    graph.vs['gdsc_mut'] = [{} for _ in graph.vs]
    graph.vs['gdsc_tgt'] = [{} for _ in graph.vs]
    for s in sens:
        if condition(s):
            target_up = set([it for sl in \
                [mapper.map_name(g, 'genesymbol', 'uniprot') for g in s['target'] \
                    if g not in drna] \
                for it in sl]) & in_netw
            mutated_up = set([it for sl in \
                [mapper.map_name(g, 'genesymbol', 'uniprot') for g in s['mutated']] \
                for it in sl]) & in_netw
            if len(mutated_up) > 0 and len(target_up) > 0:
                for tgt in target_up:
                    v = name2vid(graph, tgt)
                    graph.vs[v]['gdsc_tgt'][s['id']] = (s['deltapIC50'], s['drug'], 
                        s['pvalue'], mutated_up, s['chembls'])
                for mut in mutated_up:
                    v = name2vid(graph, mut)
                    graph.vs[v]['gdsc_mut'][s['id']] = (s['deltapIC50'], s['drug'], 
                        s['pvalue'], target_up, s['chembls'])
                for mut in mutated_up:
                    for tgt in target_up:
                        if s['id'] not in pairs:
                            pairs[s['id']] = []
                        sys.stdout.write('.')
                        sys.stdout.flush()
                        pairs[s['id']].append((mut, tgt, s['drug'], 
                            s['deltapIC50'], s['chembls']))
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
    textcol = 'black'):
    fig, ax = plt.subplots()
    sns.set(font = 'Helvetica Neue LT Std')
    # ax = sns.regplot(x, y, fit_reg = False, 
    #    color = embl_colors, saturation = 0.66)
    for i, xi in enumerate(x):
        ax.scatter(xi, y[i], color = colors[i])
    ax.set_xlabel(xlab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    ax.set_ylabel(ylab, weight = 'light', fontsize = 12, 
        variant = 'normal', color = textcol, stretch = 'normal')
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(9)
        tick.label.set_color(textcol)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(9)
        tick.label.set_color(textcol)
    fig.savefig(fname)


mysql_gelati = (None,'mapping_gelati')
mysql_chembl = (None,'chembl_ebi')

# 888888888888888888888888888 #

net = pypath.Pypath(9606, mysql=mysql_gelati, name="demo")

net.init_network(pfile = 'cache/plus_phospho.pickle')
#net.init_network()
#net.load_resources(lst={'mimp': good['mimp']})
#net.load_resources(lst={'pnetworks': good['pnetworks']})
#net.load_resources(lst={'psite_noref': good['psite_noref']})
#net.save_network(pfile = 'cache/plus_phospho.pickle')

pairs = apply_sensitivity(net.graph, net.mapper)

len([p for p in pairs.values()[0] if len(p[4]) == 0])

net.get_directed(conv_edges = True, mutual = True)
net.graph = net.dgraph

# looking up all paths between targets 
# and mutated proteins up to length 3
paths = {}
prg = Progress(len(pairs), 'Looking up paths', 1)
for ps in pairs.values():
    prg.step()
    for p in ps:
        paths[p] = find_all_paths2(net.graph, p[1], p[0], maxlen = 3, silent = True)

prg.terminate()

# here, paths between random pairs of targets and mutated proteins 
# is the null model:

pairs_list = set([(it[0], it[1]) for sl in pairs.values() for it in sl])
all_targets = set([i[1] for i in pairs_list])
all_mutated = set([i[0] for i in pairs_list])
rpaths = {}
prg = Progress(len(all_targets) * len(all_mutated), 'Looking up paths', 1)
for i in xrange(len(pairs_list)):
    prg.step()
    while True:
        rtgt = all_targets[int(random.random() * len(all_targets))]
        rmut = all_mutated[int(random.random() * len(all_mutated))]
        rpair = (rmut, rtgt)
        if rpair not in pairs_list  and rpair not in rpaths:
            rpaths[rpair] = find_all_paths2(net.graph, rtgt, rmut, 
                maxlen = 3, silent = True)
            break

prg.terminate()

# ************* #
# shortest paths:

sp = net.graph.shortest_paths()

resspaths, resranspaths = shortest_paths(net.graph, pairs)
resrevspaths, resrevranspaths = shortest_paths(net.graph, pairs, rev = True)
senspaths, senranspaths = shortest_paths(net.graph, pairs, neg = True)
senrevspaths, senrevranspaths = shortest_paths(net.graph, pairs, neg = True, rev = True)

resspfreq = shortest_path_lengths(resspaths, resranspaths)
resrevspfreq = shortest_path_lengths(resrevspaths, resrevranspaths)
senspfreq = shortest_path_lengths(senspaths, senranspaths)
senrevspfreq = shortest_path_lengths(senrevspaths, senrevranspaths)

splenfreq = shortest_path_lengths(spaths, rspaths)
revsplenfreq = shortest_path_lengths(revspaths, revrspaths)

with open('results/shortest_paths_sensitivity.tab', 'w') as f:
    for i in resspfreq:
        f.write('%s\t%s\t%s' % ('resistance', 'tgt-mut', 'observed') + \
            '\t%u\t%u\n' % i[0:2])
        f.write('%s\t%s\t%s' % ('resistance', 'tgt-mut', 'random') + \
            '\t%u\t%u\n' % (i[0], i[2]))
    for i in resrevspfreq:
        f.write('%s\t%s\t%s' % ('resistance', 'mut-tgt', 'observed') + \
            '\t%u\t%u\n' % i[0:2])
        f.write('%s\t%s\t%s' % ('resistance', 'mut-tgt', 'random') + \
            '\t%u\t%u\n' % (i[0], i[2]))
    for i in senspfreq:
        f.write('%s\t%s\t%s' % ('sensitivity', 'tgt-mut', 'observed') + \
            '\t%u\t%u\n' % i[0:2])
        f.write('%s\t%s\t%s' % ('sensitivity', 'tgt-mut', 'random') + \
            '\t%u\t%u\n' % (i[0], i[2]))
    for i in senrevspfreq:
        f.write('%s\t%s\t%s' % ('sensitivity', 'mut-tgt', 'observed') + \
            '\t%u\t%u\n' % i[0:2])
        f.write('%s\t%s\t%s' % ('sensitivity', 'mut-tgt', 'random') + \
            '\t%u\t%u\n' % (i[0], i[2]))


# maxflow between drug targets and mutated genes:

sp = net.graph.shortest_paths()

resmxflws, resranmxflws = shortest_paths(net.graph, pairs, maxflow = True)
resrevmxflws, resrevranmxflws = shortest_paths(net.graph, pairs, rev = True, maxflow = True)
senmxflws, senranmxflws = shortest_paths(net.graph, pairs, neg = True, maxflow = True)
senrevmxflws, senrevranmxflws = shortest_paths(net.graph, pairs, neg = True, 
    rev = True, maxflow = True)

resmffreq = shortest_path_lengths(resmxflws, resranmxflws)
resrevmffreq = shortest_path_lengths(resrevmxflws, resrevranmxflws)
senmffreq = shortest_path_lengths(senmxflws, senranmxflws)
senrevmffreq = shortest_path_lengths(senrevmxflws, senrevranmxflws)

mflenfreq = shortest_path_lengths(mxflws, rmxflws)
revmflenfreq = shortest_path_lengths(revmxflws, revrmxflws)

scatterplot([[i[3] for i in sorted(resmxflws)],
    [i[3] for i in sorted(senmxflws)]],
    [[np.log10(resmxflws[i]) for i in sorted(resmxflws)],
     [np.log10(senmxflws[i]) for i in sorted(senmxflws)]],
    'Effect size (delta IC50)', 'Maximum flow (target --> mutated)',
    'effect-maxflow.pdf', colors = [i[0] for i in embl_pal[1:]], 
    fontface = 'Helvetica Neue LT Std',
    textcol = embl_pal[-1][2])

scatterplot([[i[3] for i in sorted(resrevmxflws)],
    [i[3] for i in sorted(senrevmxflws)]],
    [[np.log10(resrevmxflws[i]) for i in sorted(resrevmxflws)],
     [np.log10(senrevmxflws[i]) for i in sorted(senrevmxflws)]],
    'Effect size (delta IC50)', 'Maximum flow (mutated --> target)',
    'effect-maxflow-rev.pdf', colors = [i[0] for i in embl_pal[1:]], 
    fontface = 'Helvetica Neue LT Std',
    textcol = embl_pal[-1][2])

scatterplot([[i[3] for i in sorted(resmxflws)],
    [i[3] for i in sorted(senmxflws)],
    [i[3] for i in sorted(resrevmxflws)],
    [i[3] for i in sorted(senrevmxflws)]],
    [[resmxflws[i] for i in sorted(resmxflws)],
     [senmxflws[i] for i in sorted(senmxflws)],
     [resrevmxflws[i] for i in sorted(resrevmxflws)],
     [senrevmxflws[i] for i in sorted(senrevmxflws)]],
    'Effect size (delta IC50)', 'Maximum flow',
    'effect-maxflow-alll.pdf', colors = [i[0] for i in embl_pal], 
    fontface = 'Helvetica Neue LT Std',
    textcol = embl_pal[-1][2])

stats.ttest_ind(np.array([resmxflws[i] for i in sorted(resmxflws)]),
    np.array([senmxflws[i] for i in sorted(senmxflws)]))

stats.ttest_ind(np.array([resrevmxflws[i] for i in sorted(resrevmxflws)]),
    np.array([senrevmxflws[i] for i in sorted(senrevmxflws)]))

stats.ttest_ind(np.array([resrevranmxflws[i] for i in sorted(resrevranmxflws)]),
    np.array([resrevmxflws[i] for i in sorted(resrevmxflws)]))

# maxflow is significantly greater than in random case
# between mutated genes and drug targets:
stats.ttest_ind(np.array([senrevranmxflws[i] for i in sorted(senrevranmxflws)]),
    np.array([senrevmxflws[i] for i in sorted(senrevmxflws)]))

boxplot(data = [
        np.array([np.log10(resmxflws[i]) for i in sorted(resmxflws)]),
        np.array([np.log10(resranmxflws[i]) for i in sorted(resranmxflws)]),
        np.array([np.log10(resrevmxflws[i]) for i in sorted(resrevmxflws)]),
        np.array([np.log10(resrevranmxflws[i]) for i in sorted(resrevranmxflws)]),
        np.array([np.log10(senmxflws[i]) for i in sorted(senmxflws)]),
        np.array([np.log10(senranmxflws[i]) for i in sorted(senranmxflws)]),
        np.array([np.log10(senrevmxflws[i]) for i in sorted(senrevmxflws)]),
        np.array([np.log10(senrevranmxflws[i]) for i in sorted(senrevranmxflws)]),
    ], 
    labels = ['Resistance', '', '', '', 'Sensitivity', '', '', ''], 
    xlab = '',
    ylab = 'Maximum flow values (log)',
    fname = 'sensitivity-maxflow.pdf', fontface = 'Helvetica Neue LT Std',
    textcol = embl_pal[-1][2])

boxplot(data = [
        np.array([resspaths[i] for i in sorted(resspaths)]),
        np.array([resranspaths[i] for i in sorted(resranspaths)]),
        np.array([resrevspaths[i] for i in sorted(resrevspaths)]),
        np.array([resrevranspaths[i] for i in sorted(resrevranspaths)]),
        np.array([senspaths[i] for i in sorted(senspaths)]),
        np.array([senranspaths[i] for i in sorted(senranspaths)]),
        np.array([senrevspaths[i] for i in sorted(senrevspaths)]),
        np.array([senrevranspaths[i] for i in sorted(senrevranspaths)]),
    ], 
    labels = ['Resistance', '', '', '', 'Sensitivity', '', '', ''], 
    xlab = '',
    ylab = 'Shortest path lengths',
    fname = 'sensitivity-spaths.pdf', fontface = 'Helvetica Neue LT Std',
    textcol = embl_pal[-1][2])


with open('results/shortest_paths_sensitivity.tab', 'w') as f:
    for i in resspfreq:
        f.write('%s\t%s\t%s' % ('resistance', 'tgt-mut', 'observed') + \
            '\t%u\t%u\n' % i[0:2])
        f.write('%s\t%s\t%s' % ('resistance', 'tgt-mut', 'random') + \
            '\t%u\t%u\n' % (i[0], i[2]))
    for i in resrevspfreq:
        f.write('%s\t%s\t%s' % ('resistance', 'mut-tgt', 'observed') + \
            '\t%u\t%u\n' % i[0:2])
        f.write('%s\t%s\t%s' % ('resistance', 'mut-tgt', 'random') + \
            '\t%u\t%u\n' % (i[0], i[2]))
    for i in senspfreq:
        f.write('%s\t%s\t%s' % ('sensitivity', 'tgt-mut', 'observed') + \
            '\t%u\t%u\n' % i[0:2])
        f.write('%s\t%s\t%s' % ('sensitivity', 'tgt-mut', 'random') + \
            '\t%u\t%u\n' % (i[0], i[2]))
    for i in senrevspfreq:
        f.write('%s\t%s\t%s' % ('sensitivity', 'mut-tgt', 'observed') + \
            '\t%u\t%u\n' % i[0:2])
        f.write('%s\t%s\t%s' % ('sensitivity', 'mut-tgt', 'random') + \
            '\t%u\t%u\n' % (i[0], i[2]))
