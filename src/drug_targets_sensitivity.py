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
import copy
import igraph
import louvain

# stats and plotting modules #

import math
import ranking
import numpy as np
from numpy.random import randn
import pandas as pd
from scipy import stats
from scipy.misc import comb
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as hc
import hcluster as hc2
import matplotlib.patches as mpatches

# from bioigraph #

import bioigraph
from bioigraph import chembl
from bioigraph.common import *
import _sensitivity as sens
from bioigraph import progress

# 0 ## Constants

drugfile = '/home/denes/Dokumentumok/gdsc/drugs_all_unif_manual_mapped2.tsv'
chembl_mysql = (None, 'chembl_ebi')
mysql_gelati = (None,'mapping_gelati')

# 1 ## Reading GDSC drug info

drugid_chembl = {}
chembl_drugid = {}
drug_data = {}
nominal_tgts = {}
with open(drugfile, 'r') as f:
    l = f.readline()
    for l in f:
        l = [x.strip().replace('"', '') for x in l.split('\t')]
        drugid = l[0]
        chembls = sens.select_chembls(l)
        drugid_chembl[drugid] = chembls
        for ch in chembls:
            if ch not in chembl_drugid:
                chembl_drugid[ch] = []
            chembl_drugid[ch].append(drugid)
        drug_data[drugid] = [l[6], l[7]]
        targets = l[27].split(';')
        if '' in targets:
            targets.remove('')
        nominal_tgts[drugid] = targets

# 2 ## Initializing signaling network

net = bioigraph.BioGraph(9606, mysql=mysql_gelati, name="demo")

#net.init_network(pfile = 'cache/plus_phospho.pickle')
#net.load_resources(lst={'mimp': good['mimp']})
#net.load_resources(lst={'pnetworks': good['pnetworks']})
#net.load_resources(lst={'psite_noref': good['psite_noref']})
#net.save_network(pfile = 'cache/plus_phospho.pickle')
#net.get_directed(conv_edges = True, mutual = True)
#net.graph = net.dgraph
#net.save_network(pfile = 'cache/plus_phospho_directed.pickle')
net.init_network(pfile = 'cache/plus_phospho_directed.pickle')
net.genesymbol_labels()
net.set_receptors()
net.set_tfs()

# 3 ## Reading sensitivity results and mapping them onto the network

pairs = sens.apply_sensitivity(net.graph, net.mapper, nominal_tgts)
pairs_list = set([(it[0], it[1]) for sl in pairs.values() for it in sl])
all_targets = set([i[1] for i in pairs_list])
all_mutated = set([i[0] for i in pairs_list])
len([p for p in pairs.values()[0] if len(p[4]) == 0])
# vertex attributes created by apply_sensitivity():
# 'gdsc_tgt': {a_id: (delta_pIC50, drug_name, anova_p_value, 
#    [mutated_uniprots], chembls, drugid)}
# 'gdsc_mut': {a_id: (delta_pIC50, drug_name, anova_p_value, 
#    [target_uniprots], chembls, drugid)}

# 3 ## Querying ChEMBL

c = chembl.Chembl(chembl_mysql = chembl_mysql)
all_chembls = uniqList([ch for chs in drugid_chembl.values() for ch in chs])
c.compounds_targets(all_chembls, id_type = 'chembl', assay_types = ['B'], 
    pchembl = True)
c.compounds_by_target()
c.targets_by_compound()

# 4 ## Compound-target interactions on the network at pChEMBL >= 5.0

net.graph.vs['chembl_cpds'] = [[] if v['name'] not in c.compounds \
    else [cpd for cpd in c.compounds[v['name']] \
        if max([float(pc) for pc in cpd['pchembl']]) >= 5.0]
    for v in net.graph.vs]

# 5 ## Comparing nominal targets vs. ChEMBL targets

nominal_chembl = \
dict([(ch, 
([tgt for tgts in [nominal_tgts[drugid]] \
    for drugid in chembl_drugid[ch] \
    for tgt in tgts]
,
[] if ch not in c.targets \
    else [tgt['uniprot'] for tgt in c.targets[ch] \
        if max([float(pc) for pc in tgt['pchembl']]) >= 5.0]
)) for ch in all_chembls])

sens.scatterplot([[len(nominal_chembl[ch][1]) \
                for ch in all_chembls]],
                [[len(nominal_chembl[ch][0]) \
                for ch in all_chembls]],
            'Number of targets in ChEMBL binding assays\nhaving pChEMBL >= 5.0',
            'Number of nominal targets', 
    'nominal_vs_chembl.pdf', colors = [i[0] for i in sens.embl_pal[1:]], 
    fontface = 'Helvetica Neue LT Std',
    textcol = sens.embl_pal[-1][2])

## with chembl ids:
# having nominal target(s):
len([ch for ch, tg in nominal_chembl.iteritems() if len(tg[0]) > 0])
# having target(s) from ChEMBL:
len([ch for ch, tg in nominal_chembl.iteritems() if len(tg[1]) > 0])
# having target(s) from both:
len([ch for ch, tg in nominal_chembl.iteritems() if len(tg[0]) > 0 and len(tg[1]) > 0])
# having no targets:
len([ch for ch, tg in nominal_chembl.iteritems() if len(tg[0]) == 0 and len(tg[1]) == 0])

# this is a barplot, but also made a venn diagram in inkscape, that's better
_sens.barplot(x = ['Having nominal targets', 'Having targets from ChEMBL', 'Having targets from both', 'Having no targets'], y = [
    len([ch for ch, tg in nominal_chembl.iteritems() if len(tg[0]) > 0]),
    len([ch for ch, tg in nominal_chembl.iteritems() if len(tg[1]) > 0]),
    len([ch for ch, tg in nominal_chembl.iteritems() if len(tg[0]) > 0 and len(tg[1]) > 0]),
    len([ch for ch, tg in nominal_chembl.iteritems() if len(tg[0]) == 0 and len(tg[1]) == 0])
    ], data = None, fname = 'gdsc-comp-targets.pdf', lab_size = 11, 
    ylab = 'Number of GDSC compounds (unique ChEMBL IDs)', 
    xlab = 'Target data')

## with drug ids:
# having nominal target(s):
len(uniqList(flatList([chembl_drugid[ch] for ch, tg in nominal_chembl.iteritems() if len(tg[0]) > 0])))
# having target(s) from ChEMBL:
len(uniqList(flatList([chembl_drugid[ch] for ch, tg in nominal_chembl.iteritems() if len(tg[1]) > 0])))
# having target(s) from both:
len(uniqList(flatList([chembl_drugid[ch] for ch, tg in nominal_chembl.iteritems() if len(tg[0]) > 0 and len(tg[1]) > 0])))
# having no targets:
len(uniqList(flatList([chembl_drugid[ch] for ch, tg in nominal_chembl.iteritems() if len(tg[0]) == 0 and len(tg[1]) == 0])))

# 6 ## 
gi = net.get_giant()
spc = gi.community_spinglass(spins = 500)
ug = net.graph.as_undirected(combine_edges = 'ignore')
fgc = ug.community_fastgreedy()
fgc = fgc.as_clustering()
wtc = net.graph.community_walktrap()
wtc = wtc.as_clustering()
net.graph.vs['walktrap'] = wtc.membership
imc = net.graph.community_infomap()
net.graph.vs['infomap'] = imc.membership

def community_onto_graph(graph, community, name):
    if isinstance(community, igraph.clustering.VertexDendrogram):
        community = community.as_clustering()
    graph.vs[name] = community.membership
    graph.es[name] = [(community.membership[e.source], \
        community.membership[e.target]) for e in graph.es]

def copy_vertex_attr(graph_source, graph_target, attr, default = None, idattr = 'name',
    idattr_target = None):
    idattr_target = idattr if idattr_target is None else idattr_target
    graph_target.vs[attr] = [default for _ in graph_target.vs]
    name2vid = dict([(v[idattr_target], v.index) for v in graph_target.vs])
    for v1 in graph_source.vs:
        if v1[idattr] in name2vid:
            graph_target.vs[name2vid[v1[idattr]]][attr] = v1[attr]

def copy_edge_attr(graph_source, graph_target, attr, default = None, idattr = 'name',
    idattr_target = None):
    idattr_target = idattr if idattr_target is None else idattr_target
    graph_target.es[attr] = [default for _ in graph_target.es]
    name2vid = dict([(v[idattr_target], v.index) for v in graph_target.vs])
    for e1 in graph_source.es:
        if graph_source.vs[e1.source][idattr] in name2vid and \
            graph_source.vs[e1.target][idattr] in name2vid:
            e2 = graph_target.get_eid(name2vid[graph_source.vs[e1.source][idattr]], 
                name2vid[graph_source.vs[e1.target][idattr]], error = False)
            if e2 != -1:
                graph_target.es[e2][attr] = e1[attr]
            if not graph_source.is_directed() and graph_target.is_directed():
                e2 = graph_target.get_eid(name2vid[graph_source.vs[e1.target][idattr]], 
                    name2vid[graph_source.vs[e1.source][idattr]], error = False)
                if e2 != -1:
                    graph_target.es[e2][attr] = e1[attr]

def community_onto_original(copy, original, community, name):
    community_onto_graph(copy, community, name)
    copy_vertex_attr(copy, original, name)
    copy_edge_attr(copy, original, name)

def community_metagraph(graph, group_attr, giant = True):
    groups = sorted(list(set(graph.vs[group_attr])))
    group_members = dict([(gr, []) for gr in groups])
    for v in graph.vs:
        group_members[v[group_attr]].append(v.index)
    meta = igraph.Graph(len(groups), directed = True)
    meta.vs['name'] = groups
    adjlist = dict([((v.index,v[group_attr]), \
        [(n, net.graph.vs[n][group_attr]) \
            for n in graph.neighbors(v.index, mode = 'OUT')] \
        ) for v in graph.vs])
    edgelist = {}
    for gr1 in meta.vs:
        for gr2 in meta.vs:
            if gr1 != gr2:
                one = one.index
                two = two.index
                edgelist[(one, two)] = sum([len([m2 \
                        for m2 in adjlist[m1, one] \
                        if m2[1] == two ]) 
                    for m1 in group_members[one]])
                edgelist[(two, one)] = sum([len([m1 \
                        for m1 in adjlist[m2, two] \
                        if m1[1] == one]) 
                    for m2 in group_members[two]])
    meta.add_edges([e for e, w in edgelist.iteritems() if w != 0])
    meta.es['weight'] = [edgelist[(e.source, e.target)] for e in meta.es]
    if giant:
        comp = meta.components()
        comp_sizes = comp.sizes()
        gi_cmp_idx = comp_sizes.index(max(comp_sizes))
        to_del = [idx for idx, c in enumerate(comp.membership) if c != gi_cmp_idx]
        meta.delete_vertices(to_del)
    return meta

def group_significance(graph, group_attr):
    gr_sigf = {}
    graph_dens = graph.density()
    groups = sorted(list(set(graph.vs[group_attr])))
    prg = progress.Progress(len(groups), 'Calculating group significance', 1)
    for gr in groups:
        prg.step()
        gr_size = graph.vs[group_attr].count(gr)
        gr_ecount = len([e for e in graph.es \
            if graph.vs[e.source][group_attr] == gr \
            and graph.vs[e.target][group_attr] == gr])
        ecount_full = gr_size * (gr_size - 1) if graph.is_directed() else comb(gr_size, 2)
        pk = 0.0 if ecount_full == 0.0 else gr_ecount / float(ecount_full)
        sigf = math.e ** (-1 * ecount_full * _KL(pk, graph_dens))
        gr_sigf[gr] = 1 / sigf if sigf > 0.0 else 0.0
    prg.terminate()
    maxsig = max(gr_sigf.values())
    # gr_sigf = dict([(k, v/maxsig) for k, v in gr_sigf.iteritems()])
    return gr_sigf

def _KL(q,p):
    """ The binary Kullback-Leibler divergence. """
    KL = 0;
    if q > 0 and p > 0:
      KL += q * math.log(q/p);
    if q < 1 and p < 1:
      KL += (1-q) * math.log((1-q)/(1-p));
    return KL

def community_metagraph2(graph, group_attr, giant = True):
    groups = sorted(list(set(graph.vs[group_attr])))
    group_members = dict([(gr, []) for gr in groups])
    for v in graph.vs:
        group_members[v[group_attr]].append(v.index)
    meta = igraph.Graph(len(groups), directed = True)
    meta.vs['name'] = groups
    # group significance to vertex attribute of the meta-graph:
    groups_sigf = group_significance(graph, group_attr)
    meta.vs['sig'] = [groups_sigf[v['name']] for v in meta.vs]
    meta.vs['members'] = [[graph.vs[m]['name'] for m in group_members[v.index]] \
        for v in meta.vs]
    meta.vs['members_gname'] = [[graph.vs[m]['label'] for m in group_members[v.index]] \
        for v in meta.vs]
    meta.vs['size'] = [len(group_members[v.index]) for v in meta.vs]
    # creating adjacency list of intergroup connections:
    adjlist = [[(n, graph.vs[n][group_attr]) \
            for n in graph.neighbors(v.index, mode = 'OUT') \
                if graph.vs[n][group_attr] != v[group_attr]] \
        for v in graph.vs]
    adjlist_in = [[(n, graph.vs[n][group_attr]) \
            for n in graph.neighbors(v.index, mode = 'IN') \
                if graph.vs[n][group_attr] != v[group_attr]] \
        for v in graph.vs]
    # calculating edges of the meta-graph:
    edgelist = {}
    for gr1 in meta.vs:
        for gr2 in meta.vs:
            if gr1 != gr2:
                one = gr1.index
                two = gr2.index
                # directed weights:
                ## number of edges from group one to group two
                e12 = sum([len([t for t in adjlist[s] if t[0] in group_members[two]]) \
                    for s in group_members[one]])
                ## number of edges from group two to group one
                e21 = sum([len([t for t in adjlist[s] if t[0] in group_members[one]]) \
                    for s in group_members[two]])
                ## all edges from group one to outside
                out1 = float(sum([len(adjlist[s]) for s in group_members[one]]))
                ## all edges from group two to outside
                out2 = float(sum([len(adjlist[s]) for s in group_members[two]]))
                ## all edges from outside to group one
                in1 = float(sum([len(adjlist_in[s]) for s in group_members[one]]))
                ## all edges from outside to group two
                in2 = float(sum([len(adjlist_in[s]) for s in group_members[two]]))
                ## weight of edge in direction one > two:
                ## product of proportion of edges to group two and proportion 
                ## of edges from group one
                w12 = ((e12 / out1) if out1 > 0 else 0) * ((e12 / in2) if in2 > 0 else 0)
                ## weight of edge in direction two > one
                w21 = ((e21 / out2) if out2 > 0 else 0) * ((e21 / in1) if in1 > 0 else 0)
                edgelist[(one, two)] = w12
                edgelist[(two, one)] = w21
    meta.add_edges([e for e, w in edgelist.iteritems() if w != 0.0])
    meta.es['weight'] = [edgelist[(e.source, e.target)] for e in meta.es]
    if giant:
        comp = meta.components()
        comp_sizes = comp.sizes()
        gi_cmp_idx = comp_sizes.index(max(comp_sizes))
        to_del = [idx for idx, c in enumerate(comp.membership) if c != gi_cmp_idx]
        meta.delete_vertices(to_del)
    return meta

lsig = louvain.find_partition(net.graph, method = 'Significance')
community_onto_graph(net.graph, imc, 'infomap')
community_onto_graph(net.graph, lsig, 'lvn_sig')
community_onto_graph(net.graph, wtc, 'walktrap')
community_onto_original(gi, net.graph, spc, 'spinglass')
community_onto_original(ug, net.graph, fgc, 'fastgreedy')

im_meta = community_metagraph(net.graph, 'infomap')
im_meta = community_metagraph2(net.graph, 'infomap')

lsig_meta = community_metagraph2(net.graph, 'lvn_sig')

fig, ax = plt.subplots()
sns.set(font = 'Helvetica Neue LT Std')
ax = sns.kdeplot(np.array(lsig_meta80.es['weight']), shade = True, 
    color = '#007b7f')
fig.tight_layout()
fig.savefig('density-lsig-meta80-weights.pdf')

p90 = np.percentile(lsig_meta.es['weight'], 90)
p80 = np.percentile(lsig_meta.es['weight'], 80)
p70 = np.percentile(lsig_meta.es['weight'], 70)
lsig_meta80 = copy.deepcopy(lsig_meta)
lsig_meta80.delete_edges([e.index for e in lsig_meta80.es \
    if e['weight'] < p80])

sens.scatterplot([lsig_meta.vs['size']], [lsig_meta.vs['sig']], 'Community size', 
    'Community significance', 'size-signf-lsig.pdf', ylog = True,
    colors = [i[0] for i in sens.embl_pal[1:]])

lsig_meta80.layout = lsig_meta80.layout_fruchterman_reingold(weights = 'weight', 
        repulserad = lsig_meta80.vcount() ** 2.8, maxiter = 1000, 
        area = lsig_meta80.vcount() ** 2.3)
im_meta.vs['label'] = im_meta.vs['name']

sf = cairo.PDFSurface('lsig-metagraph-80.pdf', 600, 600)
igplot = igraph.plot(lsig_meta80, 
    target = sf,
    layout = lsig_meta80.layout,
    vertex_label = [','.join(v['members_gname'][:3]) for v in lsig_meta80.vs],
    vertex_label_size = 3,
    vertex_size = [sz/np.mean(lsig_meta80.vs['size'])*3 for sz in lsig_meta80.vs['size']],
    vertex_frame_width = 0, 
    vertex_color = '#22888888', vertex_label_color = '#88222288', 
    edge_width = lsig_meta80.es['weight'], edge_color = '#33333322', 
    edge_arrow_size = 0.1,
    edge_arrow_width = 0.1,
    edge_label_color = '#22228888')

igplot.redraw()
igplot.save()

### 

nominal_pairs = \
    name2vid_pairs(
        net.graph, 
        uniqList(
            flatList(
                    [[(p[0], p[1]) for p in pp] for pp in pairs.values()]
    )))
chembl_pairs = []
pchembl_pairs = {}
for pp in pairs.values():
    for p in pp:
        for ch in p[4]:
            if ch in c.targets:
                for tgt in c.targets[ch]:
                    pchembls = [float(pc) for pc in tgt['pchembl']]
                    if tgt['uniprot'] in net.graph.vs['name'] and max(pchembls) >= 0.5:
                        pair = (p[0], tgt['uniprot'])
                        chembl_pairs.append(pair)
                        if pair not in pchembl_pairs:
                            pchembl_pairs[pair] = {}
                        pchembl_pairs[pair][ch] = pchembls

chembl_pairs = name2vid_pairs(net.graph, uniqList(chembl_pairs))

all_pairs = uniqList(nominal_pairs + chembl_pairs)

all_targets = set([i[1] for i in all_pairs])
all_mutated = set([i[0] for i in all_pairs])
non_pairs = [(m, t) for m in all_mutated for t in all_targets \
    if m != t and (m, t) not in all_pairs]

total_pairs = [(m, t) for m in all_mutated for t in all_targets]
total_pairs += reverse_pairs(total_pairs)
total_pairs = uniqList(total_pairs)

weights = groups_topology(net.graph, total_pairs, group_affinity, 
    group_attr = 'lvn_sig', meta = lsig_meta80)

###
net.update_vname()
prg = progress.Progress(len(pairs), 'Calculating topological metric', 1)
# results to be stored in these lists:
mut_tgt = []
tgt_mut = []
ww_mut_tgt = []
ww_tgt_mut = []
ew_mut_tgt = []
ew_tgt_mut = []
diff_mut_tgt = []
# iterating all sensitivity/resistance relationships:
for aid, pp in pairs.iteritems():
    prg.step()
    chembls = pp[0][4]
    deltapic50 = pp[0][3]
    # proceed only if we have data about this compound from ChEMBL:
    if bool(len(chembls)) and bool(sum([ch in c.targets for ch in chembls])):
        targets = {}
        # iterating ChEMBLs: there might be multiple ChEMBLs, 
        # but in real it is one compound (at least we assume):
        for ch in chembls:
            # if this ChEMBL ID has targets:
            if ch in c.targets:
                # iterating its targets:
                for tgt in c.targets[ch]:
                    # only if the target in the network
                    if tgt['uniprot'] in net.nodDct:
                        # the vertex id of the target:
                        tgt_vid = net.nodDct[tgt['uniprot']]
                        # collecting pChEMBL values for the target:
                        if tgt_vid not in targets:
                            targets[tgt_vid] = []
                        targets[tgt_vid] += [float(pch) for pch in tgt['pchembl']]
        # taking the maximum of all pChEMBLs:
        targets = dict([(u, max(pchs)) for u, pchs in targets.iteritems()])
        # collect the mutated genes:
        mutated = set(uniqList([net.nodDct[p[0]] for p in pp]))
        w_mut_tgt = []
        w_tgt_mut = []
        for mut_vid in mutated:
            # weights of this relationship:
            w_mut_tgt.append(sum([rweights4[(mut_vid, tgt_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]))
            w_tgt_mut.append(sum([rweights4[(tgt_vid, mut_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]))
            diff_mut_tgt.append((sum([rweights4[(tgt_vid, mut_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]) - \
                sum([rweights4[(mut_vid, tgt_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]), deltapic50))
        ww_mut_tgt += w_mut_tgt
        ww_tgt_mut += w_tgt_mut
        ew_mut_tgt += [(w, deltapic50) for w in w_mut_tgt]
        ew_tgt_mut += [(w, deltapic50) for w in w_tgt_mut]
        # control: all other mutated:
        c_mut_tgt = []
        c_tgt_mut = []
        for mut_vid in all_mutated - mutated:
            # weights of this relationship:
            c_mut_tgt.append(sum([rweights4[(mut_vid, tgt_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]))
            c_tgt_mut.append(sum([rweights4[(tgt_vid, mut_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]))
        # get ranks:
        r_mut_tgt = [rank(x, c_mut_tgt) for x in w_mut_tgt]
        r_tgt_mut = [rank(x, c_tgt_mut) for x in w_tgt_mut]
        mut_tgt += [(x, deltapic50) for x in r_mut_tgt]
        tgt_mut += [(x, deltapic50) for x in r_tgt_mut]

prg.terminate()
###

sens.scatterplot([[x[0] for x in mut_tgt], [x[0] for x in tgt_mut]],
                [[x[1] for x in mut_tgt], [x[1] for x in tgt_mut]],
            'Rank', 'Delta pIC50', 
    'dpic50_vs_rank_rw.pdf', colors = [i[0] for i in sens.embl_pal[1:]], 
    fontface = 'Helvetica Neue LT Std',
    textcol = sens.embl_pal[-1][2])

sens.scatterplot([[x[0] for x in ew_mut_tgt], [x[0] for x in ew_tgt_mut]],
                [[x[1] for x in ew_mut_tgt], [x[1] for x in ew_tgt_mut]],
            'Weight', 'Delta pIC50', 
    'dpic50_vs_weight_rw.pdf', colors = [i[0] for i in sens.embl_pal[1:]], 
    fontface = 'Helvetica Neue LT Std',
    textcol = sens.embl_pal[-1][2])


fig, ax = plt.subplots()
sns.set(font = 'Helvetica Neue LT Std')
ax = sns.kdeplot(np.array(ww_mut_tgt), shade = True, 
    color = '#007b7f')
ax = sns.kdeplot(np.array(ww_tgt_mut), shade = True, 
    color = '#FEEA9C')
fig.tight_layout()
fig.savefig('weights-raw-rw.pdf')



###

def rank(value, lst, norm = True):
    r = ranking.Ranking(sorted(lst + [value]), reverse = True, 
        strategy = ranking.FRACTIONAL).rank(value)
    if norm:
        return r / float(len(lst))
    else:
        return r

rev_nominal_pairs = reverse_pairs(nominal_pairs)
rev_chembl_pairs = reverse_pairs(chembl_pairs)
rev_non_pairs = reverse_pairs(non_pairs)

nm_lsig = groups_topology(net.graph, nominal_pairs, group_affinity, 
    group_attr = 'lvn_sig', meta = lsig_meta80)
ch_lsig = groups_topology(net.graph, chembl_pairs, group_affinity, 
    group_attr = 'lvn_sig', meta = lsig_meta80)
no_lsig = groups_topology(net.graph, non_pairs, group_affinity, 
    group_attr = 'lvn_sig', meta = lsig_meta80)
r_nm_lsig = groups_topology(net.graph, rev_nominal_pairs, group_affinity, 
    group_attr = 'lvn_sig', meta = lsig_meta80)
r_ch_lsig = groups_topology(net.graph, rev_chembl_pairs, group_affinity, 
    group_attr = 'lvn_sig', meta = lsig_meta80)
r_no_lsig = groups_topology(net.graph, rev_non_pairs, group_affinity, 
    group_attr = 'lvn_sig', meta = lsig_meta80)

def reverse_pairs(pairs):
    return [(p[1], p[0]) for p in pairs]

def name2vid_pairs(graph, pairs):
    n2v = dict([(v['name'], v.index) for v in graph.vs])
    return [(n2v[p[0]], n2v[p[1]]) for p in pairs]

def vid2name_pairs(graph, pairs):
    v2n = dict([(v.index, v['name']) for v in graph.vs])
    return [(v2n[p[0]], v2n[p[1]]) for p in pairs]

def groups_topology(graph, pairs, topology, **kwargs):
    return dict([(p, topology(graph, p[0], p[1], **kwargs)) for p in pairs])

def group_affinity(graph, source, target, group_attr, meta):
    gr1 = graph.vs[source][group_attr]
    gr2 = graph.vs[target][group_attr]
    if gr1 in meta.vs['name'] and gr2 in meta.vs['name']:
        return _group_affinity(meta, gr1, gr2)
    else:
        return 1.0 if gr1 == gr2 else 0.0

def _group_affinity(meta, gr1, gr2):
    gr1i = meta.vs['name'].index(gr1)
    gr2i = meta.vs['name'].index(gr2)
    paths = meta.get_all_shortest_paths(gr1i, gr2i, mode = 'OUT')
    weights = []
    for path in paths:
        this_weight = 1.0
        for i in xrange(len(path) - 1):
            edge = meta.get_eid(path[i], path[i + 1])
            this_weight *= meta.es[edge]['weight']
        weights.append(this_weight)
    return 0.0 if len(paths) == 0 else np.mean(weights)

def _random_walk(adj, start, maxlen):
    visited = set([])
    current_node = start
    for step in xrange(maxlen):
        routes = adj[current_node] - visited
        if len(routes) > 0:
            current_node = random.choice(list(routes))
            visited.add(current_node)
            yield step, current_node
        else:
            break

def random_walk(adj, start, maxlen, n):
    score = 1.0 / float(n)
    result = np.zeros((len(adj), maxlen), dtype = float)
    for rep in xrange(n):
        for step, current_node in _random_walk(adj, start, maxlen):
            result[current_node, step] += score
    return result

def random_walks(graph, start, maxlen = 5, n = 1000):
    result = {}
    adj = [set([nb.index for nb in v.neighbors(mode = 'OUT')]) for v in graph.vs]
    start = start if type(start) is list else [start]
    for st in start:
        result[st] = random_walk(adj, st, maxlen, n)
    return result

def random_walk_weights(graph, pairs, maxlen = 5, n = 1000):
    walks = random_walks(graph, uniqList([p[0] for p in pairs]), maxlen, n)
    return dict((p, sum(walks[p[0]][p[1]])) for p in pairs)

rweights4 = random_walk_weights(net.graph, total_pairs, n = 10000)

###

net.update_vname()
prg = progress.Progress(len(pairs), 'Calculating topological metric', 1)
# results to be stored in these lists:
mut_tgt = []
tgt_mut = []
ww_mut_tgt = []
ww_tgt_mut = []
ew_mut_tgt = []
ew_tgt_mut = []
diff_mut_tgt = []
# iterating all sensitivity/resistance relationships:
for aid, pp in pairs.iteritems():
    prg.step()
    chembls = pp[0][4]
    deltapic50 = pp[0][3]
    # proceed only if we have data about this compound from ChEMBL:
    if bool(len(chembls)) and bool(sum([ch in c.targets for ch in chembls])):
        targets = {}
        # iterating ChEMBLs: there might be multiple ChEMBLs, 
        # but in real it is one compound (at least we assume):
        for ch in chembls:
            # if this ChEMBL ID has targets:
            if ch in c.targets:
                # iterating its targets:
                for tgt in c.targets[ch]:
                    # only if the target in the network
                    if tgt['uniprot'] in net.nodDct:
                        # the vertex id of the target:
                        tgt_vid = net.nodDct[tgt['uniprot']]
                        # collecting pChEMBL values for the target:
                        if tgt_vid not in targets:
                            targets[tgt_vid] = []
                        targets[tgt_vid] += [float(pch) for pch in tgt['pchembl']]
        # taking the maximum of all pChEMBLs:
        targets = dict([(u, max(pchs)) for u, pchs in targets.iteritems()])
        # collect the mutated genes:
        mutated = set(uniqList([net.nodDct[p[0]] for p in pp]))
        w_mut_tgt = []
        w_tgt_mut = []
        for mut_vid in mutated:
            # weights of this relationship:
            w_mut_tgt.append(sum([rweights4[(mut_vid, tgt_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]))
            w_tgt_mut.append(sum([rweights4[(tgt_vid, mut_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]))
            diff_mut_tgt.append((sum([rweights4[(tgt_vid, mut_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]) - \
                sum([rweights4[(mut_vid, tgt_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]), deltapic50))
        ww_mut_tgt += w_mut_tgt
        ww_tgt_mut += w_tgt_mut
        ew_mut_tgt += [(w, deltapic50) for w in w_mut_tgt]
        ew_tgt_mut += [(w, deltapic50) for w in w_tgt_mut]
        # control: all other mutated:
        c_mut_tgt = []
        c_tgt_mut = []
        for mut_vid in all_mutated - mutated:
            # weights of this relationship:
            c_mut_tgt.append(sum([rweights4[(mut_vid, tgt_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]))
            c_tgt_mut.append(sum([rweights4[(tgt_vid, mut_vid)] * pch \
                for tgt_vid, pch in targets.iteritems()]))
        # get ranks:
        r_mut_tgt = [rank(x, c_mut_tgt) for x in w_mut_tgt]
        r_tgt_mut = [rank(x, c_tgt_mut) for x in w_tgt_mut]
        mut_tgt += [(x, deltapic50) for x in r_mut_tgt]
        tgt_mut += [(x, deltapic50) for x in r_tgt_mut]

prg.terminate()