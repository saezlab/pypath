#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

# external modules:
import math
from igraph import *

# from this module:
import pypath.colorgen as colorgen


def inter_group(g, groups):
    # group-group adjacency matrix weighted by connecting edges,
    # in nested dicts
    allGroups = {}
    for v in g.vs:
        allGroups[v[groups]] = 0
    interGroup = allGroups.copy()
    for key in allGroups.keys():
        interGroup[key] = allGroups.copy()
    for e in g.es:
        interGroup[g.vs[e.source][groups]][g.vs[e.target][groups]] += 1
        interGroup[g.vs[e.target][groups]][g.vs[e.source][groups]] += 1
    return interGroup


def meta_graph(interGroup):
    # returns graph where each node representing one group or module
    meta = Graph(len(interGroup.keys()))
    meta.vs["name"] = [None]
    gr = 0
    for key in interGroup.keys():
        meta.vs[gr]["name"] = key
        gr += 1
    metaEdges = []
    for one, others in interGroup.iteritems():
        for two, weight in others.iteritems():
            if one < two:
                grOne = meta.vs.find(name=one)
                grTwo = meta.vs.find(name=two)
                if weight > 0:
                    metaEdges.append((grOne.index, grTwo.index))
    meta.add_edges(metaEdges)
    return meta


def meta_weight(meta, interGroup):
    meta.es["weight"] = [None]
    for one, others in interGroup.iteritems():
        for two, weight in others.iteritems():
            grOne = meta.vs.find(name=one)
            grTwo = meta.vs.find(name=two)
            e = meta.es.select(_source=grOne.index, _target=grTwo.index)
            if len(e) == 0:
                e = meta.es.select(_source=grTwo.index, _target=grOne.index)
            if len(e) > 0:
                e[0]["weight"] = weight
                i = e[0].index
    meanDeg = sum(meta.vs.degree()) / float(len(meta.vs))
    for e in meta.es:
        degA = meta.vs[e.source].degree()
        degB = meta.vs[e.target].degree()
        e["weight"] = e["weight"] / \
            (math.log(float(min(degA, degB)), 2) + 1.0) ** 2
    return meta.es["weight"]


def meta_layout(g, groups, layout, **kwargs):
    groups = groups2attribute(g, groups)
    # TODO: overlapping groups; now assume discrete groups
    interGroup = inter_group(g, groups)
    meta = meta_graph(interGroup)
    if layout in set(["fr", "fruchterman_reingold"]):
        meta.es["weight"] = meta_weight(meta, interGroup)
        meta.lo = meta.layout(
            "fr",
            weights=meta.es["weight"],
            repulserad=meta.vcount()**2,
            **kwargs)
    else:
        meta.lo = meta.layout(layout, **kwargs)
    dist = layout_distances(meta.lo)
    cl = VertexClustering(g, membership=g.vs[groups])
    subGraphs = cl.subgraphs()
    g.vs["lo"] = [None]
    for i, gr in enumerate(subGraphs):
        thisGroupVs = filter(lambda x: g.vs[x][groups] == i, g.vs.indices)
        cntr = meta.lo[i]
        lo = gr.layout("fr", repulserad=gr.vcount()**2)
        loo = layout_norm(
            lo, [float(cntr[0] - dist[i]), float(cntr[0] + dist[i])],
            [float(cntr[1] - dist[i]), float(cntr[1] + dist[i])])
        for j, v in enumerate(thisGroupVs):
            g.vs[v]["lo"] = loo[j]
    return Layout(coords=g.vs["lo"], dim=2)


def layout_norm(lo, xlim, ylim):
    x = []
    y = []
    for i in lo.coords:
        x.append(i[0])
        y.append(i[1])
    minx = min(x)
    maxx = max(x)
    rangex = maxx - minx
    miny = min(y)
    maxy = max(y)
    rangey = maxy - miny
    scalex = (xlim[1] - xlim[0]) / rangex
    scaley = (ylim[1] - ylim[0]) / rangey
    movex = (xlim[1] + xlim[0]) / 2 - (minx + maxx) / 2
    movey = (ylim[1] + ylim[0]) / 2 - (miny + maxy) / 2
    xx = []
    yy = []
    nLayout = []
    for i, c in enumerate(lo.coords):
        nLayout.append([(c[0] - minx) * scalex + xlim[0],
                        (c[1] - miny) * scaley + ylim[0]])
    return Layout(coords=nLayout, dim=2)


def layout_distances(layout):
    minDs = []
    for i in range(0, len(layout)):
        minD = 100000
        for j in range(0, len(layout)):
            if i != j:
                d = math.sqrt((layout[i][0] - layout[j][0])**2 + (layout[i][
                    1] - layout[j][1])**2)
                minD = d if d < minD else minD
        minDs.append(minD * 0.3)
    return minDs


def groups2attribute(g, groups):
    groupType = groups.__class__.__name__
    if groupType != "str":
        g.vs["groups"] = [None]
        if groupType == "VertexDendrogram":
            g.vs["groups"] = groups.as_clustering().membership
        elif groupType == "VertexClustering":
            g.vs["groups"] == groups.membership
        elif groupType == "list":
            g.vs["groups"] = groups
        return "groups"
    else:
        return groups


def layout_intergroup(g, groups, **kwargs):
    groups = groups2attribute(g, groups)
    # TODO: overlapping groups; now assume discrete groups
    g.es["layout_weight"] = [None]
    for e in g.es:
        indA = e.source
        indB = e.target
        indE = e.index
        if g.vs[indA][groups] == g.vs[indB][groups]:
            g.es[indE]["layout_weight"] = 10.0
        else:
            g.es[indE]["layout_weight"] = 0.0001
    lo = g.layout(
        "fr",
        weights=g.es["layout_weight"],
        repulserad=g.vcount()**2.4,
        **kwargs)
    return Layout(coords=lo, dim=2)


def layout_modular_circle(g, groups, **kwargs):
    return meta_layout(g, groups, "circle", **kwargs)


def layout_modular_fr(g, groups, **kwargs):
    return meta_layout(g, groups, "fr", **kwargs)


def group_colors(g, groups, edges=False):
    groups = groups2attribute(g, groups)
    gr = set(g.vs[groups])
    n = len(gr)
    palette = colorgen.gethexrgbs(n)
    palette.reverse()
    grCol = {}
    for i in gr:
        grCol[i] = palette.pop()
    g.vs[groups + "_color"] = [None]
    for v in g.vs:
        i = v.index
        g.vs[i][groups + "_color"] = grCol[g.vs[i][groups]]
    if edges:
        g.es[groups + "_color"] = [None]
        for e in g.es:
            e[groups + "_color"] = colormix(g.vs[e.source][groups + "_color"],
                                            g.vs[e.target][groups + "_color"])
    return g.vs[groups + "_color"]
