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

import igraph
import cairo
import math

from pypath.ig_drawing import DefaultGraphDrawerFFsupport
from pypath.pypath import Pypath
from pypath.data_formats import best

net = Pypath()
net.load_resources({'arn': best['arn']})
net.genesymbol_labels()

net.graph.vertex_label_family = 'HelveticaNeueLT Std Lt'
net.graph.vertex_label_dist = 30.0
net.graph.edge_label_family = 'Futura Std Book'
net.graph.vs['label_family'] = ['Garamond Premr Pro' \
    if v.index < net.graph.vcount() / 2 \
    else 'HelveticaNeueLT Std Lt' for v in net.graph.vs]

vlabelf = ['Garamond Premr Pro' \
    if v.index < net.graph.vcount() / 2 \
    else 'HelveticaNeueLT Std Lt' for v in net.graph.vs]

net.graph.es['label_family'] = ['Garamond Premr Pro' \
    if e.index < net.graph.ecount() / 2 \
    else 'HelveticaNeueLT Std Lt' for e in net.graph.es]

net.graph.es['label'] = [str(e.index) for e in net.graph.es]

sf = cairo.PDFSurface('igraph-font-family-test.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(12, 12, 1012, 1012)
plot = igraph.plot(net.graph, bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_frame_width = 0, vertex_color = '#22888888', vertex_label_color = '#88222288', 
    edge_width = 1, edge_color = '#33333366', 
    edge_label_color = '#22228888')
plot.redraw()
plot.save()

edgelist = [
    ('sári', 'dénes', 5),
    ('dénes', 'sári', 4),
    ('dénes', 'málna', 4),
    ('málna', 'dénes', 4),
    ('málna', 'niki', 3),
    ('niki', 'málna', 3),
    ('niki', 'dénes', 3),
    ('dénes', 'niki', 3),
    ('sári', 'málna', -2),
    ('málna', 'sári', -1),
    ('sári', 'niki', -3),
    ('niki', 'sári', -1),
    ('franci', 'sári', 5),
    ('sári', 'franci', 4),
    ('dénes', 'franci', 3),
    ('franci', 'dénes', 2),
    ('málna', 'franci', 1),
    ('franci', 'málna', 1),
    ('niki', 'franci', -1),
    ('franci', 'niki', -2)
]

lujza = igraph.Graph.TupleList(edgelist, directed = True, edge_attrs = ['affinity'])

lujza.es['color'] = ['#88222288' if e['affinity'] < 0 else '#22882288' \
    for e in lujza.es]
lujza.es['width'] = [abs(e['affinity'])*3 for e in lujza.es]
lujza.es['arrow_width'] = [abs(e['affinity'])/float(5) for e in lujza.es]
lujza.es['arrow_size'] = [abs(e['affinity']) for e in lujza.es]
lujza.vs['label'] = lujza.vs['name']
lujza.vs['label_family'] = ['Garamond Premr Pro'] * lujza.vcount()
lujza.vs['size'] = [sum(
    [lujza.es[lujza.get_eid(v.index, n.index)]['affinity'] \
        for n in v.neighbors(mode = 'OUT')] + \
    [lujza.es[lujza.get_eid(n.index, v.index)]['affinity'] \
        for n in v.neighbors(mode = 'IN')]
    ) * 4 for v in lujza.vs]

sf = cairo.PDFSurface('lujza_network.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)
plot = igraph.plot(lujza, bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_frame_width = 0, vertex_color = '#228888FF', 
    vertex_label_color = '#777777FF',  vertex_label_family = 'Garamond Premr Pro',
    vertex_label_size = 36,  vertex_label_dist = 2, 
    edge_curved = 0.2)
plot.redraw()
plot.save()

### 
weights_raw = {
    'dénes': {
        'franci': (3, -1),
        'málna': (4, -1),
        'niki': (3, -1),
        'sári': (4, -1)
    },
    'franci': {
        'dénes': (3, -2),
        'málna': (2, -3),
        'niki': (0, -5),
        'sári': (5, -1)
    },
    'málna': {
        'dénes': (4, -1),
        'franci': (2, -1),
        'niki': (4, -2),
        'sári': (2, -4)
    },
    'niki': {
        'dénes': (5, 0),
        'franci': (4, 0),
        'málna': (5, -1),
        'sári': (3, -2)
    },
    'sári': {
        'dénes': (4, -2),
        'franci': (3, -3),
        'málna': (0, -4),
        'niki': (0, -5)
    }
}

weights_norm1 = dict([(src, {}) for src in weights_raw.keys()])
for src, tgt in weights_raw.iteritems():
    pos_sum = float(sum([w[0] for w in tgt.values()]))
    neg_sum = float(sum([w[1] for w in tgt.values()]))
    weights_norm1[src] = dict([(n, w[0] / pos_sum * 5 - w[1] / neg_sum * 3) \
        for n, w in tgt.iteritems()])

edgelist = [it for sl in [[(src, tgt, w) for tgt, w in inner.iteritems()] \
    for src, inner in weights_norm1.iteritems()] for it in sl]

lujza = igraph.Graph.TupleList(edgelist, directed = True, edge_attrs = ['affinity'])

lujza.es['color'] = ['#88222288' if e['affinity'] < 0 else '#22882288' \
    for e in lujza.es]
lujza.es['width'] = [abs(e['affinity'])*3 for e in lujza.es]
lujza.es['arrow_width'] = [abs(e['affinity'])/float(0.5) for e in lujza.es]
lujza.es['arrow_size'] = [abs(e['affinity']) for e in lujza.es]
lujza.vs['label'] = lujza.vs['name']
lujza.vs['label_family'] = ['Garamond Premr Pro'] * lujza.vcount()
lujza.vs['size'] = [abs(sum(
    [lujza.es[lujza.get_eid(v.index, n.index)]['affinity'] \
        for n in v.neighbors(mode = 'OUT')] + \
    [lujza.es[lujza.get_eid(n.index, v.index)]['affinity'] \
        for n in v.neighbors(mode = 'IN')]
    ) * 10) for v in lujza.vs]

sf = cairo.PDFSurface('lujza_network_2.pdf', 1024, 1024)
bbox = igraph.drawing.utils.BoundingBox(124, 124, 900, 900)
plot = igraph.plot(lujza, bbox = bbox, target = sf, 
    drawer_factory = DefaultGraphDrawerFFsupport, 
    vertex_frame_width = 0, vertex_color = '#228888FF', 
    vertex_label_color = '#777777FF',  vertex_label_family = 'Garamond Premr Pro',
    vertex_label_size = 36,  vertex_label_dist = 2, 
    edge_arrow_width = 1,
    edge_curved = 0.2)
plot.redraw()
plot.save()