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

import pypath
from pypath import chembl
from pypath.common import *
from pypath.data_formats import *

import numpy as np
from numpy.random import randn
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as hc
import hcluster as hc2
import matplotlib.patches as mpatches

def clusterd(data, xlab = 'X label', ylab = 'Y label', method = 'ward', 
    fname = 'hc-test.pdf', fontface = 'sans-serif', legend = None,
    textcol = 'black', labelcols = None, orientation = 'right', **kwargs):
    # fig, ax = plt.subplots()
    link = hc2.linkage(data, method = method)
    fig = plt.figure(figsize=(11.3, data.shape[0] * 0.20))
    sns.set(font = 'Helvetica Neue LT Std')
    mpl.rcParams['lines.linewidth'] = 0.5
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

# ## ##

mysql_gelati = (None,'mapping_gelati')
mysql_chembl = (None,'chembl_ebi')

net = pypath.Pypath(9606, mysql=mysql_gelati, name="demo")

net.init_network(pfile = 'cache/plus_phospho.pickle')
net.load_resources(lst={'acsn': ugly['acsn']})

labels = sorted(uniqList([it for sl in [v['sources'] for v in net.graph.vs] for it in sl]))
labels.remove(None)

vpa = [[1 if l in v['sources'] else 0 for v in net.graph.vs] for l in labels]
epa = [[1 if l in e['sources'] else 0 for e in net.graph.es] for l in labels]

vpa = np.array(vpa)
epa = np.array(epa)

clusterd(vpa, labels = labels, orientation = 'left', ylab = 'Pathway resources', 
    method = 'ward', 
    xlab = 'Similarity: protein presence/absence',
    fname = 'resources_proteins_ward.pdf')

clusterd(vpa, labels = labels, orientation = 'left', ylab = 'Pathway resources', 
    method = 'ward', 
    xlab = 'Similarity: interactions',
    fname = 'resources_interactions_ward.pdf')