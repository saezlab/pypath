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

#
    # This can be run non interactively
    # to regenerate the article figures.
#

# generic modules #

from collections import Counter
import locale
locale.setlocale(locale.LC_ALL, 'en_GB.UTF-8')

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
from scipy.cluster import hierarchy as hc

# from pypath #

import pypath
from pypath.common import *
from pypath.data_formats import omnipath as best
from pypath import data_formats
from pypath import dataio
from pypath import plot

# functions

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

# parameters
omnipath = 'All'
prefix = "int"
lab_size = (18, 21)
axis_lab_size = 36
table2file = 'int_curation_stats_stripped.tex'
fisherFile = 'fisher_tests_int'
fontW = 'medium'

console(':: Creating new network object')
net = pypath.PyPath(9606)

# see `default_network.py` for the initialization of this default network
console(':: Loading network')
net.init_network(data_formats.interaction_htp)

# Table 2 for the article
console(':: Table 2: statistics of interaction and curation content and overlaps,'\
    'writing LaTeX table into %s'%table2file)
net.curation_tab(latex_hdr = True, fname = table2file)

# loading additional data onto the network
console(':: Loading annotations')
net.genesymbol_labels()
net.set_tfs()
net.set_receptors()
net.set_chembl_mysql('chembl_ebi')
# net.set_drugtargets()
net.set_kinases()
net.set_druggability()
net.load_disgenet()
net.load_corum()
sens.in_complex(net)
net.in_complex()
net.load_ptms()
net.read_list_file(pypath.data_formats.cgc)
data_formats.intogen_cancer.inFile = '../../a2/src/pypath/data/intogene_cancerdrivers.tsv'
net.read_list_file(pypath.data_formats.intogen_cancer)

# generating individual networks for all resorces
console(':: Separating original resources')
sep = net.separate()

# lists of all receptors/druggable proteins/kinases/tfs/disease related proteins
console(':: Loading lists of all receptors/druggable proteins/kinases/'\
    'tfs/disease related proteins')
net.lists['rec'] = uniqList(flatList([net.mapper.map_name(rec, 'genesymbol', 'uniprot') \
    for rec in dataio.get_hpmr()]))

net.lists['dgb'] = uniqList(flatList([net.mapper.map_name(dgb, 'genesymbol', 'uniprot') \
    for dgb in dataio.get_dgidb()]))

net.lists['kin'] = uniqList(flatList([net.mapper.map_name(kin, 'genesymbol', 'uniprot') \
    for kin in dataio.get_kinases()]))

net.lists['tfs'] = uniqList(flatList([net.mapper.map_name(tf, 'ensg', 'uniprot') \
    for tf in dataio.get_tfcensus()['ensg']]))

net.lists['dis'] = uniqList(flatList([\
    net.mapper.map_name(dis['genesymbol'], 'genesymbol', 'uniprot') \
    for dis in dataio.get_disgenet()]))

# defining the proteome as the set of all human swissprot ids
console(':: Loading the human proteome')
proteome = dataio.all_uniprots(swissprot = 'yes')

fi = open(fisherFile, 'w')
# Fisher's exact test for enrichment of disease related proteins
# in OmniPath compared to their ratio in the whole proteome
console(':: Fisher\'s exact test for enrichment of disease related proteins in the network'\
    'compared to their abundance in the proteome')
contDisg = np.array([[len(proteome), net.graph.vcount()], [len(net.lists['dis']), len([1 for v in net.graph.vs if len(v['dis']) > 0])]])
fi.write('Disease related proteins:\t%s\t%s\n' % stats.fisher_exact(contDisg))

# Fisher's exact test for enrichment of cancer driver proteins
# in OmniPath compared to their ratio in the whole proteome
console(':: Fisher\'s exact test for enrichment of cancer drivers in the network'\
    'compared to their abundance in the proteome')
contCanc = np.array([[len(proteome), net.graph.vcount()], [len(uniqList(net.lists['CancerGeneCensus'] + net.lists['IntOGen'])), 
    len(set(net.lists['CancerGeneCensus'] + net.lists['IntOGen']) & set(net.graph.vs['name']))]])
fi.write('Cancer driver proteins:\t%s\t%s\n' % stats.fisher_exact(contCanc))

# Fisher's exact test for enrichment of druggable proteins
# in OmniPath compared to their ratio in the whole proteome
console(':: Fisher\'s exact test for enrichment of druggable proteins in the network'\
    'compared to their abundance in the proteome')
contDgb = np.array([[len(proteome), net.graph.vcount()], [len(net.lists['dgb']), len([1 for v in net.graph.vs if v['dgb']])]])
fi.write('Druggable proteins:\t%s\t%s\n' % stats.fisher_exact(contDgb))

# Fisher's exact test for enrichment of receptors
# in OmniPath compared to their ratio in the whole proteome
console(':: Fisher\'s exact test for enrichment of receptors in the network'\
    'compared to their abundance in the proteome')
contRec = np.array([[len(proteome), net.graph.vcount()], [len(net.lists['rec']), len([1 for v in net.graph.vs if v['rec']])]])
fi.write('Receptors:\t%s\t%s\n' % stats.fisher_exact(contRec))

# Fisher's exact test for enrichment of TFs
# in OmniPath compared to their ratio in the whole proteome
console(':: Fisher\'s exact test for enrichment of TFs in the network'\
    'compared to their abundance in the proteome')
contTf = np.array([[len(proteome), net.graph.vcount()], [len(net.lists['tfs']), len([1 for v in net.graph.vs if v['tf']])]])
fi.write('Transcription factors:\t%s\t%s' % stats.fisher_exact(contTf))

fi.close()

# source-vcount barplot
d = zip(*[(s, len([v for v in net.graph.vs if s in v['sources']])) for s in net.sources] + \
    [(omnipath, net.graph.vcount())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'proteins_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Number of proteins', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Interaction databases', order = 'y',
    y_break = (0.41, 0.15))

vcount_ordr = list(bp.ordr)
bp.finish()

# source-ecount barplot
console(':: Plotting `Number of proteins` barplot')
d = zip(*[(s, len([e for e in net.graph.es if s in e['sources']])) for s in net.sources] + \
    [(omnipath, net.graph.ecount())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'interactions_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Number of interactions', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Interaction databases', order = vcount_ordr, fin = False, desc = False, 
    #y_break = (0.24, 0.24)
    )

bp.ax.yaxis.labelpad = 15
bp.finish()

# density sens.barplot
console(':: Plotting `Graph density` barplot')
d = zip(*[(s, g.density()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.density())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'density_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Graph density', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Interaction databases', order = 'y', 
    y_break = (0.15, 0.1))

# density with same order as protein count:
d = zip(*[(s, g.density()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.density())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'density_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Graph density', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Interaction databases', order = vcount_ordr, desc = False, 
    y_break = (0.15, 0.1))

# vcount - ecount scatterplot:
topdata = {
'Proteins': [len([v for v in net.graph.vs if s in v['sources']]) for s in net.sources],
'Interactions': [len([e for e in net.graph.es if s in e['sources']]) for s in net.sources],
'Density': [sep[s].density() for s in net.sources],
'Redensity': [1.0 / sep[s].density() for s in net.sources],
'Transitivity': [sep[s].transitivity_undirected() for s in net.sources],
'Diameter': [sep[s].diameter() for s in net.sources],
'Number of receptors': [len(set(sep[s].vs['name']) & set(net.lists['rec'])) for s in net.sources],
'Number of TFs': [len(set(sep[s].vs['name']) & set(net.lists['tfs'])) for s in net.sources],
'Number of kinases': [len(set(sep[s].vs['name']) & set(net.lists['kin'])) for s in net.sources],
'Number of druggable proteins': [len(set(sep[s].vs['name']) & set(net.lists['dgb'])) for s in net.sources],
'Number of disease related proteins': [len(set(sep[s].vs['name']) & set(net.lists['dis'])) for s in net.sources],
'Number of cancer drivers': [len(set(sep[s].vs['name']) & \
        (set(net.lists['IntOGen']) | set(net.lists['CancerGeneCensus']))) \
    for s in net.sources],
'Complexes': [len(sens.complexes_in_network(g)) \
    for s, g in sep.iteritems()],
'PTMs': [len([1 for e in g.es if len(e['ptm']) != 0]) \
    for s, g in sep.iteritems()],
'Database': net.sources
}
topdf = pd.DataFrame(topdata, index = net.sources)


topdf.sort_values(by = ['Proteins'], inplace = True)

#import pandas as pd
#from pypath import plot
#import cPickle as pickle

#topdf = pickle.load(open('cache/topdf.pickle', 'rb'))

console(':: Plotting `Proteins vs. number of interactions` scatterplot')
sp = plot.ScatterPlus('Proteins', 'Interactions', sizes = 'Redensity', 
        labels = 'Database', data = topdf,
        xlog = True, ylog = True, xlim = [-20.0, 5000.0], 
        ylim = [-20.0, 20000.0], 
        fname = 'vcount-ecount_int-log-4.pdf', font_family = 'Helvetica Neue LT Std', 
        font_style = 'normal', font_weight = 'medium', font_variant = 'normal',
        font_stretch = 'normal', confi = True, 
        xlab = 'Number of proteins', ylab = 'Number of interacting pairs', 
        axis_lab_size = 23.76, annotation_size = 11.8, size = 20.0, size_scaling = 0.66,
        lab_angle = 90, lab_size = (11.8, 13.86), color = '#007b7f', 
        order = False, desc = True, legend = True, legtitle = 'Density',
        size_to_value = lambda x: 1 / x,
        value_to_size = lambda x: 1 / x,
        fin = True)

console(':: Plotting `Disease related proteins vs. cancer drivers` scatterplot')
sp = plot.ScatterPlus('Number of disease related proteins', 
        'Number of cancer drivers', sizes = 'Proteins', 
        labels = 'Database', data = topdf,
        xlim = [30.0, 2300.0], ylim = [0.0, 2000.0],
        xlog = True, ylog = True, 
        fname = 'dis-cancer_int-log.pdf', font_family = 'Helvetica Neue LT Std', 
        font_style = 'normal', font_weight = 'medium', font_variant = 'normal',
        font_stretch = 'normal', confi = True, 
        xlab = 'Number of\ndisease related proteins', ylab = 'Number of cancer drivers', 
        axis_lab_size = 23.76, annotation_size = 11.8, size = 20.0, size_scaling = 0.66,
        lab_angle = 90, lab_size = (11.8, 13.86), color = '#007b7f', 
        order = False, desc = True, legend = True, legtitle = 'Total number\nof proteins',
        legstrip = (1,1),
        fin = True)

console(':: Plotting `Receptors vs. TFs` scatterplot')
sp = plot.ScatterPlus('Number of receptors', 
        'Number of TFs', sizes = 'Proteins', 
        labels = None, data = topdf,
        xlog = 'symlog', ylog = 'symlog', 
        xlim = [50.0, 1000.0], ylim = [0.0, 10000.0],
        fname = 'tf-rec_int-log.pdf', font_family = 'Helvetica Neue LT Std', 
        font_style = 'normal', font_weight = 'medium', font_variant = 'normal',
        font_stretch = 'normal', confi = True, 
        xlab = 'Number of receptors', ylab = 'Number of TFs', 
        axis_lab_size = 23.76, annotation_size = 11.8, size = 20.0, size_scaling = 0.66,
        lab_angle = 90, lab_size = (11.8, 13.86), color = '#007b7f', 
        order = False, desc = True, legend = True, legtitle = 'Total number\nof proteins',
        fin = True)

sp = plot.ScatterPlus('Number of receptors', 
        'Number of TFs', sizes = 'Proteins', 
        labels = 'Database', data = topdf,
        xlim = [-10.0, 550.0], ylim = [-5.0, 400.0],
        fname = 'tf-rec_int.pdf', font_family = 'Helvetica Neue LT Std', 
        font_style = 'normal', font_weight = 'medium', font_variant = 'normal',
        font_stretch = 'normal', confi = True, 
        xlab = 'Number of receptors', ylab = 'Number of TFs', 
        axis_lab_size = 23.76, annotation_size = 8.0, size = 20.0, size_scaling = 0.66,
        lab_angle = 90, lab_size = (11.8, 13.86), color = '#007b7f', 
        order = False, desc = True, legend = True, legtitle = 'Total number\nof proteins',
        legstrip = (1,1),
        fin = True)

sp = plot.ScatterPlus('Complexes', 
        'PTMs', sizes = 'Interactions', 
        labels = 'Database', data = topdf,
        xlog = 'symlog', ylog = 'symlog',
        #xlim = [-50.0, 680.0], ylim = [10.0, 10000.0],
        xlim = [2.0, 1800.0], ylim = [0.0, 5000.0],
        fname = 'comp-ptm_int.pdf', font_family = 'Helvetica Neue LT Std', 
        font_style = 'normal', font_weight = 'medium', font_variant = 'normal',
        font_stretch = 'normal', confi = True, 
        xlab = 'Number of complexes', ylab = 'Number of PTMs', 
        axis_lab_size = 23.76, annotation_size = 8.0, size = 20.0, size_scaling = 0.66,
        lab_angle = 90, lab_size = (11.8, 13.86), color = '#007b7f', 
        order = False, desc = True, legend = True, legtitle = 'Total number\nof interactions',
        legstrip = (1,1), legloc = 2, 
        fin = True)

# transitivity barplot
console(':: Plotting `Graph transitivity` barplot')
d = zip(*[(s, g.transitivity_undirected()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.transitivity_undirected())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'transitivity_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Graph global transitivity', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Interaction databases', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'transitivity_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Graph global transitivity', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'],
    xlab = 'Interaction databases', order = vcount_ordr, desc = False)

# diameter barplot
console(':: Plotting `Graph diameter` barplot')
d = zip(*[(s, g.diameter()) for s, g in sep.iteritems()] + \
    [(omnipath, net.graph.diameter())])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'diameter_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Graph diameter', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'diameter_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Graph diameter', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False)

# receptors barplot
console(':: Plotting `Number of receptors` barplot')
d = zip(*[(s, len([v for v in g.vs if v['rec']])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec']))])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptors_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Number of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = 'y',
    #y_break = (0.6, 0.15)
    )
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptors_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Number of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False,
    #y_break = (0.6, 0.15)
    )

# receptors prop barplot
console(':: Plotting `Proportion of receptors` barplot')
d = zip(*[(s, len([v for v in g.vs if v['rec']])/float(g.vcount())*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec'])/float(net.graph.vcount())*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorprop_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = 'y',
    #y_break = (0.60, 0.12)
    )
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorprop_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of receptors', color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False,
    #y_break = (0.60, 0.12)
    )

# receptors coverage barplot
console(':: Plotting `Coverage of receptors` barplot')
d = zip(*[(s, len([v for v in g.vs if v['rec']])/float(len(net.lists['rec']))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['rec'])/float(len(net.lists['rec']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorcov_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of receptors covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'receptorcov_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of receptors covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False)

# transcription factors sens.barplot
console(':: Plotting `Number of TFs` barplot')
d = zip(*[(s, len([v for v in g.vs if v['tf']])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf']))])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfs_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Number of TFs', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = 'y',
    y_break = False)
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfs_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = 'Number of TFs', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False, 
    y_break = False)

# transcription factor prop sens.barplot
console(':: Plotting `Proportion of TFs` barplot')
d = zip(*[(s, len([v for v in g.vs if v['tf']])/float(g.vcount())*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf'])/float(net.graph.vcount())*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfprop_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of transcription factors',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfprop_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of transcription factors',
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False)

console(':: Plotting `Coverage of TFs` barplot')
d = zip(*[(s, len([v for v in g.vs if v['tf']])/float(len(net.lists['tfs']))*100) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum(net.graph.vs['tf'])/float(len(net.lists['tfs']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfcov_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of TFs covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = 'y', 
    y_break = False)
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'tfcov_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of TFs covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False, 
    y_break = False)

# receptors & TFs coverage barplot together
console(':: Plotting `Receptors-TFs-kinases-druggable proteins` grouped barplot')
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
nul = [t.set_fontsize(lab_size[0] * 0.7) for t in leg.texts]
nul = [t.set_fontsize(lab_size[0] * 0.7) for t in ax.get_yticklabels()]
xlim = ax.set_xlim(-0.5, float(len(net.sources)) + 1.0)
ylab = ax.set_ylabel(r'% of proteins covered', fontsize = axis_lab_size * 0.45)
xlab = ax.set_xlabel('Interaction databases', fontsize = axis_lab_size * 0.45)
xticks = ax.set_xticks(x + 2 * w)
xticklabs = ax.set_xticklabels(d['Database'], rotation = 90, fontsize = lab_size[0] * 0.6)
ax.xaxis.label.set_fontweight('medium')
ax.yaxis.label.set_fontweight('medium')
ax.xaxis.grid(False)

plt.tight_layout()

fig.savefig('interesting-proteins-cov_int-by-db.pdf')

plt.close()

# disease associations:
console(':: Plotting `Disease related proteins` barplot')
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
    data = None, fname = 'discov_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of disease-gene' + '\nassociations covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = 'y', 
    y_break = False)
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'discov_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    ylab = r'% of disease-gene' + '\nassociations covered', 
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False,
    y_break = False)

# number of complexes sens.barplot
console(':: Plotting `Complexes` barplot')
d = zip(*[(s, len(sens.complexes_in_network(g))) \
    for s, g in sep.iteritems()] + \
    [(omnipath, len(sens.complexes_in_network(net.graph)))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'complexes_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of complexes',
    #\n(out of %u)'%\
        #len(sens.complexes_in_network(net.graph)), 
    xlab = 'Interaction databases', order = 'y', y_break = (0.55, 0.15))
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'complexes_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of complexes',
    #\n(out of %u)'%\
        #len(sens.complexes_in_network(net.graph)), 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False, 
    y_break = (0.55, 0.15))

# number of ptms sens.barplot
console(':: Plotting `PTMs covered` barplot')
d = zip(*[(s, sum([len(e['ptm']) for e in g.es])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, sum([len(e['ptm']) for e in net.graph.es]))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ptms_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of PTMs',
    #\n(out of %u PTMs)'%\
        #sum([len(e['ptm']) for e in net.graph.es]), 
    xlab = 'Interaction databases', order = 'y',
    #y_break = (0.65, 0.15)
    )
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ptms_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Number of PTMs',
    #\n(out of %u PTMs)'%\
        #sum([len(e['ptm']) for e in net.graph.es]), 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False, 
    #y_break = (0.65, 0.15)
    )

# number of ptms sens.barplot
console(':: Plotting `Interactions with PTM` barplot')
d = zip(*[(s, len([1 for e in g.es if len(e['ptm']) != 0])) \
    for s, g in sep.iteritems()] + \
    [(omnipath, len([1 for e in net.graph.es if len(e['ptm']) != 0]))])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'havingptm_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Interactions having PTM', 
    xlab = 'Interaction databases', order = 'y',
    #y_break = (0.49, 0.1)
    )
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'havingptm_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = 'Interactions having PTM', 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False, 
    #y_break = (0.49, 0.1)
    )

# TODO difference between COSMIC and Intogene

# cosmic cancer gene census coverage sens.barplot
console(':: Plotting `Cancer Gene Census coverage` barplot')
d = zip(*[(s, len(set(net.lists['CancerGeneCensus']) & set(g.vs['name'])) / \
        float(len(net.lists['CancerGeneCensus']))*100) \
        for s, g in sep.iteritems()] + \
    [(omnipath, len(set(net.lists['CancerGeneCensus']) & set(net.graph.vs['name'])) / \
        float(len(net.lists['CancerGeneCensus']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ccgccov_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of CGC genes', 
    xlab = 'Interaction databases', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    fname = 'ccgccov_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of CGC genes', 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False)

# intogene cancer driver coverage sens.barplot
console(':: Plotting `IntOGen cancer drivers coverage` barplot')
d = zip(*[(s, len(set(net.lists['IntOGen']) & set(g.vs['name'])) / \
        float(len(net.lists['IntOGen']))*100) \
        for s, g in sep.iteritems()] + \
    [(omnipath, len(set(net.lists['IntOGen']) & set(net.graph.vs['name'])) / \
        float(len(net.lists['IntOGen']))*100)])
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov_int-by-db.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of IntOGen genes', 
    xlab = 'Interaction databases', order = 'y')
bp = plot.Barplot(x = d[0], y = d[1], 
    data = None, fname = 'intocov_int-by-db-o.pdf', lab_size = lab_size, 
    axis_lab_size = axis_lab_size, font_weight = fontW,
    color = ['#007B7F'] * (len(d[1]) - 1) + ['#6EA945'], 
    ylab = r'% of IntOGen genes', 
    xlab = 'Interaction databases', order = vcount_ordr, desc = False)

# ## #
# PTMs in PTM resources and in the network
# ## #
console(':: Calculating PTM database statistics')

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

net.uniq_ptms()

d = zip(*[(s, len(uniqList(m)),
            len([p for e in net.graph.es for p in e['ptm'] if s.lower() in [ps.lower() for ps in p.sources]])
        ) for s, m in ptms.items()] + \
        [('All', len(uniqList(flatList(ptms.values()))), len(uniqList(flatList(net.graph.es['ptm']))))])

d = dict(zip(['Database', 'All', 'In_network'], 
    [list(dd) for dd in d]))
d = pd.DataFrame(d, index = d['Database'])

d = d.sort_values(by = 'All', ascending = False)

console(':: Plotting PTM database statistics')
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
    mpl.patches.Patch(color=c['inn'], label='Mapped on Interaction databases')]
leg = ax.legend(handles = lhandles)
xlim = ax.set_xlim(-0.5, 9.25)
ylab = ax.set_ylabel('Enzyme-substrate\ninteractions', fontsize = axis_lab_size * 0.66, fontweight = fontW)
xlab = ax.set_xlabel('PTM resources', fontsize = axis_lab_size * 0.66, fontweight = fontW)
xticks = ax.set_xticks(x + w)
xticklabs = ax.set_xticklabels(d['Database'], rotation = 90, fontsize = lab_size[0] * 0.66)
nul = [t.set_fontsize(lab_size[0]*0.77) for t in ax.get_yticklabels()]
nul = [t.set_fontsize(lab_size[0]*0.77) for t in leg.texts]
ax.xaxis.grid(False)

plt.tight_layout()

fig.savefig('ptms-by-ptmdb_int.pdf')

plt.close()

## ## ##

##

# ref ecount barplot
console(':: Calculating HTP/LTP statistics')
net.init_network(pfile = 'cache/default_network_raw.pickle')
refc = Counter(flatList((r.pmid for r in e['references']) for e in net.graph.es))

# percentage of high throughput interactions
htdata = {}
htsrcs_prev = set(net.sources)
for htlim in reversed(xrange(5, 201)):
    htrefs = set([i[0] for i in refc.most_common() if i[1] > htlim])
    htedgs = [e.index for e in net.graph.es if \
        len(set([r.pmid for r in e['references']]) - htrefs) == 0]
    htsrcs = uniqList(flatList([net.graph.es[e]['sources'] for e in htedgs]))
    htsrcs_new = set(uniqList(flatList([net.graph.es[e]['sources'] for e in htedgs])))
    diff = htsrcs_new - htsrcs_prev
    if len(diff):
        print list(diff), 'has HTP over %u' % htlim
    htsrcs_prev = htsrcs_new
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

console(':: Plotting HTP/LTP statistics')
fig, axs = plt.subplots(5, figsize = (10, 20), sharex=True)
pl = axs[0].plot(sorted(htdata.keys()), [htdata[h]['rnum'] for h in sorted(htdata.keys())], 
    '-', color = '#007B7F')
ylab = axs[0].set_ylabel('Number of references', fontsize = axis_lab_size * 0.45, fontweight = fontW)
pl = axs[1].plot(sorted(htdata.keys()), [htdata[h]['enum'] for h in sorted(htdata.keys())], 
    '-', color = '#6EA945')
ylab = axs[1].set_ylabel('Number of edges', fontsize = axis_lab_size * 0.45, fontweight = fontW)
pl = axs[2].plot(sorted(htdata.keys()), [htdata[h]['snum'] for h in sorted(htdata.keys())], 
    '-', color = '#DA0025')
ylab = axs[2].set_ylabel('Number of resources', fontsize = axis_lab_size * 0.45, fontweight = fontW)
pl = axs[3].plot(sorted(htdata.keys()), [htdata[h]['lenum'] for h in sorted(htdata.keys())], 
    '-', color = '#996A44')
ylab = axs[3].set_ylabel('LT network edge count', fontsize = axis_lab_size * 0.45, fontweight = fontW)
pl = axs[4].plot(sorted(htdata.keys()), [htdata[h]['lvnum'] for h in sorted(htdata.keys())], 
    '-', color = '#FCCC06')
xlab = axs[4].set_xlabel('HT limit', fontsize = axis_lab_size * 0.66, fontweight = fontW)
ylab = axs[4].set_ylabel('LT network node count', fontsize = axis_lab_size * 0.45, fontweight = fontW)
fig.savefig('ht-limits_int.pdf')
plt.close(fig)

console(':: Done. Bye.')
