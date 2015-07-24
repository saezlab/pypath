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

# generic modules #

import sys
import os
import re
import cPickle as pickle
import copy
import igraph
import louvain
import cairo
import heapq
import operator

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
from bioigraph.data_formats import best, good, ugly, transcription
import _sensitivity as sens
from bioigraph import progress
from bioigraph import dataio
from bioigraph import go, gsea
from bioigraph.ig_drawing import DefaultGraphDrawerFFsupport


net = bioigraph.BioGraph(9606)
net.init_network(pfile = 'cache/default_plus_acsn.pickle')
net.genesymbol_labels()

psrc = dict((s, ([v['name'] for v in net.graph.vs if s in v['sources']])) for s in net.sources)

g = gsea.GSEA(mapper = net.mapper)
# hallmarks
g.load_collection('H')
# cancer modules
g.load_collection('CM')
# curated
g.load_collection('C2')
# oncogenic
g.load_collection('C6')

enr = gsea.GSEABinaryEnrichmentSet(basic_set = net.graph.vs['name'], gsea = g)

enr_db = {}
for s, uniprots in psrc.iteritems():
    enr.new_set(uniprots)
    enr_db[s] = enr.top_genesets(length = 100)