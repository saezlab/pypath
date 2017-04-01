#!/usr/bin/env python

# Denes Turei EMBL 2017
# turei.denes@gmail.com

import itertools

import pypath
from pypath import data_formats as df

# miRNA-target mRNA interactions:
pa = pypath.PyPath()
pa.init_network(df.mirna_target)

# 6128 interactions between 3035 nodes
# from 4 resources have been loaded
#
# number of miRNAs:
len(set([v for v in pa.graph.vs['name'] if v[:5] == 'MIMAT']))
# 543 # this covers ~25% of human miRNA
#
# number of proteins:
len(set([v for v in pa.graph.vs['name'] if v[:5] != 'MIMAT']))
# 2492
#
# number of references:
len(set(itertools.chain(*pa.graph.es['references'])))
# 3889
#
# interaction-reference associations:
sum(len(e['references']) for e in pa.graph.es)
# 7902

# TF-miRNA (transcriptional regulation of miRNA):
pa = pypath.PyPath()
pa.init_network(df.tf_mirna)

# 1274 interactions between 343 nodes
# from 2 resources have been loaded

# number of miRNAs:
len(set([v for v in pa.graph.vs['name'] if v[:5] == 'MIMAT']))
# 227 # ~12% of human miRNA
#
# number of proteins:
len(set([v for v in pa.graph.vs['name'] if v[:5] != 'MIMAT']))
# 116
#
# number of references:
len(set(itertools.chain(*pa.graph.es['references'])))
# 10
#
# interaction-reference associations:
sum(len(e['references']) for e in pa.graph.es)
# 38

# lncRNA-protein interactions:
pa = pypath.PyPath()
pa.init_network(df.lncrna_protein)
# 89 interactions between 112 nodes
# from 2 resources have been loaded

# number of lncRNAs:
len([v for v in pa.graph.vs['type'] if v == 'lncrna'])
# 40 # ~0.05% of human lncRNA
#
# number of proteins:
len([v for v in pa.graph.vs['type'] if v == 'protein'])
# 72
#
# number of references:
len(set(itertools.chain(*pa.graph.es['references'])))
# 85
#
# interaction-reference associations:
sum(len(e['references']) for e in pa.graph.es)
# 123

# transcriptional regulation of proteins:
pa = pypath.PyPath()
pa.init_network(df.transcription)
# 92380 interactions between 17055 nodes
# from 6 resources have been loaded

# number of TFs:
len(set(e.source for e in pa.graph.es))
# 8536 # this is way too much, something is wrong here
#
# number of regulated proteins:
len([v for v in pa.graph.vs['type'] if v == 'protein'])
# 16780
#
# number of references:
len(set(itertools.chain(*pa.graph.es['references'])))
# 1365
#
# interaction-reference associations:
sum(len(e['references']) for e in pa.graph.es)
# 53336

# let's see the smaller and perhaps higher quality resources:
pa = pypath.PyPath()
pa.init_network(df.transcription, exclude = ['htri', 'encode_dist', 'encode_prox'])
# 4293 interactions between 2884 nodes
# from 3 resources have been loaded

# number of TFs:
len(set(e.source for e in pa.graph.es))
# 1180 # this covers ~50% of human TFs
#
# number of regulated proteins:
len([v for v in pa.graph.vs['type'] if v == 'protein'])
# 2884
#
# number of references:
len(set(itertools.chain(*pa.graph.es['references'])))
# 503
#
# interaction-reference associations:
sum(len(e['references']) for e in pa.graph.es)
# 4509

# protein-protein interactions
pa = pypath.PyPath()
pa.load_omnipath(kinase_substrate_extra = True)
# 4293 interactions between 2884 nodes
# from 3 resources have been loaded

# number of proterins:
pa.graph.vcount()
# 14640
#
# number of interactions:
pa.graph.ecount()
# 119423
#
# number of references:
len(set(itertools.chain(*pa.graph.es['references'])))
# 48742
#
# interaction-reference associations:
sum(len(e['references']) for e in pa.graph.es)
# 199322

# all together:
# this requires allocation of ~3G memory
pa = pypath.PyPath()
pa.load_omnipath(kinase_substrate_extra = True)
pa.load_resources(df.transcription)
pa.load_resources(df.mirna_target)
pa.load_resources(df.tf_mirna)
pa.load_resources(df.lncrna_protein)
# 218283 interactions between 19277 nodes
# from 46 resources have been loaded

# number of molecular entities:
pa.graph.vcount()
# 19277
#
# number of connected pairs:
pa.graph.ecount()
# 218283
#
# number of references:
len(set(itertools.chain(*pa.graph.es['references'])))
# 53961
#
# interaction-reference associations:
sum(len(e['references']) for e in pa.graph.es)
# 260679
