#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2017 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import sys
from collections import Counter, OrderedDict

# from this package:
from pypath import dataio
import pypath.progress
from pypath import enrich
from pypath.common import *


class GOAnnotation(object):
    def __init__(self, organism=9606):
        go = dataio.get_go_quick(organism=organism)
        self.c = go['terms']['C']
        self.f = go['terms']['F']
        self.p = go['terms']['P']
        self.name = go['names']
        self.term = dict([(v, k) for k, v in self.name.iteritems()])

    def get_name(self, term):
        return None if term not in self.name else self.name[term]

    def get_term(self, name):
        return None if name not in self.term else self.term[name]

    def get_annot(self, uniprot, aspect):
        dic = getattr(self, aspect.lower())
        return [] if uniprot not in dic else dic[uniprot]

    def get_annots(self, uniprot):
        result = {}
        for asp in ['C', 'F', 'P']:
            result[asp.upper()] = self.get_annot(uniprot, asp)
        return result


def load_go(graph, aspect=['C', 'F', 'P']):
    '''
    @graph : igraph.Graph
    Any igraph.Graph object with uniprot IDs in its `name` vertex attribute.
    '''
    aspect = aspect if type(aspect) is list else [aspect]
    graph.vs['go'] = [{'C': [], 'F': [], 'P': []} for _ in graph.vs]
    go = dataio.get_go_goa()
    prg = progress.Progress(graph.vcount(), 'Loading GO annotations', 9)
    for v in graph.vs:
        prg.step()
        for asp in aspect:
            if v['name'] in go[asp]:
                v['go'][asp] = go[asp][v['name']]
    prg.terminate()


class GOEnrichmentSet(enrich.EnrichmentSet):
    def __init__(self,
                 aspect,
                 organism=9606,
                 annotation=None,
                 basic_set=None,
                 alpha=0.05,
                 correction_method='hommel'):
        self.aspect = aspect
        self.organism = organism
        self.alpha = alpha
        self.correction_method = correction_method
        self.annotation = GOAnnotation(organism=self.organism) \
            if annotation is None else annotation
        self.basic_set = basic_set if basic_set is not None \
            else self.get_basic_set()
        self.counts_pop = self.count(self.basic_set)
        self.pop_size = len(self.basic_set)
        self.set_annot = None
        self.set_size = None
        self.counts_set = None
        self.top_terms = self.top_names
        self.top_accessions = self.top_ids

    def new_set(self, set_names=None, set_annot=None):
        self.set_annot = set_annot if set_annot is not None \
            else self.get_annot(set_names)
        self.set_size = len(self.set_annot)
        self.counts_set = self.count(self.set_annot)
        self.calculate()

    def calculate(self):
        data = dict([(term, (cnt, self.counts_pop[term], self.set_size,
                             self.annotation.name[term]))
                     for term, cnt in self.counts_set.iteritems()])
        enrich.EnrichmentSet.__init__(
            self,
            data,
            self.pop_size,
            alpha=self.alpha,
            correction_method=self.correction_method)

    def get_basic_set(self):
        swissprots = set(
            dataio.all_uniprots(
                organism=self.organism, swissprot='yes'))
        return dict(
            filter(lambda x: x[0] in swissprots,
                   getattr(self.annotation, self.aspect.lower()).iteritems()))

    def get_annot(self, set_names):
        return dict(
            filter(lambda x: x[0] in set_names, self.basic_set.iteritems()))

    def count(self, data):
        return Counter(flatList(list(vals) for vals in data.values()))
        # return dict((name, count/float(len(data))) for name, count in
        # cnt.iteritems())

    def __str__(self):
        if self.set_annot is None:
            resp = '\n\t:: No calculations performed yet. Please define '\
                'a set of genes with `new_set()`.\n\n'
        else:
            resp = '\n :: Top significantly enriched terms (max. 10):\n\n\t'\
                + '\n\t'.join([t[0].upper() + t[1:] for t in
                               self.top_terms(length=10, significant=True)]) + '\n'
        return resp
