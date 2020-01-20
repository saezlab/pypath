#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#


from future.utils import iteritems
from past.builtins import xrange, range

import pypath.enrich as enrich
import pypath.inputs.main as dataio


class GOEnrichmentSet(enrich.EnrichmentSet):
    
    def __init__(self,
                 aspect,
                 organism=9606,
                 annotation=None,
                 basic_set=None,
                 alpha=0.05,
                 correction_method='hommel'):
        """
        Does not work at the moment as cfisher module should be replaced
        by scipy in the ``pypath.enrich`` module.
        """
        
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
                     for term, cnt in iteritems(self.counts_set)])
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
                   iteritems(getattr(self.annotation, self.aspect.lower()))))

    def get_annot(self, set_names):
        
        return dict(
            filter(lambda x: x[0] in set_names, iteritems(self.basic_set)))

    def count(self, data):
        
        return Counter(flat_list(list(vals) for vals in data.values()))
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
