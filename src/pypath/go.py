#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from future.utils import iteritems
from past.builtins import xrange, range

import sys
import imp
from collections import Counter, OrderedDict

# from this package:
import pypath.dataio as dataio
import pypath.progress as progress
import pypath.enrich as enrich
from pypath.common import *


class GeneOntology(object):
    
    all_relations = {
        'is_a', 'part_of', 'occurs_in', 'regulates',
        'positively_regulates', 'negatively_regulates',
    }
    
    def __init__(self):
        """
        Loads data about Gene Ontology terms and their relations.
        """
        
        terms = dataio.go_terms_quickgo()
        
        self.ancestors = self._merge_aspects(
            dataio.go_ancestors_quickgo()
        )
        self.descendants = self._merge_aspects(
            dataio.go_descendants_quickgo()
        )
        
        self.name = dict(i for ii in terms.values() for i in iteritems(ii))
        self.term = dict(reversed(i) for i in iteritems(self.name))
    
    def reload(self):
        """Reloads the object from the module level."""

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def get_name(self, term):
        """
        For a GO accession number returns the name of the term.
        """
        
        return None if term not in self.name else self.name[term]

    def get_term(self, name):
        """
        For a GO term name returns its GO accession number.
        """
        
        return None if name not in self.term else self.term[name]
    
    @staticmethod
    def _merge_aspects(dct):
        
        dct['P'].update(dct['C'])
        dct['P'].update(dct['F'])
        return dct['P']
    
    def subgraph_nodes(self, direction, terms, relations = None):
        """
        Returns a set of all nodes either in the subgraph of ancestors 
        or descendants of a single term or a set of terms.
        
        :param str direction:
            Possible values: `ancestors` or `descendants`.
        """
        
        relations = relations or self.all_relations
        
        if isinstance(terms, basestring):
            
            terms = {terms}
        
        graph = getattr(self, direction)
        subgraph = set()
        
        for term in terms:
            
            for related, relation in graph[term]:
                
                if relation not in relations:
                    
                    continue
                
                if related not in subgraph:
                    
                    subgraph.update(
                        self.subgraph_nodes(direction, related, relations)
                    )
                    subgraph.add(related)
        
        return subgraph
    
    def all_ancestors(self, terms, relations = None):
        """
        Returns a set of all ancestors of a single term or a set of terms.
        """
        
        return self.subgraph_nodes('ancestors', terms, relations)
    
    def all_descendants(self, terms, relations = None):
        """
        Returns a set of all descendants of a single term or a set of terms.
        """
        
        return self.subgraph_nodes('descendants', terms, relations)


class GOAnnotation(object):
    
    def __init__(self, organism = 9606, ontology = None):
        """
        For one organism loads Gene Ontology annotations, in addition it
        accepts or creates a ``GeneOntology`` object. It 
        """
        
        self.ontology = ontology or GeneOntology()
        
        annot = dataio.go_annotations_goa(organism = organism)
        self.c = annot['C']
        self.f = annot['F']
        self.p = annot['P']

    def get_name(self, term):
        """
        For a GO accession number returns the name of the term.
        """
        
        return self.ontology.get_name(term)

    def get_term(self, name):
        """
        For a GO term name returns its GO accession number.
        """
        
        return self.ontology.get_term(name)

    def get_annot(self, uniprot, aspect):
        """
        For a UniProt ID returns its direct annotations from one aspect
        of Gene Ontology.
        Returns set.
        """
        
        annot = getattr(self, aspect.lower())
        return set() if uniprot not in annot else annot[uniprot]

    def get_annots(self, uniprot):
        """
        For a UniProt ID returns its direct annotations from all aspects
        of Gene Ontology.
        Returns set.
        """
        
        return set.union(
            self.get_annot(uniprot, asp)
            for asp in ('C', 'F', 'P')
        )
    
    def select_by_go(self):
        
        pass


def annotate(graph, organism = 9606, aspects = ('C', 'F', 'P')):
    """
    Adds Gene Ontology annotations to the nodes of a graph.
    
    :param igraph.Graph graph:
        Any ``igraph.Graph`` object with uniprot IDs
        in its ``name`` vertex attribute.
    """
    
    aspects = aspects if type(aspects) in {list, tuple} else (aspects, )
    
    graph.vs['go'] = [
        {'C': set(), 'F': set(), 'P': set()}
        for _ in xrange(graph.vcount())
    ]
    
    terms, annot = dataio.go_annotations_goa(organism = organism)
    
    prg = progress.Progress(graph.vcount(), 'Loading GO annotations', 9)
    
    for v in graph.vs:
        
        prg.step()
        
        for asp in aspects:
            
            if v['name'] in annot[asp]:
                
                v['go'][asp] = annot[asp][v['name']]
    
    prg.terminate()


# old name as synonym
load_go = annotate


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
