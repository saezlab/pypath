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
import numpy as np

import pypath.dataio as dataio
import pypath.progress as progress
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
    
    aspects = ('C', 'F', 'P')
    
    def __init__(self, organism = 9606, ontology = None):
        """
        For one organism loads Gene Ontology annotations, in addition it
        accepts or creates a ``GeneOntology`` object.
        """
        
        self.ontology = ontology or GeneOntology()
        
        annot = dataio.go_annotations_goa(organism = organism)
        self.c = annot['C']
        self.f = annot['F']
        self.p = annot['P']
        
        self._ancestors_annotate()
    
    def reload(self):
        """Reloads the object from the module level."""
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def _ancestors_annotate(self):
        
        for asp in self.aspects:
            
            setattr(
                self,
                '%s_full' % asp.lower(),
                dict(
                    (
                        uniprot,
                        self.ontology.all_ancestors(annot)
                    )
                    for uniprot, annot in
                    iteritems(getattr(self, asp.lower()))
                )
            )
    
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
            for asp in self.aspects
        )
    
    def _has_term(self, uniprot, term, aspect):
        
        annot = getattr(self, '%s_full' % aspect.lower())
        
        return uniprot in annot and term in annot[uniprot]
    
    def _has_any_term(self, uniprot, terms, aspect):
        
        annot = getattr(self, '%s_full' % aspect.lower())
        
        return uniprot in annot and terms & annot[uniprot]
    
    def has_term(self, uniprot, term):
        """
        Tells if an UniProt ID is annotated with a GO term.
        """
        
        return any(
            self._has_term(uniprot, term, aspect)
            for asp in self.aspects
        )
    
    def has_any_term(self, uniprot, terms):
        """
        Tells if an UniProt ID is annotated with any of a set of GO terms.
        """
        
        return any(
            self._has_any_term(uniprot, terms, aspect)
            for asp in self.aspects
        )
    
    def i_select_by_term(self, uniprots, term):
        """
        Accepts a list of UniProt IDs and one or more gene ontology terms
        and returns a set of indices of those UniProts which are annotated
        with any of the terms.
        
        :param str,set term:
            A single GO term or set of terms.
        """
        
        method = self.has_any_term if isinstance(term, set) else self.has_term
        
        return set(
            i
            for i, uniprot in enumerate(uniprots)
            if method(uniprot, term)
        )
    
    def select_by_term(self, uniprots, term):
        """
        Accepts a list of UniProt IDs and one or more gene ontology terms
        and returns the UniProts which are annotated with any of the terms.
        
        :param str,set term:
            A single GO term or set of terms.
        """
        
        return set(
            numpy.array(uniprots)[
                list(self.i_select_by_term(uniprots, term))
            ]
        )


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
