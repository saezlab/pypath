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
import collections


try:
    
    import pybel
    
    if hasattr(pybel, 'ob'):
        
        sys.stdout.write(
            'pypath.bel: You have the `openbabel` module installed '
            'instead of `pybel`.\n'
            'To be able to use `pybel`, create a virtual env and install '
            'it by `pip install pybel`.\n'
        )
        
        # unimport openbabel
        del sys.modules['pybel']
        pybel = None
    
except ModuleNotFoundError:
    
    sys.stdout.write(
        'pypath.bel: module `pybel` not available.\n'
        'You won\'t be able to read or write BEL models.\n'
    )
    
    pybel = None


Relationship = collections.namedtuple(
    'Relationship',
    ('subject', 'predicate', 'object', 'references'),
)


class Bel(object):
    """
    Converts pypath objects to BEL format.
    
    Parameters
    ----------
    resource : object
        Object to be converted.
        E.g. ``pypath.main.PyPath`` or
        ``pypath.ptm.PtmAggregator`` or
        ``pypath.complex.ComplexAggregator`` or
        ``pypath.network.NetworkResource``.
    only_sources : set
        Process data only from these original resources.
    
    Examples
    --------
    >>> from pypath import main, data_formats, bel
    >>> pa = main.PyPath()
    >>> pa.init_network(data_formats.pathway)
    >>> be = bel.Bel(resource = pa)
    >>> be.resource_to_relationships()
    >>> be.relationships_to_bel()
    >>> be.export_bel(fname = 'omnipath_pathways.bel')
    """
    
    def __init__(
            self,
            resource,
            only_sources = None,
        ):
        
        self.relationships = []
        self.bel = None
        self.resource = resource
        self.only_sources = only_sources
    
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        
        self.resource_to_relationships()
        self.relationships_to_bel()
    
    
    def resource_to_relationships(self):
        """
        Converts the resource object to list of BEL relationships.
        """
        
        if hasattr(self.resource, 'graph'):
            # PyPath object
            
            self.resource_to_relationships_graph()
            
        elif hasattr(self.resource, 'enz_sub'):
            # PtmAggregator object
            
            self.resource_to_relationships_enzyme_substrate()
            
        elif hasattr(self.resource, 'complexes'):
            # ComplexAggregator object
            
            self.resource_to_relationships_complex()
            
        elif hasattr(self.resource, 'network'):
            # NetworkResource object
            
            self.resource_to_relationships_network()
    
    
    def resource_to_relationships_graph(self):
        """
        Converts a PyPath igraph object into list of BEL relationships.
        """
        
        for edge in self.resource.graph.es:
            
            directions = edge['dirs']
            
            for direction in (directions.straight, directions.reverse):
                
                if not directions.dirs[direction]:
                    # this direction does not exist
                    
                    continue
                
                dir_sources = directions.get_dir(direction, sources = True)
                
                if self.only_sources and not dir_sources & self.only_sources:
                    # this direction not provided
                    # in the currently enabled set of sources
                    
                    continue
                
                predicates = set()
                
                activation, inhibition = (
                    directions.get_sign(direction, sources = True)
                )
                
                if self._check_sign(activation):
                    
                    predicates.add('directlyIncreases')
                
                if self._check_sign(inhibition):
                    
                    predicates.add('directlyDecreases')
                
                if not predicates:
                    # use `regulates` if sign is unknown
                    
                    predicates.add('regulates')
                
                references = self._references(edge, direction)
                
                for predicate in predicates:
                    
                    rel = Relationship(
                        subject    = self._protein(direction[0]),
                        predicate  = predicate,
                        object     = self._protein(direction[1]),
                        references = references,
                    )
                    
                    self.relationships.append(rel)
        
        if not self._has_direction(directions):
            # add an undirected relationship
            # if no direction available
            
            references = self._references(edge, 'undirected')
            
            rel = Relationship(
                subject    = self._protein(directions.nodes[0]),
                predicate  = 'association',
                object     = self._protein(directions.nodes[1]),
                references = references,
            )
            
            self.relationships.append(rel)
    
    
    def _references(self, edge, direction):
        
        by_dir = edge['refs_by_dir']
        references = by_dir[direction] if direction in by_dir else set()
        
        if self.only_sources:
            
            references = (
                references &
                set.union(
                    *(
                        edge['refs_by_source'][src]
                        for src in (self.only_sources & edge['sources'])
                    )
                )
            )
        
        references = set(ref.pmid for ref in references)
        
        return references
    
    
    def _check_sign(self, this_sign_sources):
        
        return (
            this_sign_sources and
            (
                not self.only_sources or
                this_sign_sources & self.only_sources
            )
        )
    
    
    def _has_direction(self, directions):
        
        if not self.only_sources:
            
            return directions.is_directed()
            
        else:
            
            return (
                (
                    directions.sources_straight() |
                    directions.sources_reverse()
                ) &
                self.only_sources
            )
    
    
    @staticmethod
    def _protein(identifier, id_type = 'uniprot'):
        
        return 'p(%s:%s)' % (id_type.upper(), identifier)
    
    
    def resource_to_relationships_enzyme_substrate(self):
        
        pass
    
    
    def resource_to_relationships_complex(self):
        
        raise NotImplementedError
    
    
    def resource_to_relationships_network(self):
        
        raise NotImplementedError
    
    
    def relationships_to_bel(self):
        """
        Converts the relationships into a ``pybel.BELGraph`` object.
        """
        
        self.bel = pybel.BELGraph()
        bel_parser = pybel.parser.BELParser(self.bel)
        
        for rel in self.relationships:
            
            for pubmed in rel.references or ('0',):
                
                bel_parser.control_parser.clear()
                bel_parser.citation = {
                    pybel.constants.CITATION_REFERENCE: pubmed,
                    pybel.constants.CITATION_TYPE:
                        pybel.constants.CITATION_TYPE_PUBMED,
                }
                #TODO: what to add here?
                bel_parser.control_parser.evidence = 'PubMed:%s' % pubmed
                
                
                bel_string = ' '.join(rel[:3])
                print(bel_string)
                bel_parser.control_parser.parseString(bel_string)
        
        #TODO: add other types of objects to this graph
    
    
    def export_relationships(self, fname):
        """
        Exports relationships into 3 columns table.
        
        fname : str
            Filename.
        """
        
        with open(fname, 'w') as fp:
            
            _ = fp.write('Subject\tPredicate\tObject\n')
            
            _ = fp.write(
                '\n'.join(
                    '\t'.join(rel)
                    for rel in self.relationships
                )
            )
    
    
    def export_bel(self, fname):
        """
        Exports the BEL model into file.
        
        fname : str
            Filename.
        """
        
        with open(fname, 'w') as fp:
            
            pybel.to_bel(self.bel, file = fp)
