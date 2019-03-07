#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Convert OmniPath components to PyBEL."""

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


import imp
import itertools as itt
import sys
from typing import Set

import pybel.constants as pc
from pybel.dsl import Protein

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


class Bel:
    """Converts pypath objects to BEL format.
    
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
    >>> be.export_bel(fname = 'omnipath_pathways.bel')
    """

    def __init__(
            self,
            resource,
            only_sources=None,
    ) -> None:
        self.bel_graph = pybel.BELGraph()
        self.resource = resource
        self.only_sources = only_sources

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def main(self):
        """Convert the resource object to list of BEL relationships."""
        if hasattr(self.resource, 'graph'):  # PyPath object
            self.resource_to_relationships_graph(self.resource.graph)

        elif hasattr(self.resource, 'enz_sub'):  # PtmAggregator object
            self.resource_to_relationships_enzyme_substrate(self.resource.enz_sub)

        elif hasattr(self.resource, 'complexes'):  # ComplexAggregator object
            self.resource_to_relationships_complexes(self.resource.complexes)

        elif hasattr(self.resource, 'network'):  # NetworkResource object
            self.resource_to_relationships_network(self.resource.network)

    def resource_to_relationships_graph(self, graph):
        """Convert a PyPath igraph object into list of BEL relationships."""
        for edge in graph.es:
            directions = edge['dirs']

            for direction in (directions.straight, directions.reverse):
                if not directions.dirs[direction]:
                    # this direction does not exist
                    continue

                dir_sources = directions.get_dir(direction, sources=True)

                if self.only_sources and not dir_sources & self.only_sources:
                    # this direction not provided
                    # in the currently enabled set of sources
                    continue

                predicates = set()

                activation, inhibition = (
                    directions.get_sign(direction, sources=True)
                )

                if self._check_sign(activation):
                    predicates.add(pc.DIRECTLY_INCREASES)

                if self._check_sign(inhibition):
                    predicates.add(pc.DIRECTLY_DECREASES)

                if not predicates:
                    # use `regulates` if sign is unknown
                    predicates.add(pc.REGULATES)

                source = self._protein(direction[0])
                target = self._protein(direction[1])
                citations = self._references(edge, direction)

                for predicate, citation in itt.product(predicates, citations):
                    self.bel_graph.add_qualified_edge(
                        source, target,
                        relation=predicate,
                        citation=citation,
                        evidence='From OmniPath',
                    )

            if not self._has_direction(directions):
                # add an undirected relationship
                # if no direction available

                citations = self._references(edge, 'undirected')
                source = self._protein(directions.nodes[0])
                target = self._protein(directions.nodes[1])

                for citation in citations:
                    self.bel_graph.add_qualified_edge(
                        source, target,
                        relation='association',
                        citation=citation,
                        evidence='From OmniPath',
                    )

    def _references(self, edge, direction) -> Set[str]:
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

        return {str(ref.pmid) for ref in references}

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

        return (
                (
                        directions.sources_straight() |
                        directions.sources_reverse()
                ) & self.only_sources
        )

    @staticmethod
    def _protein(identifier, id_type='uniprot') -> Protein:
        return Protein(namespace=id_type.upper(), name=identifier)

    def resource_to_relationships_enzyme_substrate(self, enz_sub):
        pass

    def resource_to_relationships_complexes(self, complexes):
        raise NotImplementedError

    def resource_to_relationships_network(self, network):
        raise NotImplementedError

    def export_relationships(self, fname):
        """
        Exports relationships into 3 columns table.
        
        fname : str
            Filename.
        """
        with open(fname, 'w') as fp:
            print('Subject\tPredicate\tObject', file=fp)
            for u, v, d in self.bel_graph.edges(data=True):
                print('\t'.join((u.name, d[pc.RELATION, v.name])), file=fp)

    def export_bel(self, fname, **kwargs):
        """
        Exports the BEL model into file.
        
        fname : str
            Filename.
        """
        pybel.to_bel_path(self.bel_graph, fname, **kwargs)
