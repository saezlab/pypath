#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""Converts OmniPath components to PyBEL.

This module also provides a command line interface that can either be acccessed
via: ``python -m pypath.bel`` or directly as ``bio2bel_omnipath``.

This file is part of the `pypath` python module'

Copyright
2014-2022
EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University

Authors: Dénes Türei (turei.denes@gmail.com)
                 Nicolàs Palacio
#
                 Olga Ivanova
#           Ahmet Rifaioglu

Distributed under the GPLv3 License.
See accompanying file LICENSE.txt or copy at
     http://www.gnu.org/licenses/gpl-3.0.html

Website: http://pypath.omnipathdb.org/
"""

from future.utils import iteritems

import itertools
import sys
from typing import Optional, Set, Union

import click

import pypath.share.progress as progress
import pypath.share.common as common
import pypath.share.session as session_mod

_logger = session_mod.Logger(name = 'bel')

try:
    import pybel

    if hasattr(pybel, 'ob'):
        
        _logger._log(
            'You have the `openbabel` module installed '
            'instead of `pybel`. '
            'To be able to use `pybel`, create a virtual env and install '
            'it by `pip install pybel`.'
        )

        # unimport openbabel
        del sys.modules['pybel']
        pybel = None

    else:

        import pybel.constants as pc
        import pybel.dsl
        from bio2bel.manager.bel_manager import BELManagerMixin

except ModuleNotFoundError:
    
    _logger._log(
        'Module `pybel` not available. '
        'You won\'t be able to read or write BEL models.'
    )
    pybel = None

__all__ = [
    'Bel',
    'main',
]

## if we use the types here we need to import all these modules
## I am not sure if we want this -- Denes
#Resource = Union[
    #PyPath,
    #PtmAggregator,
    #AbstractComplexResource,
    #Set[Union[PyPath, PtmAggregator, AbstractComplexResource]]
#]


class Bel(BELManagerMixin, session_mod.Logger):
    """
    Converts pypath objects to BEL format.
    
    Parameters
    ----------
    only_sources :
        Process data only from these original resources.
    
    Examples
    --------
    >>> import os
    >>> from pypath import main, data_formats, bel
    >>> pa = main.PyPath()
    >>> pa.init_network(data_formats.pathway)
    >>> be = bel.Bel(resource = pa)
    >>> be.main()
    >>> be.to_bel_json(
    ...     os.path.join(
    ...         os.path.expanduser('~'),
    ...         'Desktop',
    ...         'omnipath.bel.json'
    ...     )
    ... )
    """
    
    def __init__(
            self,
            resource,
            ## note - type removed here:
            # resource: Resource
            only_sources: Optional[Set[str]] = None,
            init: bool = False,
            graph_name = 'OmniPath',
    ) -> None:
        
        session_mod.Logger.__init__(self, name = 'bel')
        self.graph_name = graph_name
        self.reset_bel_graph()
        self.resources = (
            resource
                if isinstance(resource, (list, tuple, set)) else
            (resource,)
        )
        self.resource = None
        self.only_sources = only_sources
        self.set_name()
        
        if init:
            self.main()
    
    
    @property
    def initialized(self) -> bool:
        """
        Check if the BEL graph has been initialized.
        """
        
        return 0 < self.bel_graph.number_of_nodes()
    
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        import importlib as imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def main(self):
        """
        Convert the resource object to list of BEL relationships.
        """
        
        self._log('Building bel graph from the resource(s) provided.')
        
        for resource in self.resources:
            
            self.resource = resource
            _ = self.load_resource()
        
        return self
    
    
    def load_resource(self):
        """
        Calls the appropriate processing method for the current resource
        in order to add its content to the BEL graph.
        """
        
        if hasattr(self.resource, 'graph'):
            
            self.resource_to_relationships_graph(self.resource.graph)
        
        elif hasattr(self.resource, 'enz_sub'):
            
            self.resource_to_relationships_enzyme_substrate(
                self.resource.enz_sub
            )
        
        elif hasattr(self.resource, 'complexes'):
            self.resource_to_relationships_complexes(self.resource.complexes)
            
        # FIXME NetworkResource does not exist...
        # this will work once the new reader and network module will be ready
        elif hasattr(self.resource, 'network'):  # NetworkResource object
            self.resource_to_relationships_network(self.resource.network)
            
        else:
            
            self._log(
                'Unknown resource type, don\'t know how to convert '
                'to bel graph: `%s`.' % type(self.resource)
            )
        
        self.resource = None
        
        return self
    
    
    def add_resource(self, resource):
        """
        Adds a resource to the list of resources and calls the processing
        method to add its content to the BEL graph.
        """
        
        self.resources = tuple(self.resources) + (resource,)
        self.resource = resource
        self.load_resource()
        self.resource = None
    
    
    def __iadd__(self, resource):
        
        self.add_resource(resource)
        
        return self
    
    
    def reset_bel_graph(self):
        """
        Assigns a new, empty ``pybel.BELGraph`` instance to the ``bel_graph``
        attribute.
        """
        
        self.bel_graph = pybel.BELGraph()
    
    
    def set_name(self, name = None):
        
        self.bel_graph.name = name or self.graph_name
    
    
    def resource_to_relationships_graph(
            self,
            graph,
        ) -> None:
        """
        Convert a PyPath igraph object into list of BEL relationships.
        """
        
        self._log('Building bel graph from PyPath object (igraph graph).')
        
        edges = graph.es
        prg = progress.Progress(
            len(edges),
            'Building bel graph from PyPath object (igraph graph).',
            1,
        )
        for edge in edges:
            prg.step()
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
                evid_cits = self._references(edge, direction)
                
                for (
                    predicate, (evid, cits)
                ) in itertools.product(predicates, evid_cits):
                    
                    for cit in cits:
                        
                        self.bel_graph.add_qualified_edge(
                            source,
                            target,
                            relation = predicate,
                            citation = cit,
                            evidence = 'OmniPath',
                        )
                        self.bel_graph.add_qualified_edge(
                            source,
                            target,
                            relation = predicate,
                            citation = cit,
                            evidence = evid,
                        )

            if not self._has_direction(directions):
                # add an undirected relationship
                # if no direction available

                evid_cits = self._references(edge, 'undirected')
                source = self._protein(directions.nodes[0])
                target = self._protein(directions.nodes[1])

                for evid, cits in evid_cits:
                    
                    for cit in cits:
                        
                        self.bel_graph.add_association(
                            source, target,
                            citation = cit,
                            evidence = 'OmniPath',
                        )
                        self.bel_graph.add_association(
                            source, target,
                            citation = cit,
                            evidence = evid,
                        )
        
        prg.terminate()
        self._log('Building bel graph from PyPath object finished.')
    
    
    def _references(self, edge, direction) -> Set[str]:
        by_dir = edge['refs_by_dir']
        refs_this_dir = by_dir[direction] if direction in by_dir else set()
        
        return tuple(
            (
                src,
                tuple(str(ref.pmid) for ref in refs & refs_this_dir)
            )
            for src, refs in iteritems(edge['refs_by_source'])
            if not self.only_sources or src in self.only_sources
        )
    
    
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
    def _protein(
            identifier,
            id_type = 'uniprot',
            variants = None,
        ) -> pybel.dsl.Protein:
        
        return pybel.dsl.Protein(
            namespace = id_type.upper(),
            name = identifier,
            variants = variants,
        )

    def resource_to_relationships_enzyme_substrate(self, enz_sub):
        self._log(
            'Building bel graph from `pypath.ptm.PtmAggregator` object.'
        )
        
        for enz_sub in self.resource:
            
            de = enz_sub.ptm.typ.startswith('de')
            mod_type = enz_sub.ptm.typ[2:] if de else enz_sub.ptm.typ
            
            mod_namespace = (
                None if mod_type in common.pmod_other_to_bel else 'OmniPath'
            )
            bel_mod_type = (
                common.pmod_other_to_bel[mod_type]
                    if mod_namespace is None else
                mod_type
            )
            mod_identifier = None if mod_namespace is None else mod_type
            mod = pybel.dsl.pmod(
                name = bel_mod_type,
                position = enz_sub.ptm.residue.number,
                code = common.aminoa_1_to_3_letter[enz_sub.ptm.residue.name],
                identifier = mod_identifier,
                namespace = mod_namespace,
            )
            
            enzyme = self._protein(enz_sub.domain.protein)
            substrate = self._protein(enz_sub.ptm.protein, variants = mod)
            predicate = pc.DIRECTLY_DECREASES if de else pc.DIRECTLY_INCREASES
            
            citations = (enz_sub.refs - {'', '-'}) or {''}
            
            for evid in enz_sub.sources | {'OmniPath'}:
                
                for cit in citations:
                    
                    self.bel_graph.add_qualified_edge(
                        enzyme,
                        substrate,
                        relation = predicate,
                        citation = cit,
                        evidence = evid,
                    )
        
        self._log('Building bel graph from enzyme-substrate data finished.')

    def resource_to_relationships_complexes(self, complexes):
        
        self._log(
            'Building bel graph from `pypath.complex.ComplexAggregator` '
            'object. Sorry this is not implemented yet!'
        )
        raise NotImplementedError

    def resource_to_relationships_network(self, network):
        
        self._log(
            'Building bel graph from `pypath.network.Network` '
            'object. Sorry this is not implemented yet!'
        )
        raise NotImplementedError

    def export_relationships(self, path: str) -> None:
        """Export relationships into 3 columns table.
        
        path : str
            Filename.
        """
        
        with open(path, 'w') as fp:
            
            print('Subject\tPredicate\tObject', file = fp)
            
            for u, v, d in self.bel_graph.edges(data = True):
                
                print('\t'.join((u.name, d[pc.RELATION], v.name)), file = fp)
    
    
    def to_bel(self) -> pybel.BELGraph:
        """Export the BEL graph."""
        if not self.initialized:
            self.main()
        return self.bel_graph

    def to_bel_path(self, path: str, **kwargs):
        """Export the BEL model as a BEL script.
        
        path : str
            Path to the output file.
        """
        pybel.to_bel_path(self.to_bel(), path, **kwargs)

    def to_bel_json(self, path: str, **kwargs):
        """Export the BEL model as a node-link JSON file."""
        pybel.to_json_path(self.to_bel(), path, **kwargs)

    @classmethod
    def get_cli(cls) -> click.Group:
        """Get the command line interface main group."""

        def get_resource(dataset: str):
            
            if dataset == 'graph':
                
                from pypath import main
                from pypath import data_formats
                pa = main.PyPath()
                pa.init_network(data_formats.pathway)
                return pa
            # TODO add more options
            elif dataset == 'complexes':
                
                from pypath import complexes
                co = complexes.ComplexAggregator()
                return co
                
            elif dataset == 'ptms':
                
                from pypath import ptm
                es = ptm.PtmAggregator()
                return es
                
            elif dataset == 'annotations':
                
                from pypath import annot
                an = ptm.AnnotationTable(keep_annotators = True)
                return an
                
            else:
                
                _logger._log('Unknown resource type: `%s`.' % dataset)

        @click.group()
        @click.option(
            '-r',
            '--resource-name',
            type = click.Choice([
                'graph', 'complexes', 'ptms', 'annotations',
            ]),
            default = 'graph',
        )
        @click.pass_context
        def _main(ctx: click.Context, resource_name: str):
            """Bio2BEL OmniPath CLI."""
            resource = get_resource(resource_name)
            ctx.obj = cls(resource = resource)

        @_main.group()
        def bel():
            """BEL utilities."""

        cls._cli_add_to_bel(bel)
        cls._cli_add_upload_bel(bel)

        return _main


main = Bel.get_cli()

if __name__ == '__main__':
    main()
