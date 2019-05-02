# -*- coding: utf-8 -*-

"""Convert OmniPath components to PyBEL.

This module also installs a command line interface that can either be acccessed
via: ``python -m pypath.bel`` or directly as ``bio2bel_omnipath``.

This file is part of the `pypath` python module

Copyright
2014-2019
EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University

File author(s): Dénes Türei (turei.denes@gmail.com)
                 Nicolàs Palacio

Distributed under the GPLv3 License.
See accompanying file LICENSE.txt or copy at
     http://www.gnu.org/licenses/gpl-3.0.html

Website: http://pypath.omnipathdb.org/
"""

import itertools as itt
import sys
from typing import Optional, Set, Union

import click

from pypath import progress
import pypath.session_mod as session_mod

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

except ModuleNotFoundError:
    
    _logger._log(
        'Module `pybel` not available. '
        'You won\'t be able to read or write BEL models.'
    )
    pybel = None

import pybel.constants as pc
import pybel.dsl
from bio2bel.manager.bel_manager import BELManagerMixin

__all__ = [
    'Bel',
    'main',
]

Resource = Union[PyPath, PtmAggregator, AbstractComplexResource]


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
            resource: Resource,
            only_sources: Optional[Set[str]] = None,
            init: bool = False,
    ) -> None:
        
        session_mod.Logger.__init__(self, name = 'bel')
        self.bel_graph = pybel.BELGraph()
        self.resource = resource
        self.only_sources = only_sources

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
        import imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def main(self):
        """
        Convert the resource object to list of BEL relationships.
        """
        
        self._log('Building bel graph from the resource provided.')
        
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
        
        return self
    
    
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
                citations = self._references(edge, direction)

                for predicate, citation in itt.product(predicates, citations):
                    
                    self.bel_graph.add_qualified_edge(
                        source,
                        target,
                        relation = predicate,
                        citation = citation,
                        evidence = 'From OmniPath',
                    )

            if not self._has_direction(directions):
                # add an undirected relationship
                # if no direction available

                citations = self._references(edge, 'undirected')
                source = self._protein(directions.nodes[0])
                target = self._protein(directions.nodes[1])

                for citation in citations:
                    self.bel_graph.add_association(
                        source, target,
                        citation=citation,
                        evidence='From OmniPath',
                    )
        
        self.bel_graph.name = 'OmniPath_interactions'
        
        prg.terminate()
        self._log('Building bel graph from PyPath object finished.')

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
    def _protein(identifier, id_type='uniprot') -> pybel.dsl.Protein:
        return pybel.dsl.Protein(namespace=id_type.upper(), name=identifier)

    def resource_to_relationships_enzyme_substrate(self, enz_sub):
        self._log(
            'Building bel graph from `pypath.ptm.PtmAggregator` '
            'object. Sorry this is not implemented yet!'
        )
        raise NotImplementedError

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
