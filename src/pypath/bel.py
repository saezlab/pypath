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
from tqdm import tqdm

from pypath.complex import AbstractComplexResource
from pypath.main import PyPath
from pypath.ptm import PtmAggregator

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

import pybel.constants as pc
import pybel.dsl
from bio2bel.manager.bel_manager import BELManagerMixin

__all__ = [
    'Bel',
    'main',
]

Resource = Union[PyPath, PtmAggregator, AbstractComplexResource]


class Bel(BELManagerMixin):
    """Converts pypath objects to BEL format.
    
    Parameters
    ----------
    only_sources :
        Process data only from these original resources.
    
    Examples
    --------
    >>> import os
    >>> from pypath import PyPath, data_formats, bel
    >>> pa = PyPath()
    >>> pa.init_network(data_formats.pathway)
    >>> be = bel.Bel(resource=pa)
    >>> be.main()
    >>> be.to_bel_json(os.path.join(os.path.expanduser('~'), 'Desktop', 'omnipath.bel.json'))
    """

    def __init__(
            self,
            resource: Resource,
            only_sources: Optional[Set[str]] = None,
            init: bool = False,
    ) -> None:
        self.bel_graph = pybel.BELGraph()
        self.resource = resource
        self.only_sources = only_sources

        if init:
            self.main()

    @property
    def initialized(self) -> bool:
        """Check if the BEL graph has been initialized."""
        return 0 < self.bel_graph.number_of_nodes()

    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        import imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

    def main(self):
        """Convert the resource object to list of BEL relationships."""

        if isinstance(self.resource, PyPath):
            self.resource_to_relationships_graph(self.resource.graph)

        elif isinstance(self.resource, PtmAggregator):
            self.resource_to_relationships_enzyme_substrate(self.resource.enz_sub)

        elif isinstance(self.resource, AbstractComplexResource):
            self.resource_to_relationships_complexes(self.resource.complexes)

        # FIXME NetworkResource does not exist...
        elif hasattr(self.resource, 'network'):  # NetworkResource object
            self.resource_to_relationships_network(self.resource.network)

        return self

    def resource_to_relationships_graph(self, graph, use_tqdm: bool = False) -> None:
        """Convert a PyPath igraph object into list of BEL relationships."""
        edges = graph.es
        if use_tqdm:
            edges = tqdm(edges)
        for edge in edges:
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
                    self.bel_graph.add_association(
                        source, target,
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
    def _protein(identifier, id_type='uniprot') -> pybel.dsl.Protein:
        return pybel.dsl.Protein(namespace=id_type.upper(), name=identifier)

    def resource_to_relationships_enzyme_substrate(self, enz_sub):
        pass

    def resource_to_relationships_complexes(self, complexes):
        raise NotImplementedError

    def resource_to_relationships_network(self, network):
        raise NotImplementedError

    def export_relationships(self, path: str) -> None:
        """Export relationships into 3 columns table.
        
        path : str
            Filename.
        """
        with open(path, 'w') as fp:
            print('Subject\tPredicate\tObject', file=fp)
            for u, v, d in self.bel_graph.edges(data=True):
                print('\t'.join((u.name, d[pc.RELATION], v.name)), file=fp)

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

        def get_resource(dataset):
            raise NotImplementedError  # TODO @deeenes

        @click.group()
        @click.option('--dataset')
        @click.pass_context
        def _main(ctx: click.Context, dataset: str):
            """Bio2BEL OmniPath CLI."""
            resource = get_resource(dataset)
            ctx.obj = cls(resource=resource)

        @_main.group()
        def bel():
            """BEL utilities."""

        cls._cli_add_to_bel(bel)
        cls._cli_add_upload_bel(bel)

        return _main


main = Bel.get_cli()

if __name__ == '__main__':
    main()
