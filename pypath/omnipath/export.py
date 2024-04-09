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

import os
import importlib as imp
import re
import json
import copy
import itertools
from future.utils import iteritems

import pandas as pd

import pypath.share.progress as progress
import pypath.resources.urls as urls
import pypath.resources.data_formats as data_formats
import pypath.share.settings as settings
import pypath.core.entity as entity
import pypath.share.session as session
import pypath.resources.network as netres
import pypath.share.common as common

strip_json = re.compile(r'[\[\]{}\"]')
simple_types = {bool, int, float, type(None)}


class Export(session.Logger):

    default_header_uniquepairs = [
        'UniProt_A',
        'GeneSymbol_A',
        'UniProt_B',
        'GeneSymbol_B',
        'Databases',
        'PubMed_IDs',
        'Undirected',
        'Direction_A-B',
        'Direction_B-A',
        'Stimulatory_A-B',
        'Inhibitory_A-B',
        'Stimulatory_B-A',
        'Inhibitory_B-A',
        'Category',
    ]

    default_dtypes_uniquepairs = {
        'UniProt_A': 'category',
        'GeneSymbol_A': 'category',
        'UniProt_B': 'category',
        'GeneSymbol_B': 'category',
        'Undirected': 'category',
        'Direction_A-B': 'category',
        'Direction_B-A': 'category',
        'Stimulatory_A-B': 'category',
        'Inhibitory_A-B': 'category',
        'Stimulatory_B-A': 'category',
        'Inhibitory_B-A': 'category',
        'Category': 'category',
    }

    default_header_bydirs = [
        'source',
        'target',
        'source_genesymbol',
        'target_genesymbol',
        'is_directed',
        'is_stimulation',
        'is_inhibition',
        'consensus_direction',
        'consensus_stimulation',
        'consensus_inhibition',
        'sources',
        'references',
    ]

    default_dtypes_bydirs = {
        'source': 'category',
        'target': 'category',
        'source_genesymbol': 'category',
        'target_genesymbol': 'category',
        'is_directed': 'int8',
        'is_stimulation': 'int8',
        'is_inhibition': 'int8',
        'consensus_direction': 'int8',
        'consensus_stimulation': 'int8',
        'consensus_inhibition': 'int8',
        'sources': 'category',
        'references': 'category',
        'curation_effort': 'int16',
        'extra_attrs': 'category',
        'evidences': 'category',
        'entity_type_source': 'category',
        'entity_type_target': 'category',
    }

    def __init__(
            self,
            network = None,
            only_sources = None,
            extra_node_attrs = None,
            extra_edge_attrs = None,
            outfile = None,
            default_vertex_attr_processor = None,
            default_edge_attr_processor = None,
            pa = None,
        ):

        session.Logger.__init__(self, name = 'export')

        self._log('Export object created for network.')

        self.extra_node_attrs = extra_node_attrs or {}
        self.extra_edge_attrs = extra_edge_attrs or {}
        self.outfile = outfile
        self.network = network or pa
        self.pa      = self.network
        self._set_graph()
        self.only_sources = only_sources

        if isinstance(self.only_sources, list):

            self.only_sources = set(self.only_sources)

        self.default_vertex_attr_processor = (
            default_vertex_attr_processor or
            self.default_vertex_attr_processor
        )
        self.default_edge_attr_processor = (
            default_edge_attr_processor or
            self.default_edge_attr_processor
        )

    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def _set_graph(self):

        self.graph = (
            self.pa._get_undirected()
                if hasattr(self.pa, '_get_undirected') else
            None
        )


    def make_df(
            self,
            unique_pairs = True,
            extra_node_attrs = None,
            extra_edge_attrs = None,
        ):
        """
        Creates a data frame from the network.

        By default UniProt IDs, Gene Symbols, source databases, literature
        references, directionality and sign information and interaction type
        are included.

        Args
        -----
        :param bool unique_pairs:
            If `True` each line corresponds to a unique pair of molecules,
            all directionality and sign information are covered in other
            columns. If `False`, order of `A` and `B` IDs corresponds to
            the direction while sign covered in further columns.
        :param dict extra_node_attrs:
            Additional node attributes to be included in the exported table.
            Keys are column names used in the header while values are names
            of vertex attributes. Values also might be methods which then
            will be called then on each vertex. These should return strings
            or their result will be converted to string.
            In the header `_A` and `_B` suffixes will be appended to the
            column names so the values can be assigned to A and B side
            interaction partners.
        :param dict extra_edge_attrs:
            Additional edge attributes to be included in the exported table.
            Keys are column names used in the header while values are names
            of edge attributes or callables accepting an edge as single
            argument.
        :param str outfile:
            Name of the output file. If `None` a file name
            "netrowk-<session id>.tab" is used.
        """

        self._log('Creating data frame of type `%s`.' % (
            'unique pairs' if unique_pairs else 'by direction'
        ))

        kwargs = locals()
        _ = kwargs.pop('self')

        if self.graph:

            self._make_df_igraph(**kwargs)

        else:

            self._make_df_network(**kwargs)


    def _make_df_network(
            self,
            unique_pairs = True,
            extra_node_attrs = None,
            extra_edge_attrs = None,
        ):
        """
        See docs at method ``make_df``.
        """

        self._log('Creating data frame from `core.network.Network` object.')

        if unique_pairs:

            msg = (
                'Data frame with unique pairs from `core.network.Network` '
                'is not implemented yet, only possible to create it from '
                '`legacy.main.PyPath` object.'
            )

            self._log(msg)
            raise NotImplementedError(msg)

        self.extra_node_attrs = extra_node_attrs or self.extra_node_attrs
        self.extra_edge_attrs = extra_edge_attrs or self.extra_edge_attrs

        header = self.get_header(unique_pairs = unique_pairs)

        dtypes = (
            self.default_dtypes_uniquepairs
                if unique_pairs else
            self.default_dtypes_bydirs
        )
        dtypes = dict(i for i in dtypes.items() if i[0] in header)

        result = []

        for ia in self.network:

            result.extend(self.process_interaction(ia))

        self.df = pd.DataFrame(result, columns = header)
        self.df = self.df.astype(dtypes)


    def process_interaction(self, ia):

        result = []

        consensus = ia.consensus()

        for _dir in ('a_b', 'b_a'):

            nodes = getattr(ia, _dir)
            directed = bool(ia.direction[nodes])
            directed_rev = bool(ia.direction[tuple(reversed(nodes))])

            if (
                (
                    not directed and
                    (_dir == 'b_a' or directed_rev)
                ) or (
                    ia.is_loop() and
                    _dir == 'b_a'
                )
            ):

                continue

            positive = getattr(ia, 'positive_%s' % _dir)()
            negative = getattr(ia, 'negative_%s' % _dir)()

            resources = ';'.join(
                sorted(set(
                    '%s%s' % (
                        res,
                        ('_%s' % via) if via else '',
                    )
                    for res, via in
                    itertools.chain(
                        ia.get_resource_names_via(
                            direction = 'undirected',
                            via = None,
                        ),
                        ia.get_resource_names_via(
                            direction = nodes,
                            via = None,
                        ),
                    )
                ))
            )

            references = ';'.join(
                sorted(set(
                    '%s:%s' % (
                        ev.resource.via or ev.resource.name,
                        ref.pmid
                    )
                    for ev in itertools.chain(
                        ia.get_evidences(direction = 'undirected'),
                        ia.get_evidences(direction = nodes),
                    )
                    for ref in ev.references
                    if not ev.resource.via
                ))
            )

            this_row = [
                nodes[0].identifier,
                nodes[1].identifier,
                nodes[0].label,
                nodes[1].label,
                int(directed),
                int(positive),
                int(negative),
                self.match_consensus(
                    consensus,
                    nodes,
                ),
                self.match_consensus(
                    consensus,
                    nodes,
                    'positive',
                ),
                self.match_consensus(
                    consensus,
                    nodes,
                    'negative',
                ),
                resources,
                references,
            ]

            this_row = self.add_extra_fields(ia, this_row, nodes)

            result.append(this_row)

        return result


    @staticmethod
    def match_consensus(consensus, nodes, effect = None):

        param = list(nodes) + ['directed']

        if effect:

            param.append(effect)

        return int(
            any(
                co[:len(param)] == param
                for co in consensus
            )
        )


    def _make_df_igraph(
            self,
            unique_pairs = True,
            extra_node_attrs = None,
            extra_edge_attrs = None,
        ):
        """
        See docs at method ``make_df``.
        """

        self._log('Creating data frame from `legacy.main.PyPath` object.')

        result = []

        self.pa.genesymbol_labels()

        self.extra_node_attrs = extra_node_attrs or self.extra_node_attrs
        self.extra_edge_attrs = extra_edge_attrs or self.extra_edge_attrs

        dtypes = (
            self.default_dtypes_uniquepairs
                if unique_pairs else
            self.default_dtypes_bydirs
        )

        header = self.get_header(unique_pairs = unique_pairs)

        prg = progress.Progress(
            total = self.graph.ecount(),
            name = 'Creating table',
            interval = 31
        )

        for e in self.graph.es:

            # adding default fields
            lines = (
                self._process_edge_uniquepairs_igraph(e)
                    if unique_pairs else
                self._process_edge_bydirection_igraph(e)
            )

            result.extend(lines)

            prg.step()

        prg.terminate()

        self.df = pd.DataFrame(result, columns = header)
        self.df = self.df.astype(dtypes)


    def get_header(
            self,
            unique_pairs = True,
        ):
        """
        Creates a data frame header (list of field names) according to the
        data frame type and the extra fields.
        """

        suffix_a = 'A' if unique_pairs else 'source'
        suffix_b = 'B' if unique_pairs else 'target'

        header = copy.copy(
            self.default_header_uniquepairs
            if unique_pairs else
            self.default_header_bydirs
        )
        header += self.extra_edge_attrs.keys()
        header += [
            '%s_%s' % (x, suffix_a)
            for x in self.extra_node_attrs.keys()
        ]
        header += [
            '%s_%s' % (x, suffix_b)
            for x in self.extra_node_attrs.keys()
        ]

        return header


    def _process_edge_uniquepairs_igraph(self, e):
        """
        Returns a table row representing a network edge with covering all
        annotations in a single row: directionality represented by fields
        like `Direction_A-B` and `Direction_B-A` effect sign a similar way.

        Args
        -----
        :param igraph.Edge e:
            An edge from a pypath igraph object.
        """

        if self.only_sources and not e['sources'] & self.only_sources:

            return []

        vertex_a = self.graph.vs[e.source]
        vertex_b = self.graph.vs[e.target]

        name_a, label_a = entity.Entity.igraph_vertex_name_label(vertex_a)
        name_b, label_b = entity.Entity.igraph_vertex_name_label(vertex_b)

        return [
            list(itertools.chain(
                (
                    # uniprots, genesymbols
                        name_a,
                        label_a,
                        name_b,
                        label_b,
                    # sources, references
                        ';'.join(list(e['sources'])),
                        ';'.join(
                            sorted(
                                set(r.pmid for r in e['references']),
                                key = int
                            )
                        ),
                    # directions
                        ';'.join(e['dirs'].get_dir(
                            'undirected', sources=True)),
                        ';'.join(e['dirs'].get_dir(
                            (name_a, name_b), sources=True)),
                        ';'.join(e['dirs'].get_dir(
                            (name_b, name_a), sources=True))
                ),
                (
                    # signs
                    ';'.join(a) for a in
                    itertools.chain(
                        e['dirs'].get_sign((name_a, name_b), sources=True),
                        e['dirs'].get_sign((name_b, name_a), sources=True),
                    )
                ),
                (
                    # category
                    ';'.join(e['type']),
                )
            ))
        ]

    def _process_edge_bydirection_igraph(self, e):
        """
        Returns one or more table rows representing a network edge a way
        that opposite direction connections contained in separate rows.
        Directionality and sign information covered in 3 columns:
        `is_directed`, `is_inhibition` and `is_stimulation`.
        This is the row format used in the webservice.

        Args
        -----
        :param igraph.Edge e:
            An edge from a pypath igraph object.
        """

        lines = []

        consensus_edges = set(map(tuple, e['dirs'].consensus_edges()))
        consensus_dir = set(c[:3] for c in consensus_edges)

        for d in ['straight', 'reverse']:

            uniprots = getattr(e['dirs'], d)

            if e['dirs'].dirs[uniprots]:

                is_stimulation = int(e['dirs'].is_stimulation(uniprots))
                is_inhibition = int(e['dirs'].is_inhibition(uniprots))

                this_edge = [
                    entity.Entity.entity_name_str(uniprots[0]),
                    entity.Entity.entity_name_str(uniprots[1]),
                    self.pa.name_to_label(uniprots[0]),
                    self.pa.name_to_label(uniprots[1]),
                    1, # is_directed
                    is_stimulation,
                    is_inhibition,
                    int(
                        (uniprots[0], uniprots[1], 'directed')
                        in consensus_dir
                    ),
                    int(
                        (uniprots[0], uniprots[1], 'directed', 'positive')
                        in consensus_edges
                    ),
                    int(
                        (uniprots[0], uniprots[1], 'directed', 'negative')
                        in consensus_edges
                    ),
                ]

                dsources = (
                    e['dirs'].get_dir(uniprots, sources=True) |
                    e['dirs'].get_dir('undirected', sources=True)
                )

                if self.only_sources and not dsources & self.only_sources:

                    continue

                this_edge.extend([
                    ';'.join(sorted(dsources)),
                    ';'.join(
                        r.pmid
                        for r in itertools.chain(*(
                            rs for s, rs in iteritems(e['refs_by_source'])
                            if s in dsources
                        ))
                    )
                ])

                this_edge = self.add_extra_fields(e, this_edge, uniprots)

                lines.append(this_edge)

        if not e['dirs'].is_directed():

            this_edge = [
                entity.Entity.entity_name_str(e['dirs'].nodes[0]),
                entity.Entity.entity_name_str(e['dirs'].nodes[1]),
                self.pa.name_to_label(e['dirs'].nodes[0]),
                self.pa.name_to_label(e['dirs'].nodes[1]),
                0,
                0,
                0,
                0,
                0,
                0,
                ';'.join(sorted(e['sources'])),
                ';'.join([r.pmid for r in e['references']]),
                self._dip_urls(e),
            ]

            this_edge = (
                self.add_extra_fields_igraph(e, this_edge, 'undirected')
            )

            lines.append(this_edge)

        return lines


    def add_extra_fields(self, e, line, dr = None):
        """
        Takes one table row and using the `igraph.Edge` object and the
        direction provided adds the extra node and edge attribute fields
        as they are defined in `extra_node_attrs` and `extra_edge_attrs`.

        Returns the row with extra fields added.

        Args
        -----
        :param igraph.Edge e:
            One edge.
        :param list line:
            A table row.
        :param tuple,str dr:
            Direction key. A tuple of names (most often UniProt IDs) or
            `undirected`.
        """

        # extra edge attributes on demand of the user
        for k, v in iteritems(self.extra_edge_attrs):

            line.append(
                self.generic_attr_processor(v, e, dr)
                    if hasattr(v, '__call__') else
                self.default_edge_attr_processor(
                    e[v]
                        if (
                            self.graph and
                            v in self.graph.es.attributes()
                        ) else
                    getattr(e, v)
                        if hasattr(e, v) else
                    e.attrs[v]
                        if v in e.attrs else
                    None
                )
            )

        # extra vertex attributes

        nodes = (
            (self.graph.vs[e.source], self.graph.vs[e.target])
                if self.graph else
            dr
        )

        for vertex in nodes:

            for k, v in iteritems(self.extra_node_attrs):

                line.append(
                    self.generic_attr_processor(v, vertex, dr)
                    if hasattr(v, '__call__') else
                    self.default_vertex_attr_processor(
                        vertex[v]
                            if (
                                self.graph and
                                v in self.graph.vs.attributes()
                            ) else
                        getattr(vertex, v)
                            if hasattr(vertex, v) else
                        vertex.attrs[v]
                            if v in vertex.attrs else
                        None
                    )
                )

        return line


    @staticmethod
    def default_vertex_attr_processor(vattr):

        return (
            vattr if type(vattr) in simple_types else
            ';'.join([
                x.strip()
                for x in strip_json.sub('',
                    json.dumps(
                        common.sets_to_sorted_lists(vattr)
                    )).split(',')
            ])
        )

    @staticmethod
    def default_edge_attr_processor(eattr):

        return (
            eattr if type(eattr) in simple_types else
            ';'.join([
                x.strip()
                for x in
                strip_json.sub(
                    '',
                    json.dumps(
                        common.sets_to_sorted_lists(eattr)
                    )
                ).split(',')
            ])
        )

    @staticmethod
    def generic_attr_processor(proc, obj, dr = None):
        """
        Wraps the attribute processor to handle unknown number of arguments.

        Not knowing if the attribute processor expects one or two arguments,
        have no better way than try: if calling with 2 arguments fails with
        `TypeError` we call with one argument.
        """

        dr = (
            'undirected'
                if (
                    dr and
                    hasattr(obj, 'is_directed') and
                    not obj.is_directed()
                ) else
            dr
        )

        try:

            return proc(obj, dr)

        except TypeError as e:

            try:

                return proc(obj)

            except TypeError:

                raise e


    def write_tab(self, outfile = None, **kwargs):
        """
        Writes the data frame into a tab separated file.

        Args
        -----
        :param **kwargs:
            Forwarded to `make_df()`.
        """

        if not hasattr(self, 'df'):

            self.make_df(**kwargs)

        self.write(outfile = outfile)


    def write(self, outfile = None):

        outfile = outfile or self.outfile or os.path.join(
            self.pa.outdir, 'network-%s.tab' % self.pa.session
        )

        self.df.to_csv(outfile, sep = '\t', index = False)


    def webservice_interactions_df(self):

        datasets = (
            'omnipath',
            'kinaseextra',
            'ligrecextra',
            'pathwayextra',
            'mirnatarget',
            'dorothea',
            'collectri',
            'tf_target',
            'lncrna_mrna',
            'tf_mirna',
            'small_molecule',
        )

        def get_dataset_callback(dataset: str) -> callable:

            def has_dataset(e, d) -> bool:

                return e.has_dataset(dataset, direction = d)

            return has_dataset


        dataset_args = {
            dataset: get_dataset_callback(dataset)
            for dataset in datasets
        }

        self.make_df(
            unique_pairs = False,
            extra_node_attrs = {
                'ncbi_tax_id': 'taxon',
                'entity_type': 'entity_type',
            },
            extra_edge_attrs = {
                **dataset_args,
                **{
                    'dorothea_curated': lambda e, d: (
                        e._get_attr('DoRothEA', 'curated', d)
                    ),
                    'dorothea_chipseq': lambda e, d: (
                        e._get_attr('DoRothEA', 'chipseq', d)
                    ),
                    'dorothea_tfbs': lambda e, d: (
                        e._get_attr('DoRothEA', 'tfbs', d)
                    ),
                    'dorothea_coexp': lambda e, d: (
                        e._get_attr('DoRothEA', 'coexp', d)
                    ),
                    'dorothea_level': lambda e, d: (
                        ';'.join(e.dorothea_levels(d))
                    ),
                    'type': lambda e, d: (
                        list(e.get_interaction_types(direction = d))[0]
                    ),
                    'curation_effort': lambda e, d: (
                        e.count_curation_effort(direction = d) + (
                            e.count_curation_effort(direction = 'undirected')
                                if isinstance(d, tuple) else
                            0
                        )
                    ),
                    'extra_attrs': lambda e, d: e.serialize_attrs(d),
                    'evidences': lambda e, d: e.serialize_evidences(d),
                },
            },
        )


    def webservice_interactions_df_legacy(self):

        sources_omnipath = set(
            f.name for f in data_formats.omnipath.values()
        )
        sources_extra_directions = settings.get('network_extra_directions')
        sources_kinase_extra = set(
            f.name for f in data_formats.ptm_misc.values()
        )
        sources_ligrec_extra = set(
            f.name for f in data_formats.ligand_receptor.values()
        )
        sources_pathway_extra = set(
            f.name for f in data_formats.pathway_noref.values()
        )
        sources_mirna = set(
            f.name for f in data_formats.mirna_target.values()
        )

        self.make_df(
            unique_pairs = False,
            extra_node_attrs = {
                'ncbi_tax_id': 'ncbi_tax_id'
            },
            extra_edge_attrs = {
                'omnipath': lambda e, d: (
                    (
                        bool(e['dirs'].sources[d] & sources_omnipath) or
                        (
                            bool(
                                e['dirs'].sources['undirected'] &
                                sources_omnipath
                            ) and
                            bool(
                                e['dirs'].sources[d] &
                                sources_extra_directions
                            )
                        )
                    ) and
                    'PPI' in e['type']
                ),
                'kinaseextra': lambda e, d: (
                    bool(e['dirs'].sources[d] & sources_kinase_extra) and
                    'PPI' in e['type']
                ),
                'ligrecextra': lambda e, d: (
                    bool(e['dirs'].sources[d] & sources_ligrec_extra) and
                    'PPI' in e['type']
                ),
                'pathwayextra': lambda e, d: (
                    bool(e['dirs'].sources[d] & sources_pathway_extra) and
                    'TF' in e['type']
                ),
                'mirnatarget': lambda e, d: (
                    bool(e['dirs'].sources[d] & sources_mirna) and
                    'MTI' in e['type']
                ),
                'dorothea': lambda e, d: (
                    'TF' in e['sources_by_type'] and bool(
                        e['sources_by_type']['TF'] &
                        e['dirs'].sources[d]
                    )
                ),
                'dorothea_curated': 'dorothea_curated',
                'dorothea_chipseq': 'dorothea_chipseq',
                'dorothea_tfbs':    'dorothea_tfbs',
                'dorothea_coexp':   'dorothea_coexp',
                'dorothea_level': lambda e, d: (
                    ';'.join(sorted(e['dorothea_level'])) if
                    'dorothea_level' in e.attributes() and
                    'TF' in e['sources_by_type'] and bool(
                        e['sources_by_type']['TF'] &
                        e['dirs'].sources[d]
                    ) else ''),
                # quite wrong (taking only the first one):
                'type': lambda e: e['type'][0],
                'curation_effort': lambda e, d: (
                    e.count_curation_effort(direction = d) + (
                        e.count_curation_effort(direction = 'undirected')
                            if isinstance(d, tuple) else
                        0
                    )
                ),
            }
        )

    @classmethod
    def sources_table(
            cls,
            pa,
            only_sources = None,
            unique_pairs = True,
            extra_edge_attrs = None,
            extra_node_attrs = None,
            outfile = None,
            default_vertex_attr_processor = None,
            default_edge_attr_processor   = None,
        ):
        """
        Creates a data frame which contains a column for each source database
        with binary values showing presence-absence of interactions across
        resources.
        """

        new = cls(
            pa = pa,
            only_sources = only_sources,
            extra_edge_attrs = extra_edge_attrs or {},
            extra_node_attrs = extra_node_attrs or {},
            outfile = outfile,
            default_vertex_attr_processor = default_vertex_attr_processor,
            default_edge_attr_processor   = default_edge_attr_processor,
        )

        new.make_df(unique_pairs = unique_pairs)

        src_attr = 'Databases' if unique_pairs else 'sources'

        src_all = sorted(only_sources or pa.sources)

        src_cols = dict((src, []) for src in src_all)

        for i, row in new.df.iterrows():

            this_row_src = set(row[src_attr].split(';'))

            for src in src_all:

                src_cols[src].append(int(src in this_row_src))

        for src in src_all:

            new.df.insert(
                loc = new.df.columns.get_loc(src_attr),
                column = src,
                value = src_cols[src]
            )

        return new


    def _dip_urls(self, e):

        attrs = e.attrs if hasattr(e, 'attrs') else e.attributes

        result = []

        if 'dip_id' in attrs:

            dip_ids = sorted(common.to_set(attrs['dip_id']))

            for dip_id in dip_ids:

                try:
                    result.append(
                        urls.urls['dip']['ik'] % (
                            int(dip_id.split('-')[1][:-1])
                        )
                    )
                except:

                    self._log('Could not find DIP ID: %s' % dip_id)

        return ';'.join(result)
