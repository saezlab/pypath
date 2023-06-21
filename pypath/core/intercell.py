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

from future.utils import iteritems

import re
import importlib as imp
import collections
import itertools

import numpy as np
import pandas as pd

import pypath.share.settings as settings
import pypath.share.common as common
import pypath.share.session as session
import pypath.core.annot as annot
import pypath.core.intercell_annot as intercell_annot
import pypath.core.network as network_mod
import pypath.internals.annot_formats as af


class IntercellAnnotation(annot.CustomAnnotation):


    def __init__(
            self,
            class_definitions = None,
            excludes = None,
            excludes_extra = None,
            cellphonedb_categories = None,
            baccin_categories = None,
            hpmr_categories = None,
            surfaceome_categories = None,
            gpcrdb_categories = None,
            icellnet_categories = None,
            build = True,
            composite_resource_name = None,
            **kwargs
        ):
        """
        Builds a database about roles of proteins and complexes in
        intercellular communication. The built-in category definitions
        defining the default contents of this database can be found in the
        ``pypath.core.intercell_annot`` module.

        :param tuple class_definitions:
            A series of annotation class definitions, each represented by
            an instance of ``pypath.internals.annot_formats.AnnotDef``.
            These definitions carry the attributes and instructions to
            populate the classes.
        :param dict excludes:
            A dict with parent category names (strings) or category keys
            (tuples) as keys and sets if identifiers as values.
            The identifiers in this dict will be excluded from all the
            respective categories while building the database. E.g. if
            the UniProt ID `P00533` (EGFR) is in the set under the key of
            `adhesion` it will be excluded from the category `adhesion` and
            all it's direct children.
        :param dict excludes_extra:
            Same kind of dict as `excludes` but it will be added to the
            built-in default. The built in and the provided extra sets
            will be merged. If you want to overwrite or modify the built-in
            sets provide your custom dict as `excludes`.
        :param bool build:
            Execute the build upon instantiation or set up an empty object
            the build can be executed on later.
        """

        if not hasattr(self, '_log_name'):

            session.Logger.__init__(self, name = 'intercell')

        class_definitions = (
            class_definitions or
            intercell_annot.annot_combined_classes
        )
        excludes = (
            excludes or
            intercell_annot.excludes
        )

        locals_ = locals()
        self._resource_categories = dict(
            (
                res,
                locals_['%s_categories' % res]
                    if locals_['%s_categories' % res] is not None else
                settings.get('intercell_%s_categories' % res)
            )
            for res in (
                'baccin',
                'cellphonedb',
                'hpmr',
                'surfaceome',
                'gpcrdb',
                'icellnet',
            )
        )

        annot.CustomAnnotation.__init__(
            self,
            class_definitions = class_definitions,
            excludes = excludes,
            excludes_extra = excludes_extra,
            build = build,
            composite_resource_name = composite_resource_name,
            **kwargs
        )


    def reload(self):
        """
        Reloads the object from the module level.
        """

        imp.reload(af)
        imp.reload(annot)
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)

        for k, v in iteritems(self.classes):

            k.__class__ = getattr(af, k.__class__.__name__)
            v.__class__ = getattr(af, v.__class__.__name__)


    def set_classes(self):

        self.class_names = set(itertools.chain(
            *intercell_annot.class_types.values()
        ))
        self.class_types = dict(
            (cls, typ)
            for typ, ccls in intercell_annot.class_types.items()
            for cls in ccls
        )

        self.main_classes = {}

        for cls in set(self.classes.keys()):

            mainclass = None

            cls_split = cls.split('_')

            for j in range(len(cls_split)):

                this_part = '_'.join(cls_split[:j])

                if this_part in self.class_names:

                    mainclass = this_part

            self.main_classes[cls] = mainclass


    def add_classes_to_df(self):

        if not hasattr(self, 'df'):

            return

        self.df['mainclass'] = (
            np.array([self.main_classes[c] for c in self.df.category])
        )
        self.df['mainclass'] = self.df['mainclass'].astype('category')
        self.df['class_type'] = (
            np.array([
                (
                    self.class_types[c]
                        if c in self.class_types else
                    'sub'
                )
                for c in self.df.category
            ])
        )
        self.df['class_type'] = self.df['class_type'].astype('category')


    def collect_classes(self):

        self.class_names = set(
            itertools.chain(
                *intercell_annot.class_types.values()
            )
        )

        self.class_types = dict(
            (cls, typ)
            for typ, ccls in intercell_annot.class_types.items()
            for cls in ccls
        )

        self.children = collections.defaultdict(set)
        self.parents = {}
        self.class_labels = {}
        self.resource_labels = {}

        for cls in self.classes.keys():

            mainclass = None

            if cls in intercell_annot.class_types['misc']:

                self.parents[cls] = None

            else:

                cls_split = cls.split('_')

                for j in range(len(cls_split) + 1):

                    this_part = '_'.join(cls_split[:j])

                    if this_part in self.class_names:

                        mainclass = this_part

                self.children[mainclass].add(cls)
                self.parents[cls] = mainclass

                resource = cls_split[-1]

            if mainclass is not None and resource not in mainclass:

                self.resource_labels[cls] = (
                    intercell_annot.get_resource_label(resource)
                )

            self.class_labels[cls] = (
                intercell_annot.get_class_label(mainclass or cls)
            )


    def make_df(self):

        annot.CustomAnnotation.make_df(self)

        self.df_add_causality()
        self.df_add_locations()


    def load_from_pickle(self, pickle_file):

        annot.CustomAnnotation.load_from_pickle(
            self,
            pickle_file = pickle_file,
        )


    def df_add_causality(self):

        if not hasattr(self, 'df'):

            self.make_df()
            return

        for causality in ('transmitter', 'receiver'):

            self.df[causality] = [
                bool(getattr(self.classes[key], causality))
                for key in zip(
                    self.df.category,
                    self.df.parent,
                    self.df.database,
                )
            ]


    def df_add_locations(self, locations = None):

        if not hasattr(self, 'df'):

            self.make_df()
            return

        self._log('Adding location columns to data frame.')

        locations = (
            locations or
            (
                'secreted',
                'plasma_membrane_transmembrane',
                'plasma_membrane_peripheral',
            )
        )
        location_classes = dict(
            (
                location,
                self.select(location),
            )
            for location in locations
        )

        for location, entities in iteritems(location_classes):

            self.df[location] = [
                uniprot in entities
                for uniprot in self.df.uniprot
            ]


    def pre_build(self):

        annot.CustomAnnotation.pre_build(self)
        self.add_extra_categories()


    def add_extra_categories(self):

        self.add_cellphonedb_categories()
        self.add_baccin_categories()
        self.add_hpmr_categories()
        self.add_surfaceome_categories()
        self.add_gpcrdb_categories()
        self.add_icellnet_categories()


    def add_cellphonedb_categories(self):

        if self._resource_categories['cellphonedb']:

            self.ensure_annotdb()

            cellphonedb_categories = []

            for mainclass in ('receptor', 'secreted'):

                cpdb = self.annotdb.annots['CellPhoneDB']
                attr = '%s_class' % mainclass

                categories = cpdb.get_values(attr)

                for category in categories:

                    if category in {'secreted', 'receptor'}:

                        continue

                    parent = (
                        'receptor'
                            if mainclass == 'receptor' else
                        'ligand'
                    )

                    cellphonedb_categories.append(
                        af.AnnotDef(
                            name = category,
                            parent = parent,
                            resource = 'CellPhoneDB',
                            args = {
                                mainclass: bool,
                                attr: category,
                            },
                        )
                    )

            self._class_definitions_provided += tuple(cellphonedb_categories)


    def add_baccin_categories(self):

        if self._resource_categories['baccin']:

            self.ensure_annotdb()
            baccin_categories = []
            baccin = self.annotdb.annots['Baccin2019']
            fields = baccin.get_names()
            locations = {
                'surface': {'membrane', 'both'},
                'secreted': {'secreted', 'both', 'ecm'},
            }

            subclasses = baccin.get_values('subclass') - {'other', None}

            this_fields = fields[1:]

            for subclass in subclasses:

                receptor = 'receptor' in subclass

                args = {'subclass': subclass}

                for location in ('surface', 'secreted'):

                    if receptor and location == 'secreted':

                        continue

                    if not receptor:

                        args['location'] = location

                    members = baccin.select(**args)

                    if not members:

                        continue

                    parent = (
                        'receptor'
                        if receptor else
                        (
                            'cell_surface_ligand'
                                if location == 'surface' else
                            'ligand'
                        )
                    )

                    name = subclass.replace('_receptor', '')

                    baccin_categories.append(
                        af.AnnotDef(
                            name = name,
                            parent = parent,
                            resource = 'Baccin2019',
                            args = args,
                        )
                    )

            self._class_definitions_provided += tuple(baccin_categories)


    def add_hpmr_categories(self):

        resep = re.compile(r'[- /\(\),]+')

        hpmr_categories = []

        if self._resource_categories['hpmr']:

            self.ensure_annotdb()

            hpmr = self.annotdb['HPMR']

            fields = hpmr.get_names()

            for i in range(2, len(fields) + 1):

                combinations = {
                    a[:i]
                    for entity, annots in iteritems(hpmr.annot)
                    for a in annots
                }

                for values in combinations:

                    if not values[0]:

                        continue

                    this_fields = fields[1:i]
                    this_values = values[1:i]

                    if not this_values[-1]:

                        continue

                    args = dict(zip(
                        this_fields,
                        this_values,
                    ))

                    members = hpmr.select(**args)
                    parent = values[0].lower()

                    if not members:

                        continue

                    name_parts = (
                        this_values[1:]
                            if len(this_values) > 1 else
                        this_values
                    )

                    name = '_'.join(
                        name_part.strip('_').replace('_receptors', '')
                        for name_part in
                        (
                            resep.sub('_', val).lower()
                                if val else
                            None
                            for val in reversed(name_parts)
                        )
                        if name_part
                    )

                    hpmr_categories.append(
                        af.AnnotDef(
                            name = name,
                            resource = 'HPMR',
                            args = args,
                            parent = parent,
                        )
                    )

            self._class_definitions_provided += tuple(hpmr_categories)


    def add_gpcrdb_categories(self):

        resep = re.compile(r'[- /\(\),]+')

        gpcrdb_categories = []

        if self._resource_categories['gpcrdb']:

            self.ensure_annotdb()

            gpcrdb = self.annotdb['GPCRdb']

            fields = gpcrdb.get_names()

            for i in range(1, len(fields) + 1):

                combinations = {
                    a[:i]
                    for entity, annots in iteritems(gpcrdb.annot)
                    for a in annots
                }

                for values in combinations:

                    if not values[0]:

                        continue

                    this_fields = fields[:i]
                    this_values = values[:i]

                    args = dict(zip(
                        this_fields,
                        this_values,
                    ))

                    members = gpcrdb.select(**args)

                    if not members:

                        continue

                    name = '_'.join(
                        resep.sub('_', val).lower().strip('_')
                        for val in
                        this_values
                    )
                    name = name.replace('_receptors', '')

                    gpcrdb_categories.append(
                        af.AnnotDef(
                            name = name,
                            resource = 'GPCRdb',
                            args = args,
                            parent = 'receptor',
                        )
                    )

            self._class_definitions_provided += tuple(gpcrdb_categories)


    def add_surfaceome_categories(self):

        resep = re.compile(r'[- /\(\),\.]+')
        recls = re.compile(r'_(?:transporters|receptors|ion_channels)')

        mainclasses = {
            'Receptors': 'receptor',
            'Transporters': 'transporter',
            'Enzymes': 'surface_enzyme',
        }

        if self._resource_categories['surfaceome']:

            self.ensure_annotdb()
            surfaceome = self.annotdb['Surfaceome']
            surfaceome_categories = []

            for mainclass, parent in iteritems(mainclasses):

                subclasses = {
                    sc
                    for annots in surfaceome.annot.values()
                    for a in annots
                    for sc in (a.subclasses or ())
                    if (
                        a.mainclass == mainclass and
                        sc is not None and
                        not sc[0].isdigit()
                    )
                }

                for subclass in subclasses:

                    if subclass.startswith('Other'):

                        continue

                    name = '%s_%s' % (
                        resep.sub('_', subclass).lower().strip('_'),
                        mainclass.lower(),
                    )
                    _parent = (
                        'ion_channel'
                            if 'ion_channel' in name else
                        parent
                    )
                    name = recls.sub('', name)

                    surfaceome_categories.append(
                        af.AnnotDef(
                            name = name,
                            resource = 'Surfaceome',
                            args = {
                                'mainclass': mainclass,
                                'subclasses': subclass,
                            },
                            parent = _parent,
                        )
                    )

            self._class_definitions_provided += tuple(surfaceome_categories)


    def add_icellnet_categories(self):

        icellnet_categories = []

        if self._resource_categories['icellnet']:

            self.ensure_annotdb()

            icellnet = self.annotdb['ICELLNET']

            names = icellnet.get_names()[:3]

            combinations = {
                a[:3]
                for aa in icellnet.annot.values()
                for a in aa
            }

            for values in combinations:

                for l in (2, 3):

                    _fields = names[:l]
                    _values = values[:l]

                    if _values[-1] is None:

                        continue

                    args = dict(zip(_fields, _values))

                    members = icellnet.select(**args)

                    if not members:

                        continue

                    name = '_'.join(
                        val.lower().replace('.', '').replace(' ', '_')
                        for val in _values[1:]
                        if val is not None
                    )

                    icellnet_categories.append(
                        af.AnnotDef(
                            name = name,
                            resource = 'ICELLNET',
                            args = args,
                            parent = values[0],
                        )
                    )

            self._class_definitions_provided += tuple(icellnet_categories)


    def post_load(self):

        self.make_df()


    def __repr__(self):

        return (
            '<Intercell annotations: %s records about %s entities>' % (
                self.numof_records(),
                self.numof_entities(),
            )
        )


    @classmethod
    def filter_df(
            cls,
            annot_df,
            category = None,
            name = None,
            parent = None,
            database = None,
            scope = None,
            aspect = None,
            source = None,
            entities = None,
            entity_type = None,
            causality = None,
            topology = None,
            postfix = None,
        ):

        category = category or name
        args = locals()

        _topologies = {
            'pmtm': 'plasma_membrane_transmembrane',
            'pmp': 'plasma_membrane_peripheral',
            'sec': 'secreted',
        }

        entities = args.pop('entities')
        causality = args.pop('causality') or ()
        topology = args.pop('topology') or ()

        topology = [
            _topologies[top] if top in _topologies else top
            for top in common.to_set(topology)
        ]

        query = cls._process_query_args(
            df = annot_df,
            entities = entities,
            args = args,
            postfix = postfix,
        )

        if causality:

            query.append(cls._process_boolean_group_args(causality, postfix))

        if topology:

            query.append(cls._process_boolean_group_args(topology, postfix))

        args = cls._args_add_postfix(args, postfix)

        query = ' and '.join(query)

        return annot_df.query(query) if query else annot_df


    @staticmethod
    def _process_boolean_group_args(values, postfix):

        if postfix:

            values = {
                '%s%s' % (val, postfix)
                for val in common.to_list(values)
            }

        return ' or '.join(common.to_set(values))


    def network_df(
            self,
            annot_df = None,
            network = None,
            combined_df = None,
            network_args = None,
            annot_args = None,
            annot_args_source = None,
            annot_args_target = None,
            entities = None,
            only_directed = False,
            only_undirected = False,
            undirected_orientation = None,
            only_signed = None,
            only_effect = None,
            only_proteins = False,
            swap_undirected = True,
            entities_or = False,
            transmitter_receiver = False,
            only_generic = True,
            only_composite = True,
            only_functional = True,
            exclude_intracellular = True,
        ):
        """
        Combines the annotation data frame and a network data frame.
        Creates a ``pandas.DataFrame`` where each record is an interaction
        between a pair of molecular enitities labeled by their annotations.

        network : pypath.network.Network,pandas.DataFrame
            A ``pypath.network.Network`` object or a data frame with network
            data.
        combined_df : pandas.DataFrame
            Optional, a network data frame already combined with annotations
            for filtering only.
        resources : set,None
            Use only these network resources.
        entities : set,None
            Limit the network only to these molecular entities.
        entities_source : set,None
            Limit the source side of network connections only to these
            molecular entities.
        entities_target : set,None
            Limit the target side of network connections only to these
            molecular entities.
        annot_args : dict,None
            Parameters for filtering annotation classes; note, the defaults
            might include some filtering, provide an empty dict if you want
            no filtering at all; however this might result in huge data
            frame and consequently memory issues. Passed to the ``filtered``
            method.
        annot_args_source : dict,None
            Same as ``annot_args`` but only for the source side of the
            network connections.
        annot_args_target : dict,None
            Same as ``annot_args`` but only for the target side of the
            network connections.
        only_directed : bool
            Use only the directed interactions.
        only_undirected : bool
            Use only the undirected interactions. Specifically for retrieving
            and counting the interactions without direction information.
        undirected_orientation : str,None
            Ignore the direction at all interactions and make sure all of
            them have a uniform orientation. If `id`, all interactions will
            be oriented by the identifiers of the partenrs; if `category`,
            the interactions will be oriented by the categories of the
            partners.
        only_effect : int,None
            Use only the interactions with this effect. Either -1 or 1.
        only_signed : bool
            Use only the interactions with effect sign.
        only_proteins : bool
            Use only the interactions where each of the partners is a protein
            (i.e. not complex, miRNA, small molecule or other kind of entity).
        transmitter_receiver : bool
            On the source side only transmitters, on the target side only
            receivers.
        only_generic : bool
            Use only the generic classes. If specific classes allowed the
            size of the combined data frame might be huge.
        only_composite : bool
            Use only the composite classes. If resource_specific classes
            allowed the size of the combined data frame might be huge.
        only_functional : bool
            Use only the functional classes. Locational classes are often
            not relevant and they largely increase the size of the
            combined data frame.
        exclude_intracellular : bool
            Remove the intracellular parent class and it's children. These
            classes are not relevant in intercellular signaling and having
            them largely increases the size of the combined data frame.
        """

        annot_df = annot_df or self.get_df()

        if exclude_intracellular:

            if combined_df is None:

                annot_df = annot_df[annot_df.parent != 'intracellular']

            else:

                combined_df = combined_df.query(
                    'parent_a != "intracellular" and '
                    'parent_b != "intracellular"'
                )

        annot_args = annot_args or {}
        annot_args_source = annot_args_source or {}
        annot_args_target = annot_args_target or {}

        if only_generic:

            annot_args['scope'] = 'generic'

        if only_composite:

            annot_args['source'] = 'composite'

        if only_functional:

            annot_args['aspect'] = 'functional'

        if transmitter_receiver:

            annot_args_source['causality'] = 'transmitter'
            annot_args_target['causality'] = 'receiver'

        return annot.CustomAnnotation.network_df(
            self,
            annot_df = annot_df,
            network = network,
            combined_df = combined_df,
            network_args = network_args,
            annot_args = annot_args,
            annot_args_source = annot_args_source,
            annot_args_target = annot_args_target,
            entities = entities,
            only_directed = only_directed,
            only_undirected = only_undirected,
            only_signed = only_signed,
            only_effect = only_effect,
            only_proteins = only_proteins,
            swap_undirected = swap_undirected,
            entities_or = entities_or,
            undirected_orientation = undirected_orientation,
        )


    # this became a synonym
    filter_interclass_network = network_df


    def update_summaries(self):

        self.summaries = {}

        for key, group in iteritems(self.classes):

            if group.source == 'resource_specific':

                continue

            self.summaries[key] = {
                'name': group.name,
                'aspect': group.aspect,
                'transmitter': group.transmitter,
                'receiver': group.receiver,
                'resources': self.resources_in_category(key),
                'n_proteins': group.n_proteins,
                'n_mirnas': group.n_mirnas,
                'n_complexes': group.n_complexes,
            }

        self.summaries[('Total', 'Total', 'OmniPath')] = {
            'name': 'Total',
            'aspect': '',
            'transmitter': '',
            'receiver': '',
            'resources': self.all_resources(),
            'n_proteins': self.numof_proteins(),
            'n_mirnas': self.numof_mirnas(),
            'n_complexes': self.numof_complexes(),
        }


    def summaries_tab(self, outfile = None, return_table = False):

        columns = (
            ('name', 'Category'),
            ('aspect', 'Aspect'),
            ('transmitter', 'Transmitter'),
            ('receiver', 'Receiver'),
            ('n_proteins', 'Proteins'),
            ('n_mirnas', 'miRNAs'),
            ('n_complexes', 'Complexes'),
            ('resources', 'Resources'),
        )

        tab = []
        tab.append([f[1] for f in columns])

        tab.extend([
            [
                (
                    ', '.join(self.summaries[key][f[0]])
                        if isinstance(self.summaries[key][f[0]], list) else
                    str(self.summaries[key][f[0]])
                )
                for f in columns
            ]
            for key in sorted(
                self.summaries.keys(),
                key = lambda k: k[0] if k[0] != 'Total' else 'zzzz',
            )
        ])

        if outfile:

            with open(outfile, 'w') as fp:

                fp.write('\n'.join('\t'.join(row) for row in tab))

        if return_table:

            return tab


def init_db(**kwargs):

    globals()['db'] = IntercellAnnotation(**kwargs)


def get_db(**kwargs):

    if 'db' not in globals():

        init_db(**kwargs)

    return globals()['db']
