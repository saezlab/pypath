#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Deals with intercellular communication. Provides functionality for
#  custom set of annotations needed for intercellular communication translation
#  (based on GO). Thus helps to make a meta-database.
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

import re
import importlib as imp
import collections
import itertools

import numpy as np
import pandas as pd

import pypath.share.settings as settings
import pypath.share.common as common
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

        imp.reload(annot)
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


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

        self.setup_intercell_classes()


    def load_from_pickle(self, pickle_file):

        annot.CustomAnnotation.load_from_pickle(
            self,
            pickle_file = pickle_file,
        )

        self.df_add_locations()


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


    def add_cellphonedb_categories(self):

        if self._resource_categories['cellphonedb']:

            self.ensure_annotdb()

            cellphonedb_categories = []

            for mainclass in ('receptor', 'secreted'):

                cpdb = self.annotdb.annots['CellPhoneDB']
                attr = '%s_class' % mainclass

                categories = set(
                    itertools.chain(*(
                        getattr(c, attr)
                        for cc in cpdb.annot.values()
                        for c in cc
                    ))
                )

                for category in categories:

                    if category in {'secreted', 'receptor'}:

                        continue

                    parent = (
                        'receptor'
                            if category == 'receptor' else
                        'ligand'
                    )

                    cellphonedb_categories.append(
                        af.AnnotDef(
                            name = category,
                            parent = parent,
                            resource = 'CellPhoneDB',
                            args = {attr: category},
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


    def post_load(self):

        self.make_df()


    def __repr__(self):

        return (
            '<Intercell annotations: %s records about %s entities>' % (
                self.numof_records(),
                self.numof_entities(),
            )
        )


def init_db(**kwargs):

    globals()['db'] = IntercellAnnotation(**kwargs)


def get_db(**kwargs):

    if 'db' not in globals():

        init_db(**kwargs)

    return globals()['db']
