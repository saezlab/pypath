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

import importlib as imp
import collections
import itertools

import numpy as np
import pandas as pd

import pypath.share.settings as settings
import pypath.core.annot as annot
import pypath.core.intercell_annot as intercell_annot
import pypath.core.network as network_mod
import pypath.internals.annot_formats as af


IntercellRole = collections.namedtuple(
    'IntercellRole',
    ['source', 'role'],
)


class IntercellAnnotation(annot.CustomAnnotation):


    def __init__(
            self,
            class_definitions = None,
            cellphonedb_categories = None,
            baccin_categories = None,
            **kwargs
        ):

        class_definitions = (
            class_definitions or
            intercell_annot.annot_combined_classes
        )

        locals_ = locals()
        self._resource_categories = dict(
            (
                res,
                locals_['%s_categories' % res] or
                settings.get('intercell_%s_categories' % res)
            )
            for res in ('baccin', 'cellphonedb')
        )

        annot.CustomAnnotation.__init__(
            self,
            class_definitions = class_definitions,
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

        self.setup_intercell_classes()


    def setup_intercell_classes(self):

        self.set_classes()
        self.add_classes_to_df()
        self.collect_classes()


    def pre_build(self):

        self.add_extra_categories()


    def add_extra_categories(self):

        self.add_cellphonedb_categories()
        self.add_baccin_categories()


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

                    cellphonedb_categories.append(
                        af.AnnotDef(
                            name = '%s_cellphonedb' % category,
                            source = 'CellPhoneDB',
                            args = {attr: category},
                        )
                    )

            self._class_definitions_provided += tuple(cellphonedb_categories)


    def add_baccin_categories(self):

        if self._resource_categories['baccin']:

            self.ensure_annotdb()

            baccin_categories = []

            for attr in ('subclass', 'location'):

                baccin = self.annotdb.annots['Baccin2019']

                categories = set(
                    getattr(c, attr)
                    for cc in baccin.annot.values()
                    for c in cc
                )

                for category in categories:

                    if category in {'other', None}:

                        continue

                    names = (
                        (
                            'membrane_ligand_baccin',
                            'secreted_ligand_baccin',
                        )
                            if category == 'both' else
                        ('%s_ligand_baccin' % category,)
                            if attr == 'location' else
                        ('%s_baccin' % category,)
                    )

                    for name in names:

                        baccin_categories.append(
                            af.AnnotDef(
                                name = name,
                                resource = 'Baccin2019',
                                args = {attr: category},
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

            for values in itertools.product(*(
                hpmr.get_values(field)
                for field in fields
            )):

                for i in range(2, len(fields)):

                    this_fields = fields[:i]
                    this_values = values[:i]

                    args = dict(zip(
                        this_fields,
                        this_values
                    ))

                    members = hpmr.select(**args)

                    if not members:

                        continue

                    name = '_'.join(
                        resep.sub('_', val).lower()
                        for val in reversed(this_values[1:])
                    )
                    hpmr_categories.append(
                        af.AnnotDef(
                            name = name,
                            resource = 'HPMR',
                            args = args,
                            parent = values[0].lower(),
                        )
                    )

            self._class_definitions_provided += tuple(hpmr_categories)


    def add_surfaceome_categories(self):

        resep = re.compile(r'[- /\(\),\.]+')

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

                subclasses = set(itertools.chain(*(
                    a.subclasses
                    for aa in surfaceome.annot.values()
                    for a in aa
                    if (
                        aa.subclasses is not None and
                        aa.mainclass == mainclass and
                        not a[0].isdigit()
                )))

                for subclass in subclasses:

                    name = '%s_%s' % (
                        resep.sub('_', subclass).lower(),
                        mainclass.lower(),
                    )

                    surfaceome_categories.append(
                            af.AnnotDef(
                                name = name,
                                resource = 'Surfaceome',
                                args = {
                                    'mainclass': mainclass,
                                    'subclasses': subclass,
                                },
                                parent = parent,
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
