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

import pypath.share.session as session_mod
import pypath.core.annot as annot
import pypath.core.intercell_annot as intercell_annot
import pypath.core.network as network_mod

IntercellRole = collections.namedtuple(
    'IntercellRole',
    ['source', 'role'],
)


class IntercellAnnotation(annot.CustomAnnotation):

    def __init__(self, class_definitions=None, **kwargs):

        class_definitions = (
                class_definitions or intercell_annot.annot_combined_classes
        )

        annot.CustomAnnotation.__init__(
            self,
            class_definitions=class_definitions,
            **kwargs
        )

        self.make_df()

    def reload(self):
        """
        Reloads the object from the module level.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
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
            pickle_file=pickle_file,
        )

        self.setup_intercell_classes()

    def setup_intercell_classes(self):

        self.set_classes()
        self.add_classes_to_df()
        self.collect_classes()

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
