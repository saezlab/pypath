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

import imp
import collections
import itertools

import numpy as np
import pandas as pd

import pypath.session_mod as session_mod
import pypath.annot as annot
import pypath.intercell_annot as intercell_annot
import pypath.network as network_mod
import pypath.main as main_mod


IntercellRole = collections.namedtuple(
    'IntercellRole',
    ['source', 'role'],
)


class IntercellAnnotation(annot.CustomAnnotation):


    def __init__(self, class_definitions = None):

        class_definitions = (
            class_definitions or intercell_annot.annot_combined_classes
        )

        annot.CustomAnnotation.__init__(
            self,
            class_definitions = class_definitions,
        )

        self.make_df()
        self.set_classes()
        self.add_classes_to_df()


    def reload(self):
        """
        Reloads the object from the module level.
        """

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

        for cls in set(self.df.category):

            mainclass = None

            cls_split = cls.split('_')

            for j in range(len(cls_split)):

                this_part = '_'.join(cls_split[:j])

                if this_part in self.class_names:

                    mainclass = this_part

            self.main_classes[cls] = mainclass


    def add_classes_to_df(self):

        self.df['mainclass'] = (
            pd.Series(np.array([
                self.main_classes[c] for c in self.df.category
            ])).values
        )
        self.df['class_type'] = (
            pd.Series(np.array([
                (
                    self.class_types[c]
                        if c in self.class_types else
                    'sub'
                )
                for c in self.df.category
            ])).values
        )


class Intercell(IntercellAnnotation):


    def __init__(
            self,
            class_definitions = None,
            network = None,
            network_param = None,
            omnipath = False,
        ):

        session_mod.Logger.__init__(self, name = 'intercell')

        IntercellAnnotation.__init__(
            self,
            class_definitions = class_definitions,
        )

        self.network = network
        self.network_param = network_param or {}
        self.omnipath = omnipath


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        import imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def main(self):

        self.setup_network()


    def setup_network(self):

        if self.network is None:

            if not self.omnipath and 'lst' not in self.network_param:

                self.network_param['lst'] = data_formats.pathway

            self.network = main_mod.PyPath()

            if self.omnipath:

                self.network.load_omnipath(**self.network_param)

            else:

                self.network.init_network(**self.network_param)

        if isinstance(self.network, main_mod.PyPath):
            
            self.network = network_mod.Network.from_igraph(self.network)
        
        if isinstance(self.network, network_mod.Network):

            pass
