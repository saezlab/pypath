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

import pypath.session_mod as session_mod
import pypath.annot as annot
import pypath.intercell_annot as intercell_annot


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


class Intercell(session_mod.Logger):
    
    
    def __init__(self, network):
        
        session_mod.Logger.__init__(self, name = 'intercell')
        
        self.network = network
    
    
    def reload(self):
        """
        Reloads the object from the module level.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
