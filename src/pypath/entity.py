#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Exports tables for webservice.
#
#  Copyright (c) 2014-2019 - EMBL
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


import pypath.common as common
import pypath.session_mod as session_mod


class Entitiy(session_mod.Logger):
    
    
    def __init__(self, identifier, entity_type, id_type):
        
        self.identifier = identifier
        self.entity_type = entity_type
        self.id_type = id_type
    
    
    @staticmethod
    def entity_name_str(entity):
        
        return (
            entity
                if isinstance(entity, common.basestring) else
            str(entitiy)
        )
    
    
    @classmethod
    def igraph_vertex_name(cls, igraph_v):
        
        return cls.entity_name_str(igraph_v['name'])
    
    
    @staticmethod
    def igraph_vertex_label(igraph_v):
        
        return igraph_v['label']
    
    
    @classmethod
    def igraph_vertex_name_label(cls, igraph_v):
        
        return (
            cls.igraph_vertex_name(igraph_v),
            cls.igraph_vertex_label(igraph_v),
        )
