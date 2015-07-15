#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

# from this package:
import dataio

class GO(object):
    
    def __init__(self, organism = 9606):
        go = dataio.get_go_quick(organism = organism)
        self.c = go['terms']['C']
        self.f = go['terms']['F']
        self.p = go['terms']['P']
        self.name = go['names']
        self.term = dict([(v, k) for k, v in self.name.iteritems()])
    
    def get_name(self, term):
        return None if term not in self.name else self.name[term]
    
    def get_term(self, name):
        return None if name not in self.term else self.term[name]
    
    def get_annot(self, uniprot, aspect):
        dic = getattr(self, aspect.lower())
        return [] if uniprot not in dic else dic[uniprot]
    
    def get_annots(self, uniprot):
        result = {}
        for asp in ['C', 'F', 'P']:
            result[asp.upper()] = self.get_annot(uniprot, asp)
        return result
    