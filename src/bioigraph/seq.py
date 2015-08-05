#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
# this module provides a class to look up or match
# residues and regions in sequences of 
# proteins, by default in the canonical
# sequence, and optionally in other isoforms.
#

class Seq(object):
    
    def __init__(self, protein, sequence, isoform = 1):
        self.isof = {}
        self.protein = protein
        self.canonical = isoform
        self.add_seq(sequence, isoform)
    
    def add_seq(self, sequence, isoform):
        self.isof[isoform] = sequence
    
    def match(self, pattern, start, end = None, isoform = None):
        instance = self.get(start, end, isoform)
        pattern = pattern.upper()
        if instance == pattern:
            return True
        else:
            return False
    
    def get(self, start, end = None, isoform = None):
        isoform = isoform if isoform is not None else self.canonical
        end = end if end is not None else start
        return None if len(self.isof[isoform]) < max(start, end) or min(start, end) < 1 \
            else self.isof[isoform][start-1:end]
    
    def isoforms(self):
        return self.isof.keys()
    
    def has_isoform(self, isoform):
        return isoform in self.isof
    
    def get_region(self, residue, start = None, end = None, flanking = 7, isoform = None):
        isoform = self.canonical if isoform is None else isoform
        if residue > len(self.isof[isoform]):
            return None
        start = start if start is not None else residue - flanking
        end = end if end is not None else residue + flanking
        start = max(start, 1)
        end = min(end, len(self.isof[isoform]))
        return (start, end, self.isof[isoform][start-1:end])