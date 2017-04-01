#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2017- EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import pypath.uniprot_input as uniprot_input
import pypath.dataio as dataio


class ReferenceList(object):
    def __init__(self, nameType, typ, tax, inFile, **kwargs):
        self.infile = inFile
        self.nameType = nameType
        self.typ = typ
        self.tax = tax
        self.kwargs = kwargs
        if 'organism' not in self.kwargs:
            self.kwargs['organism'] = self.tax

    def load(self):
        if hasattr(dataio, self.infile):
            toCall = getattr(dataio, self.infile)
            lst = toCall(**self.kwargs)
        else:
            f = codecs.open(self.infile, encoding='utf-8', mode='r')
            lst = []
            for l in f:
                lst.append(l.strip())
            f.close()
        self.lst = set(lst)

    def __contains__(self, something):
        return something in self.lst


def get_reflists():
    return [ReferenceList('uniprot', 'protein', 9606, 'all_uniprots')]
