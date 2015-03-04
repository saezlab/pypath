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

import sys
__all__ = ['Progress']

class Progress(object):
    
    def __init__(self,total=None,name="Progress",interval=3000):
        self.name = name
        self.interval = interval
        self.total = total
        self.done = 0
        sys.stdout.write("\r"+" "*90)
        sys.stdout.write("\r\t:: %s: 0.00%%" % (self.name))
        sys.stdout.flush()
    
    def step(self,step=1):
        if self.done % self.interval == 0:
            sys.stdout.write("\r"+" "*90)
            sys.stdout.write(
                "\r\t:: %s: %.2f%%" % (
                    self.name,
                    float(self.done)/float(self.total)*100.0))
            sys.stdout.flush()
        self.done += step
    
    def terminate(self):
        sys.stdout.write("\r"+" "*90)
        sys.stdout.write("\r\t:: %s: 100.0%%" % (self.name))
        sys.stdout.write("\n")
        sys.stdout.flush()