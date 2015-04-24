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
    
    def __init__(self, total = None, name = "Progress", interval = 3000, percent = True):
        self.name = name
        self.interval = interval
        self.total = total
        self.done = 0
        self.percent = percent
        sys.stdout.write("\r"+" "*90)
        if self.percent:
            sys.stdout.write("\r\t:: %s: 0.00%%" % (self.name))
        else:
            self.total = int(self.total)
            sys.stdout.write("\r\t:: %s: 0/%u" % (self.name, self.total))
        sys.stdout.flush()
    
    def step(self, step = 1, msg = None):
        self.done += step
        if self.done % self.interval == 0:
            sys.stdout.write("\r"+" "*90)
            if self.percent:
                sys.stdout.write(
                    "\r\t:: %s: %.2f%% %s" % (
                        self.name,
                        float(self.done)/float(self.total)*100.0,
                        '' if msg is None else '[%s]'%msg))
            else:
                sys.stdout.write(
                    "\r\t:: %s: %u/%u %s" % (
                        self.name, self.done, self.total, 
                        '' if msg is None else '[%s]'%msg))
            sys.stdout.flush()
    
    def terminate(self):
        sys.stdout.write("\r"+" "*90)
        if self.percent:
            sys.stdout.write("\r\t:: %s: 100.0%%" % (self.name))
        else:
            sys.stdout.write("\r\t:: %s: %u/%u" % (self.name, self.total, self.total))
        sys.stdout.write("\n")
        sys.stdout.flush()