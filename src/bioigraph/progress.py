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
    
    def __init__(self, total = None, name = "Progress", interval = 3000, percent = True, status = 'initializing'):
        self.set_status(status)
        self.name = name
        self.interval = interval
        self.total = total
        self.done = 0
        self.percent = percent
        sys.stdout.write("\r"+" "*90)
        if self.percent:
            sys.stdout.write("\r\t:: %s: %s 0.00%%" % (self.name, '%s,'%self.status))
        else:
            self.total = int(self.total)
            sys.stdout.write("\r\t:: %s: %s 0/%u" % (self.name, '%s,'%self.status, self.total))
        sys.stdout.flush()
    
    def step(self, step = 1, msg = None, status = ''):
        self.set_status(status)
        self.done += step
        if self.done % self.interval == 0:
            sys.stdout.write("\r"+" "*90)
            if self.percent:
                sys.stdout.write(
                    "\r\t:: %s: %s %.2f%% %s" % (
                        self.name, '%s,'%self.status, 
                        float(self.done)/float(self.total)*100.0,
                        '' if msg is None else '[%s]'%msg))
            else:
                sys.stdout.write(
                    "\r\t:: %s: %s %u/%u %s" % (
                        self.name, '%s,'%self.status, self.done, self.total, 
                        '' if msg is None else '[%s]'%msg))
            sys.stdout.flush()
    
    def set_status(self, status):
        self.status = status
    
    def terminate(self, status = 'finished'):
        sys.stdout.write("\r"+" "*90)
        self.set_status(status)
        if self.percent:
            sys.stdout.write("\r\t:: %s: %s 100.0%%" % (self.name, '%s,'%self.status))
        else:
            sys.stdout.write("\r\t:: %s: %s %u/%u" % (self.name, '%s,'%self.status, self.total, self.total))
        sys.stdout.write("\n")
        sys.stdout.flush()

