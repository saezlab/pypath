#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2016 - EMBL-EBI
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
import time

__all__ = ['Progress']

class Progress(object):
    
    def __init__(self, total = None, name = "Progress",
        interval = 3000, percent = True, status = 'initializing'):
        self.status = status
        self.name = name
        self.interval = interval
        self.set_total(total)
        self.done = 0
        self.percent = percent
        self.in_ipython()
        self.last_updated = time.time()
        self.start_time = time.time()
        self.min_update_interval = 0.0 \
            if not self.ipython or self.ipython == 'terminal' \
            else 1.0
        sys.stdout.write("\r"+" "*150)
        if self.percent:
            sys.stdout.write("\r\t:: %s: %s 0.00%%" % \
                (self.name, '%s,'%self.status))
        else:
            self.total = int(self.total)
            sys.stdout.write("\r\t:: %s: %s 0/%u" % \
                (self.name, '%s,'%self.status, self.total))
        sys.stdout.flush()
    
    def step(self, step = 1, msg = None, status = 'working on it', force = False):
        self.status = status
        self.done += step
        if force or (self.done % self.interval < 1.0 and \
            time.time() - self.last_updated > self.min_update_interval):
            sys.stdout.write("\r"+" "*150)
            if self.percent:
                sys.stdout.write(
                    "\r\t:: %s: %s %.2f%% %s" % (
                        self.name,
                        '' if self.status is None else'%s,'%self.status,
                        float(self.done)/float(self.total)*100.0,
                        '' if msg is None else '[%s]'%msg))
            else:
                sys.stdout.write(
                    "\r\t:: %s: %s %u/%u %s" % (
                        self.name, '%s,'%self.status, self.done, self.total,
                        '' if msg is None else '[%s]'%msg))
            sys.stdout.flush()
            self.last_updated = time.time()
    
    def set_total(self, total):
        self.total = total or 9999999999.0
    
    def set_done(self, done):
        self.done = done
    
    def set_status(self, status):
        #if self.ipython == 'notebook':
        #    status = '%s; progress indicator disabled in IPython notebook.' % status
        self.step(step = 0, status = status, force = True)
    
    def terminate(self, status = 'finished'):
        self.finished_time = time.time()
        self.seconds_elapsed =  self.finished_time - self.start_time
        m, s = divmod(self.seconds_elapsed, 60)
        h, m = divmod(m, 60)
        self.time_elapsed = "%d:%02d:%02d" % (h, m, s)
        sys.stdout.write("\r"+" "*150)
        self.status = status
        if self.percent:
            sys.stdout.write("\r\t:: %s: %s 100.0%% [%s elapsed]" % \
                (self.name, '%s,'%self.status, self.time_elapsed))
        else:
            sys.stdout.write("\r\t:: %s: %s %u/%u [%s elapsed]" % (\
                self.name,
                '%s,'%self.status,
                self.total,
                self.total,
                self.time_elapsed)
            )
        sys.stdout.write("\n")
        sys.stdout.flush()
        self.last_updated = time.time()
    
    def in_ipython(self):
        self.ipython = False
        if 'ipykernel' in sys.modules:
            self.ipython = 'notebook'
        elif 'Ipython' in sys.modules:
            self.ipython = 'terminal'
