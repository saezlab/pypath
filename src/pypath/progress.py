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

from __future__ import print_function

import sys
import time
import os
import imp
import shlex
import struct
import platform
import subprocess

import tqdm

__all__ = ['Progress', 'OldProgress']

class Progress(object):
    
    """
    Before I had my custom progressbar here.
    Now it is a wrapper around the great progressbar `tqdm`.
    Old implementation moved to `OldProgress` class.
    """
    
    def __init__(self, total = None, name = "Progress",
             interval = None, percent = True, status = 'initializing',
             done = 0, init = True, unit = 'it'):
        
        self.name = name
        self.interval = (
            max(int(total / 100), 1) if interval is None else interval
        )
        self.total = total
        self.done = done
        self.status = status
        self.unit = unit
        self.start_time = time.time()
        self.min_update_interval = 0.1
        self.last_printed_value = 0
        
        if init:
            self.init_tqdm()
    
    def reload(self):
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def init_tqdm(self):
        self.tqdm = tqdm.tqdm(total = self.total,
                              desc = '%s: %s' % (self.name, self.status),
                              unit_scale = True,
                              unit = self.unit)
        self.last_updated = time.time()
    
    def step(self, step = 1, msg = None, status = 'busy', force = False):
        """
        Updates the progressbar by the desired number of steps.
        
        :param int step: Number of steps or items.
        """
        self.done += step
        
        if force or (self.done % self.interval < 1.0 and \
            time.time() - self.last_updated > self.min_update_interval):
            
            self.set_status(status)
            
            this_update = max(0, self.done - self.last_printed_value)
            
            if this_update == 0:
                self.tqdm.refresh()
                self.tqdm.fp.flush()
            else:
                self.tqdm.update(int(this_update))
                
            self.last_printed_value = self.done
            self.last_updated = time.time()
    
    def terminate(self, status = 'finished'):
        """
        Terminates the progressbar and destroys the tqdm object.
        """
        self.step(self.total - self.done, force = True, status = status)
        self.tqdm.close()
    
    def set_total(self, total):
        """
        Changes the total value of the progress bar.
        """
        self.total = total
        self.tqdm.total = total
        self.step(0)
    
    def set_done(self, done):
        """
        Sets the position of the progress bar.
        """
        self.done = done
        self.tqdm.n = self.done
        self.tqdm.last_print_n = self.done
        self.step(0)
    
    def set_status(self, status):
        """
        Changes the prefix of the progressbar.
        """
        if status != self.status:
            self.status = status
            self.tqdm.set_description(self.get_desc())
            self.tqdm.refresh()
    
    def get_desc(self):
        """
        Returns a formatted string of the description, consisted of
        the name and the status. The name supposed something constant
        within the life of the progressbar, while the status is there
        to give information about the current stage of the task.
        """
        return '%s%s%s%s' % (' ' * 8,
                             self.name,
                             ' -- ' if len(self.name) else '',
                             self.status)

