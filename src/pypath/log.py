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

import time
import os
import sys
import textwrap

import pypath.settings as settings


def new_logger(name = None, logdir = None, verbosity = 0):
    """
    Returns a new logger with default settings (can be customized).

    Parameters
    ----------
    name : str
        Custom name for the log.
    logdir : str
        Path to the directoty to store log files.
    verbosity : int
        Verbosity level, lowest is 0. Messages from levels above this
        won't be written to the log..

    Returns
    -------
    ``log.Logger`` instance.
    """
    
    name = name or settings.get('module_name')
    logdir = logdir or '%s_log' % name
    
    return Logger(
        fname = '%s__%s.log' % (
            name,
            Logger.timestamp().replace(' ', '_').replace(':', '.'),
        ),
        verbosity = 0,
        logdir = logdir,
    )


class Logger(object):
    
    
    strftime = time.strftime
    
    def __init__(
            self,
            fname,
            verbosity = 0,
            console_level = -1,
            logdir = None,
            max_flush_interval = 2,
            max_width = 79,
        ):
        """
        fname : str
            Log file name.
        logdir : name
            Path to the directory containing the log files.
        verbosity : int
            Messages at and below this level will be written into the
            logfile. All other messages will be dropped.
        console_level : int
            Messages below this log level will be printed not only into
            logfile but also to the console.
        """
        
        self.last_flush = time.time()
        self.wrapper = textwrap.TextWrapper(
            width = max_width,
            subsequent_indent = ' ' * 22,
            break_long_words = False,
        )
        self.max_flush_interval = max_flush_interval
        self.logdir = self.get_logdir(logdir)
        self.fname  = os.path.join(self.logdir, fname)
        self.verbosity = verbosity
        self.console_level = console_level
        self.open_logfile()
        self.msg('Welcome!')
        self.msg('Logger started, logging into `%s`.' % self.fname)
    
    
    def msg(self, msg = '', level = 0):
        """
        Writes a message into the log file.

        Parameters
        ----------
        msg : str
            Text of the message.
        level : int
            The loglevel. Decides if the message will be written or dropped.
        """
        
        if level <= self.verbosity:
            
            msg = self.wrapper.fill(msg)
            msg = self.timestamp_message(msg)
            self.fp.write(msg)
            self.flush_hook()
        
        if level <= self.console_level:
            
            self.console(msg)
    
    
    def timestamp_message(self, msg):
        
        return '[%s] %s\n' % (self.timestamp(), msg)
    
    
    def console(self, msg = ''):
        """
        Prints a message to the console.
        
        msg : str
            Text of the message.
        """
        
        sys.stdout.write('[%s] %s\n' % (self.timestamp(), msg))
        sys.stdout.flush()
    
    
    @classmethod
    def timestamp(cls):
        """
        Returns a timestamp of the current time.
        """
        
        return cls.strftime('%Y-%m-%d %H:%M:%S')
    
    
    def __del__(self):
        
        self.msg('Logger shut down, logfile `%s` closed.' % self.fname)
        self.msg('Bye.')
        self.close_logfile()
    
    
    def get_logdir(self, dirname = None):
        """
        Returns the path to log directory.
        Also creates the directory if does not exist.
        """
        
        dirname = dirname or '%s_log' % settings.get('module_name')
        
        os.makedirs(dirname, exist_ok = True)
        
        return dirname
    
    
    def open_logfile(self):
        """
        Opens the log file.
        """
        
        self.close_logfile()
        self.fp = open(self.fname, 'w')
    
    
    def close_logfile(self):
        """
        Closes the log file.
        """
        
        if hasattr(self, 'fp') and not self.fp.closed:
            
            self.fp.close()
    
    
    def flush_hook(self):
        
        if time.time() - self.last_flush > self.max_flush_interval:
            
            self.flush()
    
    
    def flush(self):
        """
        Flushes the log file.
        """
        
        if hasattr(self, 'fp') and not self.fp.closed:
            
            self.fp.flush()
            self.last_flush = time.time()
