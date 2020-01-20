#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
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

import os
import sys
import textwrap
import time
import datetime

import timeloop
# we use this for simple little tasks only
# and don't want engage another logger
timeloop.app.logging.disable(level=9999)

import pypath.share.settings as settings


__all__ = ['new_logger', 'Logger']


_log_flush_timeloop = timeloop.Timeloop()


def new_logger(name = None, logdir = None, verbosity = None, **kwargs):
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
        **kwargs
    )


class Logger(object):

    strftime = time.strftime

    def __init__(
            self,
            fname,
            verbosity = None,
            console_level = None,
            logdir = None,
            max_width = 200,
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

        @_log_flush_timeloop.job(
            interval = datetime.timedelta(
                seconds = settings.get('log_flush_interval')
            )
        )
        def _flush():

            self.flush()

        _log_flush_timeloop.start(block = False)

        self.wrapper = textwrap.TextWrapper(
            width = max_width,
            subsequent_indent = ' ' * 22,
            break_long_words = False,
        )
        self.logdir = self.get_logdir(logdir)
        self.fname  = os.path.join(self.logdir, fname)
        self.verbosity = (
            verbosity
                if verbosity is not None else
            settings.get('log_verbosity')
        )
        self.console_level = (
            console_level
                if console_level is not None else
            settings.get('console_verbosity')
        )
        self.open_logfile()

        # sending some greetings
        self.msg('Welcome!')
        self.msg('Logger started, logging into `%s`.' % self.fname)

    def msg(self, msg='', label=None, level=0, wrap=True):
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
            msg = self.label_message(msg, label=label)
            msg = self.wrapper.fill(msg) if wrap else msg
            msg = self.timestamp_message(msg)
            self.fp.write(msg.encode('utf8', errors='replace'))

        if level <= self.console_level:
            self._console(msg)

    def label_message(self, msg, label = None):
        """
        Adds a label in front of the message.
        """

        label = '[%s] ' % label if label else ''

        return '%s%s' % (label, msg)


    def timestamp_message(self, msg):
        """
        Adds a timestamp in front of the message.
        """

        return '[%s] %s\n' % (self.timestamp(), msg)


    def _console(self, msg):

        sys.stdout.write(msg)
        sys.stdout.flush()


    def console(self, msg = '', label = None):
        """
        Prints a message to the console and also to the logfile.

        msg : str
            Text of the message.
        """

        self.msg(msg = msg, label = label, level = self.console_level)


    @classmethod
    def timestamp(cls):
        """
        Returns a timestamp of the current time.
        """

        return cls.strftime('%Y-%m-%d %H:%M:%S')


    def __del__(self):

        if hasattr(_log_flush_timeloop, 'stop'):

            _log_flush_timeloop.stop()

        self.msg('Logger shut down, logfile `%s` closed.' % self.fname)
        self.msg('Bye.')
        self.close_logfile()


    def get_logdir(self, dirname = None):
        """
        Returns the path to log directory.
        Also creates the directory if does not exist.
        """

        dirname = dirname or '%s_log' % settings.get('module_name')

        if not os.path.exists(dirname):
            os.makedirs(dirname)

        return os.path.abspath(dirname)


    def open_logfile(self):
        """
        Opens the log file.
        """

        self.close_logfile()

        self.fp = open(self.fname, 'wb')


    def close_logfile(self):
        """
        Closes the log file.
        """

        if hasattr(self, 'fp') and not self.fp.closed:

            self.fp.close()


    def flush(self):
        """
        Flushes the log file.
        """

        if hasattr(self, 'fp') and not self.fp.closed:

            self.fp.flush()
