#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2018
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from past.builtins import xrange

import os
import logging

from pypath import common


class logw(object):
    """
    Session logger object.

    :arg str session:
        Session identifier (a random alphanumeric string). See
        :py:func:`common.gen_session_id()`.
    :arg str loglevel:
        Optional, ``'INFO'`` by default. Specifies the level of the
        logger. Possible levels are: ``'DEBUG'``, ``'INFO'``,
        ``'WARNING'``, ``'ERROR'`` or ``'CRITICAL'``.

    :var str logfile:
        Path to the log file (e.g. 'log/123ab.log')
    :var logging.RootLogger logger:
        Python's built-in :py:class:`logger` object.
    :var str loglevel:
        Level of logging.
    :var str session:
        Session identifier (a random alphanumeric string).
    :var str wd:
        Path of the current working directiory.
    """

    def __init__(self, session, loglevel='INFO'):
        self.session = session
        self.loglevel = loglevel
        self.__levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

        if not os.path.exists('log'):
            os.makedirs('log')

        self.logfile = 'log/' + self.session + '.log'
        self.init_logger()
        self.wd = os.getcwd()
        open(self.logfile, 'a').close() # Just creates the log file.

        if os.name == 'posix':
            recent = 'log/recent.log'

            if os.path.islink(recent) or os.path.isfile(recent):
                os.remove(recent)

            # Creates a symbolic link to logfile (of current session) in
            # log/recent.log under the current working directory
            os.symlink(self.logfile, os.path.join(self.wd, recent))

    def msg(self, indent, message, loglevel=None):
        """
        Prints a message in the log file. If *loglevel* is ``'ERROR'``,
        the message is also printed on the standard output.

        :arg int indent:
            Indentation level for the message (each level is three
            character length).
        :arg str message:
            Message to be added in the log.
        :arg str loglevel:
            Level of the log message.
        """

        loglevel = ('INFO' if not loglevel or loglevel not in self.__levels
                    else loglevel)
        # time = datetime.datetime.today().strftime('%c')
        # time = str(datetime.datetime.today().strftime('[%Y-%m-%d %H:%M:%S]'))
        offset = ''.join('###' for i in xrange(indent))
        message = ' '.join(['###', loglevel, offset, message])
        # lfile = codecs.open(self.logfile, encoding='utf-8', mode='a')
        # lfile.write(msg)
        # lfile.close()
        l = getattr(self.logger, loglevel.lower())
        l(message)

        if loglevel == 'ERROR':
            common.console(message)

    def init_logger(self):
        """
        Initializes the logger object according to the parameters set on
        the instance creation. Includes first line with date, time and
        log file name (session ID).
        """

        logging.basicConfig(filename=self.logfile,
                            format='%(asctime)s %(message)s',
                            datefmt='[%Y-%m-%d %H:%M:%S]',
                            level=getattr(logging, self.loglevel))
        self.logger = logging.getLogger(__name__)
        self.msg(1, "Logger initialized, logging to %s" % self.logfile, 'INFO')

    def __getstate__(self):
        d = dict(self.__dict__)
        del d['logger']
        return d

    def __setstate__(self, d):
        self.__dict__ = d
        self.init_logger()
