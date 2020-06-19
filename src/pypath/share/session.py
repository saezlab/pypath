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

import random
import sys
import traceback
import itertools

import pypath.share.log as log


class Session(object):


    def __init__(self, label=None, log_verbosity=0):

        self.label = label or self.gen_session_id()
        self.log_verbosity = log_verbosity
        self.start_logger()
        self.log.msg('Session `%s` started.' % self.label)


    @staticmethod
    def gen_session_id(length=5):
        """
        Returns a 5 alphanumeric characters random identifier.
        """

        abc = '0123456789abcdefghijklmnopqrstuvwxyz'
        return ''.join(random.choice(abc) for i in range(length))


    def start_logger(self):
        """
        Creates a logger for this session which will be served to
        all modules.
        """

        self.logfile = 'pypath-%s.log' % self.label
        self.log = log.Logger(self.logfile, verbosity=self.log_verbosity)


    def finish_logger(self):

        self.log.close_logfile()
        self.log.msg('Session `%s` finished.' % self.label)


    def __repr__(self):

        return '<Session %s>' % self.label


    def __del__(self):

        if hasattr(self, 'log'):

            self.log.msg('Session `%s` finished.' % self.label)


class Logger(object):

    def __init__(self, name = None):

        self._log_name = name or self.__class__.__name__
        self._logger = get_log()

    def _log(self, msg = '', level = 0):
        """
        Writes a message into the logfile.
        """

        self._logger.msg(msg = msg, label = self._log_name, level = level)

    def _console(self, msg = ''):
        """
        Writes a message to the console and also to the logfile.
        """

        self._logger.console(msg = msg, label = self._log_name)


    def _log_traceback(self):
        """
        Includes a traceback into the log.
        """

        exc_type, exc_value, exc_traceback = sys.exc_info()

        if exc_type is not None:

            f = exc_traceback.tb_frame.f_back
            stack = traceback.extract_stack(f)

        else:

            stack = traceback.extract_stack()[:-1]

        trc = 'Traceback (most recent call last):\n'
        trc_list = list(
            itertools.chain(*(
                stack_level.strip('\n').split('\n')
                for stack_level in traceback.format_list(stack)
            ))
        )

        if exc_type is not None:

            trc_list.extend(
                ('  %s' % traceback.format_exc().lstrip(trc)).split('\n')
            )

        stack_top = 0

        for i, line in enumerate(trc_list):

            if line.strip().endswith('<module>'):

                stack_top = i

        trc_list = trc_list[stack_top:]

        self._log(trc.strip())

        for traceline in trc_list:

            self._log(traceline)


def get_session():
    """
    Creates new session or returns the one already created.
    """

    mod = sys.modules[__name__]

    if not hasattr(mod, 'session'):

        new_session()

    return sys.modules[__name__].session


def get_log():
    """
    Returns the ``log.Logger`` instance belonging to the session.
    """

    return get_session().log


def new_session(label = None, log_verbosity = 0):
    """
    Creates a new session. In case one already exists it will be deleted.

    Parameters
    ----------
    label : str
        A custom name for the session.
    log_verbosity : int
        Verbosity level passed to the logger.
    """

    mod = sys.modules[__name__]

    if hasattr(mod, 'session'):

        delattr(mod, 'session')

    setattr(mod, 'session', Session(label, log_verbosity))
