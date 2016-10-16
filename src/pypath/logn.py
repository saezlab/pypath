#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
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

import os
import logging

from pypath import common


class logw(object):
    def __init__(self, session, loglevel='INFO'):
        self.session = session
        self.loglevel = loglevel
        self.levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
        if not os.path.exists('log'):
            os.makedirs('log')
        self.logfile = 'log/' + self.session + '.log'
        self.init_logger()
        self.wd = os.getcwd()
        open(self.logfile, 'a').close()
        if os.name == 'posix':
            if os.path.islink('log/recent.log') or os.path.isfile(
                    'log/recent.log'):
                os.remove('log/recent.log')
            os.symlink(self.logfile, self.wd + '/log/recent.log')

    def msg(self, indent, message, loglevel=None):
        loglevel = 'INFO' if loglevel is None or loglevel not in self.levels else loglevel
        # time = datetime.datetime.today().strftime('%c')
        # time = str(datetime.datetime.today().strftime('[%Y-%m-%d %H:%M:%S]'))
        offset = ''.join('###' for i in range(indent))
        message = ' '.join(['###', loglevel, offset, message])
        # lfile = codecs.open(self.logfile, encoding='utf-8', mode='a')
        # lfile.write(msg)
        # lfile.close()
        l = getattr(self.logger, loglevel.lower())
        l(message)
        if loglevel == 'ERROR':
            common.console(message)

    def init_logger(self):
        logging.basicConfig(
            filename=self.logfile,
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
