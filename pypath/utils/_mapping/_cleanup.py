#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

import weakref
import datetime

import timeloop

import pypath.share.settings as _settings

"""
Unload expired (non used) mapping tables.
"""


class Cleanup:

    def __init__(self, period: int = 10):

        self._mappers = []
        self._loop = timeloop.Timeloop()
        self._loop.logger.setLevel(9999)
        self._stop_all_jobs()

        @self._loop.job(interval = datetime.timedelta(seconds = self.period))
        def _cleanup():

            self.run()

        self._loop.start(block = False)


    def _stop_all_jobs(self):

        for job in self._loop.jobs:

            if job.is_alive():

                job.stop()
                job.stopped.set()

        self._loop.jobs = []


    def _setup(self):

        self.period = _settings.get(
            'mapper_cleanup_interval',
            cleanup_period
        )


    def register(self, mapper):

        self._mappers.append(weakref.ref(mapper))


    def run(self):

        self._mappers = [
            mapper.remove_expired()
            for ref in self._mappers
            if mapper := ref()
        ]


def register(mapper):

    key = '_manager'

    if key not in globals():

        globals()[key] = Cleanup()

    globals()[key].register(mapper)
