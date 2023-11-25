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

import timeloop

_mapper_cleanup_timeloop = timeloop.Timeloop()
_mapper_cleanup_timeloop.logger.setLevel(9999)

for job in _mapper_cleanup_timeloop.jobs:

if job.is_alive():

job.stop()
job.stopped.set()

_mapper_cleanup_timeloop.jobs = []


@_mapper_cleanup_timeloop.job(
interval = datetime.timedelta(
seconds = cleanup_period
)
)
def _cleanup():

remove_expired()


_mapper_cleanup_timeloop.start(block = False)
