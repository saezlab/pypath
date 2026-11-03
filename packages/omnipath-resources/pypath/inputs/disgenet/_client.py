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

from __future__ import annotations

from typing import Callable, Generator, NamedTuple

import pypath.share.session as session
import pypath.inputs.disgenet._proc as _proc
import pypath.inputs.disgenet._request as _request


class DisgenetClient(session.Logger):

    def __init__(
            self,
            query_type: _valid.Querytype | str,
            by: _valid.By | str,
            id_type: _valid.Idtype | str | None,
            identifiers: str | Iterable[str] | None,
            evidences: bool = False,
            query_param: dict[str, Any] | None = None,
            email: str | None = None,
            password: str | None = None,
            fail_on_invalid: bool = True,
            **kwargs
        ):

        if not hasattr(self, '_log_name'):

            super().__init__(name = 'disgenet_input')

        args = locals()
        args.pop('self')
        args.pop('kwargs')

        self._request = _request.DisgenetRequest(**args, **kwargs)
        self._processor = _proc.PROCESSORS.get(query_type)


    def iterraw(self) -> Generator[dict]:
        """
        Retrieve and iterate raw records.
        """

        return self._request.__iter__()


    def __iter__(self) -> Generator[NamedTuple]:

        for rec in self.iterraw():

            yield self.process(rec)


    def process(self, record: dict) -> NamedTuple:
        """
        Process raw record.
        """

        return self._processor.process(record)
