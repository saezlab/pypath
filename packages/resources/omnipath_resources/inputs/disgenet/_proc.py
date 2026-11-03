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

from typing import Any

import collections

from pypath.inputs.disgenet._schema import SCHEMA
import pypath.inputs.disgenet._records as _records

PROCESSORS = {}


class ProcessBase:

    entity_type: str | None = None

    def __init__(
            self,
            record: type | None = None,
            entity_idx: int | None = None,
        ):

        self._record = record or getattr(_records, self.__class__.__name__)
        self._entity_idx = entity_idx


    def process(self, record: dict) -> _records.NamedTuple:

        return self._record(**{
            k: self.field(k, record)
            for k in record
        })


    def field(key: str, record: dict) -> Any:

        if key in PROCESSORS:

            return PROCESSORS[key].process(record)

        else:

            return self._field(key, record)


    def _field(key: str, record: dict) -> Any:

        proc, *key_in = SCHEMA.get(key, (lambda x: x,))
        key_in = key_in[0] if key_in else key
        key_in = key_in.format(
            entity_type = self.entity_type,
            i = self.entity_idx,
        )

        return proc(record[key_in])


    @property
    def entity_idx(self) -> int | str:

        return '' if self._entity_idx is None else self._entity_idx


class Disease(ProcessBase):
    entity_type = 'disease'


class Gene(ProcessBase):
    entity_type = 'gene'


class Variant(ProcessBase):
    entity_type = 'variant'


class DiseaseDiseaseAssociation(ProcessBase): pass


class GeneDiseaseAssociation(ProcessBase): pass


class VariantDiseaseAssociation(ProcessBase): pass


PROCESSORS.update({
    'disease': Disease(),
    'disease1': Disease(1),
    'disease2': Disease(2),
    'gene': Gene(),
    'variant': Variant(),
    'dda': DiseaseDiseaseAssociation(),
    'gda': GeneDiseaseAssociation(),
    'vda': VariantDiseaseAssociation(),
})
