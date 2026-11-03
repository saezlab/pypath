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
import abc
import functools

import pypath.share.common as common
import pypath.inputs.disgenet._records as records


class Process:
    _nn = common.not_none
    _int = common.to_int
    _float = common.to_float
    _str = _nn(lambda s: s.strip())
    _tuple = _nn(lambda s, sep: tuple(i.strip() for i in s.split(sep)))
    _passthru = lambda x: x

    @staticmethod
    def _call(
        fun: Callable,
        **kwargs,
    ) -> Callable:

        if kwargs:

            fun = functools.partial(fun, **kwargs)

        def proc(*value: Any, default: Any = None, **kwargs) -> Any:

            try:

                return fun(*value, **kwargs)

            except:

                return default

        return proc


class FieldMeta(abc.ABCMeta):


    def __new__(cls, name, bases, attrs):

        typ = attrs['_type']
        proc = (
            Process._call(typ)
                if callable(typ) else
            getattr(Process, f'_{typ}', Process._passthru)
        )
        attrs['_proc'] = proc

        return super().__new__(cls, name, bases, attrs)


class FieldBase(metaclass = FieldMeta):

    _type = 'passthru'

    def __init__(self, key: str | tuple[str]):

        self._key = common.to_tuple(key)


    def process(self, instance: dict) -> Any:

        values = [instance.get(k, None) for k in self.key]

        return self._proc(**dict(zip(self.key, values)))


class Int(FieldBase):
    _type = int


class Float(FieldBase):
    _type = float


class Str(FieldBase):
    _type = str


class Tuple(FieldBase):
    _type = tuple


class NamedTuple(FieldBase):

    @classmethod
    def _type(cls, **kwargs):

        return cls._record(*(kwargs[k] for k in cls.key))


class NamedTupleSeq(FieldBase):

    @classmethod
    def _type(cls, **kwargs) -> Any:

        return tuple(
            cls._record(*args)
            for args in
            zip(*(kwargs[k] for k in cls.key))
        )


class ProteinClass(NamedTuple):
    _type = records.ProteinClass


class DiseaseClasses(NamedTupleSeq):
    _type = records.DiseaseClass
