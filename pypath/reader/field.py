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

"""
Classes for extracting values from fields in a record (row) of raw or
preprocessed data.
"""

from future.utils import iteritems

import pypath.share.common as common
import pypath_common._constants as _const


class Field(object):
    """
    Generic field processor, base class for more specific field processors.
    The default behaviour is to return the value at a specific index from
    each record.

    :param int,list idx:
        One or more index or indices.
    :param tuple,list,dict compact:
        Special compact definitions which make it easier to describe field
        definitions in a concise way.
    :param dict mapping:
        Mapping rules for extracting values from the raw field content.
    :param str,NoneType sep:
        Within field separator (to split the field content).
    :param callable method:
        A method for processing - in case the built in ways are not
        satisfying, for greatest flexibility.
    """

    def __init__(
            self,
            idx = None,
            compact = None,
            mapping = None,
            sep = None,
            method = None,
            value_type = None,
            **kwargs
        ):

        self.param = locals()
        _ = self.param.pop('self')
        self.param.update(kwargs)
        self.setup()


    def setup(self):

        for k, v in iteritems(self.param):

            setattr(self, k, v)

        self._from_compact()


    def _from_compact(self):

        if self.compact:

            if isinstance(self.compact, bool):

                self._fix_value = self.compact

            elif isinstance(self.compact, int):

                self.idx = self.compact

            elif isinstance(self.compact, tuple):

                if len(self.compact) == 2:

                    self.idx, self.sep = self.compact

                if len(self.compact) == 3:

                    self.idx, self.sep, self.mapping = self.compact


    def process(self, record):
        """
        Processes a record and returns the value of the field.

        :param list record:
            One record (row) to be processed.

        :returns:
            The processed value of the field as a ``FieldContent`` object.
            This object either provides the processed value or iterates
            through values if possibly more than one value available from
            the field.
        """

        return FieldContent(self._process(record))


    def _process(self, record):

        if hasattr(self, '_fix_value'):

            return self._fix_value

        value = record[self.idx]

        if isinstance(self.sep, str):

            value = value.split(self.sep)

        value = (
                [self._process_one(val) for val in value]
                    if isinstance(value, list) else
                self._process_one(value)
            )

        return value


    def _process_one(self, value):

        if self.value_type:

            value = value_type[value]

        if isinstance(self.mapping, _const.LIST_LIKE):

            value = value in self.mapping

        if isinstance(self.mapping, dict):

            value = self.mapping[value] if value in self.mapping else None

        return value


class FieldContent(object):
    """
    Provides a unified interface for accessing processed field contents
    either as a single value or as an iterable.
    """


    def __init__(self, content):

        self.content = content


    @property
    def value(self):

        return self.content


    def __iter__(self):

        for value in (
            self.content
                if isinstance(self.content, _const.LIST_LIKE) else
            (self.content,)
        ):

            yield value
