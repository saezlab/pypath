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

from future.utils import iteritems

import json
import itertools

import pypath.share.common as common


class AttributeHandler(object):
    """
    Base class for other classes which carry custom attributes (data) in a
    dedicated dict under the `attrs` attribute.
    """

    __slots__ = [
        'attrs',
    ]


    def __init__(self, attrs = None, **kwargs):

        self.attrs = self._add_kwargs(attrs, **kwargs)


    def update_attrs(self, attrs = None, **kwargs):
        """
        Updates the attributes stored here. The attributes with identical
        keys are merged using the :py:func:`pypath.share.common.combine_attrs`
        function.

        The new attributes can be provided three ways: an object with an
        attribute called `attrs`; a dictionary of attributes; or the
        attributes as keyword arguments.
        """

        if hasattr(attrs, 'attrs'):

            attrs = attrs.attrs

        self._update_attrs(attrs, **kwargs)


    def _update_attrs(self, attrs = None, **kwargs):

        attrs = self._add_kwargs(attrs, **kwargs)

        for key, val in iteritems(attrs):

            if key in self.attrs:

                self.attrs[key] = common.combine_attrs((self.attrs[key], val))

            else:

                self.attrs[key] = val


    @staticmethod
    def _add_kwargs(attrs = None, **kwargs):

        attrs = attrs or {}
        attrs.update(kwargs)

        return attrs


    def __iadd__(self, other):

        self.update_attrs(other)

        return self


    def __add__(self, other):

        new = self.__copy__()
        new.update_attrs(other.attrs)

        return new


    def __copy__(self):

        return self.__class__(self.attrs.copy())


    def __iter__(self):

        return iteritems(self.attrs)


    def serialize(self, **kwargs):
        """
        Generates a JSON string with the full contents of the attributes,
        without any whitespace or line break.

        Returns
            (str): The attributes JSON serialized.
        """

        return self._serialize(self.attrs, **kwargs)


    @classmethod
    def _serialize(
            cls,
            attrs,
            top_key_prefix = False,
            prefix_sep = '_',
            **kwargs
        ):

        if not attrs:

            return ''

        param = {
            'indent': None,
            'separators': (',', ':'),
        }

        param.update(kwargs)

        default = param.pop('default', lambda x: x)

        param['default'] = (
            lambda x:
                # sets are not serializable by the json module
                # hence we always convert them to lists
                list(x)
                    if isinstance(x, set) else
                # additional defaults included here
                default
        )

        if top_key_prefix:

            attrs = dict(
                item
                for top_key, val in iteritems(attrs)
                for item in cls._add_prefix(val, top_key, sep = prefix_sep)
            )

        return json.dumps(attrs, **param)


    @staticmethod
    def _add_prefix(d, prefix, sep = '_'):

        d = d if isinstance(d, dict) else d.attrs

        for key, val in iteritems(d):

            yield (
                '%s%s%s' % (prefix, sep, key),
                val
            )


    def __str__(self):

        return self.serialize()


    def __len__(self):

        return len(self.attrs)


    def __contains__(self, other):

        return other in self.attrs


    def __getitem__(self, other):

        return self.attrs[other] if other in self.attrs else None
