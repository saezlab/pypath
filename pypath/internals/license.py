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

import os
import json
import warnings

import pypath.share.common as common


# these are all ordinal scales from the smallest to the greatest freedom
_purpose_levels = {
    'composite': 25,
    'free': 20,
    'commercial': 15,
    'for_profit': 15,
    'forprofit': 15,
    'nonprofit': 10,
    'non_profit': 10,
    'academic': 5,
    'ignore': 0,
}

_sharing_levels = {
    'composite': 20,
    'free': 25,
    'share': 20,
    'deriv': 20,
    'alike': 15,
    'noderiv': 10,
    'noshare': 5,
    'ignore': 0,
}

_attrib_levels = {
    'composite': 20,
    'free': 10,
    'noattrib': 10,
    'attrib': 5,
    'ignore': 0,
}

_int_purpose_levels = common.swap_dict(_purpose_levels)
_int_sharing_levels = common.swap_dict(_sharing_levels)
_int_attrib_levels = common.swap_dict(_attrib_levels)


class LicenseFeature(object):

    levels = {}
    name = 'feature'
    int_levels = common.swap_dict(levels)


    def __init__(self, level):

        self.level = level
        _ = self.check_level(level)


    @classmethod
    def level_to_int(cls, level):
        """
        Return the value for the key ``level`` as an integer on the ordinal
        scale of the levels.
        """

        if cls.check_level(level):

            return cls.levels[level]

        return 99


    @classmethod
    def int_to_level(cls, i):
        """
        Returns a set of labels corresponding to the level nearest to the
        integer ``i`` on the ordinal scale of levels.
        """

        if i > max(cls.int_levels.keyd()):

            i = max(cls.int_levels.keys())

        elif i not in cls.int_levels:

            for _i in sorted(cls.int_levels.keys()):

                if _i > i:

                    i = _i
                    break

        return cls.int_levels[i]


    def __str__(self):

        return self.level


    def __repr__(self):

        return '<License %s: %s>' % (self.name, self.level)


    @classmethod
    def check_level(cls, level):
        """
        Checks of ``level`` is a valid key for the
        """

        if level not in cls.levels:

            cls.unknown_warning(level)
            return False

        return True


    def to_int(self):
        """
        Returns the value of the current level of this license feature as
        an integer on the ordinal scale of levels.
        """

        return self.level_to_int(self.level)


    def __int__(self):

        return self.to_int()


    @classmethod
    def ensure_int(cls, other):

        if isinstance(other, int):

            return other

        elif isinstance(other, str):

            return cls.level_to_int(other)

        else:

            try:

                return int(other)

            except (TypeError, ValueError):

                cls.unknown_warning(other)
                return 99


    def __eq__(self, other):

        i_level = self.ensure_int(other)

        return int(self) == i_level


    def __gt__(self, other):

        i_level = self.ensure_int(other)

        return int(self) > i_level


    def __ge__(self, other):

        i_level = self.ensure_int(other)

        return int(self) >= i_level


    def __lt__(self, other):

        i_level = self.ensure_int(other)

        return int(self) < i_level


    def __le__(self, other):

        i_level = self.ensure_int(other)

        return int(self) <= i_level


    def enables(self, other):

        return self >= other


    @classmethod
    def unknown_warning(cls, value):

        warnings.warn(
            'Unknown `%s` level for licenses: `%s`. '
            'This will always disable the resource.' % (
                cls.name,
                str(value),
            )
        )


class LicensePurpose(LicenseFeature):

    levels = _purpose_levels
    name = 'purpose'
    int_levels = _int_purpose_levels


    def __init__(self, level):

        LicenseFeature.__init__(self, level)


class LicenseSharing(LicenseFeature):

    levels = _sharing_levels
    name = 'sharing'
    int_levels = _int_sharing_levels


    def __init__(self, level):

        LicenseFeature.__init__(self, level)


class LicenseAttrib(LicenseFeature):

    levels = _attrib_levels
    name = 'attrib'
    int_levels = _int_attrib_levels


    def __init__(self, level):

        LicenseFeature.__init__(self, level)


class License(object):


    def __init__(
            self,
            name,
            full_name = None,
            purpose = None,
            sharing = None,
            attrib = 'attrib',
            url = None,
            **kwargs
        ):

        self.name = name
        self.full_name = full_name or name
        self.purpose = LicensePurpose(purpose)
        self.sharing = LicenseSharing(sharing)
        self.attrib = LicenseAttrib(attrib)
        self.url = url

        for k, v in iteritems(kwargs):

            setattr(self, k, v)


    def __repr__(self):

        return '<License %s>' % self.name


    def __str__(self):

        return self.name


    def enables(self, purpose, sharing = None, attrib = None):
        """
        Checks if the license enables a particular use, according to purpose,
        sharing and attribution. For example, to check if the license enables
        academic use with redistribution under a compatible license, call
        ``License.enables(purpose = 'academic', sharing = 'alike')``.
        """

        return (
            (
                not attrib or
                self.attrib.enables(attrib)
            ) and
            (
                not sharing or
                self.sharing.enables(sharing)
            ) and
            (
                self.purpose.enables(purpose)
            )
        )


    @property
    def features(self):

        return dict(
            (
                aspect,
                getattr(self, aspect).level
            )
            for aspect in ('purpose', 'sharing', 'attrib')
        )


    @property
    def features_str(self):

        return common.dict_str(self.features)


    @classmethod
    def from_json(cls, path, **kwargs):

        with open(path, 'r') as fp:

            json_data = json.load(fp)
            json_data.update(kwargs)

        return cls(**json_data)


    # some shortcut methods
    @classmethod
    def _generate_enables_methods(cls):

        def get_enables_method(aspect, level):

            @property
            def enables_method(self):

                return getattr(self, aspect).enables(level)

            return enables_method


        for aspect in ('attrib', 'sharing', 'purpose'):

            for level in globals()['_%s_levels' % aspect]:

                if level in ('composite', 'ignore', 'free', 'attrib'):

                    continue

                method = get_enables_method(aspect, level)
                method_name = level

                setattr(cls, level, method)


License._generate_enables_methods()