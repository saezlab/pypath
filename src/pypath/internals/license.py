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

from future.utils import iteritems

import os
import json
import warnings

import pypath.share.common as common


class License(object):


    _purpose_levels = {
        'commercial': 10,
        'academic': 5,
        'ignore': 0,
    }

    _sharing_levels = {
        'ignore': 20,
        'free': 15,
        'share': 10,
        'alike': 5,
        'noderiv': 2,
        'noshare': 0,
    }

    _attrib_levels = {
        'ignore': 10,
        'noattrib': 5,
        'attrib': 0,
    }

    _synonyms = {
        'for_profit': 'commercial',
        'non_profit': 'academic',
        'free': 'ignore',
        'freedom': 'ignore',
    }


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
        self.purpose = common.to_set(purpose)
        self.sharing = common.to_set(sharing)
        self.attrib = common.to_set(attrib)

        for k, v in iteritems(kwargs):

            setattr(self, k, v)

        for aspect in ('purpose', 'sharing', 'attrib'):

            setattr(
                self,
                '_%s_labels' % aspect,
                common.swap_dict(getattr(self, '_%s_levels' % aspect))
            )


    def __repr__(self):

        return '<License %s>' % self.name


    def __str__(self):

        return self.name


    @property
    def commercial(self):
        """
        Returns `True` if the license allows commercial (for-profit) use.
        """

        return self.allowed('commercial')


    # a synonym
    for_profit = commercial


    @property
    def academic(self):
        """
        Returns `True` if the license allows academic use.
        """

        return self.allowed('academic')


    # a synonym
    non_profit = academic


    @property
    def purpose_levels(self):

        return {self._levels[f] for f in self.purpose if f in self._levels}


    @property
    def max_purpose_level(self):

        return max(self.levels)


    @classmethod
    def purpose_level(cls, use):

        use = self._synonyms[use] if use in self._synonyms else use

        if use in cls._levels:

            return cls._levels[use]

        else:

            warnings.warn(
                'Unknown purpose level for licenses: `%s`. '
                'This will always disable the purpose.' % use
            )

            return 99


    def allowed(self, use):

        return self.max_level >= self.level(use)


    @classmethod
    def from_json(cls, path):

        with open(path, 'r') as fp:

            json_data = json.load(fp)

        return cls.__new__(cls, **json_data)