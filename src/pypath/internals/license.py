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

import warnings

import pypath.share.common as common


class License(object):


    _levels = {
        'commercial': 10,
        'academic': 5,
        'ignore': 0,
    }

    _synonyms = {
        'for_profit': 'commercial',
        'non_profit': 'academic',
        'free': 'ignore',
        'freedom': 'ignore',
    }


    def __init__(self, name, full_name = None, features = None):

        self.name = name
        self.full_name = full_name or name
        self.features = common.to_set(features)


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
    def levels(self):

        return {self._levels[f] for f in self.features if f in self._levels}


    @property
    def max_level(self):

        return max(self.levels)


    @classmethod
    def level(cls, use):

        use = self._synonyms[use] if use in self._synonyms else use

        if use in cls._levels:

            return cls._levels[use]

        else:

            warnings.warn(
                'Unknown usage level for licenses: `%s`. '
                'This will always disable the usage.' % use
            )

            return 99


    def allowed(self, use):

        return self.max_level >= self.level(use)