#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  (Planned for) centrally handling cache for all databases/resources.
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

import os

import pypath.internals.license as license_mod
import pypath.share.session as session
import pypath.share.settings as settings
import pypath.share.common as common


class Licenses(session.Logger):


    def __init__(self, license_dir = None):

        session.Logger.__init__(self, name = 'licenses')

        self.license_dir = license_dir or settings.get('license_dir')
        self.licenses = {}
        self.populate()


    def add_license(self, license):

        if isinstance(license, common.basestring) and os.path.exists(license):

            license = license_mod.License.from_json(
                path = license,
                source = license,
            )

        if isinstance(license, license_mod.License):

            self.licenses[license.name] = license
            self._log(
                'Added license `%s` (%s).' % (
                    license.name,
                    license.features_str,
                )
            )

        else:

            self._log('Could not process license: `%s`.' % str(license))


    def populate(self, license_dir = None):

        license_dir = license_dir or self.license_dir

        if os.path.isdir(license_dir):

            self._log(
                'Populating license database from '
                'directory `%s`.' % license_dir
            )

            for fname in os.listdir(license_dir):

                if fname.endswith('.json'):

                    self.add_license(os.path.join(license_dir, fname))

        else:

            self._log('License directory `%s` doesn\'t exist.' % license_dir)


    def __len__(self):

        return len(self.licenses)


    def __repr__(self):

        return '<License database (%u licenses)>' % len(self)


    def __getitem__(self, key):

        if isinstance(key, license_mod.License):

            if key.name not in self:

                self.licenses[key.name] = key

            return key

        elif key in self:

            return self.licenses[key]

        else:

            self._log('Missing license: `%s`.' % str(key))
            return key


    def __contains__(self, key):

        return key in self.licenses