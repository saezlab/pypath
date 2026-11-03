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

import tqdm
import lmdbm

POPULATE_BATCH_SIZE = 5000000


class OneToMany(lmdbm.Lmdb):


    def _pre_key(self, value):

        return value.encode('utf-8')


    def _post_key(self, value):

        return value.decode('utf-8')


    def _pre_value(self, value):

        return ';'.join(value).encode('utf-8')


    def _post_value(self, value):

        return value.decode('utf-8').split(';')


    def populate(self, iterator, process = lambda x: x):

        write_cache = {}
        self.n_inserted = 0

        for it in iterator:

            k, v = process(it)
            self.n_inserted += 1

            if k not in write_cache:

                if k in self:

                    write_cache[k] = self[k]

                else:

                    write_cache[k] = []

            write_cache[k].append(v)

            if len(write_cache) > POPULATE_BATCH_SIZE:

                self.update(write_cache)
                write_cache = {}

        if write_cache:

            self.update(write_cache)


    def populate_from_file(self, fileobj, process = None):

        process = process or (lambda x: tuple(x.split())[:2])

        t = tqdm.tqdm()

        while True:

            self.populate(
                fileobj.readlines(POPULATE_BATCH_SIZE),
                process = process,
            )

            t.update(self.n_inserted)

            if self.n_inserted == 0:

                break
