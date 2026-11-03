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
import rocksdb

POPULATE_BATCH_SIZE = 1000000


class ListAppender(rocksdb.interfaces.AssociativeMergeOperator):

    def merge(self, key, existing_value, value):

        if existing_value:

            return (True, existing_value + b';' + value)

        return (True, value)

    def name(self):

        return b'ListAppender'


class OneToMany2(object):


    def __init__(self, path):

        self._path = path


    def open(self):

        opt = rocksdb.Options()
        opt.create_if_missing = True
        opt.error_if_exists = False
        opt.paranoid_checks = False
        opt.write_buffer_size = 67108864
        opt.max_open_files = 300000
        opt.max_write_buffer_number = 3
        opt.target_file_size_base = 67108864
        opt.merge_operator = ListAppender()

        self.db = rocksdb.DB(self._path, opt)


    def __getitem__(self, key):

        return self.db.get(key.encode('utf-8'))


    def __setitem__(self, key, value):

        self.db.put(key.encode('utf-8'), value.encode('utf-8'))


    def populate(self, iterator, process):

        buffer = rocksdb.WriteBatch()
        self.n_inserted = 0

        for it in iterator:

            k, v = process(it)
            self.n_inserted += 1
            buffer.merge(k.encode('utf-8'), v.encode('utf-8'))

        self.db.write(buffer)


    def populate_from_file(self, fileobj, process = None):

        process = process or (lambda x: tuple(x.split())[:2])

        t = tqdm.tqdm()

        while True:

            self.populate(
                fileobj.readlines(ROCKSDB_POPULATE_BATCH_SIZE),
                process = process,
            )

            t.update(self.n_inserted)

            if self.n_inserted == 0:

                break


    def __enter__(self):

        self.open()

        return self


    def __exit__(self, *args):

        self.close()


    def __del__(self):

        self.close()


    def close(self):

        pass
