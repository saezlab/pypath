#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2021
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
import sys
import types
import itertools
import collections
import sqlite3

import pypath.share.common as common

QUERIES = {
    'list_tables':
        '''
        SELECT name
        FROM sqlite_master
        WHERE type = 'table';
        ''',
    'create_table':
        '''
        CREATE TABLE %s (
            a VARCHAR NOT NULL,
            b VARCHAR NOT NULL,
            PRIMARY KEY (a, b)
        );
        ''',
    'lookup_pair':
        '''
        SELECT a, b
        FROM %s
        WHERE a = ? AND b = ?;
        ''',
    'lookup_one':
        '''
        SELECT %s
        FROM %s
        WHERE %s = ?;
        ''',
    'lookup_many':
        '''
        SELECT %s, %s
        FROM %s
        WHERE %s IN (%s);
        ''',
    'contains_many':
        '''
        SELECT %s
        FROM %s
        WHERE %s IN (%s);
        ''',
    'insert_one':
        '''
        INSERT OR IGNORE INTO %s (a, b)
        VALUES (?, ?);
        ''',
    'insert_many':
        '''
        INSERT OR IGNORE INTO %s (a, b)
        VALUES %s;
        ''',
    'delete':
        '''
        DROP TABLE %s;
        ''',
    'peek':
        '''
        SELECT * FROM %s LIMIT ?;
        ''',
    'remove_one':
        '''
        DELETE FROM %s
        WHERE %s = ?;
        ''',
    'remove_one_both':
        '''
        DELETE FROM %s
        WHERE a = ? OR b = ?;
        ''',
    'remove_many':
        '''
        DELETE FROM %s
        WHERE %s IN (%s);
        ''',
    'remove_many_both':
        '''
        DELETE FROM %s
        WHERE a IN (%s) OR b IN (%s);
        ''',
    'count':
        '''
        SELECT COUNT(*) FROM %s;
        ''',
    'select_all':
        '''
        SELECT * FROM %s;
        ''',
}

BULK_INSERT_LENGTH = 490
FETCHMANY_BATCH_SIZE = 10000
COLUMNS = {'a': 0, 'b': 1}
UNION_SIZE_LIMIT = 10000000


class LookupTable(object):


    def __init__(self, path):

        self._path = path
        self._table = 'pypath_lookup'

        self.open()


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        import importlib as imp
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def open(self):

        self.con = sqlite3.connect(self._path)
        self.cur = self.con.cursor()
        self.ensure_table()


    def _query(self, key, *args, values = None):

        values = (
            list(values)
                if isinstance(values, (set, types.GeneratorType)) else
            values
        )
        values = (values,) if values else ()

        args = tuple(self._process_query_param(arg) for arg in args)

        query = QUERIES.get(key, key) % args

        return self._execute(query, *values)


    @staticmethod
    def _process_query_param(value):

        return (
            value
                if common.is_str(value) else
            ','.join(value)
                if isinstance(value, common.list_like) else
            str(value)
        )


    def _execute(self, query, *args):

        execute = (
            self.cur.executemany
                if args and isinstance(args[0], list) else
            self.cur.execute
        )

        return execute(query, *args)


    def table_exists(self):
        """
        Tells if the table associated to this lookup table exists in the
        sqlite database.
        """

        return (self._table,) in self._query('list_tables').fetchall()


    def create_table(self):
        """
        Creates a table in the sqlite database to store the data associated
        to this lookup table.
        """

        self._query('create_table', self._table)


    def ensure_table(self):
        """
        Makes sure the table associated to this lookup table exists in the
        sqlite database.
        """

        if not self.table_exists():

            self.create_table()


    def close(self):

        if hasattr(self, 'con') and hasattr(self.con, 'close'):

            self.con.close()


    def lookup_pair(self, a, b):

        result = self._query('lookup_pair', self._table, values = (a, b))

        return bool(result.fetchone())


    def __contains__(self, value):

        if isinstance(value, tuple):

            return self.lookup_pair(*value[:2])

        else:

            return self.lookup_one_a(value) or self.lookup_one_b(value)


    def lookup(self, key):

        return (
            self.lookup_many_a(key) or
            self.lookup_many_b(key)
                if isinstance(key, common.list_like) else
            self.lookup_one_a(key) or
            self.lookup_one_b(key)
        )


    def lookup_one(self, value, col):

        col_other = 'b' if col == 'a' else 'a'

        result = self._query(
            'lookup_one',
            col_other,
            self._table,
            col,
            values = (value,),
        )

        return {r[0] for r in result.fetchall()}


    def lookup_one_a(self, value):

        return self.lookup_one(value, 'a')


    def lookup_one_b(self, value):

        return self.lookup_one(value, 'b')


    def lookup_many(self, values, col):

        col_other = 'b' if col == 'a' else 'a'
        values = self._ensure_tuple(values)
        question_marks = self._question_marks(values)

        result = self._query(
            'lookup_many',
            col,
            col_other,
            self._table,
            col,
            question_marks,
            values = values,
        )

        output = collections.defaultdict(set)

        for a, b in result.fetchall():

            output[a].add(b)

        return dict(output)


    @staticmethod
    def _ensure_tuple(value):

        return (
            value
                if isinstance(value, tuple) else
            tuple(value)
                if isinstance(value, common.list_like) else
            (value,)
        )


    def lookup_many_a(self, values):

        return self.lookup_many(values, 'a')


    def lookup_many_b(self, values):

        return self.lookup_many(values, 'b')


    def insert_one(self, a, b):

        self._query('insert_one', self._table, values = (a, b))


    def insert_many(self, values_ = None, **kwargs):

        values_ = (
            values_.items()
                if hasattr(values_, 'items') else
            values_
                if isinstance(
                    values_,
                    common.list_like + (
                        types.GeneratorType,
                        collections.abc.ItemsView,
                    )
                ) else
            ()
        )
        values_ = (
            self._ensure_tuple(i)
            for i in itertools.chain(
                values_,
                kwargs.items()
            )
        )

        if not self.con.in_transaction:

            self._query('BEGIN;')

        try:

            while True:

                values = list(itertools.islice(values_, BULK_INSERT_LENGTH))

                if not values:

                    break

                question_marks = self._question_marks(values)

                self.cur.execute(
                    QUERIES['insert_many'] % (
                        self._table,
                        question_marks,
                    ),
                    tuple(itertools.chain(*values)),
                )

            self._query('COMMIT;')

        except:

            self._query('ROLLBACK;')
            raise


    # synonym
    update = insert_many


    @classmethod
    def _question_marks(cls, values):

        if isinstance(values, tuple):

            return ','.join('?' * len(values))

        elif isinstance(values, list):

            qm_item = cls._question_marks(values[0])
            return '(%s)' % (
                '), ('.join(
                    qm_item for _ in range(len(values))
                )
            )

        else:

            raise ValueError


    def wipe(self):

        self._delete_table()
        self.ensure_table()


    def reset(self):

        self.close()
        self.open()


    def _delete_table(self):

        if self.table_exists():

            self._query('delete', self._table)


    def delete_database(self):

        self.close()
        os.remove(self._path)


    def peek(self, limit = 10, **kwargs):

        if not self.table_exists():

            sys.stdout.write('The table does not exist.\n')
            sys.stdout.flush()
            return

        result = self._query(
            QUERIES['peek'],
            self._table,
            values = (limit,),
        ).fetchall()

        if not result:

            sys.stdout.write('The table is empty.\n')
            sys.stdout.flush()

        else:

            common.print_table(
                dict(
                    zip(
                        ('a', 'b'),
                        zip(*result)
                    )
                ),
                **kwargs
            )

            total = len(self)

            if total > limit:

                sys.stdout.write('And %u more records.\n' % (total - limit))
                sys.stdout.flush()


    def column_index(self, col):

        return COLUMNS[col]


    def remove_one(self, value, col = 'a'):

        self._query('remove_one', self._table, 'a', values = (value,))


    def remove_one_a(self, value):

        self.remove_one(value, 'a')


    def remove_one_b(self, value):

        self.remove_one(value, 'b')


    def remove_one_both(self, value):

        self._query('remove_one_both', self._table, values = (value, value))


    def remove_many(self, values, col = 'a'):

        values = self._ensure_tuple(values)
        question_marks = self._question_marks(values)

        if values:

            self._query(
                'remove_many',
                self._table,
                col,
                question_marks,
                values = values,
            )


    def remove_many_a(self, values):

        self.remove_many(values, 'a')


    def remove_many_b(self, values):

        self.remove_many(values, 'b')


    def remove_many_both(self, values):

        values = self._ensure_tuple(values)
        question_marks = self._question_marks(values)
        values = values * 2

        if values:

            self._query(
                'remove_many_both',
                self._table,
                question_marks,
                question_marks,
                values = values,
            )


    def __getitem__(self, key):

        return self.lookup(key)


    def __setitem__(self, key, value):

        self.insert_one(key, value)


    def __del__(self):

        self.close()


    def __enter__(self):

        self.open()

        return self


    def __exit__(self, *args):

        self.close()


    def __isub__(self, value):

        if common.is_str(value):

            self.remove_one_both(value)

        elif isinstance(value, common.list_like):

            self.remove_many_both(value)

        else:

            raise ValueError

        return self


    def __delitem__(self, value):

        self.remove_one_a(value)


    def __len__(self):

        return self._query('count', self._table).fetchall()[0][0]


    def __nonzero__(self):

        return self.__len__() > 0

    __bool__ = __nonzero__


    def contains_many(self, values, col = 'a'):


        if not values:

            return set()

        values = self._ensure_tuple(values)
        question_marks = self._question_marks(values)

        result = self._query(
            'contains_many',
            col,
            self._table,
            col,
            question_marks,
            values = values,
        )

        return {r[0] for r in result.fetchall()}


    def contains_many_a(self, values):

        return self.contains_many(values, col = 'a')


    def contains_many_b(self, values):

        return self.contains_many(values, col = 'b')


    def contains_many_both(self, values):

        return (
            self.contains_many(values, col = 'a') |
            self.contains_many(values, col = 'b')
        )

    __and__ = contains_many_both


    def __or__(self, values):

        self._check_size_limit()

        values = common.to_set(values)

        return self.to_set() | values


    def __xor__(self, values):

        self._check_size_limit()

        values = common.to_set(values)

        return self.to_set() ^ values


    def _check_size_limit(self):

        if self.__len__() > UNION_SIZE_LIMIT:

            raise RuntimeError(
                'Union and symmetric difference operations are not allowed '
                'on lookup tables longer than %u to avoid running out of '
                'memory. To bypass this limitation, increase the value of '
                '`%s.UNION_SIZE_LIMIT`.' % (
                    UNION_SIZE_LIMIT,
                    self.__module__,
                )
            )


    def to_set(self):

        return {i for record in self for i in record}


    @property
    def table_name(self):

        return self._table


    @property
    def path(self):

        return self._path


    def __repr__(self):

        return '<LookupTable %s>' % self._path


    def __iter__(self):

        result = self._query('select_all', self._table)

        while True:

            batch = result.fetchmany(FETCHMANY_BATCH_SIZE)

            if not batch:

                break

            for it in batch:

                yield it
