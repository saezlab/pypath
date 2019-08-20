#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2019
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

#
# This is a multi-threaded wrapper over python's MySQLdb module
#

from past.builtins import xrange, range, reduce

import sys
import codecs
import time

try:
    import pymysql as MySQLdb
    import pymysql.cursors as cursors
except:
    sys.stdout.write('\t:: No MySQL support.\n')

import hashlib
import queue
import threading

import pypath.mysql_connect as mysql_connect


class MysqlRunner(object):
    
    def __init__(self,
                 param,
                 cursor='serverside',
                 concurrent_queries=4,
                 log=None,
                 silent=False):
        '''
        param is either a tuple of the name of mysql access config file,
        and the title of the config section in it, or
        or a dict with the config itself
        '''
        self.log = log
        self.param = param
        self.silent = silent
        self.cons = {}
        self.curs = {}
        self.results = {}
        self.scheduled = {}
        self.running = {}
        self.ready = {}
        self.controls = []
        self.listeners = {}
        self.lock = threading.Lock()
        self.max_con_value = 0
        self.connect_lock = False
        self.ss_cursor = 'SSDictCursor'
        self.cs_cursor = 'DictCursor'
        self.default_cursor = self.ss_cursor \
            if cursor == 'serverside' else self.cs_cursor
        self.max_con_q = 'SHOW VARIABLES LIKE "max_connections"'
        if type(self.param) is tuple:
            self.access = mysql_connect.MysqlConnect(
                self.param[0], log=self.log)
        with open('mysql.log', 'w') as f:
            f.write('')
        if not self.test_connection():
            emsg = 'Could not connect to MySQL'
            self.send_error(emsg)
        else:
            self.max_connections()
            self.concurrent_queries = min(self.max_con_value / 2,
                                          concurrent_queries)
            self.waiting_tasks = queue.Queue()
            self.running_tasks = queue.Queue(maxsize=self.concurrent_queries)
            for i in xrange(self.concurrent_queries):
                self.add_thread()
            self.control_main = threading.Thread(
                target=self._control_main, name='main control thread')
            self.control_main.daemon = True
            self.control_main.start()

    def add_thread(self):
        cthread = threading.Thread(
            target=self._control,
            name='sub control thread #%u' % len(self.controls))
        cthread.daemon = True
        cthread.start()
        self.controls.append(cthread)

    def get_connection(self, cursor=None, **kwargs):
        try:
            return MySQLdb.connect(
                host=self.param['host'],
                user=self.param['user'],
                port=self.param['port'],
                passwd=self.param['password'],
                db=self.param['db'],
                cursorclass=getattr(cursors, cursor or self.default_cursor),
                connect_timeout=12,
                **kwargs)
        except:
            return None

    def connekt(self, cursor=None, priority=False, **kwargs):
        if type(self.param) is tuple:
            con = self.access.get_connection(
                self.param[1], cursor=cursor or self.default_cursor, **kwargs)
        else:
            if 'port' not in self.param:
                self.param['port'] = 3306
            con = self.get_connection(cursor or self.default_cursor, **kwargs)
        return con

    def test_connection(self):
        con = self.connekt(cursor=self.cs_cursor, priority=True)
        if con is not None:
            con.close()
            return True
        else:
            return False

    def send_query(self, query, cursor=None, silent=None):
        qid = self.get_qid(query)
        task = {'query': query, 'cursor': cursor, 'silent': silent, 'qid': qid}
        self.lock.acquire()
        self.scheduled[qid] = task
        self.lock.release()
        self.waiting_tasks.put(qid)

    def _control_main(self):
        while True:
            qid = self.waiting_tasks.get()
            task = self.scheduled[qid]
            self.lock.acquire()
            del self.scheduled[qid]
            self.lock.release()
            task['cursor'] = self.cs_cursor \
                if len(self.cons) + self.concurrent_queries \
                > self.max_con_value / 4 \
                else self.default_cursor \
                if task['cursor'] is None \
                else task['cursor']
            self.lock.acquire()
            self.running[qid] = task
            self.lock.release()
            self.running_tasks.put(task)
            self.waiting_tasks.task_done()

    def _control(self):
        while True:
            task = self.running_tasks.get()
            self._run_query(**task)
            self.running_tasks.task_done()

    def _run_query(self, query, cursor, qid, silent=False):
        silent = self.silent if silent is None else silent
        res = []
        if not silent:
            self.lock.acquire()
            sys.stdout.write('\t:: Waiting for MySQL...')
            sys.stdout.flush()
            self.lock.release()
        self.lock.acquire()
        with open('mysql.log', 'a') as f:
            f.write(query + '\n')
        self.lock.release()
        con = self.connekt(cursor=cursor)
        self.lock.acquire()
        self.cons[qid] = con
        self.lock.release()
        if self.cons[qid] is not None:
            try:
                self.curs[qid] = self.cons[qid].cursor()
                self.curs[qid].execute(query)
                if cursor == self.cs_cursor:
                    res = self.curs[qid].fetchall()
                    self.curs[qid].close()
                    self.cons[qid].close()
                    self.lock.acquire()
                    del self.curs[qid]
                    del self.cons[qid]
                    self.lock.release()
                else:
                    res = self._read_result(qid)
                self.lock.acquire()
                task = self.running[qid]
                self.ready[qid] = task
                del self.running[qid]
                self.results[qid] = res
                for l in self.listeners.values():
                    l.set()
                if not silent:
                    sys.stdout.write(' Done.\n')
                self.lock.release()
            except MySQLdb.Error as e:
                emsg = 'MySQL error occured. See `mysql.error` for details.'
                self.send_error(emsg)
                out  = "MySQL Error [%d]: %s\n\n" % (e.args[0], e.args[1])
                out += "Failed to execute query:\n\n"
                out += query
                self.lock.acquire()
                with codecs.open('mysql.error', 'w') as f:
                    f.write(out)
                self.lock.release()
        else:
            emsg = 'No connection to MySQL'
            self.send_error(emsg)

    def _read_result(self, qid):
        while True:
            recs = self.curs[qid].fetchmany(100)
            if len(recs) == 0:
                self.curs[qid].close()
                self.cons[qid].close()
                self.lock.acquire()
                del self.curs[qid]
                del self.cons[qid]
                del self.results[qid]
                self.lock.release()
                break
            for rec in recs:
                yield rec

    def wait_results(self, qids):
        qids = set(qids)
        self.lock.acquire()
        i = len(self.listeners)
        self.listeners[i] = threading.Event()
        self.lock.release()
        while True:
            self.lock.acquire()
            if len(qids - set(self.results.keys())) == 0:
                del self.listeners[i]
                self.lock.release()
                return None
            else:
                self.listeners[i].clear()
                self.lock.release()
                self.listeners[i].wait()

    def print_status(self):
        if len(self.cons) > 0:
            pids = ','.join(['%u' % c.thread_id() for c in self.cons.values()])
            con = self.connekt()
            cur = con.cursor()
            q = '''SELECT COUNT(DISTINCT(Id)) AS THREADS,MAX(TIME) AS TIME
                FROM information_schema.processlist 
                WHERE Id IN (%s) GROUP BY 'a';''' % pids
            cur.execute(q)
            res = cur.fetchone()
            cur.close()
            con.close()
            if res is not None:
                self.lock.acquire()
                sys.stdout.write('\r' + ' ' * 90)
                sys.stdout.write(
                    '\r\t:: MySQL: %u queries running for %u seconds.' %
                    (res['THREADS'], res['TIME']))
                sys.stdout.flush()
                self.lock.release()
                return None
        self.lock.acquire()
        sys.stdout.write('\r' + ' ' * 90)
        sys.stdout.write('\r\t:: MySQL: finished.')
        sys.stdout.flush()
        self.lock.release()

    def processlist(self):
        if len(self.cons) > 0:
            pids = ','.join(['%u' % c.thread_id() for c in self.cons.values()])
            con = self.connekt()
            cur = con.cursor()
            q = '''SELECT ID,TIME,STATE
                FROM information_schema.processlist 
                WHERE ID IN (%s);''' % pids
            cur.execute(q)
            res = cur.fetchall()
            cur.close()
            con.close()
            if len(res) > 0:
                self.lock.acquire()
                sys.stdout.write('\tID\tState\t\t\tRunning for\n')
                for p in res:
                    status = p['STATE'] if len(p[
                        'STATE']) > 0 else 'Unknown state'
                    sys.stdout.write('\t%u\t%s\t\t\t%u s\n' %
                                     (p['ID'], status, p['TIME']))
                sys.stdout.flush()
                self.lock.release()
                return None
        self.lock.acquire()
        sys.stdout.write('No queries running.\n')
        sys.stdout.write('\n')
        sys.stdout.write('\t:: Open connections: %u\n' % len(self.cons))
        sys.stdout.write('\t:: Unread resultsets: %u\n' % len(self.results))
        sys.stdout.flush()
        self.lock.release()

    def send_error(self, error_message):
        if self.log is not None:
            self.lock.acquire()
            self.log.msg(1, error_message, 'ERROR')
            self.lock.release()
        else:
            self.lock.acquire()
            sys.stdout.write('\n\t:: ' + error_message + '\n\n')
            sys.stdout.flush()
            self.lock.release()

    def get_result(self, qid):
        result = [] if qid not in self.results else self.results[qid]
        return result

    def clean(self, qid):
        if qid in self.results:
            self.lock.acquire()
            if qid in self.curs:
                self.curs[qid].close()
                del self.curs[qid]
            if qid in self.cons:
                self.cons[qid].close()
                del self.cons[qid]
            del self.results[qid]
            del self.ready[qid]
            self.lock.release()

    def max_connections(self):
        con = self.connekt(cursor='DictCursor', priority=True)
        cur = con.cursor()
        cur.execute(self.max_con_q)
        res = cur.fetchone()
        cur.close()
        con.close()
        self.max_con_value = int(res['Value']) - 1

    def num_of_connections(self):
        q = 'SHOW PROCESSLIST;'
        con = self.access.get_connection(self.param[1], self.cs_cursor)
        cur = con.cursor()
        cur.execute(q)
        res = cur.fetchall()
        cur.close()
        con.close()
        return len(res)

    def message(self, msg):
        self.lock.acquire()
        sys.stdout.write('\r%s\r' % (' ' * 90))
        sys.stdout.write('\t:: %s' % msg)
        sys.stdout.flush()
        self.lock.release()

    def get_qid(self, query):
        '''
        Returns the 32 byte md5sum of a string:
        this serves as a unique identifier of queries,
        referred as `qid` in this module.

        @query : str
            MySQL query or any other string.
        '''
        if hasattr(query, 'encode'):
            query = query.encode('utf-8')
        return hashlib.md5(query).hexdigest()
