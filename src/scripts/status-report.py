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
import time
import datetime
import pkgutil
import importlib
import inspect
import shutil
import traceback
import json
import argparse

import bs4

# Note: modules from pypath are imported in the `init_pypath` method,
# because the source we load the module from depends on the parameters,
# we might even fetch it from git before importing.

EXCLUDE = {'common'}

ARGS = {}

HTML_TEMPLATE = (
    '''
        <!DOCTYPE html>
        <html>
            <head>
            <title>Pypath inputs status report</title>
            <style>
                table {
                    border: 1px solid #333333;
                    border-collapse: collapse;
                }
                table tr:nth-child(even){
                    background-color: #f2f2f2;
                }
                table tr:hover {
                    background-color: #ddd;
                }
                th {
                    color: white;
                    font-weight: bold;
                    text-align: left;
                    background-color: #333333;
                }
                th, td {
                    padding: 8px;
                }
                tr td:nth-of-type(1),
                tr td:nth-of-type(5),
                tr td:nth-of-type(6) {
                    font-family: monospace;
                }
                tr td:nth-of-type(6), tr td:nth-of-type(8) {
                    font-size: xx-small;
                }
            </style>
            </head>
            <body>
            <h1>Pypath inputs status report</h1>
            <p id="date"></p>
            <table>
                <tbody valign="top">
                    <tr>
                        <th>Function</th>
                        <th>Started</th>
                        <th>Finished</th>
                        <th>Elapsed (s)</th>
                        <th>Result type</th>
                        <th>Result repr</th>
                        <th>Result size</th>
                        <th>Error</th>
                    </tr>
                </tbody>
            </table>
            </body>
        </html>
    '''
)

FIELDS = (
    'function',
    'start_time',
    'end_time',
    'elapsed',
    'value_type',
    'value_repr',
    'value_size',
    'error',
)

REPORT_TIME_F = '%Y-%m-%d %H:%M:%S'

PYPATH_GIT_URL = 'https://github.com/saezlab/pypath'


class StatusReport(object):
    """
    Calls all functions defined in the submodules of the `pypath.inputs`
    module and reports about their runtime, result and errors.

    Args:
        cachedir (str): The cache directory. By default a new empty cache
            directory is created at each run.
        first (int): For testing purposes, run only the first some functions.
    """

    def __init__(
            self,
            maindir = None,
            cachedir = None,
            first = None,
            from_git = None,
        ):

        self.parse_args()
        self.maindir = maindir or self.clargs.dir
        self.cachedir = cachedir or self.clargs.cachedir
        self.first = first or self.clargs.first
        self.from_git = (
            from_git
                if from_git is not None else
            not self.clargs.nogit
        )


    def main(self):

        self.start()
        self.test_inputs()
        self.finish()


    def parse_args(self):
        """
        Parses command line arguments if present.
        """

        self.clargs = argparse.ArgumentParser()
        self.clargs.add_argument(
            '-f', '--first',
            help = 'Run only the first N functions',
            type = int,
        )
        self.clargs.add_argument(
            '-c', '--cachedir',
            help = 'Cache directory path',
            type = str,
        )
        self.clargs.add_argument(
            '-d', '--dir',
            help = 'Main directory path',
            type = str,
        )
        self.clargs.add_argument(
            '-g', '--nogit',
            help = 'Do not fetch the latest pypath from git, use the one '
                'in the current directory or installed in system paths',
            action = 'store_true',
        )
        self.clargs = self.clargs.parse_args()

    def start(self):

        self.result = []
        self.set_timestamp()
        self.set_dirs()
        _log('Started generating pypath inputs status report.')


    def finish(self):

        self.set_timestamp(end = True)
        self.save_results()
        self.compile_html()
        self.copy_log()
        _log('Finished generating pypath inputs status report.')


    def test_inputs(self):
        """
        Iterates over all functions defined in `pypath.inputs`, calls
        them one by one, saves data about the returned values and errors.
        """

        _log('Testing functions from `pypath.inputs`.')

        for i, fun in enumerate(self.functions):

            if isinstance(self.first, int) and i >= self.first:

                _log('First %u functions have been executed, exiting.' % i)
                break

            self.test_input(fun)


    @property
    def module_names(self):
        """
        Iterate module names in the `pypath.inputs` module.
        """

        for mod in pkgutil.iter_modules(inputs.__path__):

            if mod.name not in EXCLUDE:

                yield 'pypath.inputs.%s' % mod.name


    @property
    def modules(self):
        """
        Iterate submodules in the `pypath.inputs` module.
        """

        for mod_name in self.module_names:

            _log('Importing module `%s`.' % mod_name)
            mod = importlib.import_module(mod_name)

            yield mod


    @property
    def functions(self):
        """
        Iterate functions defined in the
        """

        for mod in self.modules:

            for fun_name, fun in inspect.getmembers(mod, inspect.isfunction):

                if fun.__module__ == mod.__name__:

                    _log('Selected function: `%s`.' % self.function_name(fun))

                    yield fun


    def set_timestamp(self, end = False):

        setattr(
            self,
            '%s_time' % ('end' if end else 'start'),
            time.localtime(),
        )


    def set_dirs(self):
        """
        Sets up the paths to the reporting and cache directories, creates
        the directories.
        """

        self.maindir = (
            self.maindir or
            'pypath_inputs_status__%s' % time.strftime(
                '%Y%m%d-%H%M%S',
                self.start_time,
            )
        )
        self.maindir = os.path.abspath(self.maindir)
        self.cachedir = self.cachedir or os.path.join(self.maindir, 'cache')
        self.reportdir = os.path.join(self.maindir, 'report')

        os.makedirs(self.reportdir, exist_ok = True)
        os.makedirs(self.cachedir, exist_ok = True)

        os.chdir(self.maindir)

        self.init_pypath()

        settings.setup(cachedir = self.cachedir)


    def init_pypath(self):
        """
        Attempts to find the pypath module in the current working directory
        or parent directory, or clones the latest version from github
        according to the settings. Then imports the required pypath modules.
        """

        if self.from_git:

            os.system('git clone --depth 1 %s pypath_git' % PYPATH_GIT_URL)

            os.symlink(
                os.path.join('pypath_git', 'src', 'pypath'),
                'pypath',
                target_is_directory = True,
            )

        else:

            parentdir_pypath = os.path.join('..', 'pypath')

            if (
                not os.path.exists('pypath') and
                os.path.exists(parentdir_pypath)
            ):

                os.symlink(
                    parentdir_pypath,
                    'pypath',
                    target_is_directory = True,
                )

        sys.path.insert(0, os.getcwd())

        import pypath.inputs as inputs
        import pypath.share.session as session
        import pypath.share.settings as settings

        globals()['inputs'] = inputs
        globals()['session'] = session
        globals()['settings'] = settings

        _logger = session.Logger(name = 'status_report')
        globals()['_logger'] = _logger
        globals()['_log'] = _logger._log
        _log('Working directory: `%s`.' % self.maindir)
        _log('Cache directory: `%s`.' % self.cachedir)
        _log('Reporting directory: `%s`.' % self.reportdir)
        _log('Pypath from git: `%s`.' % self.from_git)
        _log('Pypath local path: `%s`.' % os.path.dirname(inputs.__path__[0]))


    @staticmethod
    def timestamp():
        """
        Current time as a string timestamp.
        """

        return time.strftime('%Y%m%d-%H%M%S')


    @staticmethod
    def function_name(fun):
        """
        Returns the name of a function as a string, together with its
        module path.
        """

        return '%s.%s' % (fun.__module__, fun.__name__)


    def test_input(self, fun):

        fun_name = self.function_name(fun)
        result = {'function': fun_name}

        try:

            _args, _kwargs = ARGS.get(fun_name, ((), {}))

            if (
                sum(
                    a.default is inspect._empty
                    for a in inspect.signature(fun).parameters.values()
                ) <=
                len(_args) + len(_kwargs)
            ):

                t0 = time.localtime()
                result['start_time'] = time.strftime(REPORT_TIME_F, t0)
                value = fun(*_args, **_kwargs)
                t1 = time.localtime()
                result['end_time'] = time.strftime(REPORT_TIME_F, t1)
                result['elapsed'] = (
                    datetime.datetime(*t1[:6]) -
                    datetime.datetime(*t0[:6])
                ).total_seconds()

                result.update(
                    value_type = type(value).__name__,
                    value_repr = self.to_str(value.__repr__(), 300),
                    value_size = (
                        len(value)
                            if hasattr(value, '__len__') else
                        None
                    ),
                )

            else:

                msg = 'Not calling `%s`, not enough arguments.' % fun_name
                _log(msg)
                result['error'] = msg

        except Exception as e:

            exc = sys.exc_info()
            result['error'] = traceback.format_exception(*exc)
            _log('Error in function `%s`:' % fun_name)
            _logger._log_traceback()

        self.result.append(result)


    def save_results(self):
        """
        Saves the results as JSON.
        """

        path = os.path.join(self.reportdir, 'result.json')

        with open(path, 'w') as fp:

            json.dump(self.result, fp)


    def compile_html(self):
        """
        Compiles a HTML report, saves it to the reporting directory.
        """

        last_commit = self.pypath_git_hash()

        soup = bs4.BeautifulSoup(HTML_TEMPLATE, 'html.parser')
        soup.find(id = 'date').append(
            bs4.BeautifulSoup(
                'Compiled between <em>%s</em> and <em>%s</em>; '
                'pypath version: %s%s' % (
                    time.strftime(REPORT_TIME_F, self.start_time),
                    time.strftime(REPORT_TIME_F, self.end_time),
                    sys.modules['pypath']._version.__version__,
                    ' (<a href="%s/tree/%s">%s</a>)' % (
                        PYPATH_GIT_URL,
                        last_commit,
                        last_commit,
                    )
                        if last_commit else
                    '',
                ),
                'html.parser',
            )
        )

        path = os.path.join(self.reportdir, 'report.html')

        for r in self.result:

            row = soup.new_tag('tr')
            error = bool(r.get('error', False))
            row['class'] = (
                row.get('class', []) +
                ['error' if error else 'success']
            )

            for field in FIELDS:

                cell = soup.new_tag('td')

                if field in r:

                    if field == 'error':

                        pre = soup.new_tag('pre')
                        pre.string = ''.join(r[field])
                        cell.append(pre)

                    else:

                        cell.string = self.to_str(r[field], maxlen = 10000)

                row.append(cell)

            soup.find('tbody').append(row)

        with open(path, 'w') as fp:

            fp.write(soup.prettify())


    @classmethod
    def to_str(cls, value, maxlen = 76):

        return (
            '{:,}'.format(value)
                if isinstance(value, (int, float)) else
            cls.truncate(str(value), maxlen = maxlen)
        )


    @staticmethod
    def truncate(value, maxlen = 76):

        if len(value) > maxlen + 14:

            value = '%s...(truncated)' % value[:maxlen]

        return value


    def pypath_git_hash(self):

        if os.path.exists('pypath_git'):

            head = os.path.join('pypath_git', '.git', 'logs', 'HEAD')

            with open(head, 'r') as fp:

                return fp.readline().split(maxsplit = 1)[-1][:7]


    def copy_log(self):
        """
        Copies the logfile to the reporting directory.
        """

        shutil.copy2(
            _logger._logger.fname,
            os.path.join(self.reportdir, 'pypath.log')
        )


if __name__ == '__main__':

    s = StatusReport()
    s.main()
