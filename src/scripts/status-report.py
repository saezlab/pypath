#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2022
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Olga Ivanova
#           Sebastian Lobentanzer
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

"""
Here we take a snapshot from the current tip of the master branch of pypath,
collect all functions from the pypath.inputs module and run them with an
empty cache directory. Basic metrics about the returned values and errors
are presented. This simple and automatized procedure is able to reveal the
most common errors: if a server is down, if the URL, format or access
restrictions changed, or pypath got broken in the recent developments.
Note: this report is only about pypath, not about the data in the OmniPath
web service; before we build the OmniPath database, we aim to fix the errors
which could have an impact on its content.
"""

import os
import re
import sys
import time
import pkgutil
import importlib
import inspect
import shutil
import traceback
import json
import argparse
import weakref
import types
import csv

import bs4

# Note: modules from pypath are imported in the `init_pypath` method,
# because the source we load the module from depends on the parameters,
# we might even fetch it from git before importing.

EXCLUDE = {
    'biogps_annotations', # calls biomart_microarrays
    'biomart_microarrays', # takes too long
    'common',
    'get_3did_dmi', # takes too long
    'get_csa', # takes too long
    'go_annotations_quickgo',
    'process_3did_dmi', # calls get_3did_dmi
    'threedcomplex_complexes',
    '_uniprot_deleted',
    'uniprot_deleted',
    'swissprot_deleted',
    'trembl_deleted',
}

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
                tr td:nth-of-type(2),
                tr td:nth-of-type(6),
                tr td:nth-of-type(7) {
                    font-family: monospace;
                }
                tr td:nth-of-type(7), tr td:nth-of-type(9) {
                    font-size: xx-small;
                }
                tr td:nth-of-type(9), td.zero {
                    color: deeppink;
                }
                a {
                    display: inline-block;
                }
                .pilcrow {
                    font-size: x-large;
                }
                a.pilcrow:link, a.pilcrow:active, a.pilcrow:visited {
                    color: #aaaaaa;
                    text-decoration: none;
                }
                a.pilcrow:hover {
                    color: deeppink;
                    text-decoration: none;
                }
            </style>
            </head>
            <body>
            <h1>Pypath inputs status report</h1>

            <p id="desc"></p>
            <p id="comp"></p>
            <!-- Input testing summary table -->
            <table id="report-summary">
                <tr>
                    <th>Modules collected:</th>
                    <td id="n_modules"></td>
                </tr>
                <tr>
                    <th>Modules failed to import:</th>
                    <td id="n_import_errors"></td>
                </tr>
                <tr>
                    <th>Functions collected:</th>
                    <td id="n_functions"></td>
                </tr>
                <tr>
                    <th>Functions run without error:</th>
                    <td id="n_success"></td>
                </tr>
                <tr>
                    <th>Functions returned empty value:</th>
                    <td id="n_empty"></td>
                </tr>
                <tr>
                    <th>Functions skipped due to lack of arguments:</th>
                    <td id="n_noargs"></td>
                </tr>
                <tr>
                    <th>Functions run with error:</th>
                    <td id="n_errors"></td>
                </tr>
            </table>
            <br>
            <!-- Input functions table -->
            <table>
                <tbody id="report" valign="top">
                    <tr>
                        <th></th>
                        <th>Function</th>
                        <th>Started</th>
                        <th>Finished</th>
                        <th>Elapsed (s)</th>
                        <th>Result type</th>
                        <th>Result repr</th>
                        <th>Result size</th>
                        <th>Error</th>
                        <th>Change since last time</th>
                        <th>Last succeeded</th>
                    </tr>
                </tbody>
            </table>
            <!-- Footer -->
            <p><em>
                <a href="https://omnipathdb.org/">The OmniPath Team</a>&nbsp;
                &#x2022;
                <a href="https://saezlab.org/">Saez Lab</a>&nbsp;
                &#x2022;
                <span id="date"></span>
            </em></p>
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
    'diff',
    'last_succeeded',
)

COUNTERS = (
    'modules',
    'import_errors',
    'functions',
    'success',
    'empty',
    'noargs',
    'errors',
)

REPORT_TIME_F = '%Y-%m-%d %H:%M:%S'

PYPATH_GIT_URL = 'https://github.com/saezlab/pypath'

CAPTURE_GENERATORS = 300

MAINDIR_PREFIX = 'pypath_inputs_status'

RESULT_JSON_PATH = ('report', 'result.json')


def _log(*args, **kwargs):
    """
    Deliver a log message if the logger is available,
    otherwise write it to stdout.
    """

    if '__log' in globals():

        __log(*args, **kwargs)

    elif args:

        sys.stdout.write(str(args[0]))
        sys.stdout.flush()


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
            pickle_dir = None,
            build_dir = None,
            first = None,
            from_git = None,
            nobuild = None,
            prev_run = None,
        ):

        self.parse_args()
        self.maindir = maindir or self.clargs.dir
        self.cachedir = cachedir or self.clargs.cachedir
        self.pickle_dir = pickle_dir or self.clargs.pickle_dir
        self.build_dir = build_dir or self.clargs.build_dir
        self.first = first or self.clargs.first
        self.nobuild = nobuild or self.clargs.nobuild
        self.prev_dir = prev_run or self.clargs.prev_run
        self.from_git = (
            from_git
                if from_git is not None else
            not self.clargs.nogit
        )
        self.reset_counters()
        self._finalizer = weakref.finalize(
            self,
            lambda obj: obj.finish(),
            self,
        )
        self.finished = False
        self.prev_result = {}


    def main(self):

        self.start()
        self.test_inputs()
        self.finish()
        self.build()


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
            '-p', '--pickle_dir',
            help = 'Database pickle dumps directory path',
            type = str,
        )
        self.clargs.add_argument(
            '-b', '--build_dir',
            help = 'Database build directory path',
            type = str,
        )
        self.clargs.add_argument(
            '-o', '--nobuild',
            help = 'Disable database build',
            action = 'store_true',
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
        self.clargs.add_argument(
            '-r', '--prev_run',
            help = 'Compare to a previous run (path to directory). If not '
                'provided, the script will try to find the directory by '
                'itself. To disable comparison to a previous run, provide '
                '"none" to this argument.',
            type = str,
        )
        self.clargs = self.clargs.parse_args()

    def start(self):

        self.reset_counters()
        self.set_timestamp()
        self.set_dirs() # calls init_pypath
        self.reset_result()
        self.read_prev_result()
        _log('Started generating pypath inputs status report.')


    def finish(self):

        if not self.finished:

            self.set_timestamp(end = True)
            self.save_results()
            self.compile_html()
            self.copy_log()
            self.reset_result()
            self.finished = True

            _log('Finished generating pypath inputs status report.')


    def reset_counters(self):

        for cntr in COUNTERS:

            setattr(self, 'n_%s' % cntr, 0)

        self.finished = False
        self.n_modules = -1


    def reset_result(self):

        self.result = []

        _log('Resetting results.')


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

            if mod.name in EXCLUDE:

                _log('Ignoring module `pypath.inputs.%s`.' % mod.name)

            else:

                yield 'pypath.inputs.%s' % mod.name


    @property
    def modules(self):
        """
        Iterate submodules in the `pypath.inputs` module.
        """

        for mod_name in self.module_names:

            _log('Importing module `%s`.' % mod_name)

            try:

                mod = importlib.import_module(mod_name)
                _log('Imported module `%s`.' % mod_name)

                yield mod

            except Exception as e:

                exc = sys.exc_info()
                _log('Failed to import module `%s`:' % mod_name)
                _logger._log_traceback()
                self.n_import_errors += 1


    @property
    def functions(self):
        """
        Iterate functions defined in the
        """

        for mod in self.modules:

            _log('Collecting functions from module `%s`.' % mod.__name__)
            self.n_modules += 1

            for fun_name, fun in inspect.getmembers(mod, inspect.isfunction):

                _log('Next function `%s`.' % fun_name)
                mod_fun_name = self.function_name(fun)

                if fun.__module__ == mod.__name__:

                    if fun_name in EXCLUDE or mod_fun_name in EXCLUDE:

                        _log('Ignoring function `%s`.' % mod_fun_name)
                        continue

                    _log('Selected function: `%s`.' % mod_fun_name)

                    yield fun

                else:

                    _log(
                        'Ignoring function `%s`: it is from '
                        'another module' % mod_fun_name
                    )


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
            '%s__%s' % (
                MAINDIR_PREFIX,
                time.strftime('%Y%m%d-%H%M%S', self.start_time),
            )
        )
        self.maindir = os.path.abspath(self.maindir)
        self.cachedir = self.cachedir or os.path.join(self.maindir, 'cache')
        self.pickle_dir = (
            self.pickle_dir or
            os.path.join(self.maindir, 'pickles')
        )
        self.build_dir = self.build_dir or os.path.join(self.maindir, 'build')
        self.reportdir = os.path.join(self.maindir, 'report')

        os.makedirs(self.reportdir, exist_ok = True)
        os.makedirs(self.cachedir, exist_ok = True)
        os.makedirs(self.pickle_dir, exist_ok = True)
        os.makedirs(self.build_dir, exist_ok = True)

        os.chdir(self.maindir)

        self.init_pypath()

        settings.setup(cachedir = self.cachedir)
        settings.setup(pickle_dir = self.pickle_dir)

        self.set_prev_dir()


    def set_prev_dir(self):
        """
        Sets up the path to a directory with a previous run of this script.
        Having a previous run makes it possible to find the newly failing
        items and the ones which were failing and just started to work again.
        """

        if self.prev_dir == 'none':

            self.prev_dir = None
            _log('Comparing to previous run is disabled.')
            return

        if isinstance(self.prev_dir, str):

            if (
                not os.path.exists(self.prev_dir) or
                not self.is_status_report_dir(self.prev_dir)
            ):

                _log(
                    'Previous run directory does not '
                    'exist or does not contain output: `%s`.' % self.prev_dir
                )
                self.prev_dir = None

        if self.prev_dir is None:

            self._find_prev_dir()

        if (
            isinstance(self.prev_dir, str) and
            os.path.exists(self.prev_dir) and
            self.is_status_report_dir(self.prev_dir)
        ):

            _log(
                'Will compare to previous run '
                'in directory `%s`.' % self.prev_dir
            )

        else:

            _log(
                'Could not find directory with previous run, '
                'have nothing to compare against the current run.'
            )


    def _find_prev_dir(self):
        """
        This is how we try to find automatically the directory with the
        previous run.
        """

        redir = re.compile(r'%s__(\d{8}-\d{6})' % MAINDIR_PREFIX)

        parent, _ = os.path.split(self.maindir)
        latest = ''

        for d in next(os.walk(parent))[1]:

            m = redir.match(d)

            if m:

                timestamp = m.group(1)
                this_path = os.path.join(parent, d)

                if (
                    timestamp > latest and
                    self.is_status_report_dir(this_path)
                ):

                    latest = timestamp
                    self.prev_dir = this_path


    @staticmethod
    def is_status_report_dir(path):
        """
        Tells if a directory seems to contain a previous run of this script.
        It checks for the most important output file, `report/result.json`.
        """

        path = '' if path is None else path
        result_json = os.path.join(path, *RESULT_JSON_PATH)

        return os.path.exists(result_json)


    def read_prev_result(self):
        """
        Reads the results of a previous run and stores it under the
        `prev_result` attribute.
        """

        self.prev_result = {}

        if self.is_status_report_dir(self.prev_dir):

            json_path = os.path.join(self.prev_dir, *RESULT_JSON_PATH)

            _log('Reading results of previous run from `%s`.' % json_path)

            with open(json_path, 'r') as fp:

                self.prev_result = json.load(fp)

            self.prev_result = dict(
                (i['function'], i)
                for i in self.prev_result
            )


    def init_pypath(self):
        """
        Attempts to find the pypath module in the current working directory
        or parent directory, or clones the latest version from github
        according to the settings. Then imports the required pypath modules.
        """

        self.pypath_from = 'from system directory'

        if self.from_git:

            os.system('git clone --depth 1 %s pypath_git' % PYPATH_GIT_URL)

            os.symlink(
                os.path.join('pypath_git', 'src', 'pypath'),
                'pypath',
                target_is_directory = True,
            )

            self.pypath_from = 'git'

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

            if os.path.exists('pypath'):

                self.pypath_from = 'local directory'

        sys.path.insert(0, os.getcwd())

        import pypath.inputs as inputs
        import pypath.share.session as session
        import pypath.share.settings as settings
        import pypath.omnipath.server.build as build
        import pypath.omnipath as omnipath

        globals()['inputs'] = inputs
        globals()['session'] = session
        globals()['settings'] = settings
        globals()['build'] = build
        globals()['omnipath'] = omnipath

        _logger = session.Logger(name = 'status_report')
        globals()['_logger'] = _logger
        globals()['__log'] = _logger._log
        _log('Working directory: `%s`.' % self.maindir)
        _log('Cache directory: `%s`.' % self.cachedir)
        _log('Pickle directory: `%s`.' % self.pickle_dir)
        _log('Build directory: `%s`.' % self.build_dir)
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
        """
        Calls one function and captures its result or the errors raised.
        """

        fun_name = self.function_name(fun)
        result = {'function': fun_name}
        self.n_functions += 1

        try:

            _args, _kwargs = ARGS.get(fun_name, ((), {}))

            if (
                sum(
                    a.default is inspect._empty
                    for a in inspect.signature(fun).parameters.values()
                ) <=
                len(_args) + len(_kwargs)
            ):

                size = None

                t0 = time.time()
                result['start_time'] = self.strftime(t0)

                _log('Calling function `%s`.' % fun_name)
                value = fun(*_args, **_kwargs)
                _log('Function `%s` returned.' % fun_name)

                if isinstance(
                    value,
                    (
                        types.GeneratorType,
                        filter,
                        map,
                        csv.DictReader,
                    )
                ):

                    _log(
                        'The function returned a generator, '
                        'we fully consume it but capture only '
                        'the first %u items.' % CAPTURE_GENERATORS
                    )
                    value_gen, value = value, []

                    for i, rec in enumerate(value_gen):

                        if i < CAPTURE_GENERATORS:

                            value.append(rec)

                    size = i + 1

                _log('Collecting information about `%s`.' % fun_name)

                t1 = time.time()
                result['end_time'] = self.strftime(t1)
                result['elapsed'] = t1 - t0

                result.update(
                    value_type = type(value).__name__,
                    value_repr = self.to_str(value.__repr__(), 300),
                    value_size = (
                        size
                            if size is not None else
                        len(value)
                            if hasattr(value, '__len__') else
                        None
                    ),
                )

                self.n_success += 1

                if not result['value_size']:

                    self.n_empty += 1

                _log('Collected information about `%s`.' % fun_name)

            else:

                msg = 'Not calling `%s`, not enough arguments.' % fun_name
                _log(msg)
                result['error'] = msg
                self.n_noargs += 1

        except Exception as e:

            _log('Error in function `%s`:' % fun_name)
            _log('Collecting information about `%s`.' % fun_name)
            t1 = time.time()
            result['end_time'] = self.strftime(t1)
            result['elapsed'] = t1 - t0
            exc = sys.exc_info()
            result['error'] = traceback.format_exception(*exc)

            _logger._log_traceback()
            self.n_errors += 1
            _log('Collected information about `%s`.' % fun_name)

        self.compare_to_prev(result)
        self.result.append(result)

        _log('Finished testing `%s`.' % fun_name)
        _log('Result length: %u' % len(self.result))


    def compare_to_prev(self, result):
        """
        Compares the current outcome of a procedure to the one captured at
        a previous run. The differences found will be stored in the result
        dict passed, under the key `diff`.
        """

        diff = {}
        result['last_succeeded'] = result['start_time']

        if result['function'] in self.prev_result:

            prev = self.prev_result[result['function']]

            if 'error' in prev and 'error' not in result:

                diff['fixed'] = True

            if (
                'value_size' in prev and
                'value_size' in result and
                prev['value_size'] != result['value_size']
            ):

                diff['size'] = result['value_size'] - prev['value_size']

            if (
                'value_type' in prev and
                'value_type' in result and
                prev['value_type'] != result['value_type']
            ):

                diff['type'] = {
                    'previous': prev['value_type'],
                    'current': result['value_type'],
                }

            if 'error' in result and 'error' not in prev:

                diff['broke'] = True

            if 'error' in result and 'last_succeeded' in prev:

                result['last_succeeded'] = prev['last_succeeded']

        else:

            diff['first'] = True

        result['diff'] = diff


    @staticmethod
    def strftime(t):
        """
        Converts time to string
        """

        return time.strftime(REPORT_TIME_F, time.localtime(t))


    def save_results(self):
        """
        Saves the results as JSON.
        """

        _log(
            'Exporting results to JSON, result items: %u.' % len(self.result)
        )

        path = os.path.join(self.reportdir, 'result.json')

        with open(path, 'w') as fp:

            json.dump(self.result, fp, indent = 2)

        _log('JSON has been exported to `%s`.' % path)


    def compile_html(self):
        """
        Compiles a HTML report, saves it to the reporting directory.
        """

        _log('Compiling HTML report, result items: %u.' % len(self.result))

        last_commit = self.pypath_git_hash()

        soup = bs4.BeautifulSoup(HTML_TEMPLATE, 'html.parser')
        soup.find(id = 'comp').append(
            bs4.BeautifulSoup(
                'Compiled between <em>%s</em> and <em>%s;</em> '
                'pypath version: %s (from %s%s' % (
                    time.strftime(REPORT_TIME_F, self.start_time),
                    time.strftime(REPORT_TIME_F, self.end_time),
                    sys.modules['pypath']._version.__version__,
                    self.pypath_from,
                    '; <a href="%s/tree/%s">%s</a>)' % (
                        PYPATH_GIT_URL,
                        last_commit,
                        last_commit,
                    )
                        if last_commit else
                    ')',
                ),
                'html.parser',
            )
        )
        soup.find(id = 'desc').string = __doc__

        for attr in COUNTERS:

            soup.find(id = 'n_%s' % attr).string = (
                str(getattr(self, 'n_%s' % attr))
            )

        soup.find(id = 'date').string = time.strftime('%Y-%m-%d')

        path = os.path.join(self.reportdir, 'report.html')

        for r in self.result:

            row = soup.new_tag('tr')
            error = bool(r.get('error', False))
            row['class'] = (
                row.get('class', []) +
                ['error' if error else 'success']
            )

            # anchor for each function
            label = '-'.join(r['function'].split('.')[-2:])
            cell = soup.new_tag('td')
            a = soup.new_tag(
                'a',
                attrs = {
                    'id': label,
                    'class': 'pilcrow',
                    'href': '#%s' % label,
                },
            )
            a.string = '\xB6'
            cell.append(a)
            row.append(cell)

            for field in FIELDS:

                cell = soup.new_tag('td')

                if field in r:

                    if field == 'error':

                        pre = soup.new_tag('pre')
                        pre.string = ''.join(r[field])
                        cell.append(pre)

                    else:

                        cell.string = self.to_str(r[field], maxlen = 10000)

                        if field == 'size' and cell.string == '0':

                            cell['class'] = cell.get('class', []) + ['zero']

                row.append(cell)

            soup.find(id = 'report').append(row)

        with open(path, 'w') as fp:

            fp.write(soup.prettify())

        _log('HTML report has been exported to `%s`.' % path)


    @classmethod
    def to_str(cls, value, maxlen = 76):

        return (
            '{:,.2f}'.format(value)
                if isinstance(value, float) else
            '{:,}'.format(value)
                if isinstance(value, int) else
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

        path = os.path.join(self.reportdir, 'pypath.log')
        _log('Copying log file to `%s`.' % path)
        time.sleep(5)

        shutil.copy2(
            _logger._logger.fname,
            path,
        )


    def build(self):
        """
        Builds the OmniPath databases.
        """

        if self.nobuild:

            _log('Database build disabled.')
            return

        self.builder = build.WebserviceTables(
            build_dir = self.build_dir,
        )

        databases = (
            'interactions',
            'complexes',
            'enz_sub',
            'annotations',
            'intercell',
        )

        for db in databases:

            try:

                _log('Starting to build database `%s`.' % db)
                getattr(self.builder, db)()
                _log('Finished building database `%s`.' % db)
                _ = [
                    omnipath.db.remove_db(db, ncbi_tax_id = ncbi_tax_id)
                    for ncbi_tax_id in (9606, 10090, 10116)
                ]

            except Exception as e:

                exc = sys.exc_info()
                _log('Failed to build database `%s`:' % db)
                _logger._log_traceback()


if __name__ == '__main__':

    s = StatusReport()
    s.main()
