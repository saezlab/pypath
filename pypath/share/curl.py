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
from past.builtins import xrange, range

import importlib as imp
import sys
import os
import io
import shutil
import struct

import pypath.share.settings as settings
import pypath.share.session as session_mod
import pypath.share.cache as cache_mod
_logger = session_mod.log()

import pycurl
try:
    from cStringIO import StringIO
    BytesIO = StringIO
except:
    try:
        from StringIO import StringIO
        from StringIO import StringIO as BytesIO
    except:
        from io import BytesIO
        from io import StringIO

try:
    import cPickle as pickle
except:
    import pickle

import urllib

try:
    import urllib2
except ImportError:
    # this works seemless in Py3:
    import urllib.request
    urllib2 = urllib.request

try:
    import urlparse
except:
    # this works seemless in Py3:
    import urllib.parse
    urlparse = urllib.parse

if not hasattr(urllib, 'quote'):
    import urllib.parse
    _urllib = urllib
    urllib = _urllib.parse

try:
    import pysftp
except:
    _logger.msg(
        'Module `pysftp` not available. '
        'Only downloading of a small number of resources '
        'relies on this module. '
        'Please install by PIP if it is necessary for you.',
        'curl',
        -1,
    )
import codecs
import gzip
import zipfile
import tarfile
import hashlib
import re

from contextlib import closing

import pypath.share.progress as progress
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.share.settings as settings

try:
    basestring
except NameError:
    basestring = str

if 'long' not in __builtins__:
    long = int

if 'unicode' not in __builtins__:
    unicode = str

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

# global contexts for modifying Curl() behviour
CACHE = None
CACHEDEL = False
CACHEPRINT = False
DRYRUN = False
PRESERVE = False
DEBUG = False

LASTCURL = None

show_cache = False

_re_url = re.compile(r'^(?:http|https|ftp)://')


def is_url(url):
    return bool(_re_url.match(url))



class _global_context(object):
    """
    This is a metaclass for context handlers working by
    setting a module level variable to certain value.
    """

    def __init__(self, name, on_off):
        """
        :param str name: Name of the module level variable.
        :param on_off: Value of the module level variable in the context.
        """
        self.name = name
        self.module = sys.modules[__name__]
        self.on_off = on_off

    def __enter__(self):
        self._store_value = getattr(self.module, self.name)
        setattr(self.module, self.name, self.on_off)

    def __exit__(self, exception_type, exception_value, traceback):
        if exception_type is not None:
            sys.stdout.write(
                '%s, %s, %s\n' %
                (str(exception_type), str(exception_value), str(traceback)))
            sys.stdout.flush()
        setattr(self.module, self.name, self._store_value)


class _global_context_on(_global_context):
    """
    This is a metaclass for context handlers working by
    setting a module level variable to `True`.
    """

    def __init__(self, name):
        """
        :param str name: Name of the module level variable.
        """
        super(_global_context_on, self).__init__(name, True)


class _global_context_off(_global_context):
    """
    This is a metaclass for context handlers working by
    setting a module level variable to `False`.
    """

    def __init__(self, name):
        """
        :param str name: Name of the module level variable.
        """
        super(_global_context_off, self).__init__(name, False)


class cache_on(_global_context_on):
    """
    This is a context handler to turn on pypath.curl.Curl() cache.
    As most of the methods use cache as their default behaviour,
    probably it won't change anything.

    Behind the scenes it sets the value of the `pypath.curl.CACHE`
    module level variable to `True` (by default it is `None`).

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        print('`curl.CACHE` is ', curl.CACHE)

        with curl.cache_on():
            print('`curl.CACHE` is ', curl.CACHE)
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(cache_on, self).__init__('CACHE')


class cache_off(_global_context_off):
    """
    This is a context handler to turn off pypath.curl.Curl() cache.
    Data will be downloaded even if it exists in cache.

    Behind the scenes it sets the value of the `pypath.curl.CACHE`
    module level variable to `False` (by default it is `None`).

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        print('`curl.CACHE` is ', curl.CACHE)

        with curl.cache_on():
            print('`curl.CACHE` is ', curl.CACHE)
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(cache_off, self).__init__('CACHE')


class cache_print_on(_global_context_on):
    """
    This is a context handler which makes pypath.curl.Curl() print
    verbose messages about its cache.

    Behind the scenes it sets the value of the `pypath.curl.CACHEPRINT`
    module level variable to `True` (by default it is `False`).

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_print_on():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(cache_print_on, self).__init__('CACHEPRINT')


class cache_print_off(_global_context_off):
    """
    This is a context handler which stops pypath.curl.Curl() to print
    verbose messages about its cache.

    Behind the scenes it sets the value of the `pypath.curl.CACHEPRINT`
    module level variable to `False`. As by default it is `False`, this
    context won't modify the default behaviour.

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_print_off():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(cache_print_off, self).__init__('CACHEPRINT')


class cache_delete_on(_global_context_on):
    """
    This is a context handler which results pypath.curl.Curl() deleting the
    cache files instead of reading it. Then it downloads the data again,
    or does nothing if the `DRYRUN` context is turned on. Upon deleting
    cache files console messages will let you know which files have been
    deleted.

    Behind the scenes it sets the value of the `pypath.curl.CACHEDEL`
    module level variable to `True` (by default it is `False`).

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_delete_on():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(cache_delete_on, self).__init__('CACHEDEL')


class cache_delete_off(_global_context_off):
    """
    This is a context handler which stops pypath.curl.Curl() deleting the
    cache files. This is the default behaviour, so this context won't
    change anything by default.

    Behind the scenes it sets the value of the `pypath.curl.CACHEDEL`
    module level variable to `False`.

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_delete_off():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(cache_delete_off, self).__init__('CACHEDEL')


class dryrun_on(_global_context_on):
    """
    This is a context handler which results pypath.curl.Curl() to do all
    setup steps, but do not perform download or cache read.

    Behind the scenes it sets the value of the `pypath.curl.DRYRUN`
    module level variable to `True` (by default it is `False`).

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_dryrun_on():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(dryrun_on, self).__init__('DRYRUN')


class dryrun_off(_global_context_off):
    """
    This is a context handler which results pypath.curl.Curl() to
    perform download or cache read. This is the default behaviour,
    so applying this context restores the default.

    Behind the scenes it sets the value of the `pypath.curl.DRYRUN`
    module level variable to `False`.

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_dryrun_off():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(dryrun_off, self).__init__('DRYRUN')


class preserve_on(_global_context_on):
    """
    This is a context handler which results pypath.curl.Curl() to make
    a reference to itself in the module level variable `LASTCURL`. This
    is useful if you have some issue with `Curl`, and you want to access
    the instance for debugging.

    Behind the scenes it sets the value of the `pypath.curl.PRESERVE`
    module level variable to `True` (by default it is `False`).

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_preserve_on():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(preserve_on, self).__init__('PRESERVE')


class preserve_off(_global_context_off):
    """
    This is a context handler which avoids pypath.curl.Curl() to make
    a reference to itself in the module level variable `LASTCURL`. By
    default it does not do this, so this context only restores the
    default.

    Behind the scenes it sets the value of the `pypath.curl.PRESERVE`
    module level variable to `False`.

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_preserve_off():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(preserve_off, self).__init__('PRESERVE')

class debug_on(_global_context_on):
    """
    This is a context handler which results pypath.curl.Curl() to print
    debug information.
    This is useful if you have some issue with `Curl`, and you want to
    see what`s going on.

    Behind the scenes it sets the value of the `pypath.curl.DEBUG`
    module level variable to `True` (by default it is `False`).

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_debug_on():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(debug_on, self).__init__('DEBUG')


class debug_off(_global_context_off):
    """
    This is a context handler which avoids pypath.curl.Curl() to print
    debug information.
    By default it does not do this, so this context only restores the
    default.

    Behind the scenes it sets the value of the `pypath.curl.DEBUG`
    module level variable to `False`.

    Example: ::

        import pypath
        from pypath import curl, data_formats

        pa = pypath.PyPath()

        with curl.cache_debug_off():
            pa.load_resources({'signor': data_formats.pathway['signor']})
    """

    def __init__(self):
        super(debug_off, self).__init__('DEBUG')


class RemoteFile(object):
    def __init__(self,
                 filename,
                 user,
                 host,
                 passwd,
                 port = 22,
                 sep = '\t',
                 header = True,
                 rownames = True):

        for key, val in iteritems(locals()):
            setattr(self, key, val)
        self.conf = fabric.config.Config()
        env.keepalive = 60
        env.connection_attempts = 5
        env.password = self.passwd

    def wcl(self):
        with closing(
                connect(self.user, self.host, self.port, HostConnectionCache(
                ))) as ssh:
            stdin, stdout, stderr = ssh.exec_command('wc -l %s' %
                                                     self.filename)
            return int(stdout.readlines()[0].split()[0]) - (1 if self.header
                                                            else 0)

    def rowns(self):
        with closing(
                connect(self.user, self.host, self.port, HostConnectionCache(
                ))) as ssh:
            stdin, stdout, stderr = ssh.exec_command(
                'awk \'BEGIN{FS = "%s"}{print $1}\' %s%s' %
                (self.sep, self.filename, ''
                 if not self.header else ' | tail -n +2'))
            return [x.strip() for x in stdout.readlines()]

    def open(self, return_header = True):
        with closing(
                connect(self.user, self.host, self.port, HostConnectionCache(
                ))) as ssh:
            with closing(ssh.open_sftp()) as sftp:
                with closing(sftp.open(self.filename)) as f:
                    if not return_header:
                        line = f.readline()
                    for line in f:
                        yield line


class FileOpener(session_mod.Logger):
    """
    This class opens a file, extracts it in case it is a
    gzip, tar.gz, tar.bz2 or zip archive, selects the requested
    files if you only need certain files from a multifile archive,
    reads the data from the file, or returns the file pointer,
    as you request. It examines the file type and size.
    """

    FORBIDDEN_CHARS = re.compile(r'[/\\<>:"\?\*\|]')

    def __init__(
            self,
            file_param,
            compr = None,
            extract = True,
            _open = True,
            set_fileobj = True,
            files_needed = None,
            large = True,
            default_mode = 'r',
            encoding = 'utf-8',
        ):

        if not hasattr(self, '_logger'):

            session_mod.Logger.__init__(self, name = 'file')

        if not hasattr(self, 'encoding') or not self.encoding:

            self.encoding = encoding

        if not hasattr(self, 'default_mode'):

            self.default_mode = default_mode

        if not hasattr(self, 'compr'):
            self.compr = compr
        if not hasattr(self, 'files_needed'):
            self.files_needed = files_needed
        if not hasattr(self, 'large'):
            self.large = large
        self.fname = file_param \
            if type(file_param) in _const.CHAR_TYPES else file_param.name
        self.fileobj = None \
            if type(file_param) in _const.CHAR_TYPES else file_param
        if not hasattr(self, 'type'):
            self.get_type()
        if _open:
            self.open()
        if extract:
            self.extract()


    def open(self):
        """
        Opens the file if exists.
        """

        if self.fileobj is None and os.path.exists(self.fname):

            if self.encoding and self.type == 'plain':

                self.fileobj = open(
                    self.fname,
                    self.default_mode,
                    encoding = (
                        None if self.default_mode == 'rb' else self.encoding
                    ),
                )

            else:

                self.fileobj = open(self.fname, 'rb')


    def extract(self):
        """
        Calls the extracting method for compressed files.
        """

        getattr(self, 'open_%s' % self.type)()


    def open_tgz(self):
        """
        Extracts files from tar gz.
        """

        self._log('Opening tar.gz file `%s`.' % self.fileobj.name)

        self.files_multipart = {}
        self.sizes = {}
        self.tarfile = tarfile.open(fileobj = self.fileobj, mode = 'r:gz')
        self.members = self.tarfile.getmembers()

        for m in self.members:

            if (self.files_needed is None or m.name in self.files_needed) \
                    and m.size != 0:
                # m.size is 0 for dierctories
                this_file = self.tarfile.extractfile(m)
                self.sizes[m.name] = m.size
                if self.large:
                    self.files_multipart[m.name] = this_file
                else:
                    self._log(
                        'Reading contents of file '
                        'from archive: `%s`.' % m.name
                    )
                    self.files_multipart[m.name] = this_file.read()
                    this_file.close()

        if not self.large:
            self.tarfile.close()
            self._log('File closed: `%s`.' % self.fileobj.name)

        self.result = self.files_multipart


    def open_gz(self):

        self._log('Opening gzip file `%s`.' % self.fileobj.name)

        self.fileobj.seek(-4, 2)
        self.size = struct.unpack('I', self.fileobj.read(4))[0]
        self.fileobj.seek(0)
        self.gzfile = gzip.GzipFile(fileobj = self.fileobj)

        if self.large:

            io.DEFAULT_BUFFER_SIZE = 4096
            self._gzfile_mode_r = io.TextIOWrapper(
                self.gzfile,
                encoding = self.encoding,
            )
            self.result = self.iterfile(
                self.gzfile
                    if self.default_mode == 'rb' else
                self._gzfile_mode_r
            )
            self._log(
                'Result is an iterator over the '
                'lines of `%s`.' % self.fileobj.name
            )

        else:

            self.result = self.gzfile.read()
            self.gzfile.close()
            self._log(
                'Data has been read from gzip file `%s`. '
                'The file has been closed' % self.fileobj.name
            )


    def open_zip(self):

        self._log('Opening zip file `%s`.' % self.fileobj.name)

        self.files_multipart = {}
        self.sizes = {}
        self.fileobj.seek(0)
        self.zipfile = zipfile.ZipFile(self.fileobj, 'r')
        self.members = self.zipfile.namelist()
        for i, m in enumerate(self.members):
            self.sizes[m] = self.zipfile.filelist[i].file_size
            if self.files_needed is None or m in self.files_needed:
                this_file = self.zipfile.open(m)
                if self.large:

                    if self.default_mode == 'rb':

                        # keeping it in binary mode
                        self.files_multipart[m] = this_file

                    else:

                        # wrapping the file for decoding
                        self.files_multipart[m] = io.TextIOWrapper(
                            this_file, encoding=self.encoding
                        )
                else:
                    self.files_multipart[m] = this_file.read()
                    this_file.close()

        if not self.large:

            self.zipfile.close()
            self._log(
                'Data has been read from zip file `%s`.'
                'File has been closed' % self.fileobj.name
            )

        self.result = self.files_multipart

    def open_plain(self):

        self._log('Opening plain text file `%s`.' % self.fileobj.name)

        self.size = os.path.getsize(self.fileobj.name)

        if self.large:

            self.result = self.iterfile(self.fileobj)

        else:

            self.result = self.fileobj.read()
            self.fileobj.close()
            self._log(
                'Contents of `%s` has been read '
                'and the file has been closed.' % self.fileobj.name
            )

    def get_type(self):

        self.multifile = False
        if self.fname[-3:].lower() == 'zip' or self.compr == 'zip':
            self.type = 'zip'
            self.multifile = True
        elif self.fname[-3:].lower() == 'tgz' or \
                self.fname[-6:].lower() == 'tar.gz' or \
                self.compr == 'tgz' or self.compr == 'tar.gz':
            self.type = 'tgz'
            self.multifile = True
        elif self.fname[-2:].lower() == 'gz' or self.compr == 'gz':
            self.type = 'gz'
        else:
            self.type = 'plain'


    @staticmethod
    def iterfile(fileobj):

        for line in fileobj:

            yield line


class Curl(FileOpener):
    """
    This class is a wrapper around pycurl.
    You can set a vast amount of parameters.
    In addition it has a caching functionality: using this downloads
    of databases/resources is performed only once.
    It handles HTTP, FTP, cookies, headers, GET and POST params,
    multipart/form data, URL quoting, redirects, timeouts, retries,
    encodings, debugging.
    It returns either downloaded data, file pointer, files extracted
    from archives (gzip, tar.gz, zip).
    It is able to show a progress and status indicator on the console.
    """

    def __init__(
            self,
            url,
            silent = True,
            get = None,
            post = None,
            req_headers = None,
            cache = True,
            debug = False,
            outf = None,
            compr = None,
            encoding = None,
            files_needed = None,
            connect_timeout = None,
            timeout = None,
            ignore_content_length = False,
            init_url = None,
            init_fun = 'get_jsessionid',
            init_use_cache  =  False,
            follow = True,
            large = False,
            default_mode = 'r',
            override_post = False,
            init_headers = False,
            return_headers = False,
            compressed = False,
            binary_data = None,
            write_cache = True,
            force_quote = False,
            sftp_user = None,
            sftp_passwd = None,
            sftp_passwd_file = '.secrets',
            sftp_port = 22,
            sftp_host = None,
            sftp_ask = None,
            setup = True,
            call = True,
            process = True,
            retries = None,
            cache_dir = None,
            bypass_url_encoding = False,
            empty_attempt_again = True,
            keep_failed = False,
            alpn = True,
            slow = False,
            http2 = True,
        ):

        if not hasattr(self, '_logger'):

            session_mod.Logger.__init__(self, name = 'curl')

        self.result = None
        self.download_failed = False
        self.status = 0
        self.get = get
        self.large = large
        self.default_mode = default_mode
        self.silent = silent
        self.debug = debug or DEBUG
        self.url = url
        self.local_file = os.path.exists(self.url)
        self.get = get
        self.force_quote = force_quote
        self.bypass_url_encoding = bypass_url_encoding
        self.empty_attempt_again = empty_attempt_again
        self.keep_failed = keep_failed
        self.alpn = alpn
        self.http2 = http2

        self._log(
            'Creating Curl object to retrieve '
            'data from `%s`' % self.url[:200] # we just don't flood the log
                                              # with super long URLs
        )

        if not self.local_file:

            self.process_url()
            self.url_fix()
            self.set_get()

        else:

            self._log('The URL is a local file path.')
            self.filename = os.path.split(self.url)[-1]

        self.compr = compr
        # self.get_type()
        self.progress = None

        self.encoding = encoding
        self.files_needed = files_needed

        self.follow_http_redirect = follow
        self.timeout = (
            settings.get('curl_extended_timeout')
                if slow else
            settings.get('curl_timeout')
        )
        self.connect_timeout = settings.get('curl_connect_timeout')
        self.ignore_content_length = ignore_content_length
        self.override_post = override_post
        self.retries = retries or settings.get('curl_retries') or 3
        self.req_headers = req_headers
        self._req_headers_list()
        self.post = post
        self.get = get
        self.binary_data = binary_data

        self.cache_dir = cache_dir
        self.cache = cache
        self.init_cache()

        if self.local_file:

            self.cache_file_name = self.url
            self.use_cache = True

        self.write_cache = write_cache
        self.outfile = outf

        self.init_url = init_url
        self.init_fun = init_fun
        self.init_use_cache = init_use_cache

        self.sftp_host = sftp_host
        self.sftp_ask = sftp_ask
        self.sftp_port = sftp_port
        self.sftp_passwd = sftp_passwd
        self.sftp_user = sftp_user
        self.sftp_passwd_file = sftp_passwd_file

        if CACHEPRINT:

            self.show_cache()

        if CACHEDEL:

            self.delete_cache_file()
            self.init_cache()

        if not self.use_cache and not DRYRUN:

            self.title = None
            self.set_title()
            if self.sftp_host is not None:
                self.sftp_url()
                self.sftp_call()
            else:
                self.progress_setup()
                if setup:
                    self.curl_setup()
                if call:
                    self.curl_call()

        elif not self.silent:

            if self.local_file:

                self._log('Loading data from local file `%s`' % self.url)

            else:

                self._log(
                    'Loading data from cache '
                    'previously downloaded from `%s`' % self.domain
                )

        if process and not self.download_failed and not DRYRUN:
            self.process_file()

        if DRYRUN:
            self.print_debug_info('INFO', 'DRYRUN PERFORMED, RETURNING NONE')

        if PRESERVE:
            self.print_debug_info('INFO', 'PRESERVING Curl() INSTANCE '
                                  'IN pypath.curl.LASTCURL')
            setattr(sys.modules[__name__], 'LASTCURL', self)


    def __del__(self):

        fileattrs = ['tarfile', 'zipfile', 'gzfile', 'fileobj']

        for fattr in fileattrs:

            if hasattr(self, fattr):

                f = getattr(self, fattr)

                if hasattr(f, 'close') and sys.getrefcount(f) <= 1:

                    f.close()


    def reload(self):

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def print_debug_info(self, typ, msg):

        msg = self._bytes_to_unicode(msg)
        self._log('CURL DEBUG INFO: %s' % typ)
        self._log(msg)


    def process_url(self):

        self.domain = (
            self.url.replace('https://', '').replace('http://', '').\
            replace('ftp://', '').split('/')[0]
        )
        self.filename = self.url.split('/')[-1].split('?')[0]


    def is_quoted(self, string):
        """
        From http://stackoverflow.com/questions/
        1637762/test-if-string-is-url-encoded-in-php
        """

        test = string
        while (urllib.unquote(test) != test):
            test = urllib.unquote(test)
        return urllib.quote(test,
                            '/%') == string or urllib.quote(test) == string


    def is_quoted_plus(self, string):

        test = string
        while (urllib.unquote_plus(test) != test):
            test = urllib.unquote_plus(test)
        return urllib.quote_plus(
            test, '& = ') == string or urllib.quote_plus(test) == string


    def url_fix(self, charset = 'utf-8'):
        """
        From http://stackoverflow.com/a/121017/854988
        """

        if self.bypass_url_encoding:

            return

        url_raw = self.url

        if type(self.url) is bytes:

            self.url = self._bytes_to_unicode(self.url, encoding = charset)

        scheme, netloc, path, qs, anchor = urlparse.urlsplit(self.url)

        if self.force_quote or not self.is_quoted(path):

            path = urllib.quote(path, '/%')

        if self.force_quote or not self.is_quoted_plus(qs):

            qs = urllib.quote_plus(qs, '& = ')

        self.url = urlparse.urlunsplit((scheme, netloc, path, qs, anchor))

        if self.url != url_raw:

            self._log('Quoted URL: `%s`.' % self.url[:200])


    def set_title(self):

        if self.title is None:
            self.title = 'Downloading `%s` from %s' % \
                (self.filename, self.domain)


    def set_post(self):

        if type(self.post) is dict:

            self.postfields = urllib.urlencode(self.post)
            self.curl.setopt(self.curl.POSTFIELDS, self.postfields)
            self.curl.setopt(self.curl.POST, 1)
            self._log('POST parameters set: %s' % self.postfields[:100])

        else:

            self.postfields = None


    def set_get(self):

        if self.get is not None:

            if isinstance(self.get, dict):
                self.qs = '&'.join(
                    map(lambda param: '%s=%s' % (param[0], param[1]),
                        map(lambda param: (
                                urllib.quote_plus(param[0]),
                                urllib.quote_plus(param[1])
                            ),
                            iteritems(self.get)
                        )
                    )
                )

            elif isinstance(self.get, list):
                self.qs = '&'.join(['='.join(param) for param in
                [[urllib.quote_plus(arg1), urllib.quote_plus(arg2)] for arg1, arg2 in
                [item.split('=') for item in self.get]]])

            self.url = '%s%s%s' % (
                self.url,
                '&' if '?' in self.url else '?',
                self.qs
            )

            self._log(
                'GET parameters added to the URL: `%s`' % self.qs[:100]
            )


    def construct_binary_data(self):
        """
        The binary data content of a `form/multipart` type request
        can be constructed from a list of tuples (<field name>, <field value>),
        where field name and value are both type of bytes.
        """

        bdr = b'---------------------------%s' % \
            common.random_string(28).encode('ascii')
        self.binary_data_param = self.binary_data
        self.binary_data = b'\r\n'.join(
            map(lambda i: b'--%s\r\nContent-Disposition: form-data;'
                b' name="%s"\r\n\r\n%s' % (bdr, i[0], i[1]),
                self.binary_data_param))
        self.binary_data = b'%s\r\n--%s--\r\n' % (self.binary_data, bdr)
        self.req_headers.append(
            'Content-Type: multipart/form-data; boundary=%s' %
            bdr.decode('ascii'))
        self.req_headers.append('Content-Length: %u' % len(self.binary_data))


    def set_binary_data(self):
        """
        Set binary data to be transmitted attached to POST request.

        `binary_data` is either a bytes string, or a filename, or
        a list of key-value pairs of a multipart form.
        """

        if self.binary_data:

            if type(self.binary_data) is list:
                self.construct_binary_data()
            if type(self.binary_data) is bytes:
                self.binary_data_size = len(self.binary_data)
                self.binary_data_file = BytesIO()
                self.binary_data_file.write(self.binary_data)
                self.binary_data_file.seek(0)
            elif os.path.exists(self.binary_data):
                self.binary_data_size = os.path.getsize(self.binary_data)
                self.binary_data_file = open(self.binary_data, 'rb')

            self.curl.setopt(pycurl.POST, 1)
            self.curl.setopt(pycurl.POSTFIELDSIZE, self.binary_data_size)
            self.curl.setopt(pycurl.READFUNCTION, self.binary_data_file.read)
            self.curl.setopt(pycurl.CUSTOMREQUEST, 'POST')
            self.curl.setopt(pycurl.POSTREDIR, 3)
            self._log('Binary data added to query (not showing).')


    def curl_init(self, url = False):

        self.curl = pycurl.Curl()
        self.set_url(url = url)
        self.curl.setopt(self.curl.SSL_VERIFYPEER, False)

        if DEBUG:

            self._log(
                'Following HTTP redirects: %s' % (
                    str(self.follow_http_redirect)
                )
            )

        self.curl.setopt(self.curl.FOLLOWLOCATION, self.follow_http_redirect)
        self.curl.setopt(self.curl.CONNECTTIMEOUT, self.connect_timeout)
        self.curl.setopt(self.curl.TIMEOUT, self.timeout)
        self.curl.setopt(self.curl.TCP_KEEPALIVE, 1)
        self.curl.setopt(self.curl.TCP_KEEPIDLE, 2)
        self.curl.setopt(self.curl.SSL_ENABLE_ALPN, self.alpn)

        if not self.http2:

            self.curl.setopt(
                self.curl.HTTP_VERSION,
                pycurl.CURL_HTTP_VERSION_1_1,
            )

        if self.ignore_content_length:

            self.curl.setopt(self.curl.IGNORE_CONTENT_LENGTH, 136)


    def set_url(self, url = False):

        url = url or self.url

        if isinstance(url, basestring):
            url = url.encode('utf-8')

        self.curl.setopt(self.curl.URL, url)


    def set_target(self):

        target_path = (
            # by default we write into the cache,
            # and later copy to `outfile` on demand;
            # see the `copy_file` method
            self.cache_file_name
                if self.write_cache else
            # otherwise we do not write to the cache
            # but only to the outfile if it is set
            self.outfile
                if self.outfile is not None else
            # if both are disabled, we discard the downloaded data
            os.devnull
        )

        self.target = open(target_path, 'wb')

        self.curl.setopt(self.curl.WRITEFUNCTION, self.target.write)


    def _req_headers_list(self):

        self.req_headers = self.req_headers or []

        if isinstance(self.req_headers, dict):

            self.req_headers = [
                '%s: %s' % hdr
                for hdr in self.req_headers.items()
            ]


    def set_req_headers(self):

        self.init_request()

        if self.override_post:
            self._log('Overriding HTTP method.')

            self.req_headers.append('X-HTTP-Method-Override: GET')
        self.curl.setopt(
            self.curl.HTTPHEADER,
            [h.encode('ascii') for h in self.req_headers]
        )

    def set_resp_headers(self):

        self.resp_headers = []
        self.curl.setopt(self.curl.HEADERFUNCTION, self.resp_headers.append)

    def set_debug(self):

        if self.debug:
            self.curl.setopt(pycurl.VERBOSE, 1)
            self.curl.setopt(pycurl.DEBUGFUNCTION, self.print_debug_info)


    def set_compressed(self):

        if self.compressed:
            self.curl.setopt(pycurl.ENCODING, 'gzip, deflate')


    def curl_setup(self, url = False):

        self.curl_init(url = url)
        self.curl_progress_setup()
        self.set_target()
        self.set_debug()
        self.set_post()
        self.set_binary_data()
        self.set_req_headers()
        self.set_resp_headers()


    def curl_call(self):

        self._log('Setting up and calling pycurl.')

        for attempt in xrange(self.retries):

            try:

                if self.debug:

                    self.print_debug_info(
                        'INFO',
                        'pypath.curl.Curl().curl_call() :: attempt #%u'
                        % attempt
                    )

                if attempt > 0:

                    # apparently we have to set it again
                    # before each perform
                    self.set_binary_data()

                self.curl.perform()

                self.target.flush()

                if (
                    self.target.name != os.devnull and
                    os.path.exists(self.target.name) and
                    os.stat(self.target.name).st_size == 0 and
                    self.empty_attempt_again
                ):

                    self._log(
                        'Empty file retrieved, attempting downlad again'
                    )

                    continue

                if self.url.startswith('http'):
                    self.status = self.curl.getinfo(pycurl.HTTP_CODE)
                    if self.status == 200:
                        self.terminate_progress()
                        break

                if self.url.startswith('ftp'):
                    self.status == 500
                    for h in self.resp_headers:
                        if h[:3] == b'226':
                            self.status = 200
                            self.terminate_progress()
                            break
                    if self.status == 200:
                        break

            except pycurl.error as e:

                err_str = str(e.args)
                err_ssl_missing_eof = (
                    "(56, 'OpenSSL SSL_read: SSL_ERROR_SYSCALL, errno 0')"
                )

                self.print_debug_info('ERROR', f'PycURL error: {err_str}')

                if err_str == err_ssl_missing_eof:

                    self._log(
                        'Ignoring OpenSSL error about missing EOF '
                        '(truncated transmission); in some cases this '
                        'error is triggered even at successful downloads.'
                    )
                    self.status = 200

                else:

                    self.status = 500
                    if self.progress is not None:
                        self.progress.terminate(status = 'failed')
                        self.progress = None

        self.curl.close()
        self.target.close()

        if self.status != 200:
            self.download_failed = True
            self._log('Download error: HTTP %u' % self.status)

        if (
            self.target.name != os.devnull and
            os.path.exists(self.target.name) and
            os.stat(self.target.name).st_size == 0 and
            self.status != 302
        ):

            self.status = 500
            self.download_failed = True
            self._log('Download error: empty file retrieved.')

        if (
            (
                self.status >= 400 or
                self.download_failed
            )
        ):

            with open(self.target.name, 'rb') as fp:

                contents = fp.read(5000)

            try:

                contents = contents.decode('utf8')

            except UnicodeDecodeError:

                contents = str(contents)

            self._log('First 5000 bytes of response: %s' % contents)

            if not self.keep_failed:

                self._log('Download failed, removing the resulted file.')
                self.remove_target()


    def remove_target(self):

            self._log('Removing file: `%s`' % self.target.name)
            self.target.close()

            if os.path.exists(self.target.name):

                try:

                    os.remove(self.target.name)

                except PermissionError:

                    self._log(
                        'Could not remove `%s`, permission denied.' %
                        self.target.name
                    )


    def progress_setup(self):
        if not self.silent and self.progress is None and not self.debug:
            self.progress = progress.Progress(
                name = self.title, interval = 1, status = 'initializing curl')


    def curl_progress_setup(self):
        if self.progress is not None:
            self.curl.setopt(pycurl.NOPROGRESS, 0)
            if hasattr(pycurl, 'XFERINFOFUNCTION'):
                self.curl.setopt(pycurl.XFERINFOFUNCTION, self.update_progress)
            elif hasattr(pycurl, 'PROGRESSFUNCTION'):
                self.curl.setopt(pycurl.PROGRESSFUNCTION, self.update_progress)

    def _bytes_to_unicode(self, string, encoding = None):

        if type(string) is unicode:
            return string
        if encoding is not None:
            return string.decode(encoding)
        else:
            try:
                return string.decode('utf-8')

            except UnicodeDecodeError:
                try:
                    return string.decode('iso-8859-1')
                except:
                    self.print_debug_info('ERROR', 'String decoding error')
                    return u''

    def unicode2bytes(self, string, encoding = None):
        if type(string) is bytes:
            return string
        if encoding is not None:
            return string.encode(encoding)
        else:
            try:
                return string.encode('ascii')
            except UnicodeEncodeError:
                try:
                    return string.encode('utf-8')
                except:
                    self.print_debug_info('ERROR', 'String encoding error')
                    return b''

    def bytes_prefix(self, b):
        if b > 1000000000:
            return (b / 1000000000.0, u'GB')
        elif b > 1000000:
            return (b / 1000000.0, u'MB')
        elif b > 1000:
            return (b / 1000.0, u'kB')
        else:
            return (float(b), u'B')

    def get_headers(self):

        self.resp_headers_dict = {}

        if hasattr(self, 'resp_headers'):

            for header_line in self.resp_headers:

                header_line = self._bytes_to_unicode(header_line)

                if ':' not in header_line:

                    continue

                name, value = header_line.split(':', 1)
                name = name.strip()
                value = value.strip()
                name = name.lower()
                self.resp_headers_dict[name] = value

    def guess_encoding(self):

        if self.encoding is None:

            if not self.use_cache:

                if 'content-type' in self.resp_headers:

                    content_type = self.resp_headers['content-type'].lower()
                    match = re.search(r'charset = (\S+)', content_type)
                    if match:
                        self.encoding = match.group(1)

        if self.encoding is None:

            self.encoding = 'utf-8'

    def get_type(self):
        self.multifile = False
        if self.filename[-3:].lower() == 'zip' or self.compr == 'zip':
            self.type = 'zip'
            self.multifile = True
        elif self.filename[-3:].lower() == 'tgz' or \
                self.filename[-6:].lower() == 'tar.gz' or \
                self.compr == 'tgz' or self.compr == 'tar.gz':
            self.type = 'tgz'
            self.multifile = True
        elif self.filename[-2:].lower() == 'gz' or self.compr == 'gz':
            self.type = 'gz'
        else:
            self.type = 'plain'

    def get_jsessionid(self):

        self.jsessionid = [u'']
        rejsess = re.compile(r'.*(JSESSIONID\s?=\s?\w*)')

        for hdr in self.resp_headers:

            jsess = rejsess.findall(hdr.decode('utf-8'))

            if len(jsess) > 0:

                self.jsessionid = [u'Cookie: %s' % jsess[0]]

        return self.jsessionid

    def update_progress(self, download_total, downloaded, upload_total,
                        uploaded):
        if self.progress is not None:
            self.total = self.bytes_prefix(download_total)
            self.done = self.bytes_prefix(downloaded)
            msg = u'%.02f%s/%.02f%s' % \
                (self.done[0], self.done[1], self.total[0], self.total[1])
            self.progress.set_total(float(download_total))
            self.progress.set_done(float(downloaded))
            self.progress.step(step = 0, msg = msg, status = 'downloading')

    def terminate_progress(self):
        if self.progress is not None:
            self.progress.terminate(status = '%.02f%s downloaded' %
                                    (self.total[0], self.total[1]))
            self.progress = None


    def init_request(self):

        if self.init_url is not None:

            if self.progress is not None:

                self.progress.set_status('requesting cookie')

            self.init_curl = Curl(
                self.init_url,
                silent = True,
                debug = self.debug,
                cache = self.init_use_cache,
                req_headers = self.req_headers,
                follow = False,
            )
            headers = getattr(self.init_curl, self.init_fun)()
            self.req_headers.extend(headers)


    # caching:

    def init_cache(self):

        self.get_hash()
        self.cache_dir_exists()
        self.get_cache_file_name()

        self._log('Cache file path: `%s`' % self.cache_file_name)

        self.select_cache_file()


    def get_hash(self):

        if isinstance(self.cache, str):

            return

        self.post_str = (
            ''
                if self.post is None else
            (
                '?' + '&'.join(sorted([
                    i[0] + ' = ' + i[1]
                    for i in iteritems(self.post)
                ]))
            )
        )

        if self.binary_data:
            bindata = str(self.binary_data)
        else:
            bindata = ''

        self.urlmd5 = hashlib.md5(
            self.unicode2bytes('%s%s%s' % \
                (self.url, self.post_str, bindata))).hexdigest()


    def cache_dir_exists(self):

        self.cache_dir = cache_mod.get_cachedir(self.cache_dir)


    def get_cache_file_name(self):

        self.cache_file_name = (
            self.cache
                if isinstance(self.cache, str) else
            os.path.join(
                os.getcwd(),
                self.cache_dir,
                self.replace_forbidden('%s-%s' % (self.urlmd5, self.filename))
            )
        )


    @classmethod
    def cache_path(self, **kwargs) -> str:
        """
        Returns the cache path without performing download or creating file.

        Args:
            kwargs:
                Arguments to `Curl`.
        """

        kwargs.update({'setup': False, 'call': False, 'process': False})

        return Curl(**kwargs).cache_file_name


    @classmethod
    def replace_forbidden(cls, name: str, repl: str = '_') -> str:
        """
        Replaces the characters that are forbidden in certain file systems.

        The slash is forbidden in Unix, while many other characters in Windows
        environments.
        """

        return cls.FORBIDDEN_CHARS.sub(repl, name)


    def delete_cache_file(self):

        if os.path.exists(self.cache_file_name):
            self.print_debug_info('INFO',
                                  'CACHE FILE = %s' % self.cache_file_name)
            self.print_debug_info('INFO', 'DELETING CACHE FILE')
            os.remove(self.cache_file_name)
            self.use_cache = False
        else:
            self.print_debug_info('INFO',
                                  'CACHE FILE = %s' % self.cache_file_name)
            self.print_debug_info('INFO', 'CACHE FILE DOES NOT EXIST')

    def select_cache_file(self):

        self.use_cache = False

        if type(CACHE) is bool:

            self.cache = CACHE

        if (
            self.cache and
            os.path.exists(self.cache_file_name) and
            # if the cache file is empty
            # try to download again
            os.stat(self.cache_file_name).st_size > 0
        ):

            self._log('Cache file found, no need for download.')

            self.use_cache = True

    def show_cache(self):

        self.print_debug_info('INFO', 'URL = %s' % self.url)
        self.print_debug_info('INFO', 'CACHE FILE = %s' % self.cache_file_name)
        self.print_debug_info(
            'INFO', 'Using cache: %s; cache file exists: %s' %
            (self.cache, os.path.exists(self.cache_file_name)))


    # open files:

    def transcode(self):

        if not self.use_cache and self.type == 'plain':

            self.guess_encoding()

            if self.encoding is not None and self.encoding != 'utf-8':

                self._log(
                    'Converting encoding from `%s` '
                    'to `utf-8`.' % self.encoding
                )

                tmp_file_name = os.path.join(os.getcwd(), self.cache_dir,
                                             'transcoding.tmp.txt')
                os.rename(self.cache_file_name, tmp_file_name)
                if self.progress is not None:
                    self.print_status('Converting %s encoded data to utf-8' %
                                      self.encoding)
                with open(tmp_file_name, 'rb') as tmp_file:
                    with open(self.cache_file_name, 'wb') as cache_file:
                        for line in tmp_file:
                            cache_file.write(
                                line.decode(self.encoding or 'utf-8').encode(
                                    'utf-8'))
                os.remove(tmp_file_name)
                self.encoding = 'utf-8'


    def copy_file(self):

        self.transcode()

        if self.outfile is not None and self.outfile != self.cache_file_name:

            if self.write_cache:

                self._log(
                    'Copying file `%s` to `%s`.' % (
                        self.cache_file_name,
                        self.outfile,
                    )
                )
                shutil.copy(self.cache_file_name, self.outfile)

            else:
                self._log(
                    'Moving file `%s` to `%s`.' % (
                        self.cache_file_name,
                        self.outfile,
                    )
                )
                os.rename(self.cache_file_name, self.outfile)

        else:

            self.outfile = self.cache_file_name


    def process_file(self):

        self.guess_encoding()
        self.get_type()
        self.copy_file()
        self.open_file()
        self.extract_file()
        self.decode_result()
        self.report_ready()


    def open_file(self):

        if not self.silent:

            self.print_status('Opening file `%s`' % self.outfile)

        super(Curl, self).__init__(self.outfile, extract = False)


    def close(self):
        """
        Closes all file objects.
        """

        if type(self.result) is dict:
            for fp in self.result.values():
                if hasattr(fp, 'close'):
                    fp.close()
        self.fileobj.close()


    def extract_file(self):

        if not self.silent:
            self._log('Extracting data from file type `%s`' % self.type)
        self.extract()


    def decode_result(self):

        if self.progress is not None:

            self._log(
                'Decoding `%s` encoded data' % (self.encoding or 'utf-8')
            )

        def _decode_result(content):

            try:
                if isinstance(content, str):

                    return content

                else:

                    return content.decode(self.encoding or 'utf-8')

            except:

                self.print_debug_info(
                    'WARNING',
                    'Failed '
                    'decoding downloaded bytes content with encoding %s. '
                    'Result might be of type bytes' %
                    (self.encoding or 'utf-8')
                )
                return content

        if not self.large:

            if type(self.result) is dict:

                for name, content in iteritems(self.result):

                    self.result[name] = _decode_result(content)

            else:

                self.result = _decode_result(self.result)


    def get_result_type(self):

        if type(self.result) is dict:
            if len(self.result):
                self.result_type = 'dict of %s' % (
                    'byte arrays'
                    if type(next(iter(self.result.values()))) is bytes else
                    'unicode strings'
                    if type(next(iter(self.result.values()))) is unicode else
                    'file objects')
            else:
                self.result_type = 'empty dict'
        else:
            self.result_type = '%s' % (
                'byte array'
                if type(self.result) is bytes else 'unicode string'
                if type(self.result) is unicode else 'file object')


    def report_ready(self):

        self.get_result_type()

        if not self.silent:

            self._log(
                'File at `%s` successfully retrieved. '
                'Resulted file type `%s, %s`. '
                'Local file at `%s`.' % (
                    self.url,
                    'plain text'
                        if self.type == 'plain' else
                    '%s extracted data' % self.type,
                    self.result_type,
                    self.outfile,
                )
            )


    def print_status(self, status):

        if self.progress is not None:
            self.terminate_progress()
        if self.debug:
            self.print_debug_info('INFO', status)
        elif not self.silent:
            self._log(status)


    # sftp part:

    def sftp_url(self):

        if self.sftp_host is not None:
            self.sftp_filename = self.url
            self.url = '%s%s' % (self.sftp_host, self.sftp_filename)


    def sftp_call(self):

        self.sftp_success = self.sftp_download()
        if self.sftp_success:
            self.status = 200
        else:
            self.status = 501


    def ask_passwd(self, use_passwd_file = True):

        if use_passwd_file and os.path.exists(self.sftp_passwd_file):
            with open(self.sftp_passwd_file, 'r') as f:
                self.sftp_user = f.readline().strip()
                self.sftp_passwd = f.readline().strip()
            return None
        sys.stdout.write(self.sftp_ask)
        sys.stdout.flush()
        while True:
            self.user = input('\n\tUsername: ')
            self.passwd = input(
                '\tPassword (leave empty if no password needed): ')
            correct = input('Are these details correct? '
                                'User: `%s`, password: `%s` [Y/n]\n' %
                                (self.user, self.passwd))
            if correct.lower().strip() not in ['', 'y', 'yes']:
                continue
            save = input(
                'Do you wish to save your login details unencripted\n'
                'to the following file, so you don\'t '
                'need to enter them next time? File: %s\n'
                'Save login details [Y/n]' %
                self.sftp_passwd_file
            )
            break
        if save.lower().strip() in ['', 'y', 'yes']:
            with open(self.sftp_passwd_file, 'w') as f:
                f.write('%s\n%s' % (self.user, self.passwd))


    def sftp_download(self):

        self.sftp_ask = (
            'Please enter your login details for %s\n' % self.host
                if self.sftp_ask is None else
            self.sftp_ask
        )
        self.sftp_passwd_file = (
            os.path.join('cache', '%s.login' % self.sftp_host)
                if self.sftp_passwd_file is None else
            self.sftp_passwd_file
        )
        if self.sftp_user is None:
            self.ask_passwd()
        while True:

            self.sftp_passwd = self.sftp_passwd or None
            cnopts = pysftp.CnOpts()
            cnopts.hostkeys = None

            with pysftp.Connection(
                host = self.sftp_host,
                username = self.sftp_user,
                password = self.sftp_passwd,
                port = self.sftp_port,
                cnopts = cnopts
            ) as con:

                try:

                    con.get(self.sftp_filename, self.cache_file_name)
                    break

                except IOError:

                    msg = 'Failed to get %s from %s\n'\
                        'Try again (1) || Enter new login details (2) '\
                        '|| Cancel (3) ?\n' % (
                            self.sftp_filename, self.sftp_host)
                    whattodo = input(msg)
                    if '1' in whattodo:
                        continue
                    if '2' in whattodo:
                        self.ask_passwd(use_passwd_file = False)
                        continue
                    if '3' in whattodo:
                        return False
        return True
