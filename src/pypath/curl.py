#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

#
# this module makes possible
# dynamic data integration, downloads
# files from various resources, in standard
# or non-standard text based and xml formats,
# processes them, sometimes parses html
#

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import imp
import sys
import os
import shutil
import struct

import pycurl
try:
    from cStringIO import StringIO
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
    urllib2 = urllib.request

import httplib2
try:
    import urlparse
except:
    # this works seemless in Py3:
    urlparse = urllib.parse

if not hasattr(urllib, 'quote'):
    _urllib = urllib
    urllib = _urllib.parse

try:
    import pysftp
except:
    sys.stdout.write('''\n\t:: Module `pyfstp` not available. 
        Only downloading of a small number of resources 
        relies on this module.
        Please install by PIP if it is necessary for you.
    ''')
import codecs
import gzip
import zipfile
import tarfile
import hashlib

try:
    from fabric.network import connect, HostConnectionCache
    from fabric.state import env
except:
    sys.stdout.write('No `fabric` available.\n')
    sys.stdout.flush()

from contextlib import closing

import pypath.progress as progress

if 'long' not in __builtins__:
    long = int

if 'unicode' not in __builtins__:
    unicode = str

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'
CACHE = None

show_cache = False

class cache_on(object):
    
    def __init__(self):
        pass
    
    def __enter__(self):
        global CACHE
        self._store_cache = globals()['CACHE']
        CACHE = True
    
    def __exit__(self, exception_type, exception_value, traceback):
        global CACHE
        if exception_type is not None:
            sys.stdout.write('%s, %s, %s\n' % \
                (str(exception_type), str(exception_value), str(traceback)))
            sys.stdout.flush()
        CACHE = self._store_cache

class cache_off(object):
    
    def __init__(self):
        pass
    
    def __enter__(self):
        global CACHE
        self._store_cache = globals()['CACHE']
        CACHE = False
    
    def __exit__(self, exception_type, exception_value, traceback):
        global CACHE
        if exception_type is not None:
            sys.stdout.write('%s, %s, %s\n' % \
                (str(exception_type), str(exception_value), str(traceback)))
            sys.stdout.flush()
        CACHE = self._store_cache

class RemoteFile(object):
    
    def __init__(self, filename, user, host, passwd, port = 22, sep = '\t', 
        header = True, rownames = True):
        for key, val in iteritems(locals()):
            setattr(self, key, val)
        env.keepalive = 60
        env.connection_attempts = 5
        env.password = self.passwd
    
    def wcl(self):
        with closing(connect(self.user, self.host, self.port, \
            HostConnectionCache())) as ssh:
            stdin, stdout, stderr = ssh.exec_command('wc -l %s'%self.filename)
            return int(stdout.readlines()[0].split()[0]) - (1 if self.header else 0)
    
    def rowns(self):
        with closing(connect(self.user, self.host, self.port, \
            HostConnectionCache())) as ssh:
            stdin, stdout, stderr = ssh.exec_command(
                'awk \'BEGIN{FS="%s"}{print $1}\' %s%s' % \
                (self.sep, self.filename, '' if not self.header else ' | tail -n +2'))
            return [x.strip() for x in stdout.readlines()]
    
    def open(self, return_header = True):
        with closing(connect(self.user, self.host, self.port, \
            HostConnectionCache())) as ssh:
            with closing(ssh.open_sftp()) as sftp:
                with closing(sftp.open(self.filename)) as f:
                    if not return_header:
                        line = f.readline()
                    for line in f:
                        yield line

class Curl(object):
    
    def __init__(self,
        url, silent = True, get = None,
        post = None, req_headers = None, cache = True,
        debug = False, outf = None, compr = None, encoding = None,
        files_needed = None, timeout = 300, init_url = None, 
        init_fun = 'get_jsessionid', follow = True, large = False,
        override_post = False, init_headers = False,
        return_headers = False, binary_data = None,
        write_cache = True, force_quote = False,
        sftp_user = None, sftp_passwd = None, sftp_passwd_file = '.secrets',
        sftp_port = 22, sftp_host = None, sftp_ask = None,
        setup = True, call = True, process = True,
        retries = 3, cache_dir = 'cache'):
        
        self.result = None
        self.download_failed = False
        self.status = 0
        self.large = large
        self.silent = silent
        self.debug = debug
        self.url = url
        self.get = get
        self.force_quote = force_quote
        self.process_url()
        self.url_fix()
        self.set_get()
        self.compr = compr
        self.get_type()
        self.progress = None
        
        self.encoding = encoding
        self.files_needed = files_needed
        
        self.follow_http_redirect = follow
        self.timeout = timeout
        self.override_post = override_post
        self.retries = retries
        self.req_headers = req_headers or []
        self.post = post
        self.get = get
        self.binary_data = binary_data
        
        self.cache_dir = cache_dir
        self.cache = cache
        self.init_cache()
        self.write_cache = write_cache
        self.outfile = outf
        
        self.init_url = init_url
        self.init_fun = init_fun
        
        self.sftp_host = sftp_host
        self.sftp_ask = sftp_ask
        self.sftp_port = sftp_port
        self.sftp_passwd = sftp_passwd
        self.sftp_user = sftp_user
        self.sftp_passwd_file = sftp_passwd_file
        
        if not self.use_cache:
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
            sys.stdout.write('\t:: Loading data from cache '\
                'previously downloaded from %s\n' % self.domain)
            sys.stdout.flush()
        if process and not self.download_failed:
            self.process_file()
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def print_debug_info(self, msg):
        msg = self.bytes2unicode(msg)
        sys.stdout.write('\n\t%s\n' % msg)
        sys.stdout.flush()
    
    def process_url(self):
        self.domain = self.url.replace('https://', '').replace('http://', '').\
            replace('ftp://', '').split('/')[0]
        self.filename = self.url.split('/')[-1].split('?')[0]
    
    def is_quoted(self, string):
        '''
        From http://stackoverflow.com/questions/1637762/test-if-string-is-url-encoded-in-php
        '''
        test = string
        while(urllib.unquote(test) != test):
            test = urllib.unquote(test)
        return urllib.quote(test, '/%') == string or urllib.quote(test) == string

    def is_quoted_plus(self, string):
        test = string
        while(urllib.unquote_plus(test) != test):
            test = urllib.unquote_plus(test)
        return urllib.quote_plus(test, '&=') == string or urllib.quote_plus(test) == string

    def url_fix(self, charset = 'utf-8'):
        """
        From http://stackoverflow.com/a/121017/854988
        """
        if type(self.url) is bytes:
            self.url = self.bytes2unicode(self.url, encoding = charset)
        scheme, netloc, path, qs, anchor = urlparse.urlsplit(self.url)
        if self.force_quote or not self.is_quoted(path):
            path = urllib.quote(path, '/%')
        if self.force_quote or not self.is_quoted_plus(qs):
            qs = urllib.quote_plus(qs, '&=')
        self.url= urlparse.urlunsplit((scheme, netloc, path, qs, anchor))
    
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
    
    def set_title(self):
        if self.title is None:
            self.title = 'Downloading `%s` from %s' % \
                (self.filename, self.domain)
    
    def set_post(self):
        if type(self.post) is dict:
            self.postfields = urllib.urlencode(self.post)
            self.curl.setopt(self.curl.POSTFIELDS, self.postfields)
            self.curl.setopt(self.curl.POST, 1)
        else:
            self.postfields = None
    
    def set_get(self):
        if self.get is not None:
            self.qs = '&'.join(
                map(
                    lambda param:
                        '%s=%s' % (param[0], param[1]),
                    map(
                        lambda param:
                            (
                                urllib.quote_plus(param[0]),
                                urllib.quote_plus(param[1])
                            ),
                        iteritems(self.get)
                    )
                )
            )
            self.url = '%s%s%s' % (
                self.url,
                '&' if '?' in self.url else '?',
                self.qs
            )
    
    def set_binary_data(self):
        if self.binary_data:
            self.binary_data_size = os.path.getsize(self.binary_data)
            self.binary_data_file = open(self.binary_data, 'rb')
            self.curl.setopt(c.POST, 1)
            filesize = os.path.getsize(self.binary_data)
            self.curl.setopt(pycurl.POSTFIELDSIZE, filesize)
            self.curl.setopt(pycurl.READFUNCTION, self.binary_data_file.read)
            self.curl.setopt(pycurl.CUSTOMREQUEST, 'POST')
            self.curl.setopt(pycurl.POSTREDIR, 3)
    
    def curl_init(self, url = False):
        self.curl = pycurl.Curl()
        self.set_url(url = url)
        self.curl.setopt(self.curl.FOLLOWLOCATION, self.follow_http_redirect)
        self.curl.setopt(self.curl.CONNECTTIMEOUT, self.timeout)
    
    def set_url(self, url = False):
        self.curl.setopt(self.curl.URL, url or self.url)
    
    def set_target(self):
        self.target = open(self.cache_file_name, 'wb')
        self.curl.setopt(self.curl.WRITEFUNCTION, self.target.write)
        
    def set_req_headers(self):
        if self.override_post:
            self.req_headers.append('X-HTTP-Method-Override: GET')
        self.curl.setopt(self.curl.HTTPHEADER, self.req_headers)
    
    def set_resp_headers(self):
        self.resp_headers = []
        self.curl.setopt(self.curl.HEADERFUNCTION, self.resp_headers.append)
    
    def set_debug(self):
        if self.debug:
            self.curl.setopt(pycurl.VERBOSE, 1)
            self.curl.setopt(pycurl.DEBUGFUNCTION, self.print_debug_info)
    
    def curl_setup(self, url = False):
        self.curl_init(url = url)
        self.curl_progress_setup()
        self.set_target()
        self.set_req_headers()
        self.set_resp_headers()
        self.set_debug()
        self.set_post()
        self.set_binary_data()
    
    def curl_call(self):
        for attempt in xrange(self.retries):
            try:
                if self.debug:
                    self.print_debug_info(
                        'pypath.curl.Curl().curl_call() :: attempt #%u' % i)
                self.curl.perform()
                if self.url.startswith('http'):
                    self.status = self.curl.getinfo(pycurl.HTTP_CODE)
                    if self.status == 200:
                        self.terminate_progress()
                        break
                if self.url.startswith('ftp'):
                    self.status == 500
                    for h in self.resp_headers:
                        if h.startswith('226'):
                            self.status = 200
                            self.terminate_progress()
                            break
            except pycurl.error as e:
                self.status = 500
                if self.progress is not None:
                    self.progress.terminate(status = 'failed')
                    self.progress = None
                self.print_debug_info('PycURL error: %u, %s' % e.args)
        if self.status != 200:
            self.download_failed = True
        self.curl.close()
        self.target.close()
    
    def progress_setup(self):
        if not self.silent and self.progress is None and not self.debug:
            self.progress = progress.Progress(name = self.title, interval = 1,
                status = 'initializing curl')
    
    def curl_progress_setup(self):
        if self.progress is not None:
            self.curl.setopt(pycurl.NOPROGRESS, 0)
            if hasattr(pycurl, 'XFERINFOFUNCTION'):
                self.curl.setopt(pycurl.XFERINFOFUNCTION, self.update_progress)
            elif hasattr(pycurl, 'PROGRESSFUNCTION'):
                self.curl.setopt(pycurl.PROGRESSFUNCTION, self.update_progress)
    
    def bytes2unicode(self, string, encoding = None):
        if type(string) is unicode:
            return string
        if encoding is not None:
            return string.decode(encoding)
        else:
            try:
                return string.decode('ascii')
            except UnicodeDecodeError:
                try:
                    return string.decode('utf-8')
                except:
                    self.print_debug_info('String decoding error')
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
                    self.print_debug_info('String encoding error')
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
        for header_line in self.resp_headers:
            header_line = self.bytes2unicode(header_line)
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
                    match = re.search('charset=(\S+)', content_type)
                    if match:
                        self.encoding = match.group(1)
    
    def get_jsessionid(self):
        self.jsessionid = [u'']
        rejsess = re.compile(r'.*(JSESSIONID=[A-Z0-9]*)')
        for hdr in self.resp_headers_dict.values():
            jsess = rejsess.findall(hdr)
            if len(jsess) > 0:
                self.jsessionid = [u'Cookie: %s' % jsess[0]]
        return self.jsessionid
    
    def update_progress(self,
        download_total, downloaded,
        upload_total, uploaded):
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
            self.progress.terminate(status = '%.02f%s downloaded' % \
                (self.total[0], self.total[1]))
            self.progress = None
    
    def init_request(self):
        if self.init_url is not None:
            if self.progress is not None:
                self.progress.set_status('requesting cookie')
            self.init_curl = Curl(self.init_url,
                silent = True, debug = self.debug)
            headers = getattr(self.init_curl, self.init_fun)()
            self.req_headers.extend(headers)
    
    # caching:
    
    def init_cache(self):
        self.get_hash()
        self.cache_dir_exists()
        self.get_cache_file_name()
        self.select_cache_file()
    
    def get_hash(self):
        self.post_str = '' if self.post is None else \
            '?' + '&'.join(sorted([i[0]+'='+i[1] \
            for i in iteritems(self.post)]))
        self.urlmd5 = hashlib.md5(self.unicode2bytes(
            '%s%s' % \
            (self.url, self.post_str)
        )).hexdigest()
    
    def cache_dir_exists(self):
        if not os.path.exists(os.path.join(os.getcwd(), self.cache_dir)):
            os.mkdir(os.path.join(os.getcwd(), self.cache_dir))
    
    def get_cache_file_name(self):
        self.cache_file_name = os.path.join(
            os.getcwd(),
            self.cache_dir,
            '%s-%s' % (self.urlmd5, self.filename)
        )
    
    def delete_cache_file(self):
        if os.path.exists(self.cache_file_name):
            os.remove(self.cache_file_name)
    
    def select_cache_file(self):
        self.use_cache = False
        if type(CACHE) is bool:
            self.cache = CACHE
        if self.cache and os.path.exists(self.cache_file_name):
            self.use_cache = True
    
    def show_cache(self):
        self.print_debug_info('URL = %s' % self.url)
        self.print_debug_info('CACHE FILE = %s' % self.cache_file_name)
        self.print_debug_info('Using cache: %s; cache file exists: %s' % \
            (self.cache, os.path.exists(self.cache_file_name)))
    
    # open files:
    
    def transcode(self):
        if not self.use_cache and self.type == 'plain':
            self.guess_encoding()
            if self.encoding is not None and self.encoding != 'utf-8':
                tmp_file_name = os.path.join(os.getcwd(),
                    self.cache_dir, 'transcoding.tmp.txt')
                os.rename(self.cache_file_name, tmp_file_name)
                if self.progress is not None:
                    self.print_status('Converting %s encoded data to utf-8' % \
                        self.encoding)
                with open(tmp_file_name, 'rb') as tmp_file:
                    with open(self.cache_file_name, 'wb') as cache_file:
                        for line in tmp_file:
                            cache_file.write(
                                line.decode(self.encoding or 'utf-8').encode('utf-8')
                            )
                os.remove(tmp_file_name)
                self.encoding = 'utf-8'
    
    def copy_file(self):
        self.transcode()
        if self.outfile is not None and self.outfile != self.cache_file_name:
            if self.write_cache:
                shutil.copy(self.cache_file_name, self.outfile)
            else:
                os.rename(self.cache_file_name, self.outfile)
        else:
            self.outfile = self.cache_file_name
    
    def process_file(self):
        self.copy_file()
        self.open_file()
        self.extract_file()
        self.decode_result()
        self.report_ready()
    
    def open_file(self):
        if not self.silent:
            self.print_status('Opening file `%s`' % self.outfile)
        self.fileobj = open(self.outfile, 'rb')
    
    def extract_file(self):
        if not self.silent:
            self.print_status('Extracting %s data' % self.type)
        getattr(self, 'open_%s' % self.type)()
    
    def open_tgz(self):
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
                    self.files_multipart[m.name] = this_file.read()
                    this_file.close()
        if not self.large:
            self.tarfile.close()
        self.result = self.files_multipart
    
    def open_gz(self):
        self.fileobj.seek(-4, 2)
        self.size = struct.unpack('I', self.fileobj.read(4))[0]
        self.fileobj.seek(0)
        self.gzfile = gzip.GzipFile(fileobj = self.fileobj, mode = 'rb')
        #try:
        if self.large:
            self.result = self.gzfile
        else:
            self.result = self.gzfile.read()
            self.gzfile.close()
        #except:
        #    self.print_status('Error at extracting gzip file')
    
    def open_zip(self):
        self.files_multipart = {}
        self.sizes = {}
        self.zipfile = zipfile.ZipFile(self.fileobj, 'r')
        self.members = self.zipfile.namelist()
        for i, m in enumerate(self.members):
            self.sizes[m] = self.zipfile.filelist[i].file_size
            if self.files_needed is None or m in self.files_needed:
                this_file = self.zipfile.open(m)
                if self.large:
                    self.files_multipart[m] = this_file
                else:
                    self.files_multipart[m] = this_file.read()
                    this_file.close()
        if not self.large:
            self.zipfile.close()
        self.result = self.files_multipart
    
    def open_plain(self):
        self.size = os.path.getsize(self.fileobj.name)
        if self.large:
            self.result = self.fileobj
        else:
            self.result = self.fileobj.read()
            self.fileobj.close()
    
    def decode_result(self):
        if self.progress is not None:
            self.print_status('Decoding %s encoded data' % \
                (self.encoding or 'utf-8'))
        def _decode_result(content):
            try:
                return content.decode(self.encoding or 'utf-8')
            except:
                self.print_debug_info('Failed '\
                    'decoding downloaded bytes content with encoding %s. '\
                    'Result might be of type bytes' % (self.encoding or 'utf-8'))
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
                    'byte arrays' \
                        if type(next(iter(self.result.values()))) is bytes \
                    else 'unicode strings' \
                        if type(next(iter(self.result.values()))) is unicode \
                    else 'file objects'
                )
            else:
                self.result_type = 'empty dict'
        else:
            self.result_type = '%s' % (
                'byte array' if type(self.result) is bytes else \
                'unicode string' if type(self.result) is unicode else \
                'file object'
            )
    
    def report_ready(self):
        self.get_result_type()
        if not self.silent:
            self.print_status(
                'Ready. Resulted `%s` of type %s. \n'\
                '\t:: Local file at `%s`.' % \
                (
                    'plain text' if self.type == 'plain' \
                        else '%s extracted data' % self.type,
                    self.result_type,
                    self.outfile
                )
            )
            sys.stdout.write('\n')
            sys.stdout.flush()
    
    def print_status(self, status):
        if self.progress is not None:
            self.terminate_progress()
        if self.debug:
            self.print_debug_info(status)
        elif not self.silent:
            sys.stdout.write('\r%s' % (' ' * 150))
            sys.stdout.write('\r\t:: %s' % status)
            sys.stdout.flush()
    
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
            self.user = raw_input('\n\tUsername: ')
            self.passwd = raw_input('\tPassword (leave empty if no password needed): ')
            correct = raw_input('Are these details correct? '\
                'User: `%s`, password: `%s` [Y/n]\n' % (self.user, self.passwd))
            if correct.lower().strip() not in ['', 'y', 'yes']:
                continue
            save = raw_input('Do you wish to save your login details unencripted\n'\
                'to the following file, so you don\'t need to enter them next '\
                'time? File: %s\nSave login details [Y/n]' % self.sftp_passwd_file)
            break
        if save.lower().strip() in ['', 'y', 'yes']:
            with open(self.sftp_passwd_file, 'w') as f:
                f.write('%s\n%s' % (self.user, self.passwd))

    def sftp_download(self):
        self.sftp_ask = 'Please enter your login details for %s\n' % self.host \
            if self.sftp_ask is None else self.sftp_ask
        self.sftp_passwd_file = os.path.join('cache', '%s.login' % self.sftp_host) \
            if self.sftp_passwd_file is None else self.sftp_passwd_file
        if self.sftp_user is None:
            self.ask_passwd()
        while True:
            self.sftp_passwd = None \
                if self.sftp_passwd.strip() == '' \
                else self.sftp_passwd
            with pysftp.Connection(host = self.sftp_host, username = self.sftp_user,
                password = self.sftp_passwd, port = self.sftp_port) as con:
                try:
                    con.get(self.sftp_filename, self.cache_file_name)
                    break
                except IOError:
                    msg = 'Failed to get %s from %s\n'\
                        'Try again (1) || Enter new login details (2) '\
                        '|| Cancel (3) ?\n' % (self.sftp_filename, self.sftp_host)
                    whattodo = raw_input(msg)
                    if '1' in whattodo:
                        continue
                    if '2' in whattodo:
                        self.ask_passwd(use_passwd_file = False)
                        continue
                    if '3' in whattodo:
                        return False
        return True
