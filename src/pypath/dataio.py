#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
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
# this module makes possible a
# dynamic data integration, download 
# files from various resources, in standard
# or non-standard text based and xml formats,
# process them, sometimes parse html
#

import pycurl
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO
try:
    import cPickle as pickle
except:
    import pickle
import sys
import os
import re
import itertools
from collections import Counter
import urllib
import urllib2
import httplib2
import urlparse
import codecs
import gzip
import zipfile
import tarfile
import xlrd
import bs4
import xml.etree.cElementTree as ET
from lxml import etree
import hashlib
import time
import copy
import struct
import json
import webbrowser
from bioservices import WSDLService
from contextlib import closing
from fabric.network import connect, HostConnectionCache
from fabric.state import env
from xlrd import open_workbook
from xlrd.biffh import XLRDError

# from this module

import data_formats
import progress
import common
import intera
import reaction
import residues
import mapping
import seq as se

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

show_cache = False

class RemoteFile(object):
    
    def __init__(self, filename, user, host, passwd, port = 22, sep = '\t', 
        header = True, rownames = True):
        for key, val in locals().iteritems():
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

def is_quoted(string):
    '''
    From http://stackoverflow.com/questions/1637762/test-if-string-is-url-encoded-in-php
    '''
    test = string
    while(urllib.unquote(test) != test):
        test = urllib.unquote(test)
    return urllib.quote(test, '/%') == string or urllib.quote(test) == string

def is_quoted_plus(string):
    test = string
    while(urllib.unquote_plus(test) != test):
        test = urllib.unquote_plus(test)
    return urllib.quote_plus(test, '&=') == string or urllib.quote_plus(test) == string

def url_fix(s, charset='utf-8', force = False):
    """
    From http://stackoverflow.com/a/121017/854988
    """
    if isinstance(s, unicode):
        s = s.encode(charset, 'ignore')
    scheme, netloc, path, qs, anchor = urlparse.urlsplit(s)
    if force or not is_quoted(path):
        path = urllib.quote(path, '/%')
    if force or not is_quoted_plus(qs):
        qs = urllib.quote_plus(qs, '&=')
    return urlparse.urlunsplit((scheme, netloc, path, qs, anchor))

def print_debug_info(debug_type, debug_msg, truncate = 1000):
    sys.stdout.write("debug(%d): %s\n" % (debug_type, debug_msg[:truncate]))
    sys.stdout.flush()

#class Dataio(object):
#    
#    __init__(self,mapper=None):

def get_headers(header_list):
    headers = {}
    for header_line in header_list:
        if ':' not in header_line:
            continue
        name, value = header_line.split(':', 1)
        name = name.strip()
        value = value.strip()
        name = name.lower()
        headers[name] = value
    return headers

def get_jsessionid(headers):
    rejsess = re.compile(r'.*(JSESSIONID=[A-Z0-9]*)')
    for hdr in headers:
        jsess = rejsess.findall(hdr)
        if len(jsess) > 0:
            return ['Cookie: %s'%jsess[0]]

def get_xsessionid(headers):
    pass

def curl(url, silent = True, post = None, req_headers = None, cache = True, 
        debug = False, outf = None, compr = None, encoding = None, 
        files_needed = None, timeout = 300, init_url = None, 
        init_fun = 'get_jsessionid', follow = True, large = False,
        override_post = False, init_headers = False, 
        write_cache = True, force_quote = False):
    url = url_fix(url, force = force_quote)
    if init_url is not None:
        init_url = url_fix(init_url, force = force_quote)
    # either from cache or from download, we load the data into StringIO:
    multifile = False
    domain = url.replace('https://', '').replace('http://','').\
        replace('ftp://','').split('/')[0]
    # first try to find file in cache:
    if cache or write_cache:
        # outf param is to give a unique name to data
        # downloaded previously by post requests
        outf = outf if outf is not None else url.split('/')[-1].split('?')[0]
        poststr = '' if post is None else \
            '?' + '&'.join(sorted([i[0]+'='+i[1] for i in post.items()]))
        try:
            urlmd5 = hashlib.md5(url+poststr).hexdigest()
        except UnicodeEncodeError:
            urlmd5 = hashlib.md5(('%s%s' % (url, poststr)).encode('utf-8')).hexdigest()
        if not os.path.exists(os.path.join(os.getcwd(),'cache')):
            os.mkdir(os.path.join(os.getcwd(),'cache'))
        cachefile = os.path.join(os.getcwd(),'cache',urlmd5+'-'+outf)
        if show_cache:
            sys.stdout.write('\tFor URL %s\n' % url)
            sys.stdout.write('\tChache file is %s' % cachefile)
        usecache = True if os.path.exists(cachefile) and cache else False
        # load from cache:
        if usecache:
            if not silent:
                sys.stdout.write('\t:: Loading %s from cache, previously '\
                    'downloaded from %s\n'%(outf,domain))
                sys.stdout.flush()
            if large:
                result = open(cachefile, 'rb')
            else:
                with open(cachefile,'rb') as f:
                    result = StringIO()
                    result.write(f.read())
    else:
        usecache = False
    # if not found in cache, download with curl:
    if not usecache:
        headers = []
        if not init_url and large:
            result = open(cachefile, 'wb')
        else:
            result = StringIO()
        c = pycurl.Curl()
        if init_url:
            c.setopt(c.URL, init_url)
        else:
            try:
                c.setopt(c.URL, url)
            except:
                return url
        c.setopt(c.FOLLOWLOCATION, follow)
        c.setopt(c.CONNECTTIMEOUT, 15)
        c.setopt(c.TIMEOUT, timeout)
        if override_post:
            if req_headers is None: req_headers = []
            req_headers.append('X-HTTP-Method-Override: GET')
        if type(req_headers) is list:
            c.setopt(c.HTTPHEADER, req_headers)
        c.setopt(c.WRITEFUNCTION, result.write)
        c.setopt(c.HEADERFUNCTION, headers.append)
        # if debug is necessary:
        if debug:
            c.setopt(pycurl.VERBOSE, 1)
            c.setopt(pycurl.DEBUGFUNCTION, print_debug_info)
        if type(post) is dict:
            postfields = urllib.urlencode(post)
            c.setopt(c.POSTFIELDS, postfields)
            c.setopt(c.POST, 1)
        if not silent:
            sys.stdout.write('\t:: Downloading data from %s. Waiting for reply...' % \
                domain)
            sys.stdout.flush()
        for i in xrange(3):
            try:
                if debug:
                    sys.stdout.write('\t:: pypath.dataio.curl() :: attempt #%u\n' % i)
                    sys.stdout.flush()
                c.perform()
                if url.startswith('http'):
                    status = c.getinfo(pycurl.HTTP_CODE)
                    if status == 200:
                        break
                if url.startswith('ftp'):
                    status = 500
                    for h in headers:
                        if h.startswith('226'):
                            status = 200
                            break
            except pycurl.error as (errno, strerror):
                status = 500
                sys.stdout.write('\tPycURL error: %u, %s\n' % (errno, strerror))
                sys.stdout.flush()
        c.close()
    # sometimes authentication or cookies are needed to access the target url:
    if init_url and not usecache:
        if not silent:
            sys.stdout.write('\b'*20 + ' '*20 + '\b'*20 + 'Success.\n')
            sys.stdout.flush()
        # here, you may define a custom function to fetch 
        # the authentication data from cookies/headers, 
        # and return with headers for the main request:
        req_headers = globals()[init_fun](headers)
        if init_headers: return req_headers
        return curl(url = url, req_headers = req_headers, silent = silent, 
            debug = debug, outf = outf, compr = compr, encoding = encoding, 
            files_needed = files_needed, timeout = timeout, large = large,
            write_cache = write_cache)
    # get the data from the file downloaded/loaded from cache:
    if usecache or status == 200:
        if type(result) is file:
            fname = result.name
            result.close()
            result = open(fname, 'r')
        # find out the encoding:
        if encoding is None:
            if not usecache:
                headers = get_headers(headers)
                encoding = None
                if 'content-type' in headers:
                    content_type = headers['content-type'].lower()
                    match = re.search('charset=(\S+)', content_type)
                    if match:
                        encoding = match.group(1)
                if encoding is None:
                    if url.startswith('ftp'):
                        encoding = 'utf-8'
                    else:
                        encoding = 'iso-8859-1'
            else:
                # in case of using the cache:
                encoding = 'utf-8'
        if not silent and not usecache:
            sys.stdout.write('\b'*20 + ' '*20 + '\b'*20 + 'Success.\n')
            sys.stdout.flush()
        result.seek(0)
        if url.endswith('tar.gz') or url.endswith('tgz') or compr == 'tgz':
            multifile = True
            results = {}
            res = tarfile.open(fileobj = result, mode = 'r:gz')
            membs = res.getmembers()
            for m in membs:
                if (files_needed is None or m.name in files_needed) \
                    and m.size != 0:
                    # m.size is 0 for dierctories
                    this_file = res.extractfile(m)
                    if large:
                        results[m.name] = this_file
                    else:
                        results[m.name] = this_file.read()
                        this_file.close()
            if not large:
                res.close()
        elif url.endswith('gz') or compr == 'gz':
            res = gzip.GzipFile(fileobj=result, mode='rb')
            if not large:
                res = res.read()
                try:
                    res = res.decode(encoding)
                    res = res.encode('utf-8')
                except:
                    # better to proceed even if there is some trouble with encodings...
                    pass
        elif url.endswith('zip') or compr == 'zip':
            multifile = True
            results = {}
            res = zipfile.ZipFile(result,'r')
            membs = res.namelist()
            for m in membs:
                if files_needed is None or m in files_needed:
                    this_file = res.open(m)
                    if large:
                        results[m] = this_file
                    else:
                        results[m] = this_file.read()
                        this_file.close()
            res.close()
        else:
            if large:
                res = result
            else:
                res = result.getvalue()
        if not multifile:
            results = {'one': res}
        if not large:
            for k in results.keys():
                # handle files with CR line endings:
                if '\r' in results[k] and '\n' not in results[k]:
                    results[k] = results[k].replace('\r','\n')
                else:
                    results[k] = results[k].replace('\r','')
                if 'encoding' != 'utf-8':
                    try:
                        results[k] = results[k].decode(encoding).encode('utf-8')
                    except:
                        pass
        if (cache or write_cache) and not usecache and not large:
            for k in results.keys():
                if not multifile and not url.endswith('gz'):
                # write the decoded data back to StringIO
                    result.truncate(0)
                    result.write(results[k])
                # if cache is turned on, but data is not from cache,
                # place it there to make available next time:
                result.seek(0)
                with open(cachefile,'wb') as f:
                    f.write(result.getvalue())
        res = results if multifile else results['one']
    else:
        # download error:
        if not silent:
            sys.stdout.write('\b'*20 + ' '*20 + '\b'*20 + \
                'Failed. (Status: %u)\n'%status)
            if status > 200:
                sys.stdout.write('\t# URL: %s\n\t# POST: %s\n' % \
                    (url, '' if type(post) is not dict else urllib.urlencode(post)))
            sys.stdout.flush()
        res = None
    # returns raw data, dict of file names and raw data in case of 
    # multiple file archives, or file object in case of large files:
    return res

#
# thanks for http://stackoverflow.com/a/3239248/854988
#

def read_xls(xls_file, sheet = '', csv_file = None, return_table = True):
    try:
        book = open_workbook(xls_file, on_demand = True)
        try:
            sheet = book.sheet_by_name(sheet)
        except XLRDError:
            sheet = book.sheet_by_index(0)
            table = [[str(c.value) for c in sheet.row(i)] for i in xrange(sheet.nrows)]
            if csv_file:
                with open(csv_file, 'w') as csv:
                    csv.write('\n'.join(['\t'.join(r) for r in table]))
            if return_table:
                return table
    except IOError:
        sys.stdout.write('No such file: %s\n' % xls_file)
    sys.stdout.flush()

def read_table(cols, fileObject = None, data = None, sep = '\t', sep2 = None, rem = [], hdr = None):
    '''
    Generic function to read data tables.
    
    fileObject : file-like
        Any file like object: file opened for read, or StringIO buffer
    cols : dict
        Dictionary of columns to read. Keys identifying fields are returned 
        in the result. Values are column numbers.
    sepLevel1 : str
        Field separator of the file.
    sepLevel2 : dict
        Subfield separators and prefixes. 
        E.g. {2: ',', 3: '|'}
    hdr : int
        Number of header lines. If None, no headers assumed.
    rem : list
        Strings to remove. For each line these elements will be replaced with ''.
    '''
    if data is None:
        if 'readline' not in fileObject.__class__.__dict__:
            funname = sys._getframe().f_code.co_name
            sys.stdout.write('\tERROR: %s() expects file like object (file opened for read'\
                ', or StringIO buffer, etc)\n'%funname)
        fileObject.seek(0)
        if hdr:
            for h in xrange(0,hdr):
                null = fileObject.readline()
                del null
        data = fileObject
    else:
        data = [l.strip() for l in data.split('\n') if len(l) > 0][hdr:]
    res = []
    for l in data:
        for r in rem:
            l = l.replace(r,'')
        l = [f.strip() for f in l.split(sep)]
        if len(l) > max(cols.values()):
            dic = {}
            for name,col in cols.iteritems():
                field = l[col].strip()
                if sep2 is not None:
                    field = [sf.strip() for sf in field.split(sep2) if len(sf) > 0]
                dic[name] = field
            res.append(dic)
    if fileObject is not None: fileObject.close()
    return res

def all_uniprots(organism = 9606, swissprot = None):
    result = []
    rev = '' if swissprot is None else ' AND reviewed:%s'%swissprot
    url = data_formats.urls['uniprot_basic']['url']
    post = {'query': 'organism:%s%s' % (str(organism), rev), 
        'format': 'tab', 'columns': 'id'}
    data = curl(url, post = post, silent = False)
    data = data.split('\n')
    del data[0]
    for l in data:
        result.append(l.strip())
    return filter(lambda x: len(x) > 0, result)

def swissprot_seq(organism = 9606, isoforms = False):
    taxids = {
        9606: 'Homo sapiens'
    }
    result = {}
    url = data_formats.urls['uniprot_basic']['url']
    post = {'query': 'organism:%s AND reviewed:yes'%str(organism), 
        'format': 'tab', 'columns': 'id,sequence'}
    data = curl(url, post = post, silent = False)
    data = data.split('\n')
    del data[0]
    for l in data:
        l = l.strip().split('\t')
        if len(l) == 2:
            result[l[0]] = se.Seq(l[0], l[1])
    if isoforms:
        data = get_isoforms()
        for unip, isoforms in data.iteritems():
            for isof, seq in isoforms.iteritems():
                if unip in result:
                    result[unip].add_seq(seq, isof)
    return result

def get_pdb():
    data = curl(data_formats.urls['uniprot_pdb']['url'], silent=False)
    if data is None:
        return None, None
    data = data.split('\n')
    u_pdb = {}
    pdb_u = {}
    pdb = None
    pdb_re = re.compile(r'[0-9A-Z]{4}')
    for l in data:
        l = re.split('[ ]{2,}',re.sub('[ ]+,[ ]+',',',re.sub('[ ]*\(','(',l)))
        if len(l[0]) == 4 and pdb_re.match(l[0]):
            pdb = l[0].lower()
            res = None if l[2] == '-' else float(l[2].replace(' A',''))
            met = l[1]
        if pdb is not None and len(l) > 1:
            uniprots = l[1] if len(l) < 4 else l[3]
            uniprots = [u.split('(')[1].replace(')','') 
                        for u in uniprots.split(',') if '(' in u]
            pdb_u[pdb] = uniprots
            for u in uniprots:
                if u not in u_pdb:
                    u_pdb[u] = []
                u_pdb[u].append((pdb,met,res))
    return u_pdb, pdb_u

def get_pfam(uniprots=None,organism=None):
    if uniprots is None and organism is None:
        return None,None
    u_pfam = {}
    pfam_u = {}
    if uniprots is not None:
        prg = progress.Progress(len(uniprots)/30,'Downloading data from UniProt',1)
        data_all = []
        for i in xrange(0,len(uniprots),30):
            to = i + 30
            thisPart = uniprots[i:to]
            thisPart = ' OR '.join(['accession:%s'%u for u in thisPart])
            post = {'query': thisPart, 'format': 'tab', 'columns': 'id,database(Pfam)'}
            for j in xrange(3):
                data = curl(data_formats.urls['uniprot_basic']['url'],post=post)
                if data is not None:
                    break
            if data is None:
                return None,None
            data = data.split('\n')
            del data[0]
            del data[-1]
            data_all += data
            prg.step()
        prg.terminate()
    else:
        if type(organism) is not int:
            try:
                organism = int(organism)
            except:
                return None,None
        organismQuery = 'organism:%u AND reviewed:yes'%organism
        post = {'query': organismQuery, 'format': 'tab', 'columns': 'id,database(Pfam)'}
        for j in xrange(3):
            data_all = curl(data_formats.urls['uniprot_basic']['url'],post=post,
                            silent=False,outf='uniprot-pfam-%u.tab'%organism)
            if data_all is not None:
                break
        if data_all is None:
            return None
        data_all = data_all.split('\n')
        del data_all[0]
    for l in data_all:
        l = l.split('\t')
        pfams = [] if len(l) < 2 else re.sub(';$','',l[1]).split(';')
        if l[0] not in u_pfam:
            u_pfam[l[0]] = []
        u_pfam[l[0]] += pfams
        for pfam in pfams:
            if pfam not in pfam_u:
                pfam_u[pfam] = []
            pfam_u[pfam].append(l[0])
    return u_pfam, pfam_u

def get_pfam_regions(uniprots = [], pfams = [], keepfile = False, dicts = 'both'):
    url = data_formats.urls['pfam_up']['url']
    outf = url.split('/')[-1]
    urlmd5 = hashlib.md5(url).hexdigest()
    if not os.path.exists(os.path.join(os.getcwd(),'cache')):
        os.mkdir(os.path.join(os.getcwd(),'cache'))
    cachefile = os.path.join(os.getcwd(),'cache',urlmd5+'-'+outf)
    u_pfam = {}
    pfam_u = {}
    uniprots = set(uniprots)
    pfams = set(pfams)
    if not os.path.exists(cachefile):
        sys.stdout.write('\t:: Downloading data from %s' % \
            url.replace('http://', '').replace('ftp://', '').split('/')[0])
        sys.stdout.flush()
        urllib.urlretrieve(url, cachefile)
        sys.stdout.write('\n')
    with open(cachefile, 'rb') as f:
        f.seek(-4, 2)
        gzsize = struct.unpack('<I', f.read())[0]
        prg = progress.Progress(gzsize, 'Processing Pfam domains', 11)
    with gzip.open(cachefile, 'r') as f:
        for l in f:
            prg.step(len(l))
            l = l.strip().split()
            if l[0] in uniprots or l[4] in pfams:
                if dicts in ['uniprot', 'both']:
                    if l[0] not in u_pfam:
                        u_pfam[l[0]] = {}
                    if l[4] not in u_pfam[l[0]]:
                        u_pfam[l[0]][l[4]] = []
                    u_pfam[l[0]][l[4]].append({
                        'isoform': int(l[1]),
                        'start': int(l[5]),
                        'end': int(l[6])
                    })
                if dicts in ['pfam', 'both']:
                    if l[4] not in pfam_u:
                        pfam_u[l[4]] = {}
                    if l[0] not in pfam_u[l[4]]:
                        pfam_u[l[4]][l[0]] = []
                    pfam_u[l[4]][l[0]].append({
                        'isoform': int(l[1]),
                        'start': int(l[5]),
                        'end': int(l[6])
                    })
    prg.terminate()
    if not keepfile:
        os.remove(cachefile)
    if dicts == 'uniprot':
        return u_pfam
    elif dicts == 'pfam':
        return pfam_u
    else:
        return u_pfam, pfam_u

def get_pfam_names():
    data = curl(data_formats.urls['pfam_pdb']['url'],silent=False)
    if data is None:
        return None, None
    dname_pfam = {}
    pfam_dname = {}
    data = data.replace('\r','').split('\n')
    del data[0]
    for l in data:
        l = l.split('\t')
        if len(l) > 5:
            pfam = l[4].split('.')[0]
            name = l[5]
            if pfam not in pfam_dname:
                pfam_dname[pfam] = []
            if name not in dname_pfam:
                dname_pfam[name] = []
            pfam_dname[pfam].append(name)
            dname_pfam[name].append(pfam)
    for k,v in pfam_dname.iteritems():
        pfam_dname[k] = list(set(v))
    for k,v in dname_pfam.iteritems():
        dname_pfam[k] = list(set(v))
    return dname_pfam, pfam_dname

def get_pfam_pdb():
    non_digit = re.compile(r'[^\d.-]+')
    data = curl(data_formats.urls['pfam_pdb']['url'],silent=False)
    if data is None:
        return None, None
    pdb_pfam = {}
    pfam_pdb = {}
    data = data.replace('\r','').split('\n')
    del data[0]
    for l in data:
        l = l.split('\t')
        if len(l) > 4:
            pfam = l[4].split('.')[0]
            pdb = l[0].lower()
            chain = l[1]
            start = int(non_digit.sub('',l[2]))
            end = int(non_digit.sub('',l[3]))
            if pdb not in pdb_pfam:
                pdb_pfam[pdb] = {}
            if pfam not in pfam_pdb:
                pfam_pdb[pfam] = {}
            pdb_pfam[pdb][pfam] = [chain,start,end]
            pfam_pdb[pfam][pdb] = [chain,start,end]
    return pdb_pfam, pfam_pdb

def get_corum():
    complexes = {}
    members = {}
    data = curl(data_formats.urls['corum']['url'], silent=False)
    if data is None:
        return None,None
    data = data.split('\n')
    del data[0]
    prg = progress.Progress(len(data),'Processing data',9)
    for l in data:
        l = l.replace('\n','').replace('\r','').split(';')
        if len(l) < 10:
            continue
        uniprots = l[4].split(',')
        name = l[1]
        shortName = l[1].split('(')[0].strip()
        pubmeds = l[7].split(',')
        spec = l[3]
        func = l[9].replace('"','')
        dise = l[10].replace('"','')
        complexes[name] = (uniprots,shortName,pubmeds,spec,func,dise)
        for u in uniprots:
            if u not in members:
                members[u] = []
            members[u].append((name,shortName,pubmeds,spec,func,dise))
        prg.step()
    prg.terminate()
    return complexes,members

def get_complexportal(species=9606,zipped=True):
    '''
    Complex dataset from IntAct.
    See more:
    http://www.ebi.ac.uk/intact/complex/
    http://nar.oxfordjournals.org/content/early/2014/10/13/nar.gku975.full.pdf
    '''
    spec = {
        9606: 'Homo_sapiens'
    }
    species = species if type(species) is not int else spec[species]
    if zipped:
        zipurl = '/'.join([data_formats.urls['complex_portal']['url'],species+'.zip'])
        files = curl(zipurl,silent=False)
        if files is None:
            return None
    else:
        url = '/'.join([data_formats.urls['complex_portal']['url'],species,''])
        lst = curl(url,silent=False)
        if lst is None:
            return None
    if zipped:
        lst = files.keys()
    else:
        lst = [w for l in [l.split(' ') for l in lst.split('\n')] \
            for w in l if w.endswith('xml')]
    prg = progress.Progress(len(lst),'Downloading complex data',1)
    errors = []
    complexes = []
    for xmlname in lst:
        if zipped:
            xml = files[xmlname]
        else:
            url = '/'.join([data_formats.urls['complex_portal']['url'],species,xmlname])
            xml = curl(url)
        if xml is None:
            msg = 'Could not get file: \n\t\t%s' % url
            errors.append(msg)
            continue
        soup = bs4.BeautifulSoup(xml, 'html.parser')
        interactors_xml = soup.find_all('interactor')
        interactors = {}
        interactions = {}
        for i in interactors_xml:
            if i.find('primaryref').attrs['db'] == 'uniprotkb':
                interactors[i.attrs['id']] = i.find('primaryref').attrs['id']
            #'tax_id': i.find('organism').attrs['ncbitaxid']
        interactions_xml = soup.find_all('interaction')
        for i in interactions_xml:
            description = ''
            pubmeds = []
            fullname = ''
            names = {}
            pdbs = []
            uniprots = []
            for a in i.find_all('attribute'):
                if a.attrs['name'] == 'curated-complex':
                    description = a.text
            for sr in i.find_all('secondaryref'):
                if sr.attrs['db'] == 'pubmed':
                    pubmeds.append(sr.attrs['id'])
                if sr.attrs['db'] == 'wwpdb':
                    pdbs.append(sr.attrs['id'])
            for pr in i.find_all('primaryref'):
                if pr.attrs['db'] in ['wwpdb', 'rcsb pdb', 'pdbe']:
                    pdbs.append(pr.attrs['id'])
            pubmeds = list(set(pubmeds))
            pdbs = list(set(pdbs))
            fullname = '' if i.find('fullname') is None else i.find('fullname').text
            for a in i.find_all('alias'):
                names[a.attrs['type']] = a.text
            for intref in i.find_all('interactorref'):
                int_id = intref.text
                if int_id in interactors:
                    uniprots.append(interactors[int_id])
            complexes.append({
                'uniprots': uniprots,
                'pdbs': pdbs,
                'pubmeds': pubmeds,
                'fullname': fullname,
                'names': names,
                'description': description
                })
        prg.step()
    prg.terminate()
    if len(errors) > 0:
        sys.stdout.write('\t:: Failed to download %u files of total %u:\n\n' %
            (len(e),len(lst)))
        for e in errors:
            sys.stdout.write('\t'+e+'\n')
        sys.stdout.flush()
    return complexes

def read_complexes_havugimana():
    '''
    Supplement Table S3/1 from Havugimana 2012
    Cell. 150(5): 1068–1081.
    '''
    complexes = []
    # TODO: place file name into settings
    infile = os.path.join(common.ROOT,'data','complexes_havugimana2012.csv')
    with codecs.open(infile,encoding='utf-8',mode='r') as f:
        hdr = f.readline()
        del hdr
        for l in f:
            l = l.replace('\n','').replace('\r','').split('\t')
            complexes.append(l[2].split(','))
    return complexes

def get_compleat():
    url = data_formats.urls['compleat']['url']
    data = curl(url,silent=False)
    data = data.replace('\r','').split('\n')
    complexes = []
    for l in data:
        l = l.split('\t')
        if len(l) > 11:
            complexes.append({
                'source': l[6],
                'spec': [s.strip() for s in l[9].split('&')],
                'pubmeds': l[10].split(','),
                'entrez': [ee for ee in [e.strip() for e in l[11].split(' ')] \
                           if len(ee) > 0],
                'functions': l[4]
            })
    return complexes

def get_pdb_chains():
    chains = curl(data_formats.urls['pdb_chains']['url'],silent=False)
    if chains is None:
        return None,None
    chains = chains.replace('\r','').split('\n')
    del chains[0]
    del chains[0]
    pdb_u = {}
    u_pdb = {}
    non_digit = re.compile(r'[^\d.-]+')
    for l in chains:
        l = l.split('\t')
        if len(l) > 8:
            if l[0] not in pdb_u:
                pdb_u[l[0]] = {}
            pdb_u[l[0]][l[1]] = {
                'uniprot': l[2],
                'chain_beg': int(non_digit.sub('',l[3])),
                'chain_end': int(non_digit.sub('',l[4])),
                'pdb_beg': int(non_digit.sub('',l[5])),
                'pdb_end': int(non_digit.sub('',l[6])),
                'uniprot_beg': int(non_digit.sub('',l[7])),
                'uniprot_end': int(non_digit.sub('',l[8]))
                }
            if pdb_u[l[0]][l[1]]['pdb_end'] - pdb_u[l[0]][l[1]]['pdb_beg'] == \
                pdb_u[l[0]][l[1]]['uniprot_end'] - pdb_u[l[0]][l[1]]['uniprot_beg']:
                pdb_u[l[0]][l[1]]['offset'] = (pdb_u[l[0]][l[1]]['uniprot_beg'] -
                                               pdb_u[l[0]][l[1]]['pdb_beg'])
            else:
                pdb_u[l[0]][l[1]]['offset'] = None
            if l[2] not in u_pdb:
                u_pdb[l[2]] = []
            u_pdb[l[2]].append({
                    'pdb': l[0],
                    'chain': l[1],
                    'chain_beg': int(non_digit.sub('',l[3])),
                    'chain_end': int(non_digit.sub('',l[4])),
                    'pdb_beg': int(non_digit.sub('',l[5])),
                    'pdb_end': int(non_digit.sub('',l[6])),
                    'uniprot_beg': int(non_digit.sub('',l[7])),
                    'uniprot_end': int(non_digit.sub('',l[8])),
                    'offset': pdb_u[l[0]][l[1]]['offset']
                })
    return u_pdb, pdb_u

def get_3dcomplexes():
    contact = curl(data_formats.urls['3dcomplexes_contact']['url'],silent=False)
    corresp = curl(data_formats.urls['3dcomplexes_correspondancy']['url'],silent=False)
    u_pdb, pdb_u = get_pdb_chains()
    del u_pdb
    if contact is None or corresp is None or pdb_u is None:
        return None
    contact = contact.split('\n')
    corresp = corresp.split('\n')
    del contact[0]
    corr_dict = {}
    for l in corresp:
        l = l.replace('\r','').split('\t')
        if len(l) > 2:
            pdb = l[0].split('.')[0]
            if pdb not in corr_dict:
                corr_dict[pdb] = {}
            corr_dict[pdb][l[1]] = l[2]
    compl_dict = {}
    for l in contact:
        l = l.replace('\r','').split('\t')
        if len(l) > 11 and int(l[11]) == 0 and int(l[10]) == 0:
            compl = l[0]
            pdb = compl.split('_')[0]
            if pdb in corr_dict:
                if l[1] in corr_dict[pdb] and l[2] in corr_dict[pdb]:
                    ch1 = corr_dict[pdb][l[1]]
                    ch2 = corr_dict[pdb][l[2]]
                    if pdb in pdb_u and ch1 in pdb_u[pdb]:
                        up1 = pdb_u[pdb][ch1]['uniprot']
                        if pdb in pdb_u and ch2 in pdb_u[pdb]:
                            up2 = pdb_u[pdb][ch2]['uniprot']
                            if compl not in compl_dict:
                                compl_dict[compl] = {}
                            uniprots = [up1,up2]
                            uniprots.sort()
                            uniprots = tuple(uniprots)
                            if uniprots not in compl_dict[compl]:
                                compl_dict[compl][uniprots] = []
                            compl_dict[compl][uniprots].append(float(l[3]))
    return compl_dict

def get_domino_interactions():
    domino = get_domino()
    inter = []
    for l in domino:
        if len(l[0]) > 0 and len(l[1]) > 0 and len(''.join(l[5])) > 0 and \
            len(''.join([l[i] for i in \
                range(10,12) + range(14,22) + range(24,26)])) != 0 and \
            l[28] != '1':
            inter.append(l)
    return inter

def get_domino_ddi():
    domi = get_domino_ptms()
    return domi['ddi']

def get_domino_ptms():
    '''
    The table comes from dataio.get_domino(), having the following fields:
    header = ['uniprot-A', 'uniprot-B', 'isoform-A', 'isoform-B', #3
    'exp. method', 'references', 'taxon-A', 'taxon-B', #7
    'role-A', 'role-B', 'binding-site-range-A', 'binding-site-range-B', #11
    'domains-A', 'domains-B', 'ptm-residue-A', 'ptm-residue-B', #15
    'ptm-type-mi-A', 'ptm-type-mi-B', 'ptm-type-A', 'ptm-type-B', #19
    'ptm-res-name-A', 'ptm-res-name-B', 'mutations-A', 'mutations-B', #23
    'mutation-effects-A', 'mutation-effects-B', 'domains-interpro-A', #26 
    'domains-interpro-B', 'negative'] #28
    '''
    domino = get_domino()
    miont = get_ontology('MI')
    dmi = []
    ddi = []
    prg = progress.Progress(len(domino), 'Processing DOMINO', 11)
    for l in domino:
        prg.step()
        if (l[14].strip() != '' or l[15].strip() != '' or 
            (l[10] != '' and l[11] != '')) and len(l[0]) > 0 and len(l[1]) > 0:
            uniprot1 = l[0]
            uniprot2 = l[1]
            # ptms
            if '-' not in l[14] and '-' not in l[15]:
                ptmre12 = [] if len(l[14]) == 0 else \
                    [int(x) for x in l[14].split(';')]
                ptmre21 = [] if len(l[15]) == 0 else \
                    [int(x) for x in l[15].split(';')]
                ptmty12 = [None for _ in ptmre12] if len(l[16]) == 0 else \
                    l[16].split(';')
                ptmty12 = [None if x not in miont else miont[x] for x in ptmty12]
                ptmrn12 = [None for _ in ptmre12] if len(l[20]) == 0 else \
                    l[20].split(';')
                ptmrn12 = [None if x is None or x == '' or \
                        len(x) < min (ptmre12[i] - 1, 11) \
                    else x[10] if ptmre12[i] > 10 \
                    else x[ptmre12[i] - 1] for i, x in enumerate(ptmrn12)]
                ptmty21 = [None for _ in ptmre21] if len(l[17]) == 0 else \
                    l[17].split(';')
                ptmty21 = [None if x not in miont else miont[x] for x in ptmty21]
                ptmrn21 = [None for _ in ptmre21] if len(l[21]) == 0 else \
                    l[21].split(';')
                ptmrn21 = [None if x is None or x == '' or \
                        len(x) < min (ptmre21[i] - 1, 11) \
                    else x[10] if ptmre21[i] > 10 \
                    else x[ptmre21[i] - 1] for i, x in enumerate(ptmrn21)]
                for i, resnum in enumerate(ptmre12):
                    res = intera.Residue(resnum, ptmrn12[i], uniprot2)
                    ptm = intera.Ptm(uniprot2, typ = ptmty12[i], residue = res, 
                        source = 'DOMINO')
                    dom = intera.Domain(uniprot1)
                    dm = intera.DomainMotif(domain = dom, ptm = ptm, 
                        sources = 'DOMINO', refs = l[5].split(';'))
            # binding sites
            if l[10] != '' and l[11] != '':
                try:
                    bssrt1 = [int(x.split('-')[0]) for x in l[10].split(';') \
                        if x != '' and x != '0']
                    bsend1 = [int(x.split('-')[1]) for x in l[10].split(';') \
                        if x != '' and x != '0']
                    bssrt2 = [int(x.split('-')[0]) for x in l[11].split(';') \
                        if x != '' and x != '0']
                    bsend2 = [int(x.split('-')[1]) for x in l[11].split(';') \
                        if x != '' and x != '0']
                except:
                    print l
                    return None
                bs1 = []
                bs2 = []
                if l[26] != '':
                    for i, n in enumerate(bssrt1):
                        bs1.append(intera.Domain(protein = uniprot1, domain = l[26], 
                            start = bssrt1[i], end = bsend1[i], 
                            domain_id_type = 'interpro', isoform = l[2]))
                else:
                    for i, n in enumerate(bssrt1):
                        mot = intera.Motif(protein = uniprot1, start = bssrt1[i], 
                            end = bsend1[i], isoform = l[2])
                        bs1.append(intera.Ptm(protein = uniprot1, motif = mot, 
                            source = 'DOMINO', isoform = l[2]))
                if l[27] != '':
                    for i, n in enumerate(bssrt2):
                        bs2.append(intera.Domain(protein = uniprot2, domain = l[27], 
                            start = bssrt2[i], end = bsend2[i], 
                            domain_id_type = 'interpro', isoform = l[3]))
                else:
                    for i, n in enumerate(bssrt2):
                        mot = intera.Motif(protein = uniprot2, start = bssrt2[i], 
                            end = bsend2[i], isoform = l[3])
                        bs2.append(intera.Ptm(protein = uniprot2, motif = mot, 
                            source = 'DOMINO'))
                for one in bs1:
                    for two in bs2:
                        if one.__class__.__name__ == 'Domain' and \
                            two.__class__.__name__ == 'Domain':
                            dd = intera.DomainDomain(one, two, sources = 'DOMINO')
                            ddi.append(dd)
                        if one.__class__.__name__ == 'Domain' and \
                            two.__class__.__name__ == 'Ptm':
                            dm = intera.DomainMotif(domain = one, ptm = two, 
                                sources = 'DOMINO', refs = l[6].split(';'))
                            dmi.append(dm)
                        if two.__class__.__name__ == 'Domain' and \
                            one.__class__.__name__ == 'Ptm':
                            dm = intera.DomainMotif(domain = two, ptm = one, 
                                sources = 'DOMINO', refs = l[6].split(';'))
                            dmi.append(dm)
    prg.terminate()
    return {'ddi': ddi, 'dmi': dmi}

def get_3dc_ddi():
    contact = curl(data_formats.urls['3dcomplexes_contact']['url'], silent=False)
    corresp = curl(data_formats.urls['3dcomplexes_correspondancy']['url'], silent=False)
    u_pdb, pdb_u = get_pdb_chains()
    del u_pdb
    if contact is None or corresp is None or pdb_u is None:
        return None
    contact = contact.split('\n')
    corresp = corresp.split('\n')
    del contact[0]
    corr_dict = {}
    ddi = []
    uniprots = []
    for l in corresp:
        l = l.replace('\r','').split('\t')
        if len(l) > 2:
            pdb = l[0].split('.')[0]
            if pdb not in corr_dict:
                corr_dict[pdb] = {}
            corr_dict[pdb][l[1]] = l[2]
    prg = progress.Progress(len(contact), 'Collecting UniProts', 9)
    for l in contact:
        prg.step()
        l = l.replace('\r','').split('\t')
        if len(l) > 11 and int(l[11]) == 0 and int(l[10]) == 0:
            pdb = l[0].split('_')[0]
            if pdb in corr_dict:
                if l[1] in corr_dict[pdb] and l[2] in corr_dict[pdb]:
                    ch1 = corr_dict[pdb][l[1]]
                    ch2 = corr_dict[pdb][l[2]]
                    if pdb in pdb_u and ch1 in pdb_u[pdb]:
                        up1 = pdb_u[pdb][ch1]['uniprot']
                    if pdb in pdb_u and ch2 in pdb_u[pdb]:
                        up2 = pdb_u[pdb][ch2]['uniprot']
                    uniprots += [up1, up2]
    prg.terminate()
    uniprots = list(set(uniprots))
    u_pfam = get_pfam_regions(uniprots, dicts = 'uniprot')
    prg = progress.Progress(len(contact), 'Processing contact information', 9)
    for l in contact:
        prg.step()
        l = l.replace('\r','').split('\t')
        if len(l) > 11 and int(l[11]) == 0 and int(l[10]) == 0:
            pdb = l[0].split('_')[0]
            pfams1 = list(set([x.split('.')[0] for x in l[7].split(';')]))
            pfams2 = list(set([x.split('.')[0] for x in l[9].split(';')]))
            if pdb in corr_dict:
                if l[1] in corr_dict[pdb] and l[2] in corr_dict[pdb]:
                    ch1 = corr_dict[pdb][l[1]]
                    ch2 = corr_dict[pdb][l[2]]
                    if pdb in pdb_u and ch1 in pdb_u[pdb]:
                        up1 = pdb_u[pdb][ch1]['uniprot']
                        if pdb in pdb_u and ch2 in pdb_u[pdb]:
                            up2 = pdb_u[pdb][ch2]['uniprot']
                            for pfam1 in pfams1:
                                for pfam2 in pfams2:
                                    pfam1_details = [{'start': None, 'end': None, 
                                        'isoform': 1}]
                                    pfam2_details = [{'start': None, 'end': None, 
                                        'isoform': 1}]
                                    if up1 in u_pfam and pfam1 in u_pfam[up1]:
                                        pfam1_details = u_pfam[up1][pfam1]
                                    if up2 in u_pfam and pfam2 in u_pfam[up2]:
                                        pfam2_details = u_pfam[up2][pfam2]
                                    for pfam1_d in pfam1_details:
                                        for pfam2_d in pfam2_details:
                                            dom1 = intera.Domain(
                                                protein = up1, 
                                                domain = pfam1, 
                                                start = pfam1_d['start'],
                                                end = pfam1_d['end'],
                                                isoform = pfam1_d['isoform'],
                                                chains = {pdb: ch1}
                                            )
                                            dom2 = intera.Domain(
                                                protein = up2, 
                                                domain = pfam2, 
                                                start = pfam2_d['start'],
                                                end = pfam2_d['end'],
                                                isoform = pfam2_d['isoform'],
                                                chains = {pdb: ch2}
                                            )
                                            dd = intera.DomainDomain(dom1, dom2, 
                                                pdbs = pdb, sources = '3DComplex',
                                                contact_residues = float(l[3]))
                                            ddi.append(dd)
    prg.terminate()
    return ddi

def pisa_bonds(lst, chains):
    non_digit = re.compile(r'[^\d.-]+')
    bonds = []
    for bond in lst.find_all('bond'):
        seqnum1 = int(non_digit.sub('',bond.find('seqnum-1').text))
        seqnum2 = int(non_digit.sub('',bond.find('seqnum-2').text))
        res1 = bond.find('res-1').text
        res1 = res1 if res1 not in common.aaletters else common.aaletters[res1]
        res2 = bond.find('res-2').text
        res2 = res2 if res2 not in common.aaletters else common.aaletters[res2]
        chain1 = bond.find('chain-1').text
        chain2 = bond.find('chain-2').text
        uniprot1 = None if chain1 not in chains else chains[chain1]
        uniprot2 = None if chain2 not in chains else chains[chain2]
        if uniprot1 is not None and uniprot2 is not None:
            bonds.append({
                'chain_1': chain1,
                'uniprot_1': uniprot1,
                'res_1':res1,
                'seqnum_1': seqnum1,
                'chain_2': chain2,
                'uniprot_2': uniprot2,
                'res_2': res2,
                'seqnum_2': seqnum2
            })
    return bonds

def get_pisa(pdblist):
    non_digit = re.compile(r'[^\d.-]+')
    bond_types = {'hbonds':'h-bonds', 'sbridges': 'salt-bridges',
                   'covbonds': 'cov-bonds', 'ssbonds': 'ss-bonds'}
    interfaces = {}
    cachefile = os.path.join('cache','pisa.pickle')
    u_pdb, pdb_u = get_pdb_chains()
    if os.path.exists(cachefile):
        try:
            interfaces = pickle.load(open(cachefile,'rb'))
        except:
            pass
    errors = []
    p = 5
    pdblist = list(set(pdblist) - set(interfaces.keys()))
    prg = progress.Progress(len(pdblist)/p,'Downloading data from PDBe PISA',1)
    for i in xrange(0,len(pdblist),p):
        to = i + p
        thisPart = pdblist[i:to]
        url = data_formats.urls['pisa_interfaces']['url'] + ','.join(thisPart)
        data = curl(url,cache=False)
        if data is None:
            msg = 'Could not download: \n\t\t%s' % url
            errors.append(msg)
            continue
        soup = bs4.BeautifulSoup(data, 'html.parser')
        unmapped_residues = []
        for pdb in soup.find_all('pdb_entry'):
            pdb_id = pdb.find('pdb_code').text.lower()
            interfaces[pdb_id] = {}
            chains = {}
            resconv = ResidueMapper()
            if pdb_id in pdb_u:
                for chain, chain_data in pdb_u[pdb_id].iteritems():
                    chains[chain] = chain_data['uniprot']
                for interface in pdb.find_all('interface'):
                    for b, t in bond_types.iteritems():
                        lst = interface.find(t)
                        if lst is not None:
                            bonds = pisa_bonds(lst,chains)
                            for bond in bonds:
                                uniprots = (bond['uniprot_1'],bond['uniprot_2'])
                                if uniprots not in interfaces[pdb_id]:
                                    css = non_digit.sub('',interface.find('css').text)
                                    css = None if len(css) == 0 else float(css)
                                    area = non_digit.sub('',interface.find('int_area').text)
                                    area = None if len(area) == 0 else float(area)
                                    solv_en = non_digit.sub('', 
                                                interface.find('int_solv_en').text)
                                    solv_en = None if len(solv_en) == 0 else float(solv_en)
                                    stab_en = non_digit.sub('', 
                                                interface.find('stab_en').text)
                                    stab_en = None if len(stab_en) == 0 else float(stab_en)
                                    interfaces[pdb_id][uniprots] = \
                                        intera.Interface(uniprots[0], 
                                            uniprots[1], source = 'PISA', pdb = pdb_id, 
                                            css = css, 
                                            solv_en = solv_en, 
                                            area = area, 
                                            stab_en = stab_en)
                                res1 = resconv.get_residue(pdb_id, bond['seqnum_1'])
                                res2 = resconv.get_residue(pdb_id, bond['seqnum_2'])
                                if res1 is not None and res2 is not None and \
                                    res1['uniprot'] == uniprots[0] and \
                                    res2['uniprot'] == uniprots[1]:
                                    interfaces[pdb_id][uniprots].add_residues(
                                        (res1['resnum'],bond['res_1'], uniprots[0]),
                                        (res2['resnum'],bond['res_2'], uniprots[1]),
                                        typ = b)
                                else:
                                    unmapped_residues.append((pdb_id, bond['seqnum_1'], 
                                        bond['seqnum_2'], uniprots[0], uniprots[1]))
        pickle.dump(interfaces, open(cachefile, 'wb'),2)
        prg.step()
    prg.terminate()
    if len(errors) > 0:
        sys.stdout.write('\t:: Failed to download %u files of total %u:\n\n' %
            (len(errors),len(lst)))
        for e in errors:
            sys.stdout.write('\t'+e+'\n')
        sys.stdout.flush()
    return interfaces, unmapped_residues

def get_3did_ddi(residues=False,ddi_flat=None,organism=9606):
    if ddi_flat is None:
        data = curl(data_formats.urls['3did_ddi']['url'],silent=False)
        tmpfile = '3did_flat_tmp'
        if data is None:
            return None
        with codecs.open(tmpfile,encoding='utf-8',mode='w') as f:
            f.write(data)
        lnum = data.count('\n')
        del data
    else:
        tmpfile = ddi_flat
    u_pfam, pfam_u = get_pfam(organism=organism)
    u_pdb, pdb_u = get_pdb_chains()
    if pfam_u is None or pdb_u is None:
        return None
    ddi = {}
    interfaces = {}
    pdblist = {}
    ddi_collect = False
    con_collect = False
    non_digit = re.compile(r'[^\d.-]+')
    with codecs.open(tmpfile,encoding='utf-8',mode='r') as f:
        prg = progress.Progress(lnum,'Reading data',33)
        for l in f:
            prg.step()
            if l.startswith('#=') and con_collect:
                interfaces[(uniprot1,uniprot2,pdb)].append(this_interface)
                con_collect = False
            if l.startswith('#=ID'):
                # new domain pair: attach previous to results:
                if ddi_collect:
                    for u1 in uniprots1:
                        for u2 in uniprots2:
                            if u1 != u2 and len(pdblist) > 0:
                                if (u1,u2) not in ddi:
                                    ddi[(u1,u2)] = {}
                                if (pfam1,pfam2) not in ddi[(u1,u2)]:
                                    ddi[(u1,u2)][(pfam1,pfam2)] = {
                                        'pdbs': pdblist}
                    ddi_collect = False
                pdblist = {}
                l = l.split('\t')
                pfam1 = l[3].split('(')[1].split('.')[0]
                pfam2 = l[4].split('.')[0]
                uniprots1 = [] if pfam1 not in pfam_u else pfam_u[pfam1]
                uniprots2 = [] if pfam2 not in pfam_u else pfam_u[pfam2]
                if len(set(uniprots1 + uniprots2)) > 1:
                    ddi_collect = True
            elif l.startswith('#=3D'):
                l = l.split('\t')
                pdb = l[1]
                chain1 = l[2].split(':')[0]
                chain2 = l[3].split(':')[0]
                if pdb in pdb_u and \
                    chain1 in pdb_u[pdb] and \
                    chain2 in pdb_u[pdb]:
                    uniprot1 = pdb_u[pdb][chain1]['uniprot']
                    uniprot2 = pdb_u[pdb][chain2]['uniprot']
                    if uniprot1 != uniprot2:
                        if pdb not in pdblist:
                            pdblist[pdb] = []
                        pdblist[pdb] = common.addToList(
                            pdblist[pdb],(uniprot1,uniprot2))
                    if residues:
                        #res1 = [int(i) for i in l[2].split(':')[1].split('-')]
                        #res2 = [int(i) for i in l[3].split(':')[1].split('-')]
                        if chain1 != chain2:
                            if pdb_u[pdb][chain1]['offset'] is not None and \
                                    pdb_u[pdb][chain2]['offset'] is not None and \
                                    pdb_u[pdb][chain1]['uniprot'] != \
                                    pdb_u[pdb][chain2]['uniprot']:
                                con_collect = True
                                offset1 = pdb_u[pdb][chain1]['offset']
                                offset2 = pdb_u[pdb][chain2]['offset']
                                this_interface = common.Interface(uniprot1,uniprot2, 
                                                            source = '3DID', pdb = pdb)
                                if (uniprot1,uniprot2,pdb) not in interfaces:
                                    interfaces[(uniprot1,uniprot2,pdb)] = []
                            else:
                                con_collect = False
            elif not residues or not con_collect:
                continue
            else:
                l = l.split('\t')
                if len(l) > 3:
                    rnum1 = int(non_digit.sub('',l[2])) + offset1
                    rnum2 = int(non_digit.sub('',l[3])) + offset2
                    this_interface.add_residues((rnum1,l[0],uniprot1),
                                                (rnum2,l[1],uniprot2))
        prg.terminate()
        prg = progress.Progress(len(ddi),'Processing interfaces',99)
        if residues:
            for u,v1 in ddi.iteritems():
                prg.step()
                for d,v2 in v1.iteritems():
                    for p in v2['pdbs'].keys():
                        if (u[0],u[1],p) in interfaces:
                            ddi[u][d]['interfaces'] = interfaces[(u[0],u[1],p)]
        prg.terminate()
    if ddi_flat is None:
        os.remove(tmpfile)
    if residues:
        return ddi, interfaces
    else:
        return ddi

def get_3did(ddi_flat = None, res = True, organism = 9606, pickl = True):
    resultfile = os.path.join('cache', '3did_ddi.pickle')
    if pickl and os.path.exists(resultfile):
        result = pickle.load(open(resultfile, 'rb'))
        if len(result) == 1: return result
        else: return result[0], result[1]
    if ddi_flat is None:
        data = curl(data_formats.urls['3did_ddi']['url'],silent=False)
        tmpfile = '3did_flat_tmp'
        if data is None:
            return None
        with codecs.open(tmpfile, encoding='utf-8', mode='w') as f:
            f.write(data)
        lnum = data.count('\n')
        del data
    elif os.path.exists(ddi_flat):
        tmpfile = ddi_flat
    else:
        return None
    u_pdb, pdb_u = get_pdb_chains()
    all_unip = set(all_uniprots(organism = organism))
    if all_unip is None or pdb_u is None:
        return None
    ddi = []
    interfaces = []
    pdb = pdb_prev = intf = None
    skip = True
    non_digit = re.compile(r'[^\d.-]+')
    rmap = residues.ResidueMapper()
    with codecs.open(tmpfile, encoding='utf-8', mode='r') as f:
        prg = progress.Progress(lnum, 'Processing 3DID domain-domain interactions', 33)
        for l in f:
            prg.step()
            l = l.split('\t')
            if l[0].startswith('#=ID'):
                pfam1 = l[3].split('.')[0][2:]
                pfam2 = l[4].split('.')[0]
            elif l[0].startswith('#=3D'):
                pdb_prev = pdb
                skip = True
                pdb = l[1]
                chain1 = l[2][0]
                chain2 = l[3][0]
                uniprot1 = uniprot2 = None
                if pdb != pdb_prev:
                    rmap.clean()
                if pdb in pdb_u:
                    if chain1 in pdb_u[pdb]:
                        uniprot1 = pdb_u[pdb][chain1]['uniprot']
                    if chain2 in pdb_u[pdb]:
                        uniprot2 = pdb_u[pdb][chain2]['uniprot']
                if uniprot1 is not None and uniprot2 is not None and \
                    uniprot1 in all_unip and uniprot2 in all_unip and \
                    uniprot1 != uniprot2:
                    skip = False
                    if intf is not None:
                        interfaces.append(intf)
                    intf = intera.Interface(uniprot1, uniprot2, '3DID', pdb)
                    u1start = u1end = u2start = u2end = {}
                    if l[2].count('-') == 1:
                        start1 = int(non_digit.sub('', l[2][2:].split('-')[0]))
                        end1 = int(non_digit.sub('', l[2][2:].split('-')[1]))
                        u1start = rmap.pdb2uniprot(pdb, start1, chains = chain1)
                        u1end = rmap.pdb2uniprot(pdb, end1, chains = chain1)
                    if l[3].count('-') == 1:
                        start2 = int(non_digit.sub('', l[3][2:].split('-')[0]))
                        end2 = int(non_digit.sub('', l[3][2:].split('-')[1]))
                        u2start = rmap.pdb2uniprot(pdb, start2, chains = chain2)
                        u2end = rmap.pdb2uniprot(pdb, end2, chains = chain2)
                    u1start = None if len (u1start) == 0 else \
                        u1start[chain1]['resnum']
                    u1end = None if len (u1end) == 0 else \
                        u1end[chain1]['resnum']
                    u2start = None if len (u2start) == 0 else \
                        u2start[chain2]['resnum']
                    u2end = None if len (u2end) == 0 else \
                        u2end[chain2]['resnum']
                    dom1 = intera.Domain(uniprot1, domain = pfam1, 
                        start = u1start, end = u1end, isoform = 1)
                    dom2 = intera.Domain(uniprot2, domain = pfam2, 
                        start = u2start, end = u2end, isoform = 1)
                    dd = intera.DomainDomain(dom1, dom2, [pdb], '3DID')
                    ddi.append(dd)
            elif not skip and res and not l[0].startswith('//'):
                conv1 = rmap.pdb2uniprot(pdb, int(non_digit.sub('', l[2])), 
                    chains = chain1)
                conv2 = rmap.pdb2uniprot(pdb, int(non_digit.sub('', l[3])), 
                    chains = chain2)
                if len(conv1) > 0 and len(conv2) > 0:
                    intf.add_residues((conv1[chain1]['resnum'], l[0], uniprot1), 
                        (conv2[chain2]['resnum'], l[1], uniprot2))
        interfaces.append(intf)
        prg.terminate()
    if ddi_flat is None:
        os.remove(tmpfile)
    if res:
        pickle.dump([ddi, interfaces], open(resultfile, 'wb'))
        return ddi, interfaces
    else:
        pickle.dump([ddi], open(resultfile, 'wb'))
        return ddi

def get_3did_dmi(dmi_flat = None):
    resultfile = os.path.join('cache', '3did_dmi.pickle')
    if os.path.exists(resultfile):
        return pickle.load(open(resultfile, 'rb'))
    if dmi_flat is None:
        data = curl(data_formats.urls['3did_dmi']['url'], silent = False)
        tmpfile = '3did_dmi_flat_tmp'
        if data is None:
            return None
        with codecs.open(tmpfile,encoding='utf-8',mode='w') as f:
            f.write(data)
        lnum = data.count('\n')
        del data
    elif os.path.exists(dmi_flat):
        tmpfile = dmi_flat
    else:
        return None
    u_pdb, pdb_u = get_pdb_chains()
    if pdb_u is None:
        return None
    dmi = {}
    non_digit = re.compile(r'[^\d.-]+')
    rmap = residues.ResidueMapper()
    with codecs.open(tmpfile,encoding='utf-8',mode='r') as f:
        prg = progress.Progress(lnum, 'Processing 3DID domain-motif interactions', 1)
        for l in f:
            prg.step()
            l = l.strip().split()
            if l[0].startswith('#=ID'):
                domain = l[3]
            if l[0].startswith('#=PT'):
                regex = l[1]
            if l[0].startswith('#=3D'):
                pdb = l[1]
                chain1 = l[2].split(':')[0]
                chain2 = l[3].split(':')[0]
                if l[2].count('-') == 1 and l[3].count('-') == 1:
                    pdb_region1 = [int(non_digit.sub('', x)) \
                        for x in l[2].split(':')[1].split('-')]
                    pdb_region2 = [int(non_digit.sub('', x)) \
                        for x in l[3].split(':')[1].split('-')]
                    u1start = rmap.pdb2uniprot(pdb, pdb_region1[0], chains = chain1)
                    u1end = rmap.pdb2uniprot(pdb, pdb_region1[1], chains = chain1)
                    u2start = rmap.pdb2uniprot(pdb, pdb_region2[0], chains = chain2)
                    u2end = rmap.pdb2uniprot(pdb, pdb_region2[1], chains = chain2)
                    if len(u1start) != 0 and len(u2start) != 0 and \
                        len(u1end) != 0 and len(u2end) != 0:
                        uniprot_key = (u1start[chain1]['uniprot'], 
                                        u2start[chain2]['uniprot'])
                        residue_key = (u1start[chain1]['resnum'], 
                                        u1end[chain1]['resnum'], 
                                        u2start[chain2]['resnum'], 
                                        u2end[chain2]['resnum'])
                        if uniprot_key not in dmi:
                            dmi[uniprot_key] = {}
                        if residue_key not in dmi[uniprot_key]:
                            dmi[uniprot_key][residue_key] = []
                        dmi[uniprot_key][residue_key].append({
                            'pdb': pdb,
                            'regex': regex,
                            'instance': l[4],
                            'domain': domain,
                            'contacts': int(non_digit.sub('', l[5])),
                            'topology': int(non_digit.sub('', l[6]))
                        })
        prg.terminate()
    if dmi_flat is None:
        os.remove(tmpfile)
    if len(rmap.download_errors) > 0:
        sys.stdout.write('Failed to download PDB-UniProt mappings for:\n'\
            '%s\n' % ', '.join(rmap.download_errors))
    pickle.dump(dmi, open(resultfile, 'wb'))
    return dmi

def process_3did_dmi():
    dmi = get_3did_dmi()
    if dmi is None:
        return None
    dname_pfam, pfam_dname = get_pfam_names()
    dname_re = re.compile(r'(.*)(_[A-Z]{3}_)(.*)')
    dmi2 = {}
    prg = progress.Progress(len(dmi), 'Processing data', 11)
    for uniprots, dmis in dmi.iteritems():
        prg.step()
        if uniprots not in dmi2:
            dmi2[uniprots] = []
        for regions, dmi_list in dmis.iteritems():
            new = True
            for dm in dmi_list:
                if new:
                    pfam = None
                    dname = None
                    mname = None
                    name_match = dname_re.match(dm['domain'])
                    if name_match:
                        dname = name_match.groups(0)[0]
                        mname = ''.join(name_match.groups(0)[1:])[1:]
                    if dname in dname_pfam:
                        pfam = dname_pfam[dname][0]
                    domain = pfam if pfam is not None else dname
                    domain_name = 'pfam' if pfam is not None else 'domain_name'
                    dom = intera.Domain(uniprots[0], domain = domain, 
                        domain_id_type = domain_name, start = regions[0], 
                        end = regions[1])
                    mot = intera.Motif(uniprots[1], regions[2], regions[3], 
                        instance = dm['instance'], regex = dm['regex'], 
                        motif_name = mname)
                    ptm = intera.Ptm(uniprots[1], motif = mot, source = '3DID')
                    dommot = intera.DomainMotif(dom, ptm, sources = '3DID')
                    new = False
                dommot.add_pdbs(dm['pdb'])
            dmi2[uniprots].append(dommot)
    prg.terminate()
    return dmi2

def get_instruct():
    '''
    Instruct contains residue numbers in UniProt sequences, it means
    no further calculations of offsets in chains of PDB structures needed.
    Chains are not given, only a set of PDB structures supporting the 
    domain-domain // protein-protein interaction.
    '''
    non_digit = re.compile(r'[^\d.-]+')
    data = curl(data_formats.urls['instruct_human']['url'],silent=False)
    if data is None:
        return None
    data = data.replace('\r','').split('\n')
    del data[0]
    instruct = []
    for l in data:
        l = l.split('\t')
        if len(l) > 12:
            domain1 = l[6]
            domain2 = l[7]
            pdb = l[12].split(';')
            uniprot1 = l[0]
            uniprot2 = l[1]
            seq1 = [[non_digit.sub('',n) for n in s.split(',')] for s in l[10].split(';')]
            seq2 = [[non_digit.sub('',n) for n in s.split(',')] for s in l[11].split(';')]
            instruct.append({
                    uniprot1: {
                        'pfam': domain1,
                        'chain': None,
                        'seq': seq1
                        },
                    uniprot2: {
                        'pfam': domain2,
                        'chain': None,
                        'seq': seq2
                    },
                    'uniprots': [uniprot1,uniprot2],
                    'source': 'Instruct',
                    'pdb': pdb,
                    'references': l[13].split(';')
                    })
    return instruct

def get_instruct_offsets():
    '''
    These offsets should be understood as from UniProt to PDB.
    '''
    non_digit = re.compile(r'[^\d.-]+')
    data = curl(data_formats.urls['instruct_offsets']['url'],silent=False)
    if data is None:
        return None
    data = data.replace('\r','').split('\n')
    del data[0]
    offsets = {}
    for l in data:
        l = l.split('\t')
        if len(l) > 2:
            pdb = l[0].lower()
            uniprot = l[1]
            try:
                offset = int(non_digit.sub('',l[2]))
                offsets[(pdb,uniprot)] = offset
            except:
                print l[2]
                print l
    return offsets

def get_i3d():
    '''
    Interaction3D contains residue numbers in given chains in 
    given PDB stuctures, so we need to add an offset to get the residue 
    numbers valid for UniProt sequences. Offsets can be obtained from
    Instruct, or from the Pfam PDB-chain-UniProt mapping table.
    '''
    dname_pfam, pfam_dname = get_pfam_names()
    if dname_pfam is None:
        sys.stdout.write('\n\t:: Could not get Pfam domain names\n\n')
    non_digit = re.compile(r'[^\d.-]+')
    data = curl(data_formats.urls['i3d_human']['url'],silent=False)
    if data is None:
        return None
    data = data.replace('\r','').split('\n')
    del data[0]
    i3d = []
    prg = progress.Progress(len(data), 'Processing domain-domain interactions', 11)
    for l in data:
        prg.step()
        l = l.split('\t')
        if len(l) > 20:
            domain1 = None if l[13] not in dname_pfam else dname_pfam[l[13]]
            domain2 = None if l[20] not in dname_pfam else dname_pfam[l[20]]
            pdb = l[5]
            uniprot1 = l[0]
            uniprot2 = l[1]
            chain1 = l[7]
            seq1 = [[int(non_digit.sub('',l[11])),int(non_digit.sub('',l[12]))]]
            chain2 = l[14]
            seq2 = [[int(non_digit.sub('',l[18])),int(non_digit.sub('',l[19]))]]
            i3d.append({
                    uniprot1: {
                        'pfam': domain1,
                        'chain': chain1,
                        'seq': seq1
                        },
                    uniprot2: {
                        'pfam': domain2,
                        'chain': chain2,
                        'seq': seq2
                    },
                    'uniprots': [uniprot1,uniprot2],
                    'source': 'I3D',
                    'pdb': [pdb],
                    'references': []
                })
    prg.terminate()
    return i3d

def get_switches_elm():
    '''
    switches.elm is a resource containing functional switches in molecular regulation,
    in domain-motif level resolution, classified into categories according to their
    mechanism.
    '''
    residue = re.compile(r'(^[A-Z])([0-9]+)')
    url = data.formats.urls['switches.elm']['url']
    data = curl(url,silent=False)
    if data is None:
        return None
    buff = StringIO()
    buff.write(data)
    cols = {
        'intramol': 3, 
        'bindingsite_a': 5,
        'bs_a_start': 6,
        'bs_a_end': 7,
        'uniprot_a': 4,
        'uniprot_b': 8,
        'bindingsite_b': 9,
        'bs_b_start': 10,
        'bs_b_end': 11,
        'affected': 12,
        'type': 13,
        'subtype': 14,
        'mechanism': 15,
        'reversible': 16,
        'outcome': 17,
        'outcomedir': 18,
        'modification': 19,
        'modsites': 20,
        'modifiers': 21,
        'effectors': 22,
        'references': 26
        }
    table = read_table(cols=cols,fileObject = buff, sep2 = subf, hdr = 1)
    mod_ont = get_ontology('MOD')
    for l in table:
        if l['modification'].startswith('MOD'):
            if l['modification'] in mod_ont:
                l['modification'] = mod_ont[l['modification']]
        l['references'] = [x.replace('PMID:','').strip() for x in l['references']]
        l['modsites'] = [(m.group(2),m.group(1)) for m in \
                         [residue.match(s.strip()) for s in l['modsites'].split(';')]]
        l['intramol'] = True if l['intramol'].strip() == 'TRUE' else False
        l['bs_a_start'] = [x.split(';') for x in l['bs_a_start'].strip()]
        l['bs_b_start'] = [x.split(';') for x in l['bs_b_start'].strip()]
        l['bs_a_end'] = [x.split(';') for x in l['bs_a_end'].strip()]
        l['bs_b_end'] = [x.split(';') for x in l['bs_b_end'].strip()]
        l['bindingsite_a'] = [x.split(';') for x in l['bindingsite_a'].strip()]
        l['bindingsite_b'] = [x.split(';') for x in l['bindingsite_b'].strip()]
        l['modifiers'] = [x.split(':') for x in l['modifiers'].strip().split(';')]
        bs_a_ids = {}
        bs_b_ids = {}
        mod_ids = {}
        for bs in l['bindingsite_a'].split(';'):
            if ':' in bs:
                bs = bs.split(':')
                if bs[0].lower() not in bs_a_ids:
                    bs_a_ids[bs[0].lower()] = []
                bs_a_ids[bs[0].lower()].append(bs[1])
        for bs in l['bindingsite_b'].split(';'):
            if ':' in bs:
                bs = bs.split(':')
                if bs[0].lower() not in bs_b_ids:
                    bs_b_ids[bs[0].lower()] = []
                bs_b_ids[bs[0].lower()].append(bs[1])
        for mod in l['modifiers'].split(';'):
            if ':' in mod:
                mod = mod.split(':')
                if mod[0].lower() not in mod_ids:
                    mod_ids[mod[0].lower()] = []
                mod_ids[mod[0].lower()].append(mod[1])
        l['bindingsite_a'] = bs_a_ids
        l['bindingsite_b'] = bs_b_ids
        l['modifiers'] = mod_ids
    return table

def get_csa(uniprots=None):
    url = data_formats.urls['catalytic_sites']['url']
    data = curl(url,silent=False)
    if data is None:
        return None
    u_pdb, pdb_u = get_pdb_chains()
    buff = StringIO()
    buff.write(data)
    cols = {
        'pdb': 0, 
        'id': 1,
        'resname': 2,
        'chain': 3,
        'resnum': 4,
        'chem_fun': 5,
        'evidence': 6,
        }
    table = read_table(cols = cols, fileObject = buff, sep = ',', hdr = 1)
    css = {}
    prg = progress.Progress(len(table),'Processing catalytic sites',11)
    for l in table:
        if l['pdb'] in pdb_u:
            if l['chain'] in pdb_u[l['pdb']]:
                uniprot = pdb_u[l['pdb']][l['chain']]['uniprot']
                if uniprots is None or uniprot in uniprots:
                    offset = pdb_u[l['pdb']][l['chain']]['offset']
                    if offset is not None:
                        l['resnum'] = int(l['resnum']) + offset
                    else:
                        this_res = residue_pdb(l['pdb'],l['chain'],l['resnum'])
                        if len(this_res) > 0:
                            l['resnum'] = int(this_res['UPCOUNT'])
                        else:
                            l['resnum'] = None
                    if l['resnum'] is not None:
                        if uniprot not in css:
                            css[uniprot] = {}
                        if l['pdb'] not in css[uniprot]:
                            css[uniprot][l['pdb']] = {}
                        if l['id'] not in css[uniprot][l['pdb']]:
                            css[uniprot][l['pdb']][l['id']] = []
                        css[uniprot][l['pdb']][l['id']].append(
                            intera.Residue(l['resname'],l['resnum'],uniprot))
        prg.step()
    prg.terminate()
    return css

def get_ontology(ontology):
    ols = WSDLService("OLS", data_formats.urls['ols']['url'])
    ont = dict((x.key,x.value) for x in ols.serv.getAllTermsFromOntology(ontology).item)
    return ont

def get_listof_ontologies():
    ols = WSDLService("OLS", data_formats.urls['ols']['url'])
    olist = dict((x.key,x.value) for x in ols.serv.getOntologyNames().item)
    return olist

def residue_pdb(pdb,chain,residue):
    url = data_formats.urls['pdbsws']['url']
    params = urllib.urlencode({'plain': 1, 'qtype': 'pdb', 
                                'id': pdb, 'chain': chain, 'res': residue})
    data = urllib2.urlopen(url + "?%s"%params)
    result = {}
    for l in data:
        if not l.startswith('//'):
            l = [x.strip() for x in l.split(':')]
            result[l[0]] = l[1]
    return result

class ResidueMapper(object):
    
    """
    This class stores and serves the PDB --> UniProt 
    residue level mapping. Attempts to download the 
    mapping, and stores it for further use. Converts 
    PDB residue numbers to the corresponding UniProt ones.
    """
    def __init__(self):
        self.clean()
    
    def load_mapping(self,pdb):
        non_digit = re.compile(r'[^\d.-]+')
        pdb = pdb.lower()
        url = data_formats.urls['pdb_align']['url'] + pdb
        data = urllib2.urlopen(url)
        mapper = {}
        soup = bs4.BeautifulSoup(data.read(), 'html.parser')
        for block in soup.find_all('block'):
            seg = block.find_all('segment')
            chain = seg[0]['intobjectid'].split('.')[1]
            uniprot = seg[1]['intobjectid']
            pdbstart = int(non_digit.sub('', seg[0]['start']))
            pdbend = int(non_digit.sub('', seg[0]['end']))
            uniprotstart = int(non_digit.sub('', seg[1]['start']))
            uniprotend = int(non_digit.sub('', seg[1]['end']))
            if chain not in mapper:
                mapper[chain] = {}
            mapper[chain][pdbend] = {
                'uniprot': uniprot, 
                'pdbstart': pdbstart,
                'uniprotstart': uniprotstart,
                'uniprotend': uniprotend}
        self.mappers[pdb] = mapper
    
    def get_residue(self,pdb,resnum,chain=None):
        pdb = pdb.lower()
        if pdb not in self.mappers:
            self.load_mapping(pdb)
        if pdb in self.mappers:
            for chain,data in self.mappers[pdb].iteritems():
                pdbends = data.keys()
                if resnum <= max(pdbends):
                    pdbend = min([x for x in [
                            e - resnum for e in pdbends
                        ] if x >= 0] ) + resnum
                    seg = data[pdbend]
                    if seg['pdbstart'] <= resnum:
                        offset = seg['uniprotstart'] - seg['pdbstart']
                        residue = {
                            'resnum': resnum + offset,
                            'offset': offset,
                            'uniprot': seg['uniprot'],
                            'chain': chain
                        }
                        return residue
        return None
    
    def clean(self):
        '''
        Removes cached mappings, freeing up memory.
        '''
        self.mappers = {}

def get_comppi():
    url = data_formats.urls['comppi']['url']
    post = {
        'fDlSet': 'comp',
        'fDlSpec': '0',
        'fDlMLoc': 'all',
        'fDlSubmit': 'Download'
    }
    data = curl(url, post = post, silent = False, compr = 'gz')
    cols = {
        'uniprot1': 0,
        'uniprot2': 8,
        'loc1': 2,
        'loc2': 10,
        'loc_score': 16
    }
    buff = StringIO()
    buff.write(data)
    data = read_table(cols = cols, fileObject = buff, hdr = 1, sep = '\t')
    return data

def get_psite_phos(raw = True, organism = 'human'):
    url = data_formats.urls['psite_kin']['url']
    data = curl(url, silent = False, compr = 'gz', encoding = 'iso-8859-1')
    cols = {
        'kinase': 1,
        'kinase_org': 3,
        'substrate': 6,
        'substrate_org': 8,
        'residue': 9,
        'motif': 11
    }
    buff = StringIO()
    buff.write(data)
    data = read_table(cols = cols, fileObject = buff, sep = '\t', hdr = 4)
    result = []
    non_digit = re.compile(r'[^\d.-]+')
    motre = re.compile(r'(_*)([A-Za-z]+)(_*)')
    for r in data:
        if organism is None or \
            (r['kinase_org'] == organism and r['substrate_org'] == organism):
            r['resaa'] = r['residue'][0]
            r['resnum'] = int(non_digit.sub('', r['residue'][1:]))
            mot = motre.match(r['motif'])
            isoform = 1 if '-' not in r['substrate'] else r['substrate'].split('-')[1]
            r['substrate'] = r['substrate'].split('-')[0]
            if mot:
                r['start'] = r['resnum'] -7 + len(mot.groups()[0])
                r['end'] = r['resnum'] + 7 - len(mot.groups()[2])
                r['instance'] = r['motif'].replace('_', '').upper()
            else:
                r['start'] = None
                r['end'] = None
                r['instance'] = None
            if raw:
                result.append(r)
            else:
                res = intera.Residue(r['resnum'], r['resaa'], r['substrate'])
                mot = intera.Motif(r['substrate'], r['start'], r['end'], 
                instance = r['instance'])
                ptm = intera.Ptm(protein = r['substrate'], residue = res, motif = mot, 
                    typ = 'phosphorylation', source = 'PhosphoSite')
                dom = intera.Domain(protein = r['kinase'])
                dommot = intera.DomainMotif(domain = dom, ptm = ptm, 
                    sources = ['PhosphoSite'])
                result.append(dommot)
    return result

def get_psite_p(organism = 'human'):
    result = []
    url = data_formats.urls['psite_p']['url']
    data = curl(url, silent = False)
    data = [r.split('\t') for r in data.split('\n')[4:]]
    nondigit = re.compile(r'[^\d]+')
    remot = re.compile(r'(_*)([A-Za-z]+)(_*)')
    for r in data:
        if len(r) > 9 and (organism is None or r[6] == organism):
            uniprot = r[1]
            isoform = 1 if '-' not in uniprot else int(uniprot.split('-')[1])
            uniprot = uniprot.split('-')[0]
            typ = r[3].lower()
            if len(typ) == 0:
                typ = r[4].split('-')[1] if '-' in r[4] else None
            aa = r[4][0]
            num = int(nondigit.sub('', r[4]))
            motif = remot.match(r[9])
            if motif:
                start = num -7 + len(motif.groups()[0])
                end = num + 7 - len(motif.groups()[2])
                instance = r[9].replace('_', '').upper()
            else:
                start = None
                end = None
                instance = None
            res = intera.Residue(num, aa, uniprot, isoform = isoform)
            mot = intera.Motif(uniprot, start, end, instance = instance, isoform = isoform)
            ptm = intera.Ptm(uniprot, typ = typ, motif = mot, residue = res, 
                source = 'PhosphoSite', isoform = isoform)
            result.append(ptm)
    return result

def get_psite_reg():
    url = data_formats.urls['psite_reg']['url']
    data = curl(url, silent = False, compr = 'gz', encoding = 'iso-8859-1')
    cols = {
        'uniprot': 2,
        'organism': 6,
        'mod': 7,
        'on_function': 11,
        'on_process': 12,
        'on_interact': 13,
        'pmids': 15
    }
    buff = StringIO()
    buff.write(data)
    data = read_table(cols = cols, fileObject = buff, sep = '\t', hdr = 4)
    regsites = {}
    for r in data:
        interact = [[y.replace(')', '').strip() for y in x.split('(')] \
            for x in r['on_interact'].strip().split(';') if len(x) > 0]
        induces = [x[0] for x in interact if x[1] == 'INDUCES']
        disrupts = [x[0] for x in interact if x[1] == 'DISRUPTS']
        mod = r['mod']
        modt = r['mod'].split('-')
        mod = list(modt[0])
        aa = mod.pop(0)
        modt = modt[1]
        res = ''.join(mod)
        if r['uniprot'] not in regsites:
            regsites[r['uniprot']] = []
        regsites[r['uniprot']].append({
            'aa': aa,
            'res': res,
            'modt': modt,
            'organism': r['organism'],
            'pmids': [x.strip() for x in r['pmids'].split(';')],
            'induces': induces,
            'disrupts': disrupts
        })
    return regsites

def regsites_tab(regsites, mapper, outfile = None):
    header = ['uniprot_a', 'isoform_a', 'a_res_aa', 'a_res_num', 
        'a_mod_type', 'effect', 'uniprot_b', 'references']
    result = []
    for uniprot, regsite in regsites.iteritems():
        isoform = '1'
        uniprot = uniprot.split('-')
        if len(uniprot) > 1:
            isoform = uniprot[1]
        uniprot = uniprot[0]
        for r in regsite:
            if r['organism'] == 'human':
                for i in r['induces']:
                    other = mapper.map_name(i,'genesymbol','uniprot')
                    for o in other:
                        if o != 'unmapped':
                            result.append([uniprot, isoform, r['aa'], 
                                r['res'], r['modt'], '+', o])
                for i in r['disrupts']:
                    other = mapper.map_name(i,'genesymbol','uniprot')
                    for o in other:
                        if o != 'unmapped':
                            result.append([uniprot, isoform, r['aa'], 
                                r['res'], r['modt'], '-', o, ';'.join(r['pmids'])])
    if outfile is not None:
        out = '\t'.join(header) + '\n'
        for r in result:
            out += '\t'.join(r) + '\n'
        with open(outfile, 'w') as f:
            f.write(out)
    return result

def get_ielm_huge(ppi, id_type = 'UniProtKB_AC', mydomains = 'HMMS', 
        maxwait = 180, cache = True, part_size = 500, headers = None):
    ranges = range(0, len(ppi), part_size)
    result = []
    done = False
    while not done:
        for r in ranges:
            this_ppi = ppi[r:r+part_size]
            sys.stdout.write('\t:: Part %u/%u: querying %u interactions.\n' % \
                (ranges.index(r) + 1, len(ranges), len(this_ppi)))
            sys.stdout.flush()
            this_res = get_ielm(this_ppi, id_type, mydomains, maxwait, 
                cache, part = True, headers = headers)
            if this_res:
                if type(this_res) is dict:
                    print type(this_res)
                    return this_res
                result += this_res
                if r == ranges[-1]:
                    done = True
            else:
                part_size = max(int(part_size * 0.8), 20)
                ranges = range(r, len(ppi[r:]), part_size)
                sys.stdout.write('\t:: One query failed. Setting part size to %u\n' % \
                    part_size)
                sys.stdout.flush()
                break
    return result

def get_ielm(ppi, id_type = 'UniProtKB_AC', mydomains = 'HMMS', 
        maxwait = 180, cache = True, part = False, part_size = 500, headers = None):
    url = data_formats.urls['proteomic_ielm']['url']
    network = ''
    from_pickle = []
    ppi_pickle = []
    ppi_query = []
    result = []
    pcache = os.path.join('cache', 'ielm.pickle')
    if not part and os.path.exists(pcache):
        from_pickle = pickle.load(open(pcache, 'rb'))
        ppi_pickle = from_pickle['ppi']
        ppi_query = list(set(ppi) - set(ppi_pickle))
        result = from_pickle['ielm']
        if len(ppi_query) == 0:
            return result
    else:
        ppi_query = ppi
    if len(ppi_query) > part_size and not part:
        this_result = get_ielm_huge(ppi_query, id_type, mydomains, maxwait, cache, 
            part_size, headers)
    for pp in ppi_query:
        network += '%s %s\r\n' % (pp[0], pp[1])
    post = {'network': network, 'databases': id_type, 'mydomains': mydomains}
    net_md5 = hashlib.md5(network).hexdigest()
    cachefile = os.path.join(os.getcwd(),'cache',net_md5 + '.ielm')
    if os.path.exists(cachefile) and cache:
        with open(cachefile, 'r') as f:
            data = f.read()
        soup = bs4.BeautifulSoup(data, 'html.parser')
        src = 'cache'
    else:
        data = curl(url, post = post, silent = False, cache = False, 
            req_headers = headers)
        soup = bs4.BeautifulSoup(data, 'html.parser')
        sessid = soup.find('input', {'name':'session_ID'})['value']
        src = 'iELM'
    if data is None:
        print ERASE_LINE + CURSOR_UP_ONE
        sys.stdout.write('\t:: Initial query failed. No data retrieved from iELM.\n')
        sys.stdout.flush()
        return None
    wait = 0
    while soup.title.text == 'iELM Wait Page' and wait < maxwait:
            #and \
            #len([i for i in soup.find_all('font', {'color': '#FF0000'}) if i.text == \
            #'File no longer available']) == 0:
        print ERASE_LINE + CURSOR_UP_ONE
        sys.stdout.write('\t:: Waiting for result. Wait time: %u sec. '\
            'Max waiting time: %u sec.' % (wait,maxwait))
        sys.stdout.flush()
        post = {
            'session_ID': sessid,
            'database': id_type, 
            'number': '', 
            'domains': mydomains}
        data = curl('http://i.elm.eu.org/wait_2/', post = post, cache = False, 
            req_headers = headers)
        if data is not None:
            soup = bs4.BeautifulSoup(data, 'html.parser')
        time.sleep(3)
        wait += 3
    if len(soup.find_all('table')) == 0:
        print ERASE_LINE + CURSOR_UP_ONE
        sys.stdout.write('\t:: No data retrieved from iELM. \n')
        sys.stdout.flush()
        soup.title.string = 'http://i.elm.eu.org/proteomic_results/%s'%sessid
        # return {'soup': soup, 'post': urllib.urlencode(post), 'netw': network}
        return None
    if cache:
        with open(cachefile, 'w') as f:
            f.write(data)
    print ERASE_LINE + CURSOR_UP_ONE
    sys.stdout.write('\t:: Data retrieved from %s in %u seconds.\n' % (src, wait))
    sys.stdout.flush()
    tbl = soup.find('table', {'id': 'example1'})
    this_result = []
    if tbl:
        url = data_formats.urls['elm_depr']['url']
        depr_list = curl(url)
        depr_list = depr_list.replace('"', '').split('\n')[5:]
        depr = [tuple(x.split('\t')) for x in depr_list if len(x) > 0]
        depr = dict(depr + [tuple([x[0].lower(), x[1]]) for x in depr])
        # redepr = re.compile(r'\b(' + '|'.join(depr.keys()) + r')\b') :(
        rows = tbl.find_all('tr')
        prg = progress.Progress(len(rows), 'Processing data (%u rows)'%(len(rows)-1), 3)
        for tr in tbl.find_all('tr'):
            thisRow = [td.text.strip() for td in tr.find_all('td')]
            if len(thisRow) > 15 and not thisRow[0].startswith('Motif'):
                # replacing deprecated ELM names:
                if thisRow[2].lower() in depr:
                    thisRow[2] = depr[thisRow[2].lower()]
                if thisRow[2].lower() in depr:
                    thisRow[2] = depr[thisRow[2].lower()]
                # thisRow[2] = redepr.sub(lambda x: depr[x.group()], thisRow[2]) :(
                this_result.append(thisRow)
            prg.step()
        prg.terminate()
    if not part:
        result = {'ppi': list(set(ppi_pickle + ppi_query)), 
            'ielm': result + this_result}
        pickle.dump(result, open(pcache, 'wb'))
    return this_result

def get_pepcyber(cache = None):
    url = data_formats.urls['pepcyber']['url']
    # this is huge, takes a few minutes! 
    data = curl(url, silent = False, timeout = 600)
    soup = bs4.BeautifulSoup(data, 'html.parser')
    rows = soup.find_all('table')[6].find_all('tr')
    result = []
    uniprots = {}
    if cache is None:
        cache = os.path.join('cache', 'pepcyber-uniprots')
    if os.path.exists(cache):
        with open(cache, 'r') as f:
            for l in f:
                l = l.split('\t')
                l += ['', '']
                uniprots[l[0].strip()] = [ l[1].strip(), l[2].strip()]
    prg = progress.Progress(len(rows), 'Retrieving and processing data', 7)
    notfound = []
    for r in rows:
        prg.step()
        thisRow = [c.text.strip() for c in r.find_all('td')]
        if len(thisRow) > 9 and thisRow[5].isdigit():
            inum = int(r.find('a')['name'])
            thisRow[9] = None if 'p' not in thisRow[4] else \
                thisRow[4][thisRow[4].index('p') + 1]
            if thisRow[2] not in uniprots or thisRow[3] not in uniprots:
                # print "Query: %s, %s, %u" % (thisRow[2], thisRow[3], inum)
                uniprots = dict(uniprots.items() + pepcyber_uniprot(inum).items())
            if thisRow[2] in uniprots and thisRow[3] in uniprots:
                thisRow += uniprots[thisRow[2]] + uniprots[thisRow[3]]
                result.append(thisRow[1:])
            else:
                notfound.append([thisRow[2], thisRow[3], inum])
    prg.terminate()
    with open(cache, 'w') as f:
        for g, u in uniprots.iteritems():
            try:
                f.write('\t'.join([g] + u) + '\n')
            except:
                pass
    return result

def pepcyber_uniprot(num):
    result = {}
    url = data_formats.urls['pepcyber_details']['url'] % num
    data = curl(url, cache = False)
    if data is None:
        return result
    soup = bs4.BeautifulSoup(data, 'html.parser')
    gname = None
    prev = ''
    for td in soup.find_all('td'):
        if prev.startswith('Gene name'):
            gname = td.text.strip().split('(')[0]
        if prev.startswith('RefSeq'):
            refseq = td.text.strip()
        if prev.startswith('SwissProt') and gname is not None:
            swprot = td.text.strip()
            if len(gname) > 0:
                # print "Result: %s, %s, %u" % (gname, swprot, num)
                result[gname] = [swprot, refseq]
            gname = None
        prev = td.text.strip()
    return result

def get_pdzbase():
    url = data_formats.urls['pdzbase']['url']
    data = curl(url, silent = False)
    soup = bs4.BeautifulSoup(data, 'html.parser')
    rows = soup.find_all('table')[3].find('table').find('table').find_all('tr')
    result = []
    for r in rows:
        thisRow = [c.text.strip() for c in r.find_all('td')]
        result.append(thisRow)
    del result[0]
    return result

def get_domino(none_values = False, outfile = None):
    result = []
    taxid = re.compile(r'taxid:(.*)\([a-zA-Z ]*\)')
    miont = re.compile(r'MI:[0-9]{4}\((.*)\)')
    binds = re.compile(r'([-0-9]*);.*')
    domai = re.compile(r'.*;.*;.*\((.*)\)')
    dipro = re.compile(r'.*;.*;.+:(IPR[0-9]*).*')
    ptmrs = re.compile(r'([-0-9]*);.*')
    ptmmi = re.compile(r'[0-9]*;(MI:[0-9]*)\(.*\);.*;.*')
    ptmrn = re.compile(r'.*sequence:[\s]*[0-9]+-[0-9]+[\s]*:[\s]*([A-Z]{10,}).*')
    ptmty = re.compile(r'[0-9]*;MI:[0-9]*\((.*)\);.*;.*')
    refrs = re.compile(r'(pubmed|doi):["]*([-0-9a-zA-Z\.\(\)/]*)["]*')
    url = data_formats.urls['domino']['url']
    data = curl(url, silent = False)
    data = data.split('\n')
    data.remove('')
    del data[0]
    header = ['uniprot-A', 'uniprot-B', 'isoform-A', 'isoform-B', 
              'exp. method', 'references', 'taxon-A', 'taxon-B', 
              'role-A', 'role-B', 'binding-site-range-A', 'binding-site-range-B', 
              'domains-A', 'domains-B', 'ptm-residue-A', 'ptm-residue-B', 
              'ptm-type-mi-A', 'ptm-type-mi-B', 'ptm-type-A', 'ptm-type-B', 
              'ptm-res-name-A', 'ptm-res-name-B', 'mutations-A', 'mutations-B', 
              'mutation-effects-A', 'mutation-effects-B', 'domains-interpro-A', 
              'domains-interpro-B', 'negative']
    for r in data:
        r = r.split('\t')
        thisRow = [
            None if ':' not in r[0] else r[0].split(':')[1].split('-')[0], 
            None if ':' not in r[1] else r[1].split(':')[1].split('-')[0], 
            '1' if '-' not in r[0] else r[0].split('-')[1], 
            '1' if '-' not in r[1] else r[1].split('-')[1], 
            None if miont.match(r[6]) is None else miont.match(r[6]).groups(1)[0],
            None if refrs.match(r[8]) is None else refrs.match(r[8]).groups(1)[1],
            None if taxid.match(r[9]) is None else taxid.match(r[9]).groups(1)[0],
            None if taxid.match(r[10]) is None else taxid.match(r[10]).groups(1)[0],
            None if miont.match(r[11]) is None else miont.match(r[11]).groups(1)[0],
            None if miont.match(r[16]) is None else miont.match(r[17]).groups(1)[0],
            ';'.join(['' if binds.match(x) is None else binds.match(x).groups(1)[0] \
                for x in r[32].split(',')]),
            ';'.join(['' if binds.match(x) is None else binds.match(x).groups(1)[0] \
                for x in r[33].split(',')]),
            ';'.join(['' if domai.match(x) is None else domai.match(x).groups(1)[0] \
                for x in r[32].split(',')]),
            ';'.join(['' if domai.match(x) is None else domai.match(x).groups(1)[0] \
                for x in r[33].split(',')]),
            ';'.join(['' if ptmrs.match(x) is None else ptmrs.match(x).groups(1)[0] \
                for x in r[34].split('|')]),
            ';'.join(['' if ptmrs.match(x) is None else ptmrs.match(x).groups(1)[0] \
                for x in r[35].split('|')]),
            ';'.join(['' if ptmmi.match(x) is None else ptmmi.match(x).groups(1)[0] \
                for x in r[34].split('|')]),
            ';'.join(['' if ptmmi.match(x) is None else ptmmi.match(x).groups(1)[0] \
                for x in r[35].split('|')]),
            ';'.join(['' if ptmty.match(x) is None else ptmty.match(x).groups(1)[0] \
                for x in r[34].split('|')]),
            ';'.join(['' if ptmty.match(x) is None else ptmty.match(x).groups(1)[0] \
                for x in r[35].split('|')]),
            ';'.join(['' if ptmrn.match(x) is None else ptmrn.match(x).groups(1)[0] \
                for x in r[34].split('|')]),
            ';'.join(['' if ptmrn.match(x) is None else ptmrn.match(x).groups(1)[0] \
                for x in r[35].split('|')]),
            ';'.join(['' if ptmrs.match(x) is None else ptmrs.match(x).groups(1)[0] \
                for x in r[36].split('|')]),
            ';'.join(['' if ptmrs.match(x) is None else ptmrs.match(x).groups(1)[0] \
                for x in r[37].split('|')]),
            ';'.join(['' if ptmty.match(x) is None else ptmty.match(x).groups(1)[0] \
                for x in r[36].split('|')]),
            ';'.join(['' if ptmty.match(x) is None else ptmty.match(x).groups(1)[0] \
                for x in r[37].split('|')]),
            '' if dipro.match(r[32]) is None else dipro.match(r[32]).groups(1)[0],
            '' if dipro.match(r[33]) is None else dipro.match(r[33]).groups(1)[0],
            '0' if r[38].strip() == '-' else '1'
        ]
        if not none_values:
            thisRow = ['' if x is None else x for x in thisRow]
        result.append(thisRow)
    if outfile:
        with open(outfile, 'w') as outf:
            outf.write('\t'.join(header) + '\n')
            for r in result:
                outf.write('\t'.join(['' if x is None else x for x in r]) + '\n')
    return result

def get_elm_domains():
    result = {}
    url = data_formats.urls['ielm_domains']['url']
    data = curl(url, silent = False)
    soup = bs4.BeautifulSoup(data, 'html.parser')
    tbl = soup.find('table').find_all('td')
    rows = [tbl[x:x+4] for x in xxrange(0, len(tbl), 4)]
    for r in rows:
        uniprot = r[1].text
        motif = r[0].text
        if uniprot not in result:
            result[uniprot] = {}
        if motif not in result[uniprot]:
            result[uniprot][motif] = []
        result[uniprot][motif].append((r[2].text,r[3].text))
    return result

def get_phosphoelm(organism = 'Homo sapiens', ltp_only = True):
    result = []
    non_digit = re.compile(r'[^\d.-]+')
    url = data_formats.urls['p_elm']['url']
    data = curl(url, silent = False)
    data = [n for d, n in data.iteritems() \
        if d.startswith(data_formats.urls['p_elm']['psites'])]
    data = data[0] if len(data) > 0 else ''
    data = [l.split('\t') for l in data.split('\n')]
    kinases = get_phelm_kinases()
    del data[0]
    for l in data:
        if len(l) == 9 and l[7] == organism and (not ltp_only or l[6] == 'LTP'):
            l[1] = 1 if '-' not in l[0] else int(l[0].split('-')[1])
            l[0] = l[0].split('-')[0]
            del l[-1]
            if len(l[5]) > 0 and l[5] in kinases:
                kinase = kinases[l[5]]
            result.append({
                'instance': None,
                'isoform': l[1],
                'resaa': l[3],
                'resnum': int(non_digit.sub('', l[2])),
                'start': None,
                'end': None,
                'substrate': l[0],
                'kinase': kinase,
                'references': l[4].split(';'),
                'experiment': l[6],
                'organism': l[7]
            })
    return result

def phelm_interactions(organism = 'Homo sapiens'):
    result = []
    data = get_phosphoelm(ltp_only = True)
    for l in data:
        result.append([l['kinase'], l['substrate'], ';'.join(l['references']),
            l['organism']])
    return result

def phelm_psites():
    result = []
    data = get_phosphoelm()
    kinases = get_phelm_kinases()
    for l in data:
        l.append('1' if '-' not in l[0] else l[0].split('-')[1])
        l[0] = l[0].split('-')[0]
        l.append('' if l[4] not in kinases else kinases[l[4]])
        result.append(l)
    return result

def get_phelm_kinases():
    result = {}
    url = data_formats.urls['p_elm_kin']['url']
    data = curl(url, silent = False)
    soup = bs4.BeautifulSoup(data, 'html.parser')
    for row in soup.find_all('table')[1].find_all('tr'):
        thisRow = [x.text for x in row.find_all('td')]
        if len(thisRow) > 2 and len(thisRow[2].strip()) > 0:
            result[thisRow[0]] = thisRow[2].strip()
    return result

def get_elm_classes():
    url = data_formats.urls['elm_class']['url']
    data = curl(url, silent = False)
    data = [x.split('\t')[:-2] for x in \
        data.replace('"', '').split('\n')[6:] if len(x) > 0]
    return dict(zip([x[1] for x in data], data))

def get_elm_instances():
    url = data_formats.urls['elm_inst']['url']
    data = curl(url, silent = False)
    data = data.replace('"', '').split('\t')
    data = data[6:]

def get_elm_interactions():
    '''
    Downlods manually curated interactions from ELM.
    This is the gold standard set of ELM.
    '''
    result = []
    url = data_formats.urls['elm_int']['url']
    data = curl(url, silent = False)
    data = data.split('\n')
    del data[0]
    for l in data:
        result.append([x.strip() for x in l.split('\t')])
    return result

def pfam_uniprot(uniprots, infile = None):
    result = {}
    url = data_formats.urls['pfam_up']['url']
    infile = infile if infile is not None \
        else os.path.join('cache', 'pfam-regions.tab.gz')
    if not os.path.exists(infile):
        sys.stdout.write('\t:: Downloading data...\n')
        sys.stdout.flush()
        urllib.urlretrieve(url, infile)
    sys.stdout.write('\t:: Processing domains...\n')
    sys.stdout.flush()
    gzfile = gzip.open(infile, mode='r')
    prg = progress.Progress(len(uniprots), 'Looking up domains', 1)
    for l in gzfile:
        l = l.split('\t')
        if l[0] in uniprots:
            prg.step()
            if l[0] not in result:
                result[l[0]] = {}
            if l[4] not in result[l[0]]:
                result[l[0]][l[4]] = []
            result[l[0]][l[4]].append([l[1], l[5],l[6]])
    prg.terminate()
    return result

def get_dbptm():
    result = []
    byre = re.compile(r'.*by\s([A-Za-z0-9\s]+)\.*')
    andre = re.compile(r',|and')
    non_digit = re.compile(r'[^\d.-]+')
    for url in data_formats.urls['dbptm']['urls']:
        extra = curl(url, silent = False)
        for k, data in extra.iteritems():
            data = [x.split('\t') for x in data.split('\n')]
            for l in data:
                if len(l) > 8:
                    resnum = int(non_digit.sub('', l[2]))
                    ptm = ({
                        'substrate': l[1],
                        'typ': l[7].lower(),
                        'resaa': l[8][6],
                        'resnum': resnum,
                        'instance': l[8],
                        'references': l[4].split(';'),
                        'source': l[5].split()[0],
                        'kinase': None if byre.match(l[3]) is None else \
                            [i.strip() for i in andre.split(byre.match(l[3]).groups(1)[0])],
                        'start': resnum - 6,
                        'end': resnum + 6
                    })
                    if ptm['kinase'] is not None:
                        if 'autocatalysis' in ptm['kinase']:
                            ptm['kinase'].append(ptm['substrate'])
                            ptm['kinase'].remove('autocatalysis')
                        ptm['kinase'] = [k.replace('host', '').strip() for k in ptm['kinase']]
                        ptm['kinase'] = [k for k in ptm['kinase'] if len(k) > 0]
                        if len(ptm['kinase']) == 0:
                            ptm['kinase'] = None
                    result.append(ptm)
    return result

def dbptm_interactions():
    result = []
    data = get_dbptm()
    for r in data:
        if r['kinase'] is not None:
            for src in r['kinase']:
                result.append([
                    src, 
                    r['substrate'],
                    ';'.join([i for i in r['references'] if i != '-'])
                ])
    return result

def get_phosphonetworks():
    result = []
    reres = re.compile(r'([A-Z])([0-9]+)')
    non_digit = re.compile(r'[^\d.-]+')
    motre = re.compile(r'(-*)([A-Za-z]+)(-*)')
    url = data_formats.urls['phosnw']['url']
    data = curl(url, silent = False)
    if data is None:
        return None
    data = data.split('\n')
    for l in data:
        if l.startswith('>'):
            substrate = l[1:].strip()
        elif len(l.split('\t')) >= 4:
            l = [x.strip() for x in l.split('\t')]
            res = reres.match(l[1])
            resnum = int(non_digit.sub('', res.groups()[1]))
            mot = motre.match(l[0])
            if mot:
                start = resnum - 7 + len(mot.groups()[0])
                end = resnum + 7 - len(mot.groups()[2])
                instance = l[0].replace('-', '').upper()
            else:
                start = None
                end = None
                instance = l[0]
            result.append({
                'instance': instance,
                'kinase': l[2],
                'resaa': res.groups()[0],
                'resnum': resnum,
                'score': float(non_digit.sub('', l[3])),
                'substrate': substrate,
                'start': start,
                'end': end
            })
    return result

def pnetworks_interactions():
    result = []
    data = get_phosphonetworks()
    for l in data:
        result.append((l['kinase'], l['substrate']))
    return [list(x) for x in list(set(result))]

def get_depod(organism = 'Homo sapiens'):
    result = []
    reunip = re.compile(r'uniprotkb:([A-Z0-9]+)')
    url = data_formats.urls['depod']['urls'][0]
    url_mitab = data_formats.urls['depod']['urls'][1]
    data = curl(url, silent = False)
    data_mitab = curl(url_mitab, silent = False)
    data = [x.split('\t') for x in data.split('\n')]
    data_mitab = [x.split('\t') for x in data_mitab.split('\n')]
    del data[0]
    del data_mitab[0]
    for i, l in enumerate(data):
        if len(l) > 6 and l[2] == 'protein substrate' and \
            l[3].strip().startswith(organism) and \
            l[4].strip() != 'N/A':
            result.append([x.strip() for y, x in enumerate(l) if y in [0, 1, 4, 6]] + \
                [reunip.findall(data_mitab[i][0]), reunip.findall(data_mitab[i][1])])
    return result

def get_mimp():
    result = []
    non_digit = re.compile(r'[^\d.-]+')
    motre = re.compile(r'(-*)([A-Za-z]+)(-*)')
    url = data_formats.urls['mimp']['url']
    data = curl(url, silent = False)
    kclass = get_kinase_class()
    if data is None:
        return None
    data = [x.split('\t') for x in data.split('\n')]
    del data[0]
    for l in data:
        if len(l) > 6 and len(l[2]) > 0:
            kinases = l[2].split(';')
            kinases_gnames = []
            for k in kinases:
                if k.endswith('GROUP'):
                    grp = k.split('_')[0]
                    if grp in kclass['groups']:
                        kinases_gnames += kclass['groups'][grp]
                    elif grp in kclass['families']:
                        kinases_gnames += kclass['families'][grp]
                    elif grp in kclass['subfamilies']:
                        kinases_gnames += kclass['subfamilies'][grp]
                else:
                    kinases_gnames.append(k)
            mot = motre.match(l[4])
            for k in kinases_gnames:
                resaa = l[4][7]
                resnum = int(non_digit.sub('', l[3]))
                if mot:
                    start = resnum - 7 + len(mot.groups()[0])
                    end = resnum + 7 - len(mot.groups()[2])
                    instance = l[4].replace('-', '').upper()
                else:
                    start = None
                    end = None
                    instance = l[4]
                result.append({
                    'instance': instance,
                    'kinase': k.upper(),
                    'resaa': resaa,
                    'resnum': resnum,
                    'npmid': int(non_digit.sub('', l[5])),
                    'substrate_refseq': l[1],
                    'substrate': l[0],
                    'start': start,
                    'end': end,
                    'databases': l[6]
                })
    return result

def mimp_interactions():
    result = []
    mimp = get_mimp()
    for m in mimp:
        result.append([m['kinase'], m['substrate']])
    return result

def phosphopoint_directions():
    directions = []
    fname = data_formats.files['phosphopoint']['data']
    with open(fname, 'r') as f:
        nul = f.readline()
        for l in f:
            l = l.split(';')
            directions.append([l[0], l[2]])
    return directions

def get_isoforms(organism = 'Homo sapiens'):
    reorg = re.compile(r'OS=([A-Z][a-z]+\s[a-z]+)')
    result = {}
    url = data_formats.urls['unip_iso']['url']
    data = curl(url, silent = False)
    data = read_fasta(data)
    for header, seq in data.iteritems():
        org = reorg.findall(header)
        if len(org) > 0 and org[0] == organism:
            prot = header.split('|')[1].split('-')
            unip = prot[0]
            isof = int(prot[1])
            if unip not in result:
                result[unip] = {}
            result[unip][isof] = seq
    return result

def read_fasta(fasta):
    result = {}
    fasta = re.split(r'\n>', fasta)
    for section in fasta:
        section = section.strip().split('\n')
        label = section.pop(0)
        seq = ''.join(section)
        result[label] = seq
    return result

def get_kinase_class():
    result = {
        'groups': {},
        'families': {},
        'subfamilies': {},
        'kinases': {}
    }
    tabs = re.compile(r'[\t]{3,}')
    reps = re.compile(r'ps[0-9]*$')
    url = data_formats.urls['kinclass']['url']
    data = curl(url, silent = False)
    data = tabs.sub('', data)
    data = [x.split('\t') for x in data.split('\n')]
    data = data[9:]
    for l in data:
        if len(l) > 4:
            kinase = reps.sub('', l[0])
            group = l[2]
            family = l[3]
            subfamily = l[4]
            if group not in result['groups']:
                result['groups'][group] = []
            result['groups'][group].append(kinase)
            if family not in result['families']:
                result['families'][family] = []
            result['families'][family].append(kinase)
            if subfamily not in result['subfamilies']:
                result['subfamilies'][subfamily] = []
            result['subfamilies'][subfamily].append(kinase)
            result['kinases'][kinase] = {
                'group': group,
                'family': family,
                'subfamily': subfamily
            }
    return result

def get_acsn():
    greek = {
        '_alpha_': 'A',
        '_beta_': 'B',
        '_gamma_': 'C',
        '_delta_': 'D',
        '_epsilon_': 'E'
    }
    regreek = re.compile(r'\b(' + '|'.join(greek.keys()) + r')\b')
    result = []
    url = data_formats.urls['acsn']['url']
    data = curl(url, silent = False)
    data = [x.split('\t') \
        for x in data.replace('\r', '').replace('*', '').strip().split('\n')]
    for l in data:
        l[0] = regreek.sub('', l[0]).split('_')[0].split('~')[0]
        l[2] = regreek.sub('', l[2]).split('_')[0].split('~')[0]
    return data

def get_abs():
    result = []
    url = data_formats.urls['abs']['url']
    data = curl(url, silent = False)
    data = [[x.replace('*', '') for x in xx.split('\t')] for xx in data.split('\n')]
    for d in data:
        if len(d) > 2:
            result.append([d[2], d[0]])
    return result

def get_pazar():
    url = data_formats.urls['pazar']['url']
    data = curl(url, silent = False)
    return [map(x.split('\t').__getitem__, (1, 4, 10)) \
        for x in ''.join(data.values()).split('\n') if len(x) > 0]

def get_htri():
    data = curl(data_formats.urls['htri']['url'], 
        init_url = data_formats.urls['htri']['init_url'], silent = False)
    return [map(x.split(';').__getitem__, (1, 3, 6)) \
        for x in data.split('\n') if len(x) > 0][1:]

def get_oreganno_old(organism = 'Homo sapiens'):
    nsep = re.compile(r'([-A-Za-z0-9]{3,})[\s/\(]*.*')
    nrem = re.compile(r'[-/]')
    result = []
    url = data_formats.urls['oreganno_old']['url']
    data = curl(url, silent = False)
    data = [[xx.strip() for xx in x.split('\t')] for x in data.split('\n') if len(x) > 0][1:]
    for l in data:
        if l[0] == organism and \
            l[10] == 'TRANSCRIPTION FACTOR BINDING SITE' and \
            l[3] == 'POSITIVE OUTCOME' and \
            not l[11].startswith('UNKNOWN') and not l[14].startswith('UNKNOWN'):
            result.append([
                l[14] if len(l[14]) < 3 else nrem.sub('', nsep.findall(l[14])[0]),
                l[11] if len(l[11]) < 3 else nrem.sub('', nsep.findall(l[11])[0]), l[18]])
    return result

def get_oreganno(organism = 'Homo sapiens'):
    nsep = re.compile(r'([-A-Za-z0-9]{3,})[\s/\(]*.*')
    nrem = re.compile(r'[-/]')
    result = []
    url = data_formats.urls['oreganno']['url']
    data = curl(url, silent = False, large = True)
    null = data.readline()
    del null
    for l in data:
        if len(l) > 0:
            l = [x.strip() for x in l.split('\t')]
            if l[1] == organism and \
                l[3] == 'TRANSCRIPTION FACTOR BINDING SITE' and \
                l[2] == 'POSITIVE OUTCOME' and \
                not l[4] == 'N/A' and not l[7] == 'N/A':
                result.append([
                    l[7] if len(l[7]) < 3 else nrem.sub('', nsep.findall(l[7])[0]),
                    l[4] if len(l[4]) < 3 else nrem.sub('', nsep.findall(l[4])[0]), 
                    l[11] if l[11] != 'N/A' else ''])
    return result

def get_cpdb_ltp():
    return get_cpdb(['HPRD', 'BioGRID', 'PhosphoPOINT', 'MINT', 'BIND', 'IntAct'])

def get_cpdb(exclude = None):
    exclude = set(exclude) if type(exclude) is list else exclude
    result = []
    url = data_formats.urls['cpdb']['url']
    data = curl(url, silent = False)
    data = [x.split('\t') for x in data.split('\n') \
        if not x.startswith('#') and len(x) > 0]
    for l in data:
        participants = l[2].split(',')
        if len(participants) == 2:
            if not exclude or len(set(l[0].split(',')) - exclude) > 0:
                result.append([participants[0], participants[1], l[0], l[1]])
    return result

def get_uniprot_sec():
    url = data_formats.urls['uniprot_sec']['url']
    data = curl(url, silent = False)
    data = [x.split() for x in data.split('\n')[30:]]
    return data

def get_go(organism = 9606, swissprot = 'yes'):
    rev = '' if swissprot is None \
            else ' AND reviewed:%s' % swissprot
    query = 'organism:%u%s' % (int(organism), rev)
    url = data_formats.urls['uniprot_basic']['url']
    post = {
        'query': query, 
        'format': 'tab', 
        'columns': 'id,go-id'
    }
    data = curl(url, post = post, silent = False)
    return dict([(x[0], [go.strip() for go in x[1].split(';')]) for x in \
        [x.split('\t') for x in data.split('\n')] if len(x) > 1])

def get_go_goa(organism = 'human'):
    result = {'P': {}, 'C': {}, 'F': {}}
    url = data_formats.urls['goa']['url'] % (organism.upper(), organism)
    data = curl(url, silent = False)
    data = [x.split('\t') for x in data.split('\n') \
        if not x.startswith('!') and len(x) > 0]
    for l in data:
        if l[1] not in result[l[8]]:
            result[l[8]][l[1]] = []
        result[l[8]][l[1]].append(l[4])
    return result

def get_go_quick(organism = 9606, slim = False):
    termuse = 'slim' if slim or slim is None else 'ancestor'
    goslim = '' if termuse == 'ancestor' \
        else '&goid=%s'%','.join(get_goslim(url = slim))
    terms = {'C': {}, 'F': {}, 'P': {}}
    names = {}
    url = data_formats.urls['quickgo']['url'] % (goslim, termuse, organism)
    data = curl(url, silent = False)
    data = [x.split('\t') for x in data.split('\n') if len(x) > 0]
    del data[0]
    for l in data:
        try:
            if l[0] not in terms[l[3][0]]:
                terms[l[3][0]][l[0]] = set([])
            terms[l[3][0]][l[0]].add(l[1])
            names[l[1]] = l[2]
        except:
            print l
    return {'terms': terms, 'names': names}

def get_goslim(url = None):
    rego = re.compile(r'GO:[0-9]{7}')
    url = url if type(url) in [str, unicode] \
        else data_formats.urls['goslim_gen']['url']
    data = curl(url, silent = False)
    result = []
    for l in data.split('\n'):
        if l.startswith('id:'):
            result += rego.findall(l)
    return result

def netpath_names():
    repwnum = re.compile(r'_([0-9]+)$')
    result = {}
    url = data_formats.urls['netpath_names']['url']
    html = curl(url, silent = False)
    soup = bs4.BeautifulSoup(html, 'html.parser')
    for a in soup.find_all('a'):
        if a.attrs['href'].startswith('pathways'):
            num = repwnum.findall(a.attrs['href'])[0]
            name = a.text
            result[num] = name
    return result

def netpath():
    result = []
    repwnum = re.compile(r'NetPath_([0-9]+)_')
    mi = '{net:sf:psidev:mi}'
    url = data_formats.urls['netpath_psimi']['url']
    data = curl(url, silent = False)
    data = dict([(k, v) for k, v in data.iteritems() if k.endswith('xml')])
    pwnames = netpath_names()
    for pwfile, rawxml in data.iteritems():
        try:
            pwnum = repwnum.findall(pwfile)[0]
        except:
            print pwfile
        pwname = pwnames[pwnum]
        root = ET.fromstring(rawxml)
        for e in root.findall(mi+'entry'):
            thisInt = ()
            db = [pr.find(mi+'primaryRef').attrib['db'] \
                for pr in e.find(mi+'source').findall(mi+'xref')]
            refs = []
            mets = []
            for ex in e.find(mi+'experimentList').findall(mi+'experimentDescription'):
                for pm in ex.find(mi+'bibref').iter(mi+'primaryRef'):
                    if pm.attrib['db'] == 'pubmed':
                        refs.append(pm.attrib['id'])
                for me in ex.find(mi+'interactionDetectionMethod').\
                    iter(mi+'shortLabel'):
                    mets.append(me.text)
            mols = {}
            for mo in e.find(mi+'interactorList').findall(mi+'interactor'):
                iid = mo.attrib['id']
                name = mo.find(mi+'names').find(mi+'shortLabel').text
                entrez = ''
                if mo.find(mi+'xref') is not None:
                    entrez = ';'.join(
                        [ac.attrib['id'] for ac in mo.find(mi+'xref')\
                            .findall(mi+'secondaryRef') 
                        if ac.attrib['db'] == 'Entrez gene'])
                mols[iid] = (name,entrez)
            theInt = e.find(mi+'interactionList').find(mi+'interaction')
            for p in theInt.find(mi+'participantList').findall(mi+'participant'):
                pid = p.find(mi+'interactorRef').text
                roles = ''
                if p.find(mi+'experimentalRoleList') is not None:
                    roles = ';'.join(
                        [rl.find(mi+'names').find(mi+'shortLabel').text 
                        for rl in p.find(mi+'experimentalRoleList')\
                            .findall(mi+'experimentalRole')])
                mols[pid] += (roles,)
            intTyp = theInt.find(mi+'interactionType').find(mi+'names')\
                .find(mi+'shortLabel').text
            for i in range(0,len(mols)-1):
                for j in range(i,len(mols)):
                    A = mols[mols.keys()[i]][0:2]
                    B = mols[mols.keys()[j]][0:2]
                    result.append(list(A) + list(B) + \
                        [';'.join(refs), ';'.join(mets), intTyp, pwname])
    return result

def get_pubmeds(pmids):
    pmids = [str(pmid) for pmid in pmids]
    url = data_formats.urls['pubmed-eutils']['url']
    cache = len(pmids) < 10
    data = {}
    prg = progress.Progress(len(pmids)/100 + 1, 'Retrieving data from NCBI e-utils', 
        1, percent = False)
    for offset in xrange(0, len(pmids), 100):
        prg.step()
        post = {'id': ','.join(pmids[offset:offset+100]), 'retmode': 'json', 'db': 'pubmed'}
        hdr = ['X-HTTP-Method-Override:GET']
        for i in xrange(3):
            try:
                res = curl(url, silent = False, cache = cache, post = post, req_headers = hdr)
                data = dict([(k,v) for k,v in json.loads(res)['result'].iteritems()] + \
                    [(k,v) for k,v in data.iteritems()])
                break
            except ValueError:
                sys.stdout.write('\t:: Error in JSON, retry %u\n'%i)
                sys.stdout.flush()
    prg.terminate()
    return data

def get_lincs_compounds():
    sys.stdout.write('\n\tReturned dict has names, brand names or company specific\n'\
        '\tIDs of compounds as keys, and tuples of PubChem, ChEMBL, ChEBI, InChi, \n'\
        '\tInChi Key, SMILES and LINCS as values.\n\n')
    sys.stdout.flush()
    return dict(\
        [(key, pair[1]) for pair in \
            [\
                (\
                    [it for sl in \
                        [filter(lambda z: len(z) > 0, y.split(';')) \
                            for y in x[1:4] \
                                if len(y) > 0
                        ] for it in sl\
                    ], \
                    (x[4], 
                        '' if len(x[7]) == 0 else 'CHEMBL%s'%x[7], 
                        '' if len(x[8]) == 0 else 'CHEBI:%s'%x[8], 
                    x[9], x[10], x[11], x[3])\
                ) \
            for x in \
                [\
                    [b.strip() for b in a.split('\t')] \
                        for a in ''.join(\
                            [\
                                s.replace(',','\t') \
                                    if i%2==0 
                                    else s.replace('\n', '') \
                                for i, s in enumerate(\
                                    curl(\
                                        data_formats.urls['lincs-compounds']['url'], 
                                        silent = False)\
                                    .split('"')\
                                )] \
                        ).split('\n')[1:] \
                            if len(a) > 0]
                ] \
        for key in pair[0]]\
    )

def get_hpmr():
    '''
    Downloads and processes the list of all human receptors from
    human receptor census (HPMR -- Human Plasma Membrane Receptome).
    Returns list of GeneSymbols.
    '''
    html = curl(data_formats.urls['hpmr']['url'], silent = False)
    soup = bs4.BeautifulSoup(html, 'html.parser')
    gnames = [gname for gname in [tr.find_all('td')[1].text \
        for tr in soup.find_all('tr', class_ = 'res_receptor_rec')] \
        if not gname.lower().startswith('similar')]
    return common.uniqList(gnames)

def get_tfcensus(classes = ['a', 'b', 'other']):
    '''
    Downloads and processes list of all human transcripton factors.
    Returns dict with lists of ENSGene IDs and HGNC Gene Names.
    '''
    ensg = []
    hgnc = []
    reensg = re.compile(r'ENSG[0-9]{11}')
    fname = os.path.join(common.ROOT, 'data', 'vaquerizas2009-s3.txt')
    with open(fname, 'r') as f:
        for l in f:
            if len(l) > 0 and l.split('\t')[0] in classes:
                ensg += reensg.findall(l)
                h = l.split('\t')[5].strip()
                if len(h) > 0:
                    hgnc.append(h)
    return {'ensg': ensg, 'hgnc': hgnc}

def get_guide2pharma(organism = 'human', endogenous = True):
    '''
    Downloads and processes Guide to Pharmacology data.
    Returns list of dicts.
    
    @organism : str
        Name of the organism, e.g. `human`.
    @endogenous : bool
        Whether to include only endogenous ligands interactions.
    '''
    positives = ['agonist', 'activator', 'potentiation', 'partial agonist', 
        'inverse antagonist', 'full agonist', 'activation', 
        'irreversible agonist', 'positive']
    negatives = ['inhibitor', 'antagonist', 'inhibition', 'irreversible inhibition',
        'inverse agonist', 'negative', 'weak inhibition', 'reversible inhibition']
    url = data_formats.urls['gtp']['url']
    data = curl(url, silent = False)
    buff = StringIO()
    buff.write(data)
    cols = {
        'receptor_uniprot': 3,
        'ligand_uniprot': 7,
        'receptor_species': 9,
        'ligand_genesymbol': 12,
        'ligand_species': 13,
        'effect0': 15,  # Agonist, Activator
                        # Inhibitor, Antagonist
        'effect1': 16,  # Potentiation, Partial agonist, inverse antagonist, Full agonist,
                        # Activator, Activation, Irreversible agonist, Positive
                        # // Inhibition, irreversible inhibition, Inverse agonist, Negative,
        'effect2': 17,  # Agonist, Activation, Activator, Partial agonist, Positive,
                        # Potentiation, //
                        # Inhibition, weak inhibition, Negative, 
                        # Reversible inhibition, Irreversible inhibition
        'endogenous': 18,
        'pubmed': 33
    }
    data = read_table(cols = cols, fileObject = buff, sep = ',', hdr = 1)
    if organism is not None:
        data = [d for d in data if d['receptor_species'].lower().strip() == organism and \
            d['ligand_species'].lower().strip() == organism]
    if endogenous:
        data = [d for d in data if d['endogenous'].strip() == 't']
    for d in data:
        d['effect'] = 0
        for n in [0, 1, 2]:
            if d['effect%u'%n].lower().strip() in positives:
                d['effect'] = 1
            if d['effect%u'%n].lower().strip() in negatives:
                d['effect'] = -1
    return data

def open_pubmed(pmid):
    '''
    Opens PubMed record in web browser.
    
    @pmid : str or int
        PubMed ID
    '''
    pmid = str(pmid)
    url = data_formats.urls['pubmed']['url'] % pmid
    webbrowser.open(url)

def only_pmids(idList, strict = True):
    '''
    Return elements unchanged which compy to PubMed ID format,
    and attempts to translate the DOIs and PMC IDs using NCBI
    E-utils.
    Returns list containing only PMIDs.
    
    @idList : list, str
        List of IDs or one single ID.
    @strict : bool
        Whether keep in the list those IDs which are not PMIDs,
        neither DOIs or PMC IDs or NIH manuscript IDs.
    '''
    if type(idList) in common.simpleTypes:
        idList = [idList]
    pmids = set([i for i in idList if i.isdigit()])
    pmcids = [i for i in idList if i.startswith('PMC')]
    dois = [i for i in idList if '/' in i]
    manuscids = [i for i in idList if i.startswith('NIHMS')]
    if not strict:
        non_pmids = set(idList) - (set(pmids) | set(dois) | set(pmcids) | set(manuscids))
        pmids = pmids | non_pmids
    if len(pmcids) > 0:
        pmids = pmids | set(pmids_list(pmcids))
    if len(dois) > 0:
        pmids = pmids | set(pmids_list(dois))
    return list(pmids)

def get_pmid(idList):
    '''
    For a list of doi or PMC IDs 
    fetches the corresponding PMIDs.
    '''
    if type(idList) in common.simpleTypes:
        idList = [idList]
    url = data_formats.urls['pubmed-eutils']['conv'] % ','.join(str(i) for i in idList)
    data = curl(url, silent = True)
    try:
        js = json.loads(data)
    except:
        js = {}
    return js

def pmids_dict(idList):
    jsn = get_pmid(idList)
    result = {
        'doi': {},
        'pmc': {}
    }
    if 'records' in jsn:
        for r in jsn['records']:
            if 'pmid' in r:
                if 'doi' in r:
                    result['doi'][r['pmid']] = r['doi']
                if 'pmcid' in r:
                    result['pmc'][r['pmid']] = r['pmcid']
    return result

def pmids_list(idList):
    jsn = get_pmid(idList)
    result = []
    if 'records' in jsn:
        for r in jsn['records']:
            if 'pmid' in r:
                result.append(r['pmid'])
    return result

def get_hprd(in_vivo = True):
    '''
    Downloads and preprocesses HPRD data.
    '''
    url = data_formats.urls['hprd_all']['url']
    files = [data_formats.urls['hprd_all']['ptm_file']]
    data = curl(url, silent = False, files_needed = files)
    if len(data) == 0:
        return []
    data = [l.split('\t') for l in data[files[0]].split('\n')][:-1]
    if in_vivo:
        data = [i for i in data if 'in vivo' in i[9].split(';')]
    return data

def hprd_interactions(in_vivo = True):
    '''
    Processes HPRD data and extracts interactions.
    Returns list of interactions.
    '''
    return [i for i in get_hprd(in_vivo = in_vivo) if i[6] != '-']

def get_hprd_ptms(in_vivo = True):
    '''
    Processes HPRD data and extracts PTMs.
    Returns list of kinase-substrate interactions.
    '''
    ptms = []
    non_digit = re.compile(r'[^\d]+')
    data = get_hprd(in_vivo = in_vivo)
    for ptm in data:
        if ptm[6] != '-':
            resnums = [int(nn) for nn in \
                [non_digit.sub('', n) for n in ptm[4].split(';')]\
                if len(nn) > 0]
            for resnum in resnums:
                ptms.append({
                    'resaa': ptm[5],
                    'resnum': resnum,
                    'typ': ptm[8].lower(),
                    'references': ptm[10].split(','),
                    'kinase': ptm[6],
                    'substrate_refseqp': ptm[3],
                    'substrate': ptm[1],
                    'start': max(resnum -7, 1),
                    'end': resnum + 7,
                    'instance': None
                })
    return ptms

def get_disgenet(dataset = 'curated'):
    '''
    Downloads and processes the list of all human disease related proteins
    from DisGeNet.
    Returns dict of dicts.
    
    @dataset : str
        Name of DisGeNet dataset to be obtained:
        `curated`, `literature`, `befree` or `all`.
    '''
    url = data_formats.urls['disgenet']['url'] % dataset
    data = curl(url, silent = False, files_needed = [url.split('/')[-1].replace('tar.gz', 'txt')])
    cols = {
        'entrez': 0,
        'genesymbol': 1,
        'umls': 3,
        'disease': 4,
        'score': 5,
        'assoc_typ': 7,
        'source': 8
    }
    data = read_table(cols = cols, data = data.values()[0], hdr = 1, sep = '\t')
    for i, d in enumerate(data):
        data[i]['score'] = float(data[i]['score'])
        data[i]['assoc_typ'] = [x.strip() for x in data[i]['assoc_typ'].split(',')]
        data[i]['source'] = [x.strip() for x in data[i]['source'].split(',')]
    return data

def load_lmpid(fname = 'LMPID_DATA_pubmed_ref.xml', organism = 9606):
    '''
    Reads and processes LMPID data from local file
    `pypath.data/LMPID_DATA_pubmed_ref.xml`.
    The file was provided by LMPID authors and is now
    redistributed with the module.
    Returns list of domain-motif interactions.
    '''
    result = []
    with open(os.path.join(common.ROOT, 'data', fname), 'r') as f:
        data = f.read()
    soup = bs4.BeautifulSoup(data, 'html.parser')
    uniprots = all_uniprots(organism = organism, swissprot = None)
    prg = progress.Progress(len(soup.find_all('record')), 'Processing data from LMPID', 21)
    for rec in soup.find_all('record'):
        prg.step()
        uniprot_bait = rec.bait_uniprot_id.text
        uniprot_prey = rec.prey_uniprot_id.text
        if uniprot_bait in uniprots and uniprot_prey in uniprots:
            result.append({
                'bait': uniprot_bait,
                'prey': uniprot_prey,
                'refs': [x.strip() for x in rec.references.text.split(',')],
                'pos': [int(x) for x in rec.sequence_position.text.split('-')],
                'inst': rec.motif_instance.text,
                'dom': rec.interacting_domain.text
            })
    prg.terminate()
    return result

def lmpid_interactions(fname = 'LMPID_DATA_pubmed_ref.xml', organism = 9606):
    '''
    Converts list of domain-motif interactions supplied by
    `pypath.dataio.load_lmpid()` to list of interactions.
    '''
    data = load_lmpid(fname = fname, organism = organism)
    return [[l['prey'], l['bait'], ';'.join(l['refs'])] for l in data]

def lmpid_dmi(fname = 'LMPID_DATA_pubmed_ref.xml', organism = 9606):
    '''
    Converts list of domain-motif interactions supplied by
    `pypath.dataio.load_lmpid()` to list of 
    `pypath.intera.DomainMotif() objects.
    '''
    data = load_lmpid(fname = fname, organism = organism)
    return [{
        'motif_protein': l['bait'],
        'domain_protein': l['prey'],
        'instance': l['inst'],
        'motif_start': l['pos'][0],
        'motif_end': l['pos'][1],
        'domain_name': l['dom'],
        'domain_name_type': 'name',
        'refs': l['refs']
    } for l in data]

def get_hsn():
    '''
    Downloads and processes HumanSignalingNetwork version 6
    (published 2014 Jan by Edwin Wang).
    Returns list of interactions.
    '''
    url = data_formats.urls['hsn']['url']
    data = curl(url, silent = False).split('\n')[1:]
    data = [r.split(',') for r in data if len(r) > 0]
    return data

def get_li2012(filename = None):
    '''
    Reads supplementary data of Li 2012 from local file.
    Returns table (list of lists).
    '''
    filename = filename if filename is not None \
        else os.path.join(common.ROOT, 'data', data_formats.urls['li2012']['file'])
    with open(filename, 'r') as f:
        data = [r.split('\t') for r in f.read().split('\n')[1:] if len(r) > 0]
    return data

def li2012_interactions():
    '''
    Converts table read by `pypath.dataio.get_li2012()` to
    list of interactions.
    '''
    result = []
    data = get_li2012()
    for l in data:
        subs_protein = l[1].split('/')[0]
        tk_protein = l[2].split()[0]
        reader_protein = l[3].split()[0]
        route = l[4]
        result.append((
            tk_protein, subs_protein, route, 'phosphorylation'
        ))
        result.append((
            subs_protein, reader_protein, route, 'phosphomotif_binding'
        ))
    return [list(l) for l in common.uniqList(result)]

def li2012_phospho():
    '''
    Converts table read by `pypath.dataio.get_li2012()` to
    list of dicts of kinase-substrate interactions.
    '''
    result = []
    non_digit = re.compile(r'[^\d]+')
    data = get_li2012()
    for l in data:
        subs_protein = l[1].split('/')[0]
        tk_protein = l[2].split()[0]
        subs_resnum = int(non_digit.sub('', l[1].split('/')[1]))
        result.append((subs_protein, tk_protein, None, 
            None, None, 'Y', subs_resnum))
    result = [dict(zip(['substrate', 'kinase', 'instance', \
        'start', 'end', 'resaa', 'resnum'], \
        list(l))) for l in common.uniqList(result)]
    return result

def li2012_dmi(mapper = None):
    '''
    Converts table read by `pypath.dataio.get_li2012()` to
    list of `pypath.intera.DomainMotif()` objects.
    Translates GeneSymbols to UniProt IDs.
    
    @mapper : pypath.mapping.Mapper()
        If not provided, a new `Mapper()` instance will be 
        initialized, reserving more memory.
    '''
    result = {}
    nondigit = re.compile(r'[^\d]+')
    se = swissprot_seq(isoforms = True)
    if type(mapper) is not mapping.Mapper: mapper = mapping.Mapper(9606)
    data = get_li2012()
    for l in data:
        subs_protein = l[1].split('/')[0]
        tk_protein = l[2].split()[0]
        reader_protein = l[3].split()[0]
        subs_uniprots = mapper.map_name(subs_protein, 'genesymbol', 'uniprot')
        tk_uniprots = mapper.map_name(tk_protein, 'genesymbol', 'uniprot')
        reader_uniprots = mapper.map_name(reader_protein, 'genesymbol', 'uniprot')
        subs_resnum = int(non_digit.sub('', l[1].split('/')[1]))
        for su in subs_uniprots:
            if su in se:
                subs_iso = None
                for iso, s in se[su].isof.iteritems():
                    if se[su].get(subs_resnum, isoform = iso) == 'Y':
                        subs_iso = iso
                        break
                if subs_iso:
                    start = min(1, subs_resnum - 7)
                    end = max(subs_resnum + 7, len(se[su].isof[subs_iso]))
                    for ku in tk_uniprots:
                        res = intera.Residue(subs_resnum, 'Y', su, isoform = subs_iso)
                        mot = intera.Motif(su, start, end, isoform = subs_iso, 
                            instance = se[su].get(start, end, isoform = subs_iso))
                        ptm = intera.Ptm(su, motif = mot, residue = res, 
                            isoform = subs_iso, source = 'Li2012')
                        dom = intera.Domain(ku)
                        dommot = intera.DomainMotif(domain = dom, ptm = ptm, 
                            sources = ['Li2012'])
                        result = {}
    return result

def take_a_trip(cachefile = 'trip.pickle'):
    '''
    Downloads TRIP data from webpage and preprocesses it.
    Saves preprocessed data into `cachefile` and next
    time loads from this file.
    
    @cachefile : str
        Filename, located in `./cache`.
        To disable cache, pass `None`.
        To download again, remove file from `./cache`.
    '''
    if cachefile is not None:
        cachefile = os.path.join('cache', 'trip.pickle')
        if os.path.exists(cachefile):
            result = pickle.load(open(cachefile, 'rb'))
            return result
    result = {
        'sc': {},
        'cc': {},
        'vvc': {},
        'vtc': {},
        'fc': {}
    }
    intrs = {}
    titles = {
        'Characterization': 'cc',
        'Screening': 'sc',
        'Validation: In vitro validation': 'vtc',
        'Validation: In vivo validation': 'vvc',
        'Functional consequence': 'fc'
    }
    interactors = {}
    base_url = data_formats.urls['trip']['base']
    show_url = data_formats.urls['trip']['show']
    # url = data_formats.urls['trip']['url']
    # json_url = data_formats.urls['trip']['json']
    # jsn = curl(json_url, silent = False)
    # data = curl(url, silent = False)
    # jsn = json.loads(jsn, encoding = 'utf-8')
    mainhtml = curl(base_url)
    mainsoup = bs4.BeautifulSoup(mainhtml, 'html.parser')
    trppages = common.flatList([[a.attrs['href'] for a in ul.find_all('a')] \
        for ul in mainsoup.find('div', id = 'trp_selector').find('ul').find_all('ul')])
    for trpp in trppages:
        trp = trpp.split('/')[-1]
        trpurl = show_url % trp
        trphtml = curl(trpurl, silent = False)
        trpsoup = bs4.BeautifulSoup(trphtml, 'html.parser')
        trp_uniprot = trip_find_uniprot(trpsoup)
        if trp_uniprot is None or len(trp_uniprot) < 6:
            sys.stdout.write('\t\tcould not find uniprot for %s\n' % trp)
            sys.stdout.flush()
        for tab in trpsoup.find_all('th', colspan = ['11', '13']):
            ttl = titles[tab.text.strip()]
            tab = tab.find_parent('table')
            trip_process_table(tab, result[ttl], intrs, trp_uniprot)
    pickle.dump(result, open(cachefile, 'wb'))
    return result

def trip_process_table(tab, result, intrs, trp_uniprot):
    '''
    Processes one HTML table downloaded from TRIP webpage.
    
    @tab : bs4.element.Tag()
        One table of interactions from TRIP webpage.
    @result : dict
        Dictionary the data should be filled in.
    @intrs : dict
        Dictionary of already converted interactor IDs.
        This serves as a cache so do not need to look up
        the same ID twice.
    @trp_uniprot : str
        UniProt ID of TRP domain containing protein.
    '''
    for row in tab.find_all('tr'):
        cells = row.find_all(['td', 'th'])
        if 'th' not in [c.name for c in cells]:
            intr = cells[2].text.strip()
            if intr not in intrs:
                intr_uniprot = trip_get_uniprot(intr)
                intrs[intr] = intr_uniprot
                if intr_uniprot is None or len(intr_uniprot) < 6:
                    print '\t\tcould not find uniprot for %s' % intr
            else:
                intr_uniprot = intrs[intr]
            if (trp_uniprot, intr_uniprot) not in result:
                result[(trp_uniprot, intr_uniprot)] = []
            result[(trp_uniprot, intr_uniprot)].append([c.text.strip() for c in cells])

def trip_get_uniprot(syn):
    '''
    Downloads table from TRIP webpage and UniProt attempts to
    look up the UniProt ID for one synonym.
    
    @syn : str
        The synonym as shown on TRIP webpage.
    '''
    url = data_formats.urls['trip']['show'] % syn
    html = curl(url)
    soup = bs4.BeautifulSoup(html, 'html.parser')
    return trip_find_uniprot(soup)

def trip_find_uniprot(soup):
    '''
    Looks up a UniProt name in table downloaded from TRIP
    webpage.
    
    @soup : bs4.BeautifulSoup
        The `BeautifulSoup` instance returned by `pypath.dataio.trip_get_uniprot()`.
    '''
    for tr in soup.find_all('div', id='tab2')[0].find_all('tr'):
        if tr.find('td') is not None and tr.find('td').text.strip() == 'Human':
            uniprot = tr.find_all('td')[2].text.strip()
            # print '\t\tuniprot found: %s' % str(uniprot)
            return uniprot
    return None

def trip_process(exclude_methods = ['Inference', 'Speculation'], 
    predictions = False, species = 'Human', strict = False):
    '''
    Downloads TRIP data by calling `pypath.dadio.take_a_trip()` and
    further provcesses it.
    Returns dict of dict with TRIP data.
    
    @exclude_methods : list
        Interaction detection methods to be discarded.
    @predictions : bool
        Whether to include predicted interactions.
    @species : str
        Organism name, e.g. `Human`.
    @strict : bool
        Whether include interactions with species not 
        used as a bait or not specified.
    '''
    nd = 'Not determined'
    spec = set([]) if strict \
        else set(['Not specified', 'Not used as a bait', ''])
    spec.add(species)
    result = {}
    data = take_a_trip()
    for uniprots in common.uniqList(common.flatList( \
        [v.keys() for v in data.values()])):
        to_process = False
        refs = set([])
        mets = set([])
        tiss = set([])
        reg = set([])
        eff = set([])
        if uniprots in data['sc']:
            for sc in data['sc'][uniprots]:
                if sc[4] in spec and sc[6] in spec and \
                    (predictions or sc[9] != 'Prediction') and \
                    sc[3] not in exclude_methods:
                    refs.add(sc[10])
                    mets.add(sc[3])
                    tiss.add(sc[7])
        if uniprots in data['vtc']:
            for vtc in data['vtc'][uniprots]:
                if vtc[4] in spec and vtc[7] in spec and \
                    vtc[3] not in exclude_methods:
                    refs.add(vtc[10])
                    mets.add(vtc[3])
        if uniprots in data['vvc']:
            for vvc in data['vvc'][uniprots]:
                if vvc[6] in spec and vvc[8] in spec and \
                    vvc[3] not in exclude_methods:
                    refs.add(vvc[10])
                    mets.add(vvc[3])
                    if len(vvc[4]) > 0:
                        tiss.add(vvc[4])
                    if len(vvc[5]) > 0:
                        tiss.add(vvc[5])
        if uniprots in data['cc']:
            for cc in data['cc'][uniprots]:
                if cc[4] in spec and cc[6] in spec and \
                    cc[3] not in exclude_methods:
                    refs.add(cc[10])
                    mets.add(cc[3])
                    if (cc[5] != nd and len(cc[5]) > 0) or \
                        (cc[7] != nd and len(cc[7]) > 0):
                        reg.add((cc[5], cc[7]))
        if uniprots in data['fc']:
            for fc in data['fc'][uniprots]:
                mets.add(fc[3])
                refs.add(fc[7])
                if len(fc[5]) > 0:
                    eff.add(fc[5])
                if len(fc[6]) > 0:
                    eff.add(fc[6])
        if len(refs) > 0:
            result[uniprots] = {
                'refs': refs,
                'methods': mets,
                'tissues': tiss,
                'effect': eff,
                'regions': reg
            }
    return result

def trip_interactions(exclude_methods = ['Inference', 'Speculation'], 
    predictions = False, species = 'Human', strict = False):
    '''
    Obtains processed TRIP interactions by calling `pypath.dataio.trip_process()`
    and returns list of interactions. All arguments are passed to 
    `trip_process()`, see their definition there.
    '''
    data = trip_process(exclude_methods, predictions, species, strict)
    def trip_effect(eff):
        pos = set(['Sensitization', 'Activation', 'Increase in plasma membrane level', 
            'Increase in lysosomal membrane level', 'New channel creation'])
        neg = set(['Desensitization', 'Decrease in plasma membrane level', 'Inhibition', 
            'Internalization from membrane by ligand', 'Retain in the endoplasmic reticulum'])
        return 'stimulation' if len(eff & pos) > 0 \
            else 'inhibition' if len(eff & neg) > 0 else 'unknown'
    return [[unipr[0], unipr[1], ';'.join(d['refs']), 
        ';'.join(d['methods']), trip_effect(d['effect'])] for unipr, d in data.iteritems()]

def load_signor_ptms(fname = 'signor_22052015.tab'):
    '''
    This function is deprecated, you should not use it.
    Loads and processes Signor PTMs from local file.
    Returns dict of dicts.
    '''
    reres = re.compile(r'([A-Za-z]{3})([0-9]+)')
    result = []
    aalet = dict((k.lower().capitalize(), v) for k, v in common.aaletters.iteritems())
    fname = os.path.join(common.ROOT, 'data', fname)
    if os.path.exists(fname):
        with open(fname, 'r') as f:
            data = [[i.strip() for i in l.split('\t')] for l in \
                [ll.strip() for ll in f.read().split('\n')[1:] \
                if len(ll.strip()) > 0]]
    for d in data:
        resm = reres.match(d[10])
        if resm is not None:
            aa = aalet[resm.groups()[0]]
            aanum = int(resm.groups()[1])
            typ = d[9]
            inst = d[11].upper()
            result.append({
                'typ': d[9],
                'resnum': aanum, 
                'instance': inst, 
                'substrate': d[6], 
                'start': aanum - 7, 
                'end': aanum + 7, 
                'kinase': d[2], 
                'resaa': aa, 
                'motif': inst
            })
    return result

def load_macrophage():
    '''
    Loads Macrophage from local file.
    Returns list of interactions.
    '''
    fname = data_formats.omnipath['macrophage'].inFile
    fname = os.path.join(common.ROOT, 'data', fname)
    with open(fname, 'r') as f:
        data = f.read()
    data = data.replace('?', '').replace('->', ',')

def kegg_pathways(mapper = None):
    '''
    Downloads and processes KEGG Pathways.
    Returns list of interactions.
    '''
    rehsa = re.compile(r'.*(hsa[0-9]+).*')
    mapper = mapper if mapper is not None else mapping.Mapper()
    hsa_list = []
    interactions = []
    htmllst = curl(data_formats.urls['kegg_pws']['list_url'], silent = True)
    lstsoup = bs4.BeautifulSoup(htmllst, 'html.parser')
    for a in lstsoup.find_all('a', href = True):
        m = rehsa.match(a['href'])
        if m:
            hsa_list.append(m.groups(0)[0])
    prg = progress.Progress(len(hsa_list), 'Processing KEGG Pathways', 1, percent = False)
    for hsa in hsa_list:
        prg.step()
        kgml = curl(data_formats.urls['kegg_pws']['kgml_url'] % hsa, silent = True)
        kgmlsoup = bs4.BeautifulSoup(kgml, 'html.parser')
        entries = {}
        for ent in kgmlsoup.find_all('entry'):
            gr = ent.find('graphics')
            if gr and 'name' in gr.attrs:
                entries[ent.attrs['id']] = [n.strip() for n in gr.attrs['name']\
                    .replace('...', '').split(',')]
        uentries = dict([(eid, common.uniqList(common.flatList(
            [mapper.map_name(gn, 'genesymbol', 'uniprot', strict = True) \
                for gn in gns]))) for eid, gns in entries.iteritems()])
        for rel in kgmlsoup.find_all('relation'):
            st = rel.find('subtype')
            if rel.attrs['entry1'] in uentries and rel.attrs['entry2'] in uentries and \
                st and 'name' in st.attrs:
                for u1 in uentries[rel.attrs['entry1']]:
                    for u2 in uentries[rel.attrs['entry2']]:
                        interactions.append((u1, u2, st.attrs['name']))
    prg.terminate()
    return common.uniqList(interactions)

def signor_urls():
    '''
    This function is deprecated.
    '''
    tsv_urls = []
    url = data_formats.urls['signor']['list_url']
    baseurl = data_formats.urls['signor']['base_url']
    html = curl(url, silent = True)
    soup = bs4.BeautifulSoup(html, 'html.parser')
    for td in soup.find_all('td', style = lambda x: x.startswith('border')):
        pw = td.text.strip().split(':')[0]
        tsv_url = baseurl % td.find('a').attrs['href']
        tsv_urls.append((pw, tsv_url))
    return tsv_urls

def signor_interactions():
    '''
    Downloads the full dataset from Signor.
    Return the file contents.
    '''
    url = data_formats.urls['signor']['all_url']
    return curl(url, silent = False, large = True)

def rolland_hi_ii_14(filename = None):
    '''
    Loads the HI-II-14 unbiased interactome from the large scale screening
    of from Rolland 2014.
    Returns list of interactions.
    '''
    filename = filename if filename is not None \
        else data_formats.urls['hiii14']['file']
    with open(filename, 'r') as f:
        null = f.readline()
        data = [l.strip().split('\t') for l in f]
    return data

def read_xls(xls_file, sheet = '', csv_file = None, return_table = True):
    '''
    Generic function to read MS Excel XLS file, and convert one sheet
    to CSV, or return as a list of lists
    '''
    try:
        book = xlrd.open_workbook(xls_file, on_demand = True)
        try:
            sheet = book.sheet_by_name(sheet)
        except XLRDError:
            sheet = book.sheet_by_index(0)
        table = [[str(c.value) for c in sheet.row(i)] for i in xrange(sheet.nrows)]
        if csv_file:
            with open(csv_file, 'w') as csv:
                csv.write('\n'.join(['\t'.join(r) for r in table]))
        if not return_table:
            table = None
        book.release_resources()
        return table
    except IOError:
        sys.stdout.write('No such file: %s\n' % xls_file)
    sys.stdout.flush()

def get_kinases():
    '''
    Downloads and processes the list of all human kinases.
    Returns a list of GeneSymbols.
    '''
    url = data_formats.urls['kinome']['url']
    xlsf = curl(url, large = True, silent = False)
    xlsname = xlsf.name
    xlsf.close()
    tbl = read_xls(xlsname)
    genesymbols = [l[23] for l in tbl[1:] if len(l[23]) > 0]
    return genesymbols

def get_dgidb():
    '''
    Downloads and processes the list of all human druggable proteins.
    Returns a list of GeneSymbols.
    '''
    genesymbols = []
    url = data_formats.urls['dgidb']['main_url']
    html = curl(url, silent = False)
    soup = bs4.BeautifulSoup(html, 'html.parser')
    cats = [o.attrs['value'] \
        for o in soup.find('select', {'id': 'gene_categories'})\
            .find_all('option')]
    for cat in cats:
        url = data_formats.urls['dgidb']['url'] % cat
        html = curl(url)
        soup = bs4.BeautifulSoup(html, 'html.parser')
        trs = soup.find('tbody').find_all('tr')
        genesymbols.extend([tr.find('td').text.strip() for tr in trs])
    return common.uniqList(genesymbols)

def reactome_sbml():
    '''
    Downloads Reactome human reactions in SBML format.
    Returns gzip.GzipFile object.
    '''
    url = data_formats.urls['reactome']['sbml']
    sbml = curl(url, silent = False, large = True)
    return sbml

def reactome_biopax(organism = 9606, cache = True):
    '''
    Downloads Reactome human reactions in SBML format.
    Returns File object.
    '''
    organisms = {
        9606: 'Homo_sapiens'
    }
    unzipped = 'cache/reactome_biopax_%s.owl' % organisms[organism]
    if not os.path.exists(unzipped) or not cache:
        url = data_formats.urls['reactome']['biopax_l3']
        bpz = curl(url, silent = False, large = True, 
            files_needed = ['%s.owl'%organisms[organism]]).values()[0]
        with open(unzipped, 'w') as _unzipped:
            while True:
                chunk = bpz.read(4096)
                if not chunk:
                    break
                _unzipped.write(chunk)
        bpz.close()
    _unzipped = open(unzipped, 'r')
    return _unzipped

def pid_biopax():
    url = data_formats.urls['nci-pid']['biopax_l3']
    return curl(url, silent = False, large = True)

def panther_biopax():
    url = data_formats.urls['panther']['biopax_l3']
    return curl(url, silent = False, large = True).values()

def acsn_biopax():
    url = data_formats.urls['acsn']['biopax_l3']
    return curl(url, silent = False, large = True)

def reactome_bs():
    sbml = reactome_sbml()
    soup = bs4.BeautifulSoup(sbml.read(), 'html.parser')
    return soup

# Process Reactome BioPAX level 3

def get_soup(elem):
    return bs4.BeautifulSoup(etree.tostring(elem), 'html.parser')

def _bp_collect_resources(elem, tag, restype = None):
    rdfpref = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    rdfres = '%sresource' % rdfpref
    return [x.get(rdfres).replace('#', '') for x in elem.iterfind(tag) \
        if rdfres in x.attrib and \
            (restype is None or x.get(rdfres).replace('#', '').startswith(restype))]

def reactions_biopax(biopax_file, organism = 9606, protein_name_type = 'UniProt', clean = True):
    '''
    Processes a BioPAX file and extracts binary interactions.
    '''
    cachefile = os.path.join('cache', '%s.processed.pickle' % os.path.split(biopax_file.name)[1])
    if os.path.exists(cachefile):
        sys.stdout.write('\t:: Loading already processed data\n')
        sys.stdout.flush()
        return pickle.load(open(cachefile, 'r'))
    # string constants
    bppref = '{http://www.biopax.org/release/biopax-level3.owl#}'
    rdfpref = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    rdfid = '%sID' % rdfpref
    rdfab = '%sabout' % rdfpref
    rdfres = '%sresource' % rdfpref
    bpprot = '%sProtein' % bppref
    bpcplx = '%sComplex' % bppref
    bpprre = '%sProteinReference' % bppref
    bpreac = '%sBiochemicalReaction' % bppref
    bpcata = '%sCatalysis' % bppref
    bpctrl = '%sControl' % bppref
    bpcoma = '%sComplexAssembly' % bppref
    bppstp = '%sPathwayStep' % bppref
    bpuxrf = '%sUnificationXref' % bppref
    bpstoi = '%sStoichiometry' % bppref
    bppubr = '%sPublicationXref' % bppref
    bppath = '%sPathway' % bppref
    bpfrfe = '%sFragmentFeature' % bppref
    bpseqi = '%sSequenceInterval' % bppref
    bpseqs = '%sSequenceSite' % bppref
    bpmodf = '%sModificationFeature' % bppref
    bpmodv = '%sSequenceModificationVocabulary' % bppref
    bpmphe = '%smemberPhysicalEntity' % bppref
    bperef = '%sentityReference' % bppref
    bpxref = '%sxref' % bppref
    bprelr = '%sRelationshipXref' % bppref
    bpcsto = '%scomponentStoichiometry' % bppref
    bpstoc = '%sstoichiometricCoefficient' % bppref
    bpphye = '%sphysicalEntity' % bppref
    bpcted = '%scontrolled' % bppref
    bpcter = '%scontroller' % bppref
    bpctyp = '%scontrolType' % bppref
    bpleft = '%sleft' % bppref
    bprgth = '%sright' % bppref
    bpsprc = '%sstepProcess' % bppref
    bpfeat = '%sfeature' % bppref
    bpfelo = '%sfeatureLocation' % bppref
    bpibeg = '%ssequenceIntervalBegin' % bppref
    bpiend = '%ssequenceIntervalEnd' % bppref
    bpseqp = '%ssequencePosition' % bppref
    bpmoty = '%smodificationType' % bppref
    bppcom = '%spathwayComponent' % bppref
    bpterm = '%sterm' % bppref
    bpdb = '%sdb' % bppref
    bpid = '%sid' % bppref
    upStr = 'UniProt'
    modvoc = data_formats.reactome_modifications
    # intermediate results
    proteins = {}
    proteinfamilies = {}
    uniprots = {}
    proteinreferences = {}
    complexes = {}
    complexvariations = {}
    stoichiometries = {}
    reactions = {}
    complexassemblies = {}
    catalyses = {}
    controls = {}
    pathways = {}
    pathwaysteps = {}
    publications = {}
    fragmentfeatures = {}
    sequenceintervals = {}
    sequencesites = {}
    modificationfeatures = {}
    modificationvocabulary = {}
    protein_name_type = protein_name_type.lower()
    # processing the XML
    bpf = reactome_biopax(organism = organism)
    bp_filesize = 0
    if type(biopax_file) is file:
        bp_filesize = os.path.getsize(biopax_file.name)
    elif type(biopax_file) is tarfile.ExFileObject:
        bp_filesize = biopax_file.size
    elif type(biopax_file) is gzip.GzipFile:
        f = open(biopax_file.name, 'rb')
        f.seek(-4, 2)
        bp_filesize = struct.unpack('<I', f.read())[0]
        f.close()
    prg = progress.Progress(bp_filesize, 'Processing BioPAX XML', 1)
    fpos = biopax_file.tell()
    bp = etree.iterparse(biopax_file, events = ('end',))
    used_elements = []
    try:
        for ev, elem in bp:
            new_fpos = biopax_file.tell()
            prg.step(new_fpos - fpos)
            fpos = new_fpos
            _id = elem.get(rdfid) if rdfid in elem.attrib else elem.get(rdfab)
            # Protein
            if elem.tag == bpprot:
                entref = elem.find(bperef)
                if entref is not None:
                    proteins[_id] = {
                        'protein': entref.get(rdfres).replace('#', ''),
                        'seqfeatures': _bp_collect_resources(elem, bpfeat),
                        'modfeatures': _bp_collect_resources(elem, bpfeat)
                    }
                else:
                    proteinfamilies[_id] = _bp_collect_resources(elem, bpmphe)
            # ProteinReference
            elif elem.tag == bpprre:
                proteinreferences[_id] = _bp_collect_resources(elem, bpxref)
            # UnificationXref
            elif elem.tag == bpuxrf or elem.tag == bprelr:
                db = elem.find(bpdb)
                if db is not None:
                    if elem.find(bpdb).text.lower().startswith(protein_name_type):
                        i = elem.find(bpid)
                        if i is not None:
                            uniprots[_id] = i.text
            # Complex
            elif elem.tag == bpcplx:
                if elem.find(bpcsto) is not None:
                    complexes[_id] = _bp_collect_resources(elem, bpcsto)
                else:
                    complexvariations[_id] = _bp_collect_resources(elem, bpmphe)
            # Stoichiometry
            elif elem.tag == bpstoi:
                stoichiometries[_id] = (elem.find(bpphye).get(rdfres).replace('#', ''), int(float(elem.find(bpstoc).text)))
            # BiochemicalReaction
            elif elem.tag == bpreac:
                reactions[_id] = {
                    'refs': _bp_collect_resources(elem, bpxref),
                    'left': _bp_collect_resources(elem, bpleft),
                    'right': _bp_collect_resources(elem, bprgth)
                }
            # ComplexAssembly
            elif elem.tag == bpcoma:
                complexassemblies[_id] = {
                    'refs': _bp_collect_resources(elem, bpxref),
                    'left': _bp_collect_resources(elem, bpleft),
                    'right': _bp_collect_resources(elem, bprgth)
                }
            # Catalysis
            elif elem.tag == bpcata:
                cter = elem.find(bpcter)
                cted = elem.find(bpcted)
                if cter is not None and cted is not None:
                    typ = elem.find(bpctyp)
                    catalyses[_id] = {
                        'controller': cter.get(rdfres).replace('#', ''), 
                        'controlled': cted.get(rdfres).replace('#', ''), 
                        'type': '' if typ is None else typ.text}
            # Control
            elif elem.tag == bpctrl:
                cter = elem.find(bpcter)
                cted = elem.find(bpcted)
                if cter is not None and cted is not None:
                    typ = elem.find(bpctyp)
                    controls[_id] = {
                        'refs': _bp_collect_resources(elem, bpxref),
                        'type': typ.text if typ is not None else '',
                        'controller': cter.get(rdfres).replace('#', ''), 
                        'controlled': cted.get(rdfres).replace('#', '')
                    }
            # PathwayStep
            elif elem.tag == bppstp:
                pathwaysteps[_id] = _bp_collect_resources(elem, bppstp)
            # PublicationXref
            elif elem.tag == bppubr:
                pmid = elem.find(bpid)
                if pmid is not None:
                    publications[_id] = pmid.text
            # FragmentFeature
            elif elem.tag == bpfrfe:
                fragmentfeatures[_id] = elem.find(bpfelo).get(rdfres).replace('#', '')
            # SequenceInterval
            elif elem.tag == bpseqi:
                beg = elem.find(bpibeg)
                end = elem.find(bpiend)
                sequenceintervals[_id] = (beg.get(rdfres).replace('#', '') if beg is not None else None,
                    elem.find(bpiend).get(rdfres).replace('#', '') if end is not None else None)
            # SequenceSite
            elif elem.tag == bpseqs:
                seqp = elem.find(bpseqp)
                if seqp is not None:
                    sequencesites[_id] = int(seqp.text)
            # ModificationFeature
            elif elem.tag == bpmodf:
                felo = elem.find(bpfelo)
                moty = elem.find(bpmoty)
                if felo is not None and moty is not None:
                    modificationfeatures[_id] = (elem.find(bpfelo).get(rdfres).replace('#', ''), 
                        elem.find(bpmoty).get(rdfres).replace('#', ''))
            # SequenceModificationVocabulary
            elif elem.tag == bpmodv:
                term = elem.find(bpterm)
                if term is not None:
                    modificationvocabulary[_id] = term.text
            # Pathway
            elif elem.tag == bppath:
                try:
                    pathways[_id] = {
                        'reactions': _bp_collect_resources(elem, bppcom),
                        'pathways': _bp_collect_resources(elem, bppcom)
                    }
                except TypeError:
                    print etree.tostring(elem)
            if clean:
                used_elements.append(elem)
                if len(used_elements) > 800:
                    for e in used_elements[:400]:
                        e.clear()
                    used_elements = used_elements[400:]
    except etree.XMLSyntaxError as e:
        sys.stdout.write('\n\tWARNING: XML processing error: %s\n' % str(e))
        sys.stdout.flush()
    prg.terminate()
    del bp
    biopax_file.close()
    # # # # # # # # # # # # # # # # # #
    # from intermediate to final results
    prg = progress.Progress(len(proteins), 'Processing proteins', 11)
    proteins_uniprots = {}
    # return proteinreferences, uniprots
    for pref, protein in proteins.iteritems():
        prg.step()
        if protein['protein'] in proteinreferences:
            for prref in proteinreferences[protein['protein']]:
                if prref in uniprots:
                    proteins_uniprots[pref] = uniprots[prref]
    prg.terminate()
    prg = progress.Progress(len(proteins), 'Processing PTMs', 11)
    proteins_modifications = {}
    for pref, protein in proteins.iteritems():
        prg.step()
        for modf in protein['modfeatures']:
            if modf in modificationfeatures:
                if modificationfeatures[modf][0] in sequencesites:
                    if modificationfeatures[modf][1] in modificationvocabulary:
                        if modificationvocabulary[modificationfeatures[modf][1]] in modvoc:
                            if pref not in proteins_modifications:
                                proteins_modifications[pref] = set([])
                            proteins_modifications[pref].add(
                                (sequencesites[modificationfeatures[modf][0]], 
                                modvoc[modificationvocabulary[modificationfeatures[modf][1]]][1],
                                modvoc[modificationvocabulary[modificationfeatures[modf][1]]][0])
                            )
    prg.terminate()
    # build a uniform dict to handle all protein based entities
    # including complexes and variations/families
    entity_uniprot = {}
    prg = progress.Progress(len(proteins_uniprots), 'Processing proteins', 11)
    for pref, protein in proteins_uniprots.iteritems():
        prg.step()
        entity_uniprot[pref] = [{
            'members': [protein],
            'ptms': {} if protein not in proteins_modifications \
                else {protein: proteins_modifications[pref]}
        }]
    prg.terminate()
    prg = progress.Progress(len(proteinfamilies), 'Processing protein families', 11)
    for pfref, prefs in proteinfamilies.iteritems():
        prg.step()
        entity_uniprot[pfref] = []
        for pref in prefs:
            if pref in proteins_uniprots:
                entity_uniprot[pfref].append({
                    'members': [proteins_uniprots[pref]],
                    'ptms': {} if pref not in proteins_modifications \
                        else {proteins_uniprots[pref]: proteins_modifications[pref]}
                })
    prg.terminate()
    # return entity_uniprot, complexes, proteins, proteinreferences, uniprots, proteinfamilies, proteins_uniprots, reactions, controls, catalyses, complexassemblies
    del proteins
    del proteinfamilies
    del proteinreferences
    prg = progress.Progress(len(complexes), 'Processing complexes', 11)
    for cref, cplex in complexes.iteritems():
        prg.step()
        if cref not in entity_uniprot:
            process_complex(0, cref, entity_uniprot, 
                complexes, complexvariations, cplex, stoichiometries)
    prg.terminate()
    del complexes
    del stoichiometries
    del proteins_uniprots
    #return entity_uniprot, proteins, proteinreferences, uniprots, complexes, stoichiometries
    # # #
    prg = progress.Progress(len(reactions) + len(complexassemblies), 'Processing reactions', 11)
    reactions_uniprots = \
        process_reactions(reactions, entity_uniprot, publications)
    complexassemblies_uniprots = \
        process_reactions(complexassemblies, entity_uniprot, publications)
    del reactions
    del complexassemblies
    # # #
    prg = progress.Progress(len(controls) + len(catalyses), 'Processing controls and catalyses', 11)
    controls_uniprots = _process_controls(dict(controls.items() + catalyses.items()), entity_uniprot, 
        dict(reactions_uniprots.items() + complexassemblies_uniprots.items()), publications)
    for caref, ca in complexassemblies_uniprots.iteritems():
        controls_uniprots[caref] = {
            'type': 'BINDING',
            'refs': [publications[r] for r in ca['refs'] if r in publications],
            'controller': None,
            'controlled': ca
        }
    del entity_uniprot
    pickle.dump(controls_uniprots, open(cachefile, 'w'))
    # return controls_uniprots, entity_uniprot, proteins, proteinreferences, uniprots, complexes, stoichiometries
    return controls_uniprots

def process_reactions(reactions, entity_uniprot, publications):
    result = {}
    for rref, rea in reactions.iteritems():
        result[rref] = {
            'refs': [publications[r] for r in rea['refs'] if r in publications],
            'left': [entity_uniprot[l] for l in rea['left'] \
                if l in entity_uniprot],
            'right': [entity_uniprot[r] for r in rea['right'] \
                if r in entity_uniprot]
        }
    return result

def _process_controls(controls, entity_uniprot, reactions_uniprots, publications):
    result = {}
    for cref, ctrl in controls.iteritems():
        result[cref] = {
            'type': ctrl['type'],
            'refs': [publications[r] for r in ctrl['refs'] if r in publications] \
                if 'refs' in ctrl else [],
            'controller': entity_uniprot[ctrl['controller']] \
                if ctrl['controller'] in entity_uniprot else None,
            'controlled': reactions_uniprots[ctrl['controlled']] \
                if ctrl['controlled'] in reactions_uniprots else None
        }
    return result

def process_complex(depth, cref, entity_uniprot, 
    complexes, complexvariations, cplex, stoichiometries):
    log = open('reactome.log', 'a')
    tabs = '\t' * (depth + 1)
    log.write('%sStarting processing %s, depth = %u\n' % (tabs[1:], cref, depth))
    this_cplex = [{'members': [], 'ptms': {}}]
    log.write('%sComplex %s have %u member entities\n' % (tabs, cref, len(cplex)))
    for stoi in cplex:
        if stoi in stoichiometries:
            ref, num = stoichiometries[stoi]
            log.write('%sNew member entity: %s, stoichiometric coeff: %u\n' % (tabs, ref, num))
            if ref.startswith('Complex') \
                and ref not in entity_uniprot:
                if ref in complexes:
                    log.write('%s%s is a complex with %u subentities, and hasn\'t been processed yet\n' % \
                        (tabs, ref, len(complexes[ref])))
                    process_complex(depth + 1, ref, entity_uniprot,
                        complexes, complexvariations, complexes[ref], stoichiometries)
                if ref in complexvariations:
                    log.write('%s%s is a complex group with %u variations, and hasn\'t been processed yet\n' % \
                        (tabs, ref, len(complexvariations[ref])))
                    entity_uniprot[ref] = []
                    for mref in complexvariations[ref]:
                        if mref not in entity_uniprot and mref in complexes:
                            log.write('%s%s is a complex with %u subentities, and hasn\'t been processed yet\n' % \
                                (tabs, mref, len(complexes[mref])))
                            process_complex(depth + 1, mref, entity_uniprot, 
                                complexes, complexvariations, complexes[mref], stoichiometries)
                        if mref in entity_uniprot:
                            log.write('%s%s is now processed, adding it as an instance of %s\n' % \
                                (tabs, mref, ref))
                            entity_uniprot[ref].extend(
                                entity_uniprot[mref]
                            )
            if ref in entity_uniprot:
                log.write('%s%s is an already processed entity, with %u variants and %u members\n' % \
                    (tabs, ref, len(entity_uniprot[ref]), 
                    len(entity_uniprot[ref][0]['members']) if len(entity_uniprot[ref]) > 0 else 0))
                log.write('%sNumber of variants after processing %s: %u x %u = %u\n' % \
                    (tabs, ref, len(this_cplex), len(entity_uniprot[ref]), 
                    len(this_cplex) * len(entity_uniprot[ref])))
                this_cplex_new = []
                for var in this_cplex:
                    i = 0
                    for new_member in entity_uniprot[ref]:
                        var_new = copy.deepcopy(var)
                        var_new['members'].extend(new_member['members'] * num)
                        for u, ptm in new_member['ptms'].iteritems():
                            if u not in var_new['ptms']:
                                var_new['ptms'][u] = set([])
                            var_new['ptms'][u] = var_new['ptms'][u] | new_member['ptms'][u]
                        this_cplex_new.append(var_new)
                        i += 1
                this_cplex = this_cplex_new
                log.write('%sNumber of variants after processing %s: %u\n' % \
                    (tabs, ref, len(this_cplex)))
                log.write('%sNumber of members in %s: %u\n' % (tabs, cref, len(this_cplex[0]['members']) if len(this_cplex) > 0 else 0))
            else:
                log.write('%sPermanently missing: %s\n' % (tabs, ref))
    log.write('%sFinished processing %s, found %u variants with %u members\n' % \
        (tabs[1:], cref, len(this_cplex), len(this_cplex[0]['members']) if len(this_cplex) > 0 else 0))
    if cref not in entity_uniprot:
        entity_uniprot[cref] = []
    entity_uniprot[cref].extend(this_cplex)

def reactome_interactions(cacheFile = None, **kwargs):
    '''
    Downloads and processes Reactome BioPAX.
    Extracts binary interactions.
    The applied criteria are very stringent, yields very few interactions.
    Requires large free memory, approx. 2G.
    '''
    cacheFile = os.path.join('cache', 'reactome.interactions.pickle') \
        if cacheFile is None else cacheFile
    if os.path.exists(cacheFile):
        interactions = pickle.load(open(cacheFile, 'rb'))
    else:
        while True:
            sys.stdout.write('\nProcessing Reactome requires huge memory.\n'\
                'Please hit `y` if you have at least 2G free memory,\n'\
                'or `n` to omit Reactome.\n'\
                'After processing once, it will be saved in \n'\
                '%s, so next time can be loaded quickly.\n\n'\
                'Process Reactome now? [y/n]\n' % cacheFile
                    )
            sys.stdout.flush()
            answer = raw_input().lower()
            if answer == 'y':
                return get_interactions('reactome', **kwargs)
            else:
                return []

def acsn_interactions(**kwargs):
    return get_interactions('acsn', **kwargs)

def pid_interactions(**kwargs):
    return get_interactions('pid', **kwargs)

def panther_interactions(**kwargs):
    return get_interactions('panther', **kwargs)

def get_interactions(source, mandatory_refs = True):
    ctrls = get_controls(source)
    return process_controls(ctrls, mandatory_refs)[0]

def get_controls(source, protein_name_type = None):
    name_types = {
        'acsn': 'HGNC',
        'reactome': 'UniProt',
        'pid': 'UniProt',
        'panther': 'UniProt'
    }
    if protein_name_type is None and source in name_types:
        protein_name_type = name_types[source]
    biopax = globals()['%s_biopax' % source]
    bpfile = biopax()
    if type(bpfile) is list:
        result = {}
        for bpf in bpfile:
            result = dict(reactions_biopax(bpf, protein_name_type = protein_name_type).items() + \
                result.items())
    else:
        result = reactions_biopax(bpfile, protein_name_type = protein_name_type)
    return result

def process_controls(controls, mandatory_refs = True):
    interactions = set([])
    ptms = []
    regulations = []
    prg = progress.Progress(len(controls), 'Processing interactions', 11)
    for c in controls.values():
        prg.step()
        if len(c['refs']) > 0 or not mandatory_refs:
            if c['controller'] is not None and len(c['controller']) > 0:
                for ctr in c['controller']:
                    if len(common.uniqList(ctr['members'])) == 1:
                        this_ctr = ctr['members'][0].split('-')[0]
                        ctd = c['controlled']
                        if ctd is not None:
                            #ctd['left'] is not None and ctd['right'] is not None:
                            for leftInst in itertools.product(*ctd['left']):
                                for rightInst in itertools.product(*ctd['right']):
                                    lr = common.uniqList(
                                        common.flatList(
                                            [l['members'] for l in leftInst] + \
                                            [r['members'] for r in rightInst]
                                        )
                                    )
                                    if len(lr) == 1:
                                        this_ctd = lr[0].split('-')[0]
                                        interactions.add(
                                            (this_ctr, this_ctd, c['type'], 
                                            ';'.join(c['refs'] \
                                                if len(c['refs']) > 0 \
                                                else ctd['refs']), 'directed')
                                        )
                                    else:
                                        modDiff = {}
                                        ptmsLeft = set([(ptms[0], ptm) \
                                            for l in leftInst \
                                                for ptms in l['ptms'].items() \
                                                    for ptm in ptms[1]])
                                        ptmsRight = set([(ptms[0], ptm) \
                                            for r in rightInst \
                                                for ptms in r['ptms'].items() \
                                                    for ptm in ptms[1]])
                                        ptmsDiff = ptmsLeft ^ ptmsRight
                                        diffUniProts = common.uniqList([ptm[0] for ptm in ptmsDiff])
                                        if len(diffUniProts) == 1:
                                            this_ctd = diffUniProts[0].split('-')[0]
                                            interactions.add(
                                                (this_ctr, this_ctd, c['type'], 
                                                ';'.join(c['refs'] \
                                                    if len(c['refs']) > 0 \
                                                    else ctd['refs']), 'directed')
                                            )
                                        else:
                                            lefts = [set(l['members']) for l in leftInst]
                                            rights = [set(r['members']) for r in rightInst]
                                            onlyLefts = [l for l in lefts if l not in rights]
                                            onlyRights = [r for r in rights if r not in lefts]
                                            diffs = []
                                            for l in onlyLefts:
                                                for r in onlyRights:
                                                    diff = l ^ r
                                                    if len(diff) == 1:
                                                        diffs.append(list(diff))
                                            diffs = common.uniqList(common.flatList(diffs))
                                            if len(diffs) == 1:
                                                this_ctd = diffs[0].split('-')[0]
                                                interactions.add(
                                                    (this_ctr, this_ctd, c['type'], 
                                                    ';'.join(c['refs'] \
                                                        if len(c['refs']) > 0 \
                                                        else ctd['refs']), 'undirected')
                                                )
            # if the controller is unknown
            # and the reaction has only 2 proteins
            # these most probably bind each other
            # to form a complex
            else:
                ctd = c['controlled']
                if ctd is not None:
                    for leftInst in itertools.product(*ctd['left']):
                        for rightInst in itertools.product(*ctd['right']):
                            lr = common.uniqList(
                                common.flatList(
                                    [l['members'] for l in leftInst] + \
                                    [r['members'] for r in rightInst]
                                )
                            )
                            if len(lr) == 2:
                                interactions.add(
                                    (lr[0].split('-')[0], lr[1].split('-')[0], c['type'], ';'.join(ctd['refs']))
                                )
    prg.terminate()
    return list(interactions), ptms, regulations

# Process Reactome SBML

def _reactome_id(obj, attr):
    return _reactome_extract_id(obj.attrs[attr])

def _reactome_extract_id(value):
    return int(value.split('_')[1])

def _reactome_res(obj):
    return _reactome_extract_res(obj.attrs['rdf:resource'])

def _reactome_extract_res(value):
    return value.split(':')[-1]

def _reactome_reactions():
    species = {}
    compartments = {}
    reactions = {}
    soup = reactome_bs()
    m = soup.find('model')
    for cp in m.find('listofcompartments').find_all('compartment'):
        compartments[_reactome_id(cp, 'id')] = cp.attrs['name']
    for sp in m.find('listofspecies').find_all('species'):
        cp = _reactome_id(sp, 'compartment')
        si = _reactome_id(sp, 'id')
        nm = sp.attrs['name']
        ids = []
        for i in sp.find('bqbiol:haspart').find_all('rdf:li'):
            ids.append(_reactome_res(i))
        ids = sorted(uniqList(ids))
        species[si] = {
            'name': nm,
            'comp': cp,
            'ids': ids
        }
    for rea in m.find('listofreactions').find_all('reaction'):
        ri = _reactome_id(rea, 'id')
        refs = []
        for r in rea.find('bqbiol:isdescribedby').find_all('rdf:li'):
            refs.append(_reactome_res(r))
        refs = sorted(uniqList(refs))
        reas = []
        for r in rea.find('listofreactants').find_all('speciesreference'):
            reas.append(_reactome_id(r, 'species'))
        reas = sorted(uniqList(reas))
        prds = []
        for p in rea.find('listofproducts').find_all('speciesreference'):
            prds.append(_reactome_id(p, 'species'))
        prds = sorted(uniqList(prds))
        note = rea.find('notes').text
        reactions[ri] = {
            'refs': refs,
            'reas': reas,
            'prds': prds,
            'note': note
        }
    return compartments, species, reactions

def _reactome_reactions_et():
    sbmlPfx = '{http://www.sbml.org/sbml/level2/version4}'
    compStr = '%scompartment' % sbmlPfx
    reacStr = '%sreaction' % sbmlPfx
    specStr = '%sspecies' % sbmlPfx
    species = {}
    compartments = {}
    reactions = {}
    sbmlfile = reactome_sbml()
    ctx = etree.iterparse(sbmlfile, events = ('end',))
    for ev, elem in ctx:
        if elem.tag == compStr:
            k, v = _reactome_compartment(elem)
            compartments[k] = v
        elif elem.tag == reacStr:
            k, v = _reactome_reaction(elem)
            reactions[k] = v
        elif elem.tag == specStr:
            k, v = _reactome_species(elem)
            species[k] = v
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    return compartments, species, reactions

def _reactome_compartment(elem):
    ci = _reactome_extract_id(elem.get('id'))
    nm = elem.get('name')
    return ci, nm

def _reactome_species(elem):
    bqBiolPfx = '{http://biomodels.net/biology-qualifiers/}'
    rdfPfx = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    hasPartStr = '%shasPart' % bqBiolPfx
    resStr = '%sresource' % rdfPfx
    si = _reactome_extract_id(elem.get('id'))
    cp = _reactome_extract_id(elem.get('compartment'))
    nm = elem.get('name')
    ids = sorted(uniqList(_reactome_collect_resources(elem, hasPartStr)))
    return si, {
        'name': nm,
        'comp': cp,
        'ids': ids
    }

def _reactome_reaction(elem):
    bqBiolPfx = '{http://biomodels.net/biology-qualifiers/}'
    rdfPfx = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    sbmlPfx = '{http://www.sbml.org/sbml/level2/version4}'
    specStr = 'species'
    spRefStr = '%sspeciesReference' % sbmlPfx
    isDescStr = '%sisDescribedBy' % bqBiolPfx
    resStr = '%sresource' % rdfPfx
    lofReaStr = '%slistOfReactants' % sbmlPfx
    lofPrdStr = '%slistOfProducts' % sbmlPfx
    ri = _reactome_extract_id(elem.get('id'))
    refs = _reactome_collect_resources(elem, isDescStr)
    reas = _reactome_collect_species(elem, lofReaStr)
    prds = _reactome_collect_species(elem, lofPrdStr)
    note = elem.find('note').text # prefix?
    return ri, {
        'refs': refs,
        'reas': reas,
        'prds': prds,
        'note': note
    }

def _reactome_collect_resources(elem, tag):
    rdfPfx = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    resStr = '%sresource' % rdfPfx
    liStr = '%sli' % rdfPfx
    res = []
    for i in  elem.find('.//%s' % tag).iterfind('.//%s' % liStr):
        res.append(_reactome_extract_res(i.get(resStr)))
    return res

def _reactome_collect_species(elem, tag):
    sbmlPfx = '{http://www.sbml.org/sbml/level2/version4}'
    spRefStr = '%sspeciesReference' % sbmlPfx
    specStr = 'species'
    res = []
    for sp in elem.find('.//%s' % tag).iterfind('.//%s' % spRefStr):
        res.apped(_reactome_extract_id(sp.get(specStr)))
    return res

def signalink_interactions():
    '''
    Reads and processes SignaLink3 interactions from local file.
    Returns list of interactions.
    '''
    repar = re.compile(r'.*\(([a-z\s]+)\)')
    notNeeded = set(['acsn', 'reactome'])
    edgesFile = os.path.join(common.ROOT, 'data', 
        data_formats.files['signalink']['edges'])
    nodesFile = os.path.join(common.ROOT, 'data', 
        data_formats.files['signalink']['nodes'])
    nodes = {}
    interactions = []
    def _get_attr(attrs, attrName):
        return _process_attr(attrs[attrName]) if attrName in attrs else ''
    def _process_attr(attr):
        m = repar.match(attr)
        if m is not None:
            return m.groups()[0]
        else:
            return attr
    with open(nodesFile, 'r') as f:
        for l in f:
            if len(l) > 0:
                l = l.split('\t')
                _id = int(l[0])
                uniprot = l[1].replace('uniprot:', '')
                pathways = [pw.split(':')[-1] for pw in l[4].split('|') \
                    if pw.split(':')[0] not in notNeeded]
                nodes[_id] = [uniprot, pathways]
    prg = progress.Progress(os.path.getsize(edgesFile), 'Reading file', 33)
    with open(edgesFile, 'r') as f:
        lPrev = None
        for l in f:
            prg.step(len(l))
            l = l.strip().split('\t')
            if lPrev is not None:
                l = lPrev + l[1:]
                lPrev = None
            if len(l) == 13:
                if l[-1] == '0':
                    dbs = [_process_attr(db.split(':')[-1]) \
                        for db in l[9].replace('"', '').split('|')]
                    dbs = list(set(dbs) - notNeeded)
                    if len(dbs) == 0:
                        continue
                    idSrc = int(l[1])
                    idTgt = int(l[2])
                    uniprotSrc = l[3].replace('uniprot:', '')
                    uniprotTgt = l[4].replace('uniprot:', '')
                    refs = [ref.split(':')[-1] for ref in l[7].split('|')]
                    attrs = dict(tuple(attr.strip().split(':', 1)) for attr in l[8].\
                        replace('"', '').split('|'))
                    interactions.append([
                        uniprotSrc,
                        uniprotTgt,
                        ';'.join(refs),
                        ';'.join(dbs),
                        _get_attr(attrs, 'effect'),
                        _get_attr(attrs, 'is_direct'),
                        _get_attr(attrs, 'is_directed'),
                        _get_attr(attrs, 'molecular_background'),
                        ';'.join(nodes[idSrc][1]),
                        ';'.join(nodes[idTgt][1])
                    ])
            else:
                lPrev = l
    prg.terminate()
    return interactions

def get_laudanna_directions():
    '''
    Downloads and processes the SignalingFlow edge attributes
    from Laudanna Lab.
    Returns list of directions.
    '''
    url = data_formats.urls['laudanna']['sigflow']
    data = curl(url, silent = False)
    data = data.split('\n')[1:]
    directions = []
    for l in data:
        if len(l) > 0:
            directions.append(l.split(' =')[0].split(' (pp) '))
    return directions

def get_laudanna_effects():
    '''
    Downloads and processes the SignalingDirection edge attributes
    from Laudanna Lab.
    Returns list of effects.
    '''
    url = data_formats.urls['laudanna']['sigdir']
    data = curl(url, silent = False)
    data = data.split('\n')[1:]
    effects = []
    for l in data:
        if len(l) > 0:
            l = l.split(' = ')
            effects.append(l[0].split(' (pp) ') + [l[1]])
    return effects

def get_acsn_effects():
    '''
    Processes ACSN data, returns list of effects.
    '''
    negatives = set(['NEGATIVE_INFLUENCE', 'UNKNOWN_NEGATIVE_INFLUENCE'])
    positives = set(['TRIGGER', 'POSITIVE_INFLUENCE', 'UNKNOWN_POSITIVE_INFLUENCE'])
    directed = set([
            'UNKNOWN_TRANSITION', 'INTERACTION_TYPE', 'KNOWN_TRANSITION_OMITTED', 
            'INHIBITION', 'UNKNOWN_POSITIVE_INFLUENCE', 'PROTEIN_INTERACTION',
            'UNKNOWN_CATALYSIS', 'POSITIVE_INFLUENCE', 'STATE_TRANSITION', 
            'TRANSLATION', 'UNKNOWN_NEGATIVE_INFLUENCE', 'NEGATIVE_INFLUENCE', 
            'MODULATION', 'TRANSCRIPTION', 'COMPLEX_EXPANSION', 'TRIGGER', 'CATALYSIS',
            'PHYSICAL_STIMULATION', 'UNKNOWN_INHIBITION', 'TRANSPORT'])
    data = acsn_ppi()
    effects = []
    for l in data:
        if len(l) == 4:
            eff = set(l[2].split(';'))
            if len(eff & negatives) > 0:
                effects.append([l[0], l[1], '-'])
            elif len(eff & positives) > 0:
                effects.append([l[0], l[1], '+'])
            elif len(eff & directed) > 0:
                effects.append([l[0], l[1], '*'])
    return effects

def get_wang_effects():
    '''
    Downloads and processes Wang Lab HumanSignalingNetwork.
    Returns list of effects.
    '''
    url = data_formats.urls['wang']['url']
    data = curl(url, silent = False)
    data = data.split('\n')
    effects = []
    nodes = {}
    reading_nodes = False
    reading_edges = False
    for l in data:
        if len(l.strip()) == 0:
            reading_nodes = False
            reading_edges = False
        l = l.split(',')
        if reading_nodes:
            nodes[l[0]] = l[1]
        if reading_edges:
            effects.append([nodes[l[0]], nodes[l[1]], l[2]])
        if l[0].startswith('Node'):
            reading_nodes = True
        if l[0].startswith('From'):
            reading_nodes = False
            reading_edges = True
    return effects

def biogrid_interactions(organism = 9606, htp_limit = 1):
    '''
    Downloads and processes BioGRID interactions.
    Keeps only the "low throughput" interactions.
    Returns list of interactions.
    
    @organism : int
        NCBI Taxonomy ID of organism.
    @htp_limit : int
        Exclude interactions only from references
        cited at more than this number of interactions.
    '''
    organism = str(organism)
    interactions = []
    refc = []
    url = data_formats.urls['biogrid']['url']
    f = curl(url, silent = False, large = True).values()[0]
    nul = f.readline()
    for l in f:
        l = l.split('\t')
        if len(l) > 17:
            if l[17].startswith('Low') and l[15] == organism and l[16] == organism:
                interactions.append([l[7], l[8], l[14]])
                refc.append(l[14])
    refc = Counter(refc)
    interactions = [i for i in interactions if refc[i[2]] <= htp_limit]
    return interactions

def acsn_ppi(keep_in_complex_interactions = True):
    '''
    Processes ACSN data from local file.
    Returns list of interactions.
    
    @keep_in_complex_interactions : bool
        Whether to include interactions from complex expansion.
    '''
    nfname = data_formats.files['acsn']['names']
    pfname = data_formats.files['acsn']['ppi']
    names = {}
    interactions = []
    with open(nfname, 'r') as f:
        for l in f:
            l = l.strip().split('\t')
            names[l[0]] = l[2:]
    with open(pfname, 'r') as f:
        nul = f.readline()
        for l in f:
            l = l.strip().split('\t')
            if l[0] in names:
                for a in names[l[0]]:
                    if l[2] in names:
                        for b in names[l[2]]:
                            if keep_in_complex_interactions:
                                if 'PROTEIN_INTERACTION' in l[1]:
                                    l[1].replace('COMPLEX_EXPANSION', 
                                        'IN_COMPLEX_INTERACTION')
                            interactions.append([a, b, l[1], l[3]])
    return interactions

def get_graphviz_attrs():
    '''
    Downloads graphviz attribute list from graphviz.org.
    Returns 3 dicts of dicts: graph_attrs, vertex_attrs and edge_attrs.
    '''
    url = data_formats.urls['graphviz']['url']
    html = curl(url)
    soup = bs4.BeautifulSoup(html)
    vertex_attrs = {}
    edge_attrs = {}
    graph_attrs = {}
    for tbl in soup.find_all('table'):
        if tbl.find('tr').text.startswith('Name'):
            for r in tbl.find_all('tr'):
                r = r.find_all('td')
                if len(r) > 0:
                    usedby = r[1].text
                    this_attr = {
                            'type': r[2].text.strip(),
                            'default': r[3].text.strip(),
                            'min': r[4].text.strip(),
                            'notes': r[5].text.strip()
                        }
                    attr_name = r[0].text.strip()
                    if 'N' in usedby:
                        vertex_attrs[attr_name] = this_attr
                    if 'E' in usedby:
                        edge_attrs[attr_name] = this_attr
                    if 'G' in usedby:
                        graph_attrs[attr_name] = this_attr
            break
    return graph_attrs, vertex_attrs, edge_attrs

def get_phosphosite(cache = True):
    '''
    Downloads curated and HTP data from Phosphosite,
    from preprocessed cache file if available.
    Processes BioPAX format.
    Returns list of interactions.
    '''
    curated_cache = data_formats.files['phosphosite']['curated']
    noref_cache = data_formats.files['phosphosite']['noref']
    if cache and os.path.exists(curated_cache) and os.path.exists(noref_cache):
        return (
            pickle.load(open(curated_cache, 'rb')),
            pickle.load(open(noref_cache, 'rb'))
        )
    result_curated = []
    result_noref = []
    url = data_formats.urls['psite_bp']['url']
    bpax = curl(url, silent = False, large = True)
    xml = ET.parse(bpax)
    xmlroot = xml.getroot()
    bpprefix = '{http://www.biopax.org/release/biopax-level3.owl#}'
    rdfprefix = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
    proteins = {}
    for p in xmlroot.iter(bpprefix+'ProteinReference'):
        psid = p.attrib[rdfprefix+'ID']
        db = p.find(bpprefix+'xref').find(bpprefix+'UnificationXref').find(bpprefix+'db').text
        up = p.find(bpprefix+'xref').find(bpprefix+'UnificationXref').find(bpprefix+'id').text
        tax = ''
        if p.find(bpprefix+'organism') is not None:
            tmp = p.find(bpprefix+'organism')
            if rdfprefix+'resource' in tmp.attrib:
                tax = tmp.attrib[rdfprefix+'resource'].split('_')[1]
        if db == 'UniProtKB':
            up = up[0:6]
        proteins[psid] = {'id': up, 'db': db, 'species': tax, 'psid': psid}
    evidences = {}
    for p in xmlroot.iter(bpprefix+'EvidenceCodeVocabulary'):
        evid = p.attrib[rdfprefix+'ID'].split('_')[1]
        evname = p.find(bpprefix+'term').text
        evidences[evid] = evname
    ev_short = {'0113': 'WB', '0427': 'MS', '0074': 'MA', '0421': 'AB'}
    nosrc = []
    notgt = []
    norefs = []
    noev = []
    noth = []
    edges = []
    for c in xmlroot.findall(bpprefix+'Catalysis'):
        if rdfprefix+'resource' in c.find(bpprefix+'controller').attrib:
            src = 'po_'+c.find(bpprefix+'controller').attrib[rdfprefix+'resource'].split('_')[1]
        else:
            srcProt =  c.find(bpprefix+'controller').find(bpprefix+'Protein')
            if srcProt is not None:
                src = 'po_'+srcProt.attrib[rdfprefix+'ID'].split('_')[1]
            else:
                nosrc.append(c)
        tgtProt = c.find(bpprefix+'controlled').iter(bpprefix+'ProteinReference')
        tgt = next(tgtProt, None)
        if tgt is not None:
            tgt = tgt.attrib[rdfprefix+'ID']
        else:
            tgtProt = c.find(bpprefix+'controlled').iter(bpprefix+'entityReference')
            tgt = next(tgtProt, None)
            if tgt is not None:
                if rdfprefix+'resource' in tgt.attrib:
                    tgt = tgt.attrib[rdfprefix+'resource'][1:]
            else:
                tgtProt = c.find(bpprefix+'controlled').iter(bpprefix+'left')
                tgt = next(tgtProt, None)
                if tgt is not None:
                    if rdfprefix+'resource' in tgt.attrib:
                        tgt = 'po_'+tgt.attrib[rdfprefix+'resource'].split('_')[1]
                else:
                    notgt.append(c)
        refs = c.iter(bpprefix+'PublicationXref')
        pmids = []
        for r in refs:
            pm = r.attrib[rdfprefix+'ID'].split('_')
            if pm[0] == 'pmid':
                pmids.append(pm[1])
        refs = c.iter(bpprefix+'evidence')
        for r in refs:
            rrefs = r.iter(bpprefix+'xref')
            for rr in rrefs:
                if rdfprefix+'resource' in rr.attrib:
                    pm = rr.attrib[rdfprefix+'resource'].split('_')
                    if pm[0] == 'pubmed':
                        pmids.append(pm[1])
        evs = []
        for e in c.iter(bpprefix+'evidenceCode'):
            if rdfprefix+'resource' in e.attrib:
                evs.append(ev_short[e.attrib[rdfprefix+'resource'].split('_')[1]])
            else:
                ev = e.find(bpprefix+'EvidenceCodeVocabulary')
                evs.append(ev_short[ev.attrib[rdfprefix+'ID'].split('_')[1]])
        for e in c.iter(bpprefix+'evidence'):
            if rdfprefix+'resource' in e.attrib:
                ev = e.attrib[rdfprefix+'resource'].split('_')
                if len(ev) == 4:
                    if len(ev[3]) == 4:
                        evs.append(ev_short[ev[3]])
        if (src is not None and tgt is not None and src in proteins and tgt in proteins and 
                proteins[src]['id'] is not None and proteins[tgt]['id'] is not None):
            edges.append({'src': proteins[src], 'tgt': proteins[tgt], 'pmids': list(set(pmids)), 
                        'evs': list(set(evs))})
            if len(evs) == 0:
                noev.append(c)
            if len(pmids) == 0:
                norefs.append(c)
            if len(evs) == 0 and len(pmids) == 0:
                noth.append(c)
    for e in edges:
        this_iaction = [e['src']['id'], e['tgt']['id'], e['src']['species'], \
            e['tgt']['species'], ';'.join(e['evs']), ';'.join(e['pmids'])]
        if len(this_iaction[-1]) > 0:
            result_curated.append(this_iaction)
        else:
            result_noref.append(this_iaction)
    pickle.dump(result_curated, open(curated_cache, 'wb'))
    pickle.dump(result_noref, open(noref_cache, 'wb'))
    return result_curated, result_noref

def get_phosphosite_curated():
    '''
    Loads literature curated PhosphoSite data,
    from preprocessed cache file if available.
    Returns list of interactions.
    '''
    curated_cache = data_formats.files['phosphosite']['curated']
    if not os.path.exists(curated_cache):
        curated, noref = get_phosphosite()
        return curated
    else:
        return pickle.load(open(cureted_cache, 'rb'))

def get_phosphosite_noref():
    '''
    Loads HTP PhosphoSite data,
    from preprocessed cache file if available.
    Returns list of interactions.
    '''
    noref_cache = data_formats.files['phosphosite']['noref']
    if not os.path.exists(noref_cache):
        curated, noref = get_phosphosite()
        return noref
    else:
        return pickle.load(open(noref_cache, 'rb'))

def phosphosite_directions(organism = 'human'):
    '''
    From curated and HTP PhosphoSite data generates a
    list of directions.
    '''
    curated, noref = get_phosphosite()
    return [i[:2] for i in curated + noref \
        if i[2] == organism and i[3] == organism]

def get_lit_bm_13():
    '''
    Downloads and processes Lit-BM-13 dataset, the high confidence
    literature curated interactions from CCSB.
    Returns list of interactions.
    '''
    url = data_formats.urls['hid']['lit-bm-13']
    data = curl(url, silent = False)
    return map(lambda l: l.strip().split('\t'), data.split('\n')[1:])

def get_ca1():
    '''
    Downloads and processes the CA1 signaling network (Ma\'ayan 2005).
    Returns list of interactions.
    '''
    url = data_formats.urls['ca1']['url']
    data = curl(url, silent = False, files_needed = ['S1.txt'])
    return filter(lambda l: len(l) == 13,
        map(lambda l: l.strip().split(), data['S1.txt'].split('\n')[1:])
    )

def get_ccmap(organism = 9606):
    '''
    Downloads and processes CancerCellMap.
    Returns list of interactions.
    
    @organism : int
        NCBI Taxonomy ID to match column #7 in nodes file.
    '''
    organism = '%u' % organism
    interactions = []
    nodes_url = data_formats.urls['ccmap']['nodes']
    edges_url = data_formats.urls['ccmap']['edges']
    nodes = curl(nodes_url, silent = False, files_needed = ['cell-map-node-attributes.txt'])
    edges = curl(edges_url, silent = False, files_needed = ['cell-map-edge-attributes.txt'])
    nodes = dict(
        map(lambda l: (l[1], l[2].split(':')),
            filter(lambda l: l[5] == 'protein' and l[6] == organism,
                filter(lambda l: len(l) == 7,
                    map(lambda l: l.strip().split('\t'), 
                        nodes['cell-map-node-attributes.txt'].split('\n')[1:])
                )
            )
        )
    )
    edges = filter(lambda l: len(l) == 7,
        map(lambda l: l.strip().split('\t'), 
            edges['cell-map-edge-attributes.txt'].split('\n')[1:])
    )
    print len(nodes), len(edges)
    for e in edges:
        if e[1] != 'IN_SAME_COMPONENT' and e[3] in nodes and e[4] in nodes:
            for src in nodes[e[3]]:
                for tgt in nodes[e[4]]:
                    interactions.append([
                            src, tgt,
                            'directed' if e[1] == 'STATE_CHANGE' else 'undirected',
                            e[6].strip(';').replace('PUBMED:', '')
                        ])
    return interactions