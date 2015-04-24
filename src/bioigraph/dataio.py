#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `bioigraph` python module
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
import urllib
import urllib2
import codecs
import gzip
import zipfile
import tarfile
import bs4
import xml.etree.ElementTree as ET
import hashlib
import time
import struct
import json
from bioservices import WSDLService
from contextlib import closing
from fabric.network import connect, HostConnectionCache
from fabric.state import env

import data_formats
import progress
import common
import intera
import residues
import seq as se

CURSOR_UP_ONE = '\x1b[1A'
ERASE_LINE = '\x1b[2K'

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

def test(debug_type, debug_msg):
    print "debug(%d): %s" % (debug_type, debug_msg)

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
        init_fun = 'get_jsessionid', follow = True, large = False):
    # either from cache or from download, we load the data into StringIO:
    multifile = False
    domain = url.replace('https://', '').replace('http://','').\
        replace('ftp://','').split('/')[0]
    # first try to find file in cache:
    if cache:
        # outf param is to give a unique name to data
        # downloaded previously by post requests
        outf = outf if outf is not None else url.split('/')[-1].split('?')[0]
        poststr = '' if post is None else \
            '?' + '&'.join(sorted([i[0]+'='+i[1] for i in post.items()]))
        urlmd5 = hashlib.md5(url+poststr).hexdigest()
        if not os.path.exists(os.path.join(os.getcwd(),'cache')):
            os.mkdir(os.path.join(os.getcwd(),'cache'))
        cachefile = os.path.join(os.getcwd(),'cache',urlmd5+'-'+outf)
        usecache = True if os.path.exists(cachefile) else False
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
            c.setopt(c.URL, url)
        c.setopt(c.FOLLOWLOCATION, follow)
        c.setopt(c.CONNECTTIMEOUT, 15)
        c.setopt(c.TIMEOUT, timeout)
        if type(req_headers) is list:
            c.setopt(c.HTTPHEADER, req_headers)
        c.setopt(c.WRITEFUNCTION, result.write)
        c.setopt(c.HEADERFUNCTION, headers.append)
        # if debug is necessary:
        if debug:
            c.setopt(pycurl.VERBOSE, 1)
            # c.setopt(pycurl.DEBUGFUNCTION, test)
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
                c.perform()
                if url.startswith('http'):
                    status = c.getinfo(pycurl.HTTP_CODE)
                if url.startswith('ftp'):
                    status = 500
                    for h in headers:
                        if h.startswith('226'):
                            status = 200
                            break
                break
            except:
                status = 500
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
        return curl(url = url, req_headers = req_headers, silent = silent, 
            debug = debug, outf = outf, compr = compr, encoding = encoding, 
            files_needed = files_needed, timeout = timeout, large = large)
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
            res.close()
        elif url.endswith('gz') or compr == 'gz':
            res = gzip.GzipFile(fileobj=result, mode='rb')
            if not large:
                res = res.read()
                res = res.decode(encoding)
                res = res.encode('utf-8')
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
                        results[k] = results[k].decode(encoding)
                        results[k] = results[k].encode('utf-8')
                    except:
                        pass
        if cache and not usecache and not large:
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

def read_table(fileObject, cols, sep = '\t', sep2 = None, rem = [], hdr = None):
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
    if 'readline' not in fileObject.__class__.__dict__:
        funname = sys._getframe().f_code.co_name
        sys.stdout.write('\tERROR: %s() expects file like object (file opened for read'\
            ', or StringIO buffer, etc)\n'%funname)
    fileObject.seek(0)
    res = []
    if hdr:
        for h in xrange(0,hdr):
            null = fileObject.readline()
            del null
    for l in fileObject:
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
    fileObject.close()
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
    return result

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
    pdb_re = re.compile('[0-9A-Z]{4}')
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
    data = curl(data_formats.urls['corum']['url'],silent=False)
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
        soup = bs4.BeautifulSoup(xml)
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
        soup = bs4.BeautifulSoup(data)
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
    table = read_table(buff,cols=cols,sepLevel2=subf,hdr=1)
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
    table = read_table(buff,cols=cols,sep=',',hdr=1)
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

def pdb_residue_mapper(pdb):
    
    return mapper

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
        soup = bs4.BeautifulSoup(data.read())
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
    data = read_table(buff, cols = cols, hdr = 1, sep = '\t')
    return data

def get_psite_phos(raw = True, organism = 'human'):
    url = data_formats.urls['psite_kin']['url']
    data = curl(url, silent = False, compr = 'gz', encoding = 'iso-8859-1')
    cols = {
        'kinase': 1,
        'kinase_org': 4,
        'substrate': 7,
        'substrate_org': 9,
        'residue': 10,
        'motif': 12
    }
    buff = StringIO()
    buff.write(data)
    data = read_table(buff, cols = cols, sep = '\t', hdr = 4)
    result = []
    non_digit = re.compile(r'[^\d.-]+')
    motre = re.compile(r'(_*)([A-Za-z]+)(_*)')
    for r in data:
        if r['kinase_org'] == organism and r['substrate_org'] == organism:
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
                if organism is None or (r['kinase_org'] == organism and \
                    r['substrate_org'] == organism):
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
    data = read_table(buff, cols = cols, sep = '\t', hdr = 4)
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
        soup = bs4.BeautifulSoup(data)
        src = 'cache'
    else:
        data = curl(url, post = post, silent = False, cache = False, 
            req_headers = headers)
        soup = bs4.BeautifulSoup(data)
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
            soup = bs4.BeautifulSoup(data)
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
    soup = bs4.BeautifulSoup(data)
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
    soup = bs4.BeautifulSoup(data)
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
    soup = bs4.BeautifulSoup(data)
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
    soup = bs4.BeautifulSoup(data)
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
    data = curl(url, silent = False, \
        files_needed = [data_formats.urls['p_elm']['psites']])
    data = data[data_formats.urls['p_elm']['psites']]
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
                'references': [non_digit.sub('', r) for r in l[4].split(';')],
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
    soup = bs4.BeautifulSoup(data)
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
    non_digit = re.compile(r'[^\d.-]+')
    for url in data_formats.urls['dbptm']['urls']:
        extra = curl(url, silent = False)
        for k, data in extra.iteritems():
            data = [x.split('\t') for x in data.split('\n')]
            for l in data:
                if len(l) > 8:
                    resnum = int(non_digit.sub('', l[2]))
                    result.append({
                        'substrate': l[1],
                        'typ': l[7].lower(),
                        'resaa': l[8][6],
                        'resnum': resnum,
                        'instance': l[8],
                        'references': l[4].split(';'),
                        'source': l[5].split()[0],
                        'kinase': None if byre.match(l[3]) is None else \
                            byre.match(l[3]).groups(1)[0].split(' and '),
                        'start': resnum - 6,
                        'end': resnum + 6
                    })
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
                    ';'.join(r['references'])
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
    reunip = re.compile(r'uniprotkb:([A-Z0-9]{6})')
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

def get_oreganno(organism = 'Homo sapiens'):
    nsep = re.compile(r'([-A-Za-z0-9]{3,})[\s/\(]*.*')
    nrem = re.compile(r'[-/]')
    result = []
    url = data_formats.urls['oreganno']['url']
    data = curl(url, silent = False)
    data = [x.split('\t') for x in data.split('\n') if len(x) > 0][1:]
    for l in data:
        if l[0] == organism and \
            l[10].startswith('TRANSCRIPTION FACTOR BINDING SITE') and \
            not l[11].startswith('UNKNOWN') and not l[14].startswith('UNKNOWN'):
            result.append([
                l[11] if len(l[11]) < 3 else nrem.sub('', nsep.findall(l[11])[0]),
                l[14] if len(l[14]) < 3 else nrem.sub('', nsep.findall(l[14])[0]), l[18]])
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

def get_go_quick(organism = 9606):
    terms = {'C': {}, 'F': {}, 'P': {}}
    names = {}
    url = data_formats.urls['quickgo']['url'] % organism
    data = curl(url, silent = False)
    data = [x.split('\t') for x in data.split('\n') if len(x) > 0]
    del data[0]
    for l in data:
        try:
            if l[0] not in terms[l[3][0]]:
                terms[l[3][0]][l[0]] = []
            terms[l[3][0]][l[0]].append(l[1])
            names[l[1]] = l[2]
        except:
            print l
    return {'terms': terms, 'names': names}

def netpath_names():
    repwnum = re.compile(r'_([0-9]+)$')
    result = {}
    url = data_formats.urls['netpath_names']['url']
    html = curl(url, silent = False)
    soup = bs4.BeautifulSoup(html)
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
                    if A[0] != B[0]:
                        result.append(list(A) + list(B) + \
                            [';'.join(refs), ';'.join(mets), intTyp, pwname])
    return result

class ProteomicsDB(object):
    
    def __init__(self):
        # BRENDA Tissue Ontology
        self.bto = get_ontology('BTO')
        # subs: uniprot, ms_level, scope
        self.exp_url = data_formats.urls['protdb_exp']['url']
        # subs: bto, swissprot_only, isoform
        self.tis_url = data_formats.urls['protdb_tis']['url']
    
    def get_expression(self, uniprot, ms_level = 1, scope = 1):
        url = self.exp_url % (uniprot, ms_level, scope)
        data = curl(url)
        try:
            data = json.loads(data)
            return data
        except:
            return {'url': url, 'response': data}
