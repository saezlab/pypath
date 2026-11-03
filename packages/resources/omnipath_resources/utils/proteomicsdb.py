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
from past.builtins import xrange, range, reduce

import json
import base64
import re
import os
import sys
import importlib as imp

try:
    import cPickle as pickle
except:
    import pickle

try:
    import ssl
except ImportError:
    sys.stdout.write("\t:: Error: no ssl support :(\n\n")
    sys.stdout.flush()

# from this module:
import pypath.share.curl as curl
import pypath.share.progress as progress
from pypath.share.common import *


class ProteomicsDB(object):
    
    def __init__(self, username, password, output_format='json'):
        '''
        This is an extensible class for downloading and processing data
        from ProteomicsDB. Now 2 of the 10 available APIs implemented here,
        but feel free to write functions for the other APIs. 
        To find out more about ProteomicsDB, take a look at Wilhelm et al. 2014, Nature:
        http://www.nature.com/nature/journal/v509/n7502/full/nature13319.html
        To read a comprehensive descritpion of the APIs, visit here:
        https://www.proteomicsdb.org/proteomicsdb/#api

        @username : str
            Registered and API enabled user for ProteomicsDB. To have such a
            user, you need first to register, AND then write an e-mail to the
            address given on the webpage. In a couple of days the admins will
            enable the API for your user.
        @password : str
            Password of the user.
        @output_format : str
            Either 'json' or 'xml'. Some functions in this module process
            JSON further and give certain objects.
        '''
        self.auth = [
            (b'Authorization: Basic %s' %
            base64.encodestring((
                "%s:%s" % (username, password)).encode('ascii')
            )).decode('ascii').rstrip('\n')
        ]
        self.port = 443
        self.output_format = output_format
        self.host = 'www.proteomicsdb.org'
        self.tissues = []
        self.expression = {}
        self.samples = {}
        self.tissues_loaded = set([])
        self.testurl = '''/proteomicsdb/logic/api/proteinpeptideresult.xsodata/InputParams(PROTEINFILTER='Q92769')/Results?$select=UNIQUE_IDENTIFIER,PROTEIN_NAME,START_POSITION,END_POSITION,PEPTIDE_SEQUENCE,PEPTIDE_MASS,Q_VALUE,RANK,SCORE,SEARCH_ENGINE&$filter=PEPTIDE_MASS%20gt%201000%20&$format=xml '''
        self.urls = {
            'tissues': 'https://www.proteomicsdb.org/proteomicsdb/logic/api/'
            'tissuelist.xsodata/CA_AVAILABLEBIOLOGICALSOURCES_API?'
            '$select=TISSUE_ID,TISSUE_NAME,TISSUE_GROUP_NAME,TISSUE_CATEGORY,'
            'SCOPE_ID,SCOPE_NAME,QUANTIFICATION_METHOD_ID,'
            'QUANTIFICATION_METHOD_NAME,MS_LEVEL&$format=%s',
            'proteinpertissue':
            'https://www.proteomicsdb.org/proteomicsdb/logic/api/'
            'proteinspertissue.xsodata/InputParams(TISSUE_ID=%%27%s%%27,'
            'CALCULATION_METHOD=%u,SWISSPROT_ONLY=%u,NO_ISOFORM=%u)/Results?'
            '$select=ENTRY_NAME,UNIQUE_IDENTIFIER,DATABASE,PROTEIN_DESCRIPTION,'
            'PEPTIDES,SAMPLE_NAME,SAMPLE_DESCRIPTION,UNNORMALIZED_EXPRESSION,'
            'NORMALIZED_EXPRESSION&$format=%s'
        }
    
    def reload(self):
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist=[modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    def query(self, api, param, silent=False, large=False):
        '''
        Retrieves data from the API. 

        @api : str
            Shold be one of the 10 API sections available.
        @param : tuple
            Tuple of the parameters according to the API.
        @large : bool
            Passed to the curl wrapper function. If True, 
            the file will be written to disk, and a file 
            object open for reading is returned; if False,
            the raw data will be returned, in case of JSON,
            converted to python object, in case of XML, as
            a string.
        '''
        url = self.urls[api] % param
        # long timeout is given, because huge files (hundreds MB) take time to
        # load
        c = curl.Curl(
            url,
            req_headers=self.auth,
            silent=silent,
            timeout=1200,
            large=large)
        
        data = c.fileobj
        self.tmp = c
        
        if self.output_format == 'json' and not large:
            self.result = self.get_json(c.result)
        else:
            self.result = c.fileobj

    def get_json(self, content):
        
        return json.loads(content)['d']['results']

    def get_tissues(self):
        '''
        Gets an annotated list of all tissues for which ProteomicsDB has 
        expression data. Result stored in `ProteomicsDB.tissues`. 
        '''
        self.query('tissues', (self.output_format, ))
        self.tissues = self.result

    def which_tissues(self, name, value):
        if len(self.tissues) == 0:
            self.get_tissues()
        value = set(value) if type(value) is list else set([value])
        match = []
        for tis in self.tissues:
            if tis[name] in value:
                match.append(tis)
        return match

    def get_proteins(self,
                     tissue_id,
                     calculation_method=0,
                     swissprot_only=1,
                     no_isoform=1):
        '''
        '''
        for i in xrange(3):
            
            self.query(
                'proteinpertissue',
                (tissue_id, calculation_method, swissprot_only, no_isoform,
                 self.output_format),
                large=True)
            
            if hasattr(self.result, 'read'):
                break

    def get_pieces(self, size=20480, delimiters=('{', '}')):
        '''
        A generator for reading huge files (hundreds of MBs).
        Reads segments of @size, searches for self-contained
        JSON objects, and returns a list of them. 

        @size : int
            Size to read at once (in Bytes).
        @delimiters : tuple
            Starting and closing delimiters. By default, these
            are curly braces, to return individual JSON objects
            of the largest possible size.
        '''
        piece = ''
        while True:
            new = self.result.read(size).decode('ascii')
            if len(new) == 0:
                yield []
                break
            buff = piece + new
            piece = ''
            open_br = 0
            results = []
            for c in buff:
                if c == delimiters[0]:
                    open_br += 1
                if open_br > 0:
                    piece += c
                if c == delimiters[1]:
                    open_br -= 1
                    if open_br == 0:
                        results.append(piece)
                        piece = ''
            yield results

    def get_expression(self, normalized=True, tissue_average=False):
        '''
        Extracts normalized or unnormalized expression data from 
        previously downloaded data, stored on disk, and 
        opened for reading in file object `ProteomicsDB.result`. 
        Optionally averages data per tissue.

        @normalized : bool
            Read normalized or unnormalized expression values.
        @tissue_average : bool
            Read and store data for each samples, or keep only the
            mean value per tissue.
        '''
        non_digit = re.compile(r'[^\d.-]+')
        #try:
        self.result.seek(0)
        if hasattr(self.result, 'read'):
            nul = self.result.read(17)
            self.current_samples = set([])
            for pp in self.get_pieces():
                for p in pp:
                    protein = json.loads(p)
                    if protein['SAMPLE_NAME'] not in self.expression:
                        self.expression[protein['SAMPLE_NAME']] = {}
                    self.current_samples.add(protein['SAMPLE_NAME'])
                    self.expression[protein['SAMPLE_NAME']]\
                        [protein['UNIQUE_IDENTIFIER']] = \
                        float(non_digit.sub('', protein['NORMALIZED_EXPRESSION']))\
                        if normalized else \
                        float(non_digit.sub('', protein[
                                'UNNORMALIZED_EXPRESSION']))
            self.result.close()
            self.result = None
        try:
            pass
        except:
            sys.stdout.write(
                'Error in pypath.proteomicsdb.py/get_expression():\n')
            sys.stdout.write('Result type: %s\n' % str(type(self.result)))
            sys.stdout.write('Result mode: %s\n' % str(self.result.mode))
            sys.stdout.write('Result name: %s\n' % str(self.result.name))
            sys.stdout.flush()

    def tissues_x_proteins(self, normalized=True, tissues=None):
        '''
        For all tissues downloads the expression of all the proteins.
        In the result, a dict of dicts will hold the expression values
        of each proteins, grouped by samples.
        '''
        self.get_tissues()
        tissues_selected = set([
            t['TISSUE_ID'] for t in self.tissues
            if tissues is None or t['TISSUE_ID'] in tissues
        ]) - self.tissues_loaded
        prg = progress.Progress(
            len(tissues_selected),
            'Downloading expression data',
            1,
            percent=False)
        for tis in tissues_selected:
            prg.step()
            sys.stdout.write('Querying tissue %s\n' % tis)
            sys.stdout.flush()
            self.get_proteins(tis)
            if not hasattr(self.result, 'read'):
                sys.stdout.write('\tFailed: %s\n' % tis)
                sys.stdout.flush()
            else:
                self.tissues_loaded.add(tis)
                self.get_expression(normalized)
                if tis not in self.samples:
                    self.samples[tis] = []
                self.samples[tis] = unique_list(self.samples[tis] + list(
                    self.current_samples))
                self.current_samples = set([])
        prg.terminate()

    def save(self, outf=None):
        self.result = None
        outf = outf if outf is not None else os.path.join(
            'cache', 'proteomicsdb.pickle')
        pickle.dump(self, open(outf, 'wb'))

    def load(self, pfile=None):
        pfile = pfile if pfile is not None else os.path.join(
            'cache', 'proteomicsdb.pickle')
        if os.path.exists(pfile):
            loaded = pickle.load(open(pfile, 'rb'))
            for k, v in iteritems(loaded.__dict__):
                if not hasattr(v, '__call__'):
                    setattr(self, k, v)
            sys.stdout.write(
                '\t:: Loaded %u samples of %u tissues from file %s\n' %
                (len(self.expression.keys()), len(self.tissues_loaded), pfile))
        else:
            sys.stdout.write('\t:: File not found: %s\n' % pfile)
        sys.stdout.flush()

    def pandas_matrix(self):
        '''
        Returns expression data in a pandas matrix. Not implemented.
        '''
        pass
