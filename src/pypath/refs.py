#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2015 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import webbrowser

import pypath.curl as curl
import pypath.common as common
import pypath.urls as urls

class Reference(object):
    
    def __init__(self, pmid):
        self.pmid = str(pmid.strip())
        
    def __eq__(self, other):
        return other.__class__.__name__ == self.__class__.__name__ \
            and self.pmid == other.pmid
    
    def __hash__(self):
        return hash(self.pmid)
    
    def open(self):
        dataio.open_pubmed(self.pmid)
    
    def __str__(self):
        return self.pmid
    
    def info(self):
        return dataio.get_pubmeds([self.pmid])


def open_pubmed(pmid):
    '''
    Opens PubMed record in web browser.
    
    @pmid : str or int
        PubMed ID
    '''
    pmid = str(pmid)
    url = urls.urls['pubmed']['url'] % pmid
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
    url = urls.urls['pubmed-eutils']['conv'] % ','.join(str(i) for i in idList)
    c = curl.Curl(url, silent = True)
    data = c.result
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
