#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2020
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

from future.utils import iteritems

import re
import bs4
import warnings

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session

_logger = session.Logger(name = 'inputs.genecards')
_log = _logger._log


_respace = re.compile(r'\s+')
_summary_sources = {
    'Gene Wiki': 'GeneWiki',
    'UniProtKB/Swiss-Prot': 'UniProt',
    'GeneCards': 'GeneCards',
}


def card(gene):
    
    url = urls.urls['genecards']['url'] % gene
    
    c = curl.Curl(url, silent = True, large = True)
    
    result = {}
    
    if c.status not in {0, 200}:
        
        _log('Failed to retrieve gene card for ID `%s`.' % gene)
        
        return result
    
    with warnings.catch_warnings():
        
        warnings.simplefilter('ignore')
        soup = bs4.BeautifulSoup(c.fileobj)
    
    summaries = soup.select_one('section#summaries')
    
    for summary in summaries.select('div.gc-subsection'):
        
        title = summary.select_one('h3').text.strip('\r\n ')
        
        if title[:7] in {'No data', 'Additio'}:
            
            continue
        
        content = _respace.sub(
            ' ',
            ' '.join(
                par.text
                for par in summary.select(':not(:nth-child(1))')
            )
        ).strip('\n\r ')
        
        for gc_name, name in iteritems(_summary_sources):
            
            if title.startswith(gc_name):
                
                title = name
                break
        
        if content:
            
            result[title] = content
    
    return result
