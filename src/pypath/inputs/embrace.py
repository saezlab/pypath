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

import os
import itertools
import collections

import pypath.inputs.common as inputs_common
import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.utils.mapping as mapping
import pypath.share.common as common
import pypath.share.session as session

_logger = session.Logger(name = 'inputs.embrace')
_log = _logger._log


def embrace_raw():
    """
    Returns Supplementary Table S11 from 10.1016/j.isci.2019.10.026
    (Sheikh et al. 2019) as a list of lists.
    """
    
    url = urls.urls['embrace']['url']
    c_nocall = curl.Curl(
        url,
        call = False,
        setup = False,
        process = False,
        silent = True,
    )
    c_nocall.get_cache_file_name()
    path = c_nocall.cache_file_name
    
    init_url = urls.urls['embrace']['article']
    req_headers = []
    
    if not os.path.exists(path):
        
        cookies = {}
        
        for step in range(3):
            
            c_init = curl.Curl(
                init_url,
                silent = True,
                large = True,
                cache = False,
                follow = False,
                req_headers = req_headers + ['user-agent: curl/7.69.1'],
                bypass_url_encoding = True,
                retries = 1,
                empty_attempt_again = False,
            )
            
            new_cookies = dict(
                tuple(
                    h.decode().split(':')[1].\
                    split(';')[0].\
                    strip().split('=', maxsplit = 1)
                )
                for h in c_init.resp_headers
                if h.lower().startswith(b'set-cookie')
            )
            cookies.update(new_cookies)
            _ = cookies.pop('__cflb', None)
            
            for h in c_init.resp_headers:
                
                if h.lower().startswith(b'location'):
                    
                    init_url = h.decode().split(':', maxsplit = 1)[1].strip()
            
            req_headers = (
                [
                    'Cookie: %s' % (
                        '; '.join(
                            '%s=%s' % cookie
                            for cookie in iteritems(cookies)
                        )
                    )
                ]
                    if cookies else
                []
            )
            
            _log(
                'HTTP %u; location: `%s`, cookies: `%s`.' % (
                    c_init.status,
                    init_url,
                    req_headers[0] if req_headers else '',
                )
            )
            
            if c_init.status != 302:
                
                break
    
    c_table = curl.Curl(
        url,
        silent = False,
        large = True,
        empty_attempt_again = False,
        req_headers = req_headers + ['user-agent: curl/7.69.1'],
    )
    path = c_table.cache_file_name
    c_table.fileobj.close()
    
    content = inputs_common.read_xls(path)
    
    EmbraceRawRecord = collections.namedtuple(
        'EmbraceRawRecord',
        content[0]
    )
    
    return [
        EmbraceRawRecord(*(line[:2] + [int(float(n)) for n in line[2:]]))
        for line in
        content[1:]
    ]


def embrace_interactions(organism = 9606):
    
    
    EmbraceInteraction = collections.namedtuple(
        'ItalkInteraction',
        [
            'ligand',
            'receptor',
            'classification',
        ]
    )
    
    
    rdata = italk_raw()
    result = []
    
    for row in rdata.itertuples():
        
        if (
            not isinstance(row[2], common.basestring) or
            not isinstance(row[4], common.basestring)
        ):
            
            continue
        
        ligands = mapping.map_name(row[2], 'genesymbol', 'uniprot')
        receptors = mapping.map_name(row[4], 'genesymbol', 'uniprot')
        cls = row[6]
        
        for ligand, receptor in itertools.product(ligands, receptors):
            
            result.append(
                ItalkInteraction(
                    ligand = ligand,
                    receptor = receptor,
                    classification = cls,
                )
            )
    
    return result


def italk_annotations():
    
    
    ItalkAnnotation = collections.namedtuple(
        'ItalkAnnotation',
        [
            'mainclass',
            'subclass',
        ]
    )
    
    
    rdata = italk_raw()
    result = collections.defaultdict(set)
    
    for row in rdata.itertuples():
        
        ligands = (
            mapping.map_name(row[2], 'genesymbol', 'uniprot')
                if isinstance(row[2], common.basestring) else
            ()
        )
        receptors = (
            mapping.map_name(row[4], 'genesymbol', 'uniprot')
                if isinstance(row[4], common.basestring) else
            ()
        )
        subclass = row[6]
        
        for mainclass, uniprot in itertools.chain(
            itertools.zip_longest((), ligands, fillvalue = 'ligand'),
            itertools.zip_longest((), receptors, fillvalue = 'receptor'),
        ):
            
            result[uniprot].add(
                ItalkAnnotation(
                    mainclass = mainclass,
                    subclass = subclass,
                )
            )
    
    return dict(result)
