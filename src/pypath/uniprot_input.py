#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright (c) 2014-2016 - EMBL-EBI
#
#  File author(s): Dénes Türei (denes@ebi.ac.uk)
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

import pypath.urls as urls
import pypath.curl as curl

def all_uniprots(organism = 9606, swissprot = None):
    rev = '' if swissprot is None else ' AND reviewed:%s'%swissprot
    url = urls.urls['uniprot_basic']['url']
    post = {'query': 'organism:%s%s' % (str(organism), rev), 
        'format': 'tab', 'columns': 'id'}
    c = curl.Curl(url, post = post, silent = False)
    data = c.result
    return list(filter(lambda x:
        len(x) > 0,
        map(lambda l:
            l.strip(),
            data.split('\n')[1:]
        )
    ))

def swissprot_seq(organism = 9606, isoforms = False):
    taxids = {
        9606: 'Homo sapiens'
    }
    result = {}
    url = urls.urls['uniprot_basic']['url']
    post = {'query': 'organism:%s AND reviewed:yes'%str(organism), 
        'format': 'tab', 'columns': 'id,sequence'}
    c = curl.Curl(url, post = post, silent = False)
    data = c.result
    data = data.split('\n')
    del data[0]
    for l in data:
        l = l.strip().split('\t')
        if len(l) == 2:
            result[l[0]] = se.Seq(l[0], l[1])
    if isoforms:
        data = get_isoforms()
        for unip, isoforms in iteritems(data):
            for isof, seq in iteritems(isoforms):
                if unip in result:
                    result[unip].add_seq(seq, isof)
    return result
