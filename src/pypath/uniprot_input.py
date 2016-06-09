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

from future.utils import iteritems

import re

import pypath.urls as urls
import pypath.curl as curl
import pypath.seq as se

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

def get_isoforms(organism = 'Homo sapiens'):
    reorg = re.compile(r'OS=([A-Z][a-z]+\s[a-z]+)')
    result = {}
    url = urls.urls['unip_iso']['url']
    c = curl.Curl(url, silent = False)
    data = c.result
    data = read_fasta(data)
    for header, seq in iteritems(data):
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
