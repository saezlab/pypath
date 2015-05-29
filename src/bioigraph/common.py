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

import os
import sys
import random
import textwrap

__all__ = ['ROOT', 'aacodes', 'aaletters', 'simpleTypes', 'uniqList', 'addToList', 
           'gen_session_id', 'sorensen_index', 'console', 'wcl', 'flatList', 
           'charTypes', 'delEmpty']

# get the location
ROOT = os.path.abspath(os.path.dirname(__file__))

default_name_type = {
    "protein": "uniprot", 
    "mirna": "mirbase", 
    "lncrna": "lncrnaname",
    "drug": "pubchem"
}

aacodes = {
    'A': 'ALA',
    'G': 'GLY',
    'M': 'MET',
    'W': 'TRP',
    'Y': 'TYR',
    'R': 'ARG',
    'C': 'CYS',
    'S': 'SER',
    'V': 'VAL',
    'I': 'ILE',
    'L': 'LEU',
    'P': 'PRO',
    'H': 'HIS',
    'T': 'THR',
    'E': 'GLU',
    'Q': 'GLN',
    'Z': 'GLX',
    'D': 'ASP',
    'N': 'ASN',
    'B': 'ASX',
    'K': 'LYS',
    'F': 'PHE',
    'J': 'XLE',
    'X': 'XAA'
}

aaletters = dict(zip(aacodes.values(),aacodes.keys()))

simpleTypes = [int, float, str, unicode]

charTypes = [str, unicode]

def uniqList(seq):
   # Not order preserving
   # from http://www.peterbe.com/plog/uniqifiers-benchmark
   keys = {}
   for e in seq:
       keys[e] = 1
   return keys.keys()

def flatList(lst):
    return [it for sl in lst for it in sl]

def delEmpty(lst):
    return [i for i in lst if len(i) > 0]

def uniqOrdList(seq, idfun = None): 
   # Order preserving
   # from http://www.peterbe.com/plog/uniqifiers-benchmark
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def addToList(lst, toadd):
    if isinstance(toadd, list):
        lst += toadd
    else:
        lst.append(toadd)
    if None in lst:
        lst.remove(None)
    return uniqList(lst)

def gen_session_id(length=5):
    abc = '0123456789abcdefghijklmnopqrstuvwxyz'
    return ''.join(random.choice(abc) for i in range(length))

def sorensen_index(a,b):
    a = set(a)
    b = set(b)
    ab = a.intersection(b)
    return float(2 * len(ab)) / float(len(a) + len(b))

def console(message):
    message = '\n\t'.join(textwrap.wrap(message,80))
    sys.stdout.write(('\n\t'+message).ljust(80))
    sys.stdout.write('\n')
    sys.stdout.flush()

def wcl(f):
    toClose = type(f) is file
    f = f if type(f) is file else open(f, 'r')
    for i, l in enumerate(f):
        pass
    if toClose:
        f.seek(0)
    else:
        f.close()
    return i + 1