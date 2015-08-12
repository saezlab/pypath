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
import math
import random
import textwrap

__all__ = ['ROOT', 'aacodes', 'aaletters', 'simpleTypes', 'uniqList', 'addToList', 
           'gen_session_id', 'sorensen_index', 'console', 'wcl', 'flatList', 
           'charTypes', 'delEmpty', '__version__', 'get_args', 
           'something', 'rotate']

# get the location
ROOT = os.path.abspath(os.path.dirname(__file__))

def _get_version():
    with open(os.path.join(ROOT, '__version__'), 'r') as v:
        return tuple([int(i) for i in v.read().strip().split('.')])

_MAJOR, _MINOR, _MICRO = _get_version()
__version__ = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
__release__ = '%d.%d' % (_MAJOR, _MINOR)


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
        try:
            keys[e] = 1
        except:
            print e
            print seq
            print keys
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

def something(anything):
    return not (anything is None or \
        (type(anything) in [list, set, dict, str, unicode] \
            and len(anything) == 0))

def gen_session_id(length=5):
    abc = '0123456789abcdefghijklmnopqrstuvwxyz'
    return ''.join(random.choice(abc) for i in range(length))

def simpson_index(a, b):
    a = set(a)
    b = set(b)
    ab = a & b
    return float(len(ab)) / float(min(len(a),len(b)))

def sorensen_index(a, b):
    a = set(a)
    b = set(b)
    ab = a & b
    return float(len(ab)) / float(len(a) + len(b))

def jaccard_index(a, b):
    a = set(a)
    b = set(b)
    ab = a & b
    return float(len(ab)) / float(len(a | b))

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

def get_args(loc_dict, remove = set([])):
    if type(remove) not in [set, list]: remove = set([remove])
    if type(remove) is list: remove = set(remove)
    remove.add('self')
    remove.add('kwargs')
    args = dict((k, v) for k, v in loc_dict.iteritems() if k not in remove)
    if 'kwargs' in loc_dict: args = dict(args.items() + loc_dict['kwargs'].items())
    return args

def rotate(point, angle, center = (0.0, 0.0)):
    """
    from http://stackoverflow.com/a/20024348/854988
    Rotates a point around center. Angle is in degrees.
    Rotation is counter-clockwise
    """
    angle = math.radians(angle)
    temp_point = point[0]-center[0] , point[1]-center[1]
    temp_point = ( temp_point[0]*math.cos(angle)-temp_point[1]*math.sin(angle), \
        temp_point[0]*math.sin(angle) + temp_point[1]*math.cos(angle))
    temp_point = temp_point[0]+center[0] , temp_point[1]+center[1]
    return temp_point