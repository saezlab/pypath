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
# Code from Janus Troelsen, see here:
# http://stackoverflow.com/a/13781114/854988
#

import colorsys
import itertools
import math
from fractions import Fraction

def zenos_dichotomy():
    """
    http://en.wikipedia.org/wiki/1/2_%2B_1/4_%2B_1/8_%2B_1/16_%2B_%C2%B7_%C2%B7_%C2%B7
    """
    for k in itertools.count():
        yield Fraction(1,2**k)

def getfracs():
    """
    [Fraction(0, 1), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4), Fraction(1, 8), Fraction(3, 8), Fraction(5, 8), Fraction(7, 8), Fraction(1, 16), Fraction(3, 16), ...]
    [0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 0.0625, 0.1875, ...]
    """
    yield 0
    for k in zenos_dichotomy():
        i = k.denominator # [1,2,4,8,16,...]
        for j in range(1,i,2):
            yield Fraction(j,i)

bias = lambda x: (math.sqrt(x/3)/Fraction(2,3)+Fraction(1,3))/Fraction(6,5) # can be used for the v in hsv to map linear values 0..1 to something that looks equidistant

def genhsv(h):
    for s in [Fraction(6,10)]: # optionally use range
        for v in [Fraction(8,10),Fraction(5,10)]: # could use range too
            yield (h, s, bias(v)) # use bias for v here if you use range

def dec2hex(d):
    h = '#'
    for i in d:
        h += "%02x" % math.floor(i*255.99)
    return h.upper()

def hex2dec(h):
    return (int(h[1:3], 16),int(h[3:5], 16),int(h[4:6], 16))

genrgb = lambda x: colorsys.hsv_to_rgb(*x)

flatten = itertools.chain.from_iterable

gethsvs = lambda x: flatten(map(genhsv,list(itertools.islice(getfracs(),x))))

getrgbs = lambda x: map(genrgb, gethsvs(x))

gethexrgbs = lambda x: map(dec2hex, getrgbs(x))

def colormix(colA,colB):
    rgbA = hex2dec(colA)
    rgbB = hex2dec(colB)
    rgbAB = ((rgbA[0]+rgbB[0])/510.0, (rgbA[1]+rgbB[1])/510.0, (rgbA[2]+rgbB[2])/510.0)
    return dec2hex(rgbAB)
