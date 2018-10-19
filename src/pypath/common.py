#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2018
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://www.ebi.ac.uk/~denes
#

from future.utils import iteritems
from past.builtins import xrange, range, reduce

import os
import sys
import math
import random
import textwrap
import hashlib

__all__ = ['ROOT', 'aacodes', 'aaletters', 'simpleTypes', 'numTypes',
           'uniqList', 'addToList', 'addToSet', 'gen_session_id',
           'sorensen_index', 'simpson_index', 'simpson_index_counts',
           'jaccard_index', 'console', 'wcl', 'flatList', 'charTypes',
           'delEmpty', 'get_args', 'something', 'rotate', 'cleanDict',
           'igraph_graphics_attrs', 'md5', 'mod_keywords', 'Namespace', 'fun',
           'taxids', 'taxa', 'phosphoelm_taxids', 'dbptm_taxids',]

# get the location
ROOT = os.path.abspath(os.path.dirname(__file__))
DATA = os.path.join(ROOT, 'data')

try:
    basestring

except NameError:
    basestring = str

default_name_type = {"protein": "uniprot",
                     "mirna": "mirbase",
                     "lncrna": "lncrnaname",
                     "drug": "pubchem"}

aacodes = {'A': 'ALA', 'G': 'GLY', 'M': 'MET', 'W': 'TRP', 'Y': 'TYR',
           'R': 'ARG', 'C': 'CYS', 'S': 'SER', 'V': 'VAL', 'I': 'ILE',
           'L': 'LEU', 'P': 'PRO', 'H': 'HIS', 'T': 'THR', 'E': 'GLU',
           'Q': 'GLN', 'Z': 'GLX', 'D': 'ASP', 'N': 'ASN', 'B': 'ASX',
           'K': 'LYS', 'F': 'PHE', 'J': 'XLE', 'X': 'XAA'}

aanames = {'alanine': 'A',
           'arginine': 'R',
           'asparagine': 'N',
           'aspartic acid': 'D',
           'cysteine': 'C',
           'glutamine': 'Q',
           'glutamic acid': 'E',
           'glycine': 'G',
           'histidine': 'H',
           'isoleucine': 'I',
           'leucine': 'L',
           'lysine': 'K',
           'methionine': 'M',
           'phenylalanine': 'F',
           'proline': 'P',
           'serine': 'S',
           'threonine': 'T',
           'tryptophan': 'W',
           'tyrosine': 'Y',
           'valine': 'V'}

mod_keywords = {'Reactome': [('phosphopantetheinylation', ['phosphopantet']),
                             ('phosphorylation', ['phospho']),
                             ('acetylneuraminylation', ['acetylneuraminyl']),
                             ('acetylation', ['acetyl']),
                             ('farnesylation', ['farnesyl']),
                             ('palmitoylation', ['palmito']),
                             ('methylation', ['methyl']),
                             ('tetradecanoylation', ['tetradecanoyl']),
                             ('decanoylation', ['decanoyl']),
                             ('palmitoleylation', ['palmytoleil']),
                             ('formylation', ['formyl']),
                             ('ubiquitination', ['ubiquitin']),
                             ('galactosylation', ['galactos']),
                             ('glutamylation', ['glutamyl']),
                             ('fucosylation', ['fucosyl']),
                             ('myristoylation', ['myristoyl']),
                             ('carboxylation', ['carboxyl']),
                             ('biotinylation', ['biotinyl']),
                             ('glycosylation', ['glycosyl']),
                             ('octanoylation', ['octanoyl']),
                             ('glycylation', ['glycyl']),
                             ('hydroxylation', ['hydroxy']),
                             ('sulfhydration', ['persulfid']),
                             ('thiolation', ['thio']),
                             ('amidation', ['amide']),
                             ('selenation', ['seleno']),
                             ('glucosylation', ['glucosyl']),
                             ('neddylation', ['neddyl']),
                             ('sumoylation', ['sumoyl']),
                             ('prenylation', ['prenyl'])],
                'ACSN': [('phosphorylation', ['phospho']),
                         ('glycylation', ['glycyl']),
                         ('ubiquitination', ['ubiquitin']),
                         ('acetylation', ['acetyl']),
                         ('myristoylation', ['myristoyl']),
                         ('prenylation', ['prenyl']),
                         ('hydroxylation', ['hydroxy'])],
                'WikiPathways': [],
                'NetPath': [('phosphorylation', ['phospho']),
                            ('glycylation', ['glycyl']),
                            ('ubiquitination', ['ubiquitin']),
                            ('acetylation', ['acetyl'])],
                'PANTHER': [('phosphorylation', ['phospho']),
                            ('acetylation', ['acetyl'])],
                'NCI-PID': [('phosphorylation', ['phospho']),
                            ('methylation', ['methyl']),
                            ('farnesylation', ['farnesyl']),
                            ('palmitoylation', ['palmito']),
                            ('myristoylation', ['myristoyl']),
                            ('glycylation', ['glycyl']),
                            ('ubiquitination', ['ubiquitin']),
                            ('acetylation', ['acetyl']),
                            ('glycosylation', ['glycosyl']),
                            ('geranylation', ['geranyl']),
                            ('hydroxylation', ['hydroxy'])],
                'KEGG': [('phosphorylation', ['hospho']),
                         ('methylation', ['methyl']),
                         ('ubiquitination', ['ubiquitin']),
                         ('acetylation', ['acetyl']),
                         ('hydroxylation', ['hydroxy']),
                         ('carboxyethylation', ['carboxyethyl']),
                         ('ribosylation', ['ribosyl']),
                         ('nitrosylation', ['nitrosyl']),
                         ('sulfoylation', ['ulfo']),
                         ('biotinylation', ['biotinyl']),
                         ('malonylation', ['malonyl']),
                         ('glutarylation', ['lutaryl'])]}

# FIXME:
# Running from file does not raise error, but if you run the lines in
# the interpreter, raises TypeError: argument of type 'module' is not iterable
# TODO: __builtins__ -> dir(__builtins__)
if 'long' not in __builtins__:
    long = int

if 'unicode' not in __builtins__:
    unicode = str

# Inverted `aacodes` dict
aaletters = dict(zip(aacodes.values(), aacodes.keys()))

# Type definitions
simpleTypes = set([int, long, float, str, unicode, bytes, bool, type(None)])
numTypes = set([int, long, float])
charTypes = set([str, unicode, bytes])


# From http://www.peterbe.com/plog/uniqifiers-benchmark
def uniqList(seq):
    """
    Takes any iterable and returns a list of unique elements on it. If
    the argument is a dictionary, returns a list of unique keys.
    **NOTE:** Does not preserve the order of the elements.

    * Arguments:
        - *seq*: Sequence to be processed, can be any iterable type.

    * Returns:
        - [list]: List of unique elements in the sequence *seq*.

    * Examples:
        >>> uniqList('aba')
        ['a', 'b']
        >>> uniqList([0, 1, 2, 1, 0])
        [0, 1, 2]
    """

    return list({}.fromkeys(seq).keys())


def uniqList1(seq): # XXX: Not used
    """
    Not order preserving
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    """

    return list(set(seq))


def uniqList2(seq): # XXX: Not used
    """
    Not order preserving
    From http://www.peterbe.com/plog/uniqifiers-benchmark
    """

    keys = {}

    for e in seq:
        keys[e] = 1

    return list(keys.keys())


def flatList(lst):
    """
    Coerces the elements of a list of iterables into a single list.

    * Arguments:
        - *lst* [list]: List to be flattened. Its elements can be also
          lists or any other iterable.

    * Returns:
        - [list]: Flattened list of *lst*.

    * Examples:
        >>> flatList([(0, 1), (1, 1), (2, 1)])
        [0, 1, 1, 1, 2, 1]
        >> flatList(['abc', 'def'])
        ['a', 'b', 'c', 'd', 'e', 'f']
    """

    return [it for sl in lst for it in sl]


def delEmpty(lst): # XXX: Only used in main.py line: 1278
    """
    Removes empty entries of a list.

    * Arguments:
        - *lst* [list]: List from which empty elements will be removed.

    * Returns:
        - [list]: Copy of *lst* without elements whose length was zero.

    * Examples:
        >>> delEmpty(['a', '', 'b', 'c'])
        ['a', 'b', 'c']
    """

    return [i for i in lst if len(i) > 0]


# Order preserving
# From http://www.peterbe.com/plog/uniqifiers-benchmark
def uniqOrdList(seq, idfun=None): # XXX: Only used in plot.py line: 510
    """
    Returns a copy of *seq* without repeated elements. Preserves the
    order.

    * Arguments:
        - *seq* [list]: Or any other iterable type. The sequence from
          which repeated elements are to be removed.
        - *idfun* [function]: Optional, ``None`` by default. Identifier
          function, for each entry of *seq*, returns a identifier of
          that entry from which uniqueness is determined. Default
          behavior is f(x) = x. See examples below.

    * Returns:
        - [list]: Copy of *seq* without the repeated elements (according
          to *idfun*). If any element is repeated, the first instance is
          kept.

    * Examples:
        >>> uniqOrdList([0, 1, 2, 1, 5])
        [0, 1, 2, 5]
        >>> uniqOrdList('abracadabra')
        ['a', 'b', 'r', 'c', 'd']
        >>> def f(x):
        ...     if x > 0:
        ...             return 0
        ...     else:
        ...             return 1
        >>> uniqOrdList([-32, -42, 1, 15, -12], idfun=f)
        [-32, 1]
        >>> def g(x): # Given a file name, return it without extension
        ...    return x.split('.')[0]
        >>> uniqOrdList(['a.png', 'a.txt', 'b.pdf'], idfun=g)
        ['a.png', 'b.pdf']
    """

    if idfun is None:

        def idfun(x):
            return x

    seen = {}
    result = []

    for item in seq:
        marker = idfun(item)

        if marker in seen:
            continue

        seen[marker] = 1
        result.append(item)

    return result


def addToList(lst, toadd):
    """
    Appends *toadd* to *lst*. Function differs from ``lst.append`` since
    is capable to handle different data types. This is, if *lst* is not
    a list, it will be converted to one. Similarly, if *toadd* is not a
    list, it will be converted and added. If *toadd* is or contains
    ``None``, these will be ommited. The returned list will only contain
    unique elements and does not necessarily preserve order.

    * Arguments:
        - *lst*: List or any other type (will be converted into a list).
          Original sequence to which *toadd* will be appended.
        - *toadd*: Element(s) to be added into *lst*.

    * Returns:
        - [list]: Contains the unique element(s) from the union of *lst*
          and *toadd*. **NOTE:** Makes use of ``common.uniqList``, does
          not preserve order of elements.

    * Examples:
        >>> addToList('ab', 'cd')
        ['ab', 'cd']
        >>> addToList('ab', ['cd', None, 'ab', 'ef'])
        ['ab', 'ef', 'cd']
        >>> addToList((0, 1, 2), 4)
        [0, 1, 2, 4]
    """

    if type(lst) is not list:

        if type(lst) in simpleTypes:
            lst = [lst]

        else:
            lst = list(lst)

    if toadd is None:
        return lst

    if type(toadd) in simpleTypes:
        lst.append(toadd)

    else:

        if type(toadd) is not list:
            toadd = list(toadd)

        lst.extend(toadd)

    if None in lst:
        lst.remove(None)

    return uniqList(lst)


def addToSet(st, toadd):
    """
    Appends *toadd* to *st*. Function is capable to handle different
    input data types. This is, if *toadd* is a list, it will be
    converted to a set and added.

    * Arguments:
        - *st* [set]: Original set to which *toadd* will be added.
        - *toadd*: Element(s) to be added into *st*.

    * Returns:
        - [set]: Contains the element(s) from the union of *st* and
          *toadd*.

    * Examples:
        >>> st = set([0, 1, 2])
        >>> addToSet(st, 3)
        set([0, 1, 2, 3])
        >>> addToSet(st, [4, 2, 5])
        set([0, 1, 2, 4, 5])
    """

    if type(toadd) in simpleTypes:
        st.add(toadd)

    if type(toadd) is list:
        toadd = set(toadd)

    if type(toadd) is set:
        st.update(toadd)

    return st


def something(anything):
    """
    Checks if *anything* is ``None`` or empty.

    * Arguments:
        - *anything*: Self-explanatory.

    * Returns:
        - [bool]: ``False`` if *anyhting* is ``None`` or any empty data
          type.

    * Examples:
        >>> something(None)
        False
        >>> something(123)
        True
        >>> something('Hello world!')
        True
        >>> something('')
        False
        >>> something([])
        False
    """

    return not (anything is None
                or (type(anything) in [list, set, dict, str, unicode]
                    and len(anything) == 0))


def gen_session_id(length=5):
    """
    Generates a random alphanumeric string.

    * Arguments:
        - *length* [int]: Optional, ``5`` by default. Specifies the
          length of the random string.

    * Returns:
        - [str]: Random alphanumeric string of the specified length.
    """

    abc = '0123456789abcdefghijklmnopqrstuvwxyz'

    return ''.join(random.choice(abc) for i in xrange(length))


# XXX: Are you sure this is the way to compute Simpson's index?
def simpson_index(a, b):
    """
    Given two sets *a* and *b*, returns the Simpson index.

    * Arguments:
        - *a* [set]: Or any iterable type (will be converted to set).
        - *b* [set]: Or any iterable type (will be converted to set).

    * Returns:
        - [float]: The Simpson index between *a* and *b*.
    """

    a = set(a)
    b = set(b)
    ab = a & b

    return float(len(ab)) / float(min(len(a), len(b)))


# XXX: Related to comment above, what is this exactly?
def simpson_index_counts(a, b, c):
    """
    * Arguments:
        - *a* []: .
        - *b* []: .
        - *c* []: .

    * Returns:
        - [float]: .
    """

    return float(c) / float(min(a, b)) if min(a, b) > 0 else 0.0


def sorensen_index(a, b):
    """
    Computes the Sorensen-Dice coefficient between two sets *a* and *b*.

    * Arguments:
        - *a* [set]: Or any iterable type (will be converted to set).
        - *b* [set]: Or any iterable type (will be converted to set).

    * Returns:
        - [float]: The Sorensen-Dice coefficient between *a* and *b*.
    """

    a = set(a)
    b = set(b)
    ab = a & b

    return float(len(ab)) / float(len(a) + len(b))


def jaccard_index(a, b):
    """
    Computes the Jaccard index between two sets *a* and *b*.

    * Arguments:
        - *a* [set]: Or any iterable type (will be converted to set).
        - *b* [set]: Or any iterable type (will be converted to set).

    * Returns:
        - [float]: The Jaccard index between *a* and *b*.
    """

    a = set(a)
    b = set(b)
    ab = a & b

    return float(len(ab)) / float(len(a | b))


def console(message):
    """
    Prints a *message* to the standard output (e.g. terminal) formatted
    to 80 characters per line plus first-level indentation.

    * Arguments:
        - *message* [str]: The message to be printed.
    """

    message = '\n\t'.join(textwrap.wrap(message, 80))
    sys.stdout.write(('\n\t' + message).ljust(80))
    sys.stdout.write('\n')
    sys.stdout.flush()


def wcl(f): # XXX: Not used (another function w/ same name defined in curl.py)
    """

    """

    toClose = type(f) is file
    f = f if type(f) is file else open(f, 'r')

    for i, l in enumerate(f):
        pass

    if toClose:
        f.seek(0)

    else:
        f.close()

    return i + 1


# XXX: Not very clear to me the purpose of this function.
def get_args(loc_dict, remove=set([])):
    """
    Given a dictionary of local variables, returns a copy of it without
    ``self``, ``kwargs`` (in the scope of a ``class``) plus any other
    specified in the keyword argument *remove*.

    * Arguments:
        - *loc_dict* [dict]: Dictionary containing the local variables
          (e.g. a call to ``locals()`` in a given scope).
        - *remove* [set]: Optional, ``set([])`` by default. Can also be
          a list. Contains the keys of the elements in *loc_dict* that
          will be removed.

    * Returns:
        - [dict]: Copy of *loc_dict* without ``self``, ``kwargs`` and
          any other element specified in *remove*.
    """

    if type(remove) not in [set, list]:
        remove = set([remove])

    if type(remove) is list:
        remove = set(remove)

    remove.add('self')
    remove.add('kwargs')
    args = dict((k, v) for k, v in iteritems(loc_dict) if k not in remove)

    if 'kwargs' in loc_dict:
        args = dict(args.items() + loc_dict['kwargs'].items())

    return args


# From http://stackoverflow.com/a/20024348/854988
def rotate(point, angle, center=(0.0, 0.0)): # XXX: Not used? Wrote the docs before checking XD
    """
    Rotates a given *point* around a *center* according to the specified
    *angle* (in degrees) in a two-dimensional space. The rotation is
    counter-clockwise.

    * Arguments:
        - *point* [tuple]: Or list. Contains the two coordinates of the
          point to be rotated.
        - *angle* [float]: Angle (in degrees) from which the point
          will be rotated with respect to *center* (counter-clockwise).
        - *center* [tuple]: Optional, ``(0.0, 0.0)`` by default. Can
          also be a list. Determines the two coordinates of the center
          relative to which the point has to be rotated.

    * Returns:
        - [tuple]: Pair of coordinates of the rotated point.
    """

    angle = math.radians(angle)
    temp_point = point[0] - center[0], point[1] - center[1]
    temp_point = (temp_point[0] * math.cos(angle)
                  - temp_point[1] * math.sin(angle),
                  temp_point[0] * math.sin(angle)
                  + temp_point[1] * math.cos(angle))
    temp_point = temp_point[0] + center[0], temp_point[1] + center[1]

    return temp_point


def cleanDict(dct):
    """
    Removes ``None`` values from  a dictionary *dct* and casts all other
    values to strings.

    * Arguments:
        - *dct* [dict]: Dictionary to be cleaned from ``None`` values.

    * Returns:
        - [dict]: Copy of *dct* without ``None`` value entries and all
          other values formatted to strings.
    """

    toDel = []

    for k, v in iteritems(dct):

        if v is None:
            toDel.append(k)

        else:
            dct[k] = str(v)

    for k in toDel:
        del dct[k]

    return dct


def md5(value):
    """
    Computes the sum of MD5 hash of a given string *value*.

    * Arguments:
        - *value* [str]: Or any other type (will be converted to
          string). Value for which the MD5 sum will be computed. Must
          follow ASCII encoding.

    * Return:
        - [str]: Hash resulting from the MD5 sum of the *value* string.
    """

    try:
        string = str(value).encode('ascii')

    except: # XXX: Bad practice to catch any exception
        string = str(value).encode('ascii') # XXX: Same as the try statement!?

    return hashlib.md5(string).hexdigest()


# XXX: Shouldn't we keep all functions and variables separated
#      (together among them)?
igraph_graphics_attrs = {'vertex': ['size', ' color', 'frame_color',
                                    'frame_width', 'shape', 'label',
                                    'label_dist', 'label_color', 'label_size',
                                    'label_angle'],
                         'edge': ['curved', 'color', 'width', 'arrow_size',
                                  'arrow_width']}

##############################################################################
#                            |   DU BIST HIER   |                            #
#                            V                  V                            #
##############################################################################


def merge_dicts(d1, d2):
    """
    Merges dicts recursively
    """

    for k2, v2 in iteritems(d2):
        t = type(v2)

        if k2 not in d1:
            d1[k2] = v2

        elif t is dict:
            d1[k2] = merge_dicts(d1[k2], v2)

        elif t is list:
            d1[k2].extend(v2)

        elif t is set:
            d1[k2].update(v2)

    return d1


def dict_set_path(d, path):
    """
    In dict of dicts ``d`` looks up the keys following ``path``,
    creates new subdicts and keys if those do not exist yet,
    and sets/merges the leaf element according to simple heuristic.
    """

    val = path[-1]
    key = path[-2]
    subd = d

    for k in path[:-2]:

        if type(subd) is dict:

            if k in subd:
                subd = subd[k]

            else:
                subd[k] = {}
                subd = subd[k]

        else:
            return d

    if key not in subd:
        subd[key] = val

    elif type(val) is dict and type(subd[key]) is dict:
        subd[key].update(val)

    elif type(subd[key]) is list:

        if type(val) is list:
            subd[key].extend(val)

        else:
            subd[key].append(val)

    elif type(subd[key]) is set:

        if type(val) is set:
            subd[key].update(val)

        else:
            subd[key].add(val)

    return d


def dict_diff(d1, d2):
    """

    """

    ldiff = {}
    rdiff = {}
    keys = set(d1.keys()) & set(d2.keys())

    for k in keys:

        if type(d1[k]) is dict and type(d2[k]) is dict:
            ldiff[k], rdiff[k] = dict_diff(d1[k], d2[k])

        elif type(d1[k]) is set and type(d2[k]) is set:
            ldiff[k], rdiff[k] = (d1[k] - d2[k], d2[k] - d1[k])

    return ldiff, rdiff


def dict_sym_diff(d1, d2):
    """

    """

    diff = {}
    keys = set(d1.keys()) & set(d2.keys())

    for k in keys:

        if type(d1[k]) is dict and type(d2[k]) is dict:
            diff[k] = dict_sym_diff(d1[k], d2[k])

        elif type(d1[k]) is set and type(d2[k]) is set:
            diff[k] = d1[k] ^ d2[k]

    return diff


def swap_dict(d):
    """
    Interchanges the keys and values of a dict.
    Results dict of sets even if the values are strings or ints
    and are unique.
    """

    _d = {}

    for key, vals in iteritems(d):
        vals = [vals] if type(vals) in simpleTypes else vals

        for val in vals:

            if val not in _d:
                _d[val] = set([])

            _d[val].add(key)

    if all(len(v) <= 1 for v in _d.values()):
        _d = dict((k, list(v)[0]) for k, v in iteritems(_d) if len(v))

    return _d


def swap_dict_simple(d):
    """
    Interchanges the keys and values of a dict.
    Assumes the values are unique and hashable,
    otherwise overwrites duplicates or raises error.
    """

    return dict((v, k) for k, v in iteritems(d))


def join_dicts(d1, d2, _from = 'keys', to = 'values'):
    """

    """

    result = {}

    if to == 'keys':
        d2 = swap_dict(d2)

    for key1, val1 in iteritems(d1):
        sources = ([key1] if _from == 'keys'
                          else [val1] if type(val1) in simpleTypes
                                      else val1)
        meds = ([key1] if _from == 'values'
                       else [val1] if type(val1) in simpleTypes
                                   else val1)
        targets = set([])

        for med in meds:

            if med in d2:

                if type(targets) is list:
                    targets.append(d2[med])

                elif type(d2[med]) in simpleTypes:
                    targets.add(d2[med])

                elif type(d2[med]) is list:
                    targets.update(set(d2[med]))

                elif type(d2[med]) is set:
                    targets.update(d2[med])

                elif d2[med].__hash__ is not None:
                    targets.add(d2[med])

                else:
                    targets = list(targets)
                    targets.append(d2[med])

        for source in sources:

            if type(targets) is list:

                if source not in result:
                    result[source] = []

                result[source].extend(targets)

            elif type(targets) is set:

                if source not in result:
                    result[source] = set([])

                result[source].update(targets)

    if all(len(x) <= 1 for x in result.values()):
        result = dict((k, list(v)[0]) for k, v in iteritems(result) if len(v))

    return result


class Namespace(object): # XXX: WHY?
    pass


def fun(): # XXX: Best name for a function
    """

    """

    print(__name__)
    print(__name__ in globals())

    for n in __name__.split('.'):
        print(n, n in globals())

    return __name__


taxids = {9606: 'human',
          10090: 'mouse',
          10116: 'rat',
          9031: 'chicken',
          9913: 'cow',
          9986: 'rabbit',
          9940: 'sheep',
          10141: 'guinea pig',
          10036: 'hamster',
          7227: 'fruit fly',
          9615: 'dog',
          9823: 'pig',
          8355: 'frog',
          9091: 'quail',
          9796: 'horse',
          9925: 'goat',
          89462: 'water buffalo',
          9598: 'monkey',
          9103: 'turkey',
          9685: 'cat',
          7604: 'starfish',
          1213717: 'torpedo',
          9669: 'ferret',
          8839: 'duck'}

taxa = dict(map(lambda i: (i[1], i[0]), taxids.items()))

phosphoelm_taxids = {9606: 'Homo sapiens',
                     10090: 'Mus musculus',
                     9913: 'Bos taurus',
                     9986: 'Oryctolagus cuniculus',
                     9615: 'Canis familiaris',
                     10029: 'Cricetulus griseus',
                     9267: 'Didelphis virginiana',
                     9031: 'Gallus gallus',
                     10036: 'Mesocricetus auratus',
                     9940: 'Ovis aries',
                     10116: 'Rattus norvegicus',
                     9823: 'Sus scrofa',
                     8355: 'Xenopus laevis'}

dbptm_taxids = {9606: 'HUMAN',
                10090: 'MOUSE',
                7227: 'DROME',
                10116: 'RAT',
                559292: 'YEAST',
                284812: 'SCHPO',
                4081: 'SOLLC',
                3702: 'ARATH',
                9940: 'SHEEP',
                9913: 'BOVIN',
                9925: 'CAPHI',
                44689: 'DICDI',
                4577: 'MAIZE',
                9823: 'PIG',
                9615: 'CANLF',
                6239: 'CAEEL',
                8455: 'XENLA',
                83333: 'ECOLI',
                1891767: 'SV40'}

psite_mod_types = [('p', 'phosphorylation'),
                   ('ac', 'acetylation'),
                   ('ga', 'galactosylation'),
                   ('gl', 'glycosylation'),
                   ('sm', 'sumoylation'),
                   ('ub', 'ubiquitination'),
                   ('me', 'methylation')]

psite_mod_types2 = [('p', 'phosphorylation'),
                    ('ac', 'acetylation'),
                    ('ga', 'galactosylation'),
                    ('gl', 'glycosylation'),
                    ('ub', 'ubiquitination'),
                    ('me', 'methylation'),
                    ('sm', 'sumoylation'),
                    ('sc', 'succinylation'),
                    ('m1', 'mono-methylation'),
                    ('m2', 'di-methylation'),
                    ('m3', 'tri-methylation'),
                    ('ad', 'adenylation'),
                    ('pa', 'palmitoylation'),
                    ('ne', 'neddylation'),
                    ('sn', 'nitrosylation'),
                    ('ca', 'caspase-cleavage')]

mirbase_taxids = {9606: 'hsa',
                  10090: 'mmu',
                  10116: 'rno',
                  7227: 'dme'}


class silent(object): # XXX: Never used
    """

    """

    def __init__(self):

        pass

    def __enter__(self):
        self.aux = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exception_type, exception_value, traceback):

        sys.stdout.close()
        sys.stdout = self.aux
