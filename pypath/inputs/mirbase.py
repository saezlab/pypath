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

import re
import collections

import pypath.inputs.uniprot as uniprot_input
import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common
import pypath_common._constants as _const
import pypath.utils.taxonomy as taxonomy


def get_mirbase_aliases(organism = 9606):
    """
    Downloads and processes mapping tables from miRBase.
    """

    if type(organism) in _const.CHAR_TYPES:
        mborganism = organism
    elif organism not in taxonomy.mirbase_taxids:

        raise ValueError(
            'Organism not known: %u. Try to pass miRBase '
            'taxon prefix as string, e.g. `hsa`.' % organism
        )

    else:

        mborganism = taxonomy.mirbase_taxids[organism]

    mat = {}
    mir = {}

    url = urls.urls['mirbase']['aliases']
    c = curl.Curl(url, silent = False, large = True)

    for l in c.result:

        l = l.strip().strip(';').split('\t')

        if l[1][:3] != mborganism:
            continue

        d = mat if l[0][:5] == 'MIMAT' else mir

        if l[0] not in d:
            d[l[0]] = set([])

        for m in l[1].split(';'):
            d[l[0]].add(m)

    return mat, mir


def mirbase_mature(organism = 9606):

    mat, mir = get_mirbase_aliases(organism)

    result = {}

    for mimat, mmats in iteritems(mat):

        for mmat in mmats:

            yield mimat, mmat


def mirbase_precursor(organism = 9606):

    mat, mir = get_mirbase_aliases(organism)

    result = {}

    for mi, mpres in iteritems(mir):

        for mpre in mpres:

            yield mi, mpre


def mirbase_precursor_to_mature(organism = 9606):

    pre = mirbase_precursor(organism)
    ids = mirbase_ids(organism)

    _ids = collections.defaultdict(set)
    _pre = collections.defaultdict(set)

    for mmat, mpre in ids:

        _ids[mpre].add(mmat)

    for preid, prename in pre:

        _pre[prename].add(preid)

    result = {}

    for prename, mpres in iteritems(_pre):

        for mpre in mpres:

            if mpre in _ids:

                for mmat in _ids[mpre]:

                    yield prename, mmat


def mirbase_ids(organism = 9606):

    reprename = re.compile(r'([-A-z]*[-]?\d+[a-z]*)(-\d*)')

    def get_pre_name(mat_name):

        return mat_name.replace(
            '*', '').replace(
            '-3p', '').replace(
            '-5p', '')

    mat, mir = get_mirbase_aliases(organism)

    mir = dict((k, set.union(set(reprename.sub(r'\1', vv) for vv in v), v))
               for k, v in iteritems(mir))

    mir = common.swap_dict(mir)

    mat = dict((k, set(get_pre_name(vv) for vv in v))
               for k, v in iteritems(mat))

    if (sum(sum(vv in mir for vv in v) for v in mat.values()) <
        sum(sum(vv.lower() in mir for vv in v) for v in mat.values())):

        mat = dict((k, set(vv.lower() for vv in v))
               for k, v in iteritems(mat))

    mat_mir = common.join_dicts(mat, mir)

    for ma, mis in iteritems(mat_mir):

        for mi in (mis if type(mis) not in _const.SIMPLE_TYPES else [mis]):

            yield ma, mi


def mirbase_mature_all(organism = 9606):

    return [i[0] for i in mirbase_ids(organism = organism)]


def mirbase_precursor_all(organism = 9606):

    return [i[1] for i in mirbase_ids(organism = organism)]
