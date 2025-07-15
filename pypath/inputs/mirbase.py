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
    NOTE: The miRBase aliases file format has changed. This function
    now returns a simpler mapping based on available data.
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

    if not c.result:
        return mat, mir

    # Handle HTML-formatted content from new miRBase site
    # New format: internal_id, db_type, external_id, name
    for line in c.result:
        # Remove HTML tags and split by <br> to get individual lines
        line = line.replace('<p>', '').replace('</p>', '')
        entries = line.split('<br>')
        
        for entry in entries:
            entry = entry.strip()
            if not entry:
                continue
                
            parts = entry.split('\t')
            if len(parts) < 4:
                continue
                
            internal_id, db_type, external_id, name = parts[:4]
            
            # Filter for human miRNAs (organism-specific names)
            # Look for names that start with organism prefix or are typical miRNA names
            is_organism_specific = (
                name.lower().startswith(mborganism.lower()) or
                name.lower().startswith('mir') or
                name.lower().startswith('let-')
            )
            
            if not is_organism_specific:
                continue
                
            # For NCBI Gene entries (db_type 5), these are likely mature miRNAs
            # For other types, treat as precursors
            if db_type == '5' and name.startswith('MIR'):
                # Mature miRNA - use external_id as key, name as alias
                if external_id not in mat:
                    mat[external_id] = set()
                mat[external_id].add(name)
                # Also add lowercase version for compatibility
                if external_id not in mat:
                    mat[external_id] = set()
                mat[external_id].add(name.lower())
            else:
                # Precursor or other - use internal_id as key, name as alias  
                if internal_id not in mir:
                    mir[internal_id] = set()
                mir[internal_id].add(name)

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
