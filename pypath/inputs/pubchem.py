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

import os

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.cache as cache
import pypath.share.lookup as lookup


def pubchem_mapping(target, source = 'cid'):
    """
    Identifier translation data from PubChem.

    Args
        target (str): The target ID type, either as it is used in the file
            names in the PubChem FTP service or as simpler, all lowercase
            strings used in this module. Possible values are parent-cid,
            component-cid, inchi, iupac, preferred-cid, sid, smiles,
            synonym.
        source (str): The source ID type. Either sid or cid.

    Returns
        (dict): A dict of sets with the source identifiers as keys and sets
            of target identifiers as values.
    """

    id_types = {
        'parent-cid': 'Parent',
        'component-cid': 'Component',
        'inchi': 'InChi',
        'iupac': 'IUPAC',
        'preferred-cid': 'Preferred',
        'pubchem-sid': 'SID',
        'sid': 'SID',
        'smiles': 'SMILES',
        'synonym': 'Synonym-unfiltered',
        'cid': 'CID',
        'pubchem-cid': 'CID',
    }

    _target = id_types.get(target, target)
    _source = id_types.get(source, source)

    if _source not in {'CID', 'SID'}:

        msg = (
            'The source identifier type must be either CID or SID, '
            'not `%s`.' % source
        )
        _log(msg)
        raise ValueError(msg)

    ftp_dir = (
        {
            'SID': 'Substance',
            'CID': 'Compound',
        }[_source]
    )

    url = urls.urls['pubchem']['ftp'] % (ftp_dir, _source, _target)
    c = curl.Curl(url, large = True, silent = False)

    db_path = os.path.join(
        cache.get_cachedir(),
        'pubchem_%s_%s.sqlite' % (_source, _target)
    )

    with lookup.ManyToMany(db_path) as result:

        result.populate(fileobj = c._gzfile_mode_r)
