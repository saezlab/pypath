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

from past.builtins import xrange, range

import collections

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.inputs.common as inputs_common
import pypath.utils.taxonomy as taxonomy


def mirtarbase_interactions(curated = True, strong = True, all = False):
    """
    Retrieves experimentally validated miRNA-target gene interactions from
    miRTarBase (mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2019/index.php).
    """

    result = []

    if curated:

        result.extend(_mirtarbase_interactions('curated'))

    if strong:

        result.extend(_mirtarbase_interactions('strong'))

    if all:

        result.extend(_mirtarbase_interactions('all'))

    return result


def _mirtarbase_interactions(dataset):

    MirtarbaseInteraction = collections.namedtuple(
        'MirtarbaseInteraction',
        (
            'mirtarbase_id',
            'mirna_name',
            'mirna_organism',
            'target_genesymbol',
            'target_entrez',
            'target_organism',
            'target_site',
            'method',
            'category',
            'pmid',
            'dataset',
        ),
    )

    url = urls.urls['mirtarbase'][dataset]
    c = curl.Curl(url, silent = False, large = True)

    if c.fileobj is None:
        return []

    # Read CSV data instead of Excel
    import csv
    
    # Handle generator result from curl
    if hasattr(c.result, '__iter__') and not isinstance(c.result, (str, bytes)):
        # Join generator lines
        data = ''.join(c.result)
    else:
        data = c.result
    
    if isinstance(data, bytes):
        data = data.decode('utf-8')
    
    reader = csv.reader(data.splitlines())
    tbl = list(reader)

    c.close()

    for i in xrange(len(tbl)):
        # Handle potential empty cells or short rows
        if len(tbl[i]) > 4 and tbl[i][4]:
            tbl[i][4] = tbl[i][4].split('.')[0]
        if len(tbl[i]) > 8 and tbl[i][8]:
            tbl[i][8] = tbl[i][8].split('.')[0]

    result = []
    for l in tbl[1:]:  # Skip header row
        if len(l) >= 9:  # Ensure minimum required columns
            try:
                interaction = MirtarbaseInteraction(
                    *l[:2],
                    taxonomy.ensure_ncbi_tax_id(l[2]),
                    *l[3:5],
                    taxonomy.ensure_ncbi_tax_id(l[5]),
                    l[6] if dataset == 'curated' else None,
                    *l[-3:],
                    dataset,
                )
                result.append(interaction)
            except (IndexError, ValueError):
                continue  # Skip malformed rows
    return result
