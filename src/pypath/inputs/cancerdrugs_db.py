#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Helps to translate from the mouse data to human data
#
#  Copyright
#  2014-2021
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Sebastian Lobentanzer
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/

"""
Web scraper for the CancerDrugs_DB database.

Todo:
    * translate chembl to CID (PubChem)
"""

import csv
import collections
import re

import pypath.share.curl as curl
import pypath.share.common as common
import pypath.resources.urls as urls
import pypath.utils.mapping as mapping
import pypath.share.session as session

_logger = session.Logger(name = 'cellcall_input')
_log = _logger._log

def cancerdrugs_db_download():
    """
    Downloads a curated set of interactions of cancer drugs licensed 
    in most parts of the world with gene targets (where available).
    From https://www.anticancerfund.org/en/cancerdrugs-db.
    This function downloads a single dataset.

    Args:
        None.

    Returns:
        A list of dicts, each is a record as it provided by the database.
    """

    url = urls.urls['cancerdrugs_db']['url']

    c = curl.Curl(url, large = True, silent = False)

    return list(csv.DictReader(c.result, delimiter = '\t'))


def cancerdrugs_db_interactions():
    """
    Returns drug-gene interactions from CancerDrugs_DB.

    Args:
        None.

    Returns:
        List of named tuples, each describing a drug-gene interaction.
        Identifiers of type PubChem and UniProt.
    """

    CancerDrugsInteraction = collections.namedtuple(
        'CancerDrugsInteraction',
        ['source', 'target'],
    )

    result = []

    data = cancerdrugs_db_download()

    for rec in data:
        
        chembl = rec['ChEMBL']

        # remove HTML link
        extr = re.search('">(.*?)</a>', chembl)

        if extr == None:
            continue
        # do we need to assert "CHEMBL" + integer structure of ID?

        chembl = extr.group(1)
        
        # targets are a number of gene symbols
        # -> translate to uniprot here?
        tars = rec['Targets']

        if tars == None:
            continue

        for tar in tars.split('; '):
            
            result.append(
                CancerDrugsInteraction(
                    source = chembl,
                    target = mapping.map_name(tar, 'genesymbol', 'uniprot'),
                )
            )

    return result


def cancerdrugs_db_annotations():
    """
    Returns drug annotations from CancerDrugs_DB.

    Args:
        None.

    Returns:
        Dict of annotations, keys are PubChem IDs, values are sets of
        annotations.
    """

    CancerDrugsAnnotation = collections.namedtuple(
        'CancerDrugsAnnotation',
        ['label', 'ATC'],
    )

    result = collections.defaultdict(set)

    data = cancerdrugs_db_download()

    for rec in data:
        chembl = rec['ChEMBL']

        # remove HTML link
        extr = re.search('">(.*?)</a>', chembl)

        if extr == None:
            continue
        # do we need to assert "CHEMBL" + integer structure of ID?

        chembl = extr.group(1)

        for _chembl in mapping.map_name(chembl, 'chembl', 'chembl'):
            result[chembl].add(CancerDrugsAnnotation(
                    label = common.upper0(rec['Product'].strip()),
                    ATC = rec['ATC'].strip(),
            ))

    return result
