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
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Sebastian Lobentanzer
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/

"""
Web scraper for the CancerDrugs_DB database.
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

    # note: drug-target interactions are from GDIdb. if they're in OnmiPath
    # already, do we skip the interactions altogether, or is it better to 
    # keep them here for consistency? the dataset is very small anyways.
    # if we keep them, do we test for consistency between the interactions?

    CancerDrugsInteraction = collections.namedtuple(
        'CancerDrugsInteraction',
        ['source', 'target'],
    )

    result = []
    unmapped_id = []
    no_tars = []
    unmapped_tars = []

    data = cancerdrugs_db_download()

    for rec in data:
        
        chembl = rec['ChEMBL']

        # remove HTML link
        extr = re.search('">(.*?)</a>', chembl)

        if extr == None:
            unmapped_id.append(rec['Product'])
            continue
        # do we need to assert "CHEMBL" + integer structure of ID?

        chembl = extr.group(1)
        pubchem = mapping.map_name(chembl, 'chembl', 'pubchem')
        
        # targets are a number of gene symbols
        # -> translate to uniprot here?
        tars = rec['Targets']

        if tars == '':
            no_tars.append(rec['Product'])
            continue

        for tar in tars.split('; '):
            tar = mapping.map_name(tar, 'genesymbol', 'uniprot')
            if tar == None:
                unmapped_tars.append(tar)
                continue
            
            result.append(
                CancerDrugsInteraction(
                    source = pubchem,
                    target = tar,
                )
            )
    
    _log(
        'Could not find CHEMBL IDs for %u '
        'CancerDrugs_DB Products.' % len(unmapped_id)
    )
    
    _log(
        '%u CancerDrugs_DB Products had no targets.'
        % len(no_tars)
    )
    
    _log(
        'Could not find UniProt IDs for %u '
        'CancerDrugs_DB Product targets.' % len(unmapped_tars)
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
        ['label', 'indications', 'last_updated'],
    )

    result = collections.defaultdict(set)
    unmapped_id = []

    data = cancerdrugs_db_download()

    for rec in data:
        chembl = rec['ChEMBL']

        # remove HTML link
        extr = re.search('">(.*?)</a>', chembl)

        if extr == None:
            unmapped_id.append(rec['Product'])
            continue
        # do we need to assert "CHEMBL" + integer structure of ID?

        chembl = extr.group(1)
        pubchem = mapping.map_name(chembl, 'chembl', 'pubchem')

        for pubchem in pubchem:
            result[pubchem].add(CancerDrugsAnnotation(
                    label = common.upper0(rec['Product'].strip()),
                    indications = rec['Indications'].strip(),
                    last_updated = rec['Last Update'].strip(),
            ))

    _log(
        'Could not find CHEMBL IDs for %u '
        'CancerDrugs_DB Products.' % len(unmapped_id)
    )

    return result
