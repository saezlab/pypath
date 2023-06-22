#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright
#  2014-2023
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: Dénes Türei (turei.denes@gmail.com)
#           Nicolàs Palacio
#           Sebastian Lobentanzer
#           Erva Ulusoy
#           Olga Ivanova
#           Ahmet Rifaioglu
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from __future__ import annotations

import collections

import pypath.share.curl as curl
import pypath.resources.urls as urls
import pypath.share.common as common


def sider_drug_names() -> dict[str, set[tuple]]:
    """
    Retrieves drug information from the SIDER database.

    Returns:
        Drug PubChem CID, name and ATC information as a list of named tuples.
    """

    SiderDrugName = collections.namedtuple(
        'SiderDrugName',
        ('name', 'atc'),
    )

    result = collections.defaultdict(set)
    attrs = {}

    for attr in ('name', 'atc'):

        url = urls.urls['sider'][f'drug_{attr}s']
        c = curl.Curl(url, large = True, silent = False)
        attrs[attr] = collections.defaultdict(list)

        for line in c.result:

            cid, value = line.strip('\n').split('\t')
            attrs[attr][cid].append(value)

    for cid in set.union(*map(set, attrs.values())):

        for atc in attrs['atc'].get(cid, (None,)):

            result[cid].add(SiderDrugName(
                name = attrs['name'].get(cid, (None,))[0],
                atc = atc,
            ))

    return dict(result)


def sider_side_effects(freq: bool = False) -> dict[str, set[tuple]]:
    """
    Retrieves side effect information from the SIDER database.

    Args:
        freq:
            Retrieve the dataset with frequency information. This is
            an independent dataset with lower coverage.

    Returns:
        Drug PubChem CID, UMLS concept ids both for label and MedDra
        and side effect name.
    """

    fields = (
        'umls_concept_on_label',
        'umls_concept_in_meddra',
        'side_effect',
    )
    fields += (('frequency',) if freq else ())
    record = collections.namedtuple(
        'SiderSideEffect%s' % ('Freq' if freq else ''),
        fields,
    )
    result = collections.defaultdict(set)

    url = urls.urls['sider']['meddra_%s' % ('freq' if freq else 'all')]

    c = curl.Curl(url, large = True, silent = False)

    # essential features' indices
    indices = (2, 8, 9, 4) if freq else (2, 4, 5)

    for line in c.result:

        line = line.strip().split('\t')

        if not line:

            continue

        result[line[0]].add(
            record(**{
                key: line[i] or None
                for key, i in zip(fields, indices)
            })
        )

    return dict(result)


def sider_side_effect_frequencies() -> list[tuple]:
    """
    Retrieves side effect information from the SIDER database.

    Returns:
        Drug CID, UMLS concept ids both for label and MedDRA,
        frequency information and side effect name.

    Attention! -> `sider_side_effects` function returns about 20k more rows
    than this dataset, but without frequency information.
    """

    return sider_side_effects(freq = True)

def sider_meddra_tsv() -> list[tuple]:
    """
    Retrieves MedDRA side effect information from the SIDER database.

    Returns:
        A list of named tuples containing the following fields:
        - cid: Drug PubChem CID
        - meddra_id: MedDRA ID for the side effect
        - side_effect_name: Name of the side effect
    """
    fields = (
        'cid',
        'kind_of_term',
        'meddra_id',
        'side_effect_name',
    )
    
    fields = common.to_list(fields)
    
    url_meddra_tsv = urls.urls['sider']['meddra_tsv']
    
    c = curl.Curl(
        url_meddra_tsv,
        large=True,
        silent=False
    )
    
    result = set()
    record = collections.namedtuple('Drug', fields[0:1] + fields[2:])
    
    for line in c.result:
        
        if not line.strip():
            continue
        
        line = line.strip().split('\t')
        line = dict(zip(fields, line))
        
        if line['kind_of_term'] == 'PT':
            result.add(
                record(
                    cid=line['cid'],
                    meddra_id=line['meddra_id'],
                    side_effect_name=line['side_effect_name']
                )
            )
    
    return list(result)
