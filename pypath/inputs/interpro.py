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

from __future__ import annotations

import json
from typing import List, Dict, Union, Literal

import collections
from lxml import etree
import gzip
import shutil
import math

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.session as session
import pypath.inputs.common as inputs_common
from pypath.inputs import uniprot

_logger = session.Logger(name = 'inputs.interpro')
_log = _logger._log


def interpro_entries() -> List[tuple]:
    """
    Downloads detailed InterPro entry information.

    Returns
        A list of named tuples, each representing information about
        one InterPro entry.
    """

    InterproEntry = collections.namedtuple(
            'InterproEntry',
            (
                'interpro_id',
                'protein_count',
                'name',
                'type',
                'publications',
                'parent_list',
                'child_list',
                'member_list'

            ),
        )

    result=[]

    url = urls.urls['interpro']['entries']
    path = curl.Curl(
        url,
        silent = False,
        large = False
    ).fileobj.name
    with gzip.open(path, 'rb') as f_in:
        with open(path.split('.gz')[0], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    parser = etree.iterparse(path.split('.gz')[0], events = ('end',), tag='interpro')

    for ev, elem in parser:

        if elem.find('pub_list') is not None:

            pubs=[]

            for pub in elem.find('pub_list'):

                pubs.append(pub.attrib['id'])

        else:

            pubs=''

        if elem.find('parent_list') is not None:

            parent_ids=[]

            for parent in elem.find('parent_list'):

                parent_ids.append(parent.attrib['ipr_ref'])

        else:

            parent_ids=''

        if elem.find('child_list') is not None:

            child_ids=[]

            for child in elem.find('child_list'):

                child_ids.append(child.attrib['ipr_ref'])
        else:

            child_ids=''

        member_ids={}
        for member in elem.find('member_list'):

            if member.attrib['db'] in member_ids:

                member_ids[member.attrib['db']].append(member.attrib['dbkey'])

            else:

                member_ids[member.attrib['db']]=[]
                member_ids[member.attrib['db']].append(member.attrib['dbkey'])

        result.append(
            InterproEntry(
                interpro_id= elem.attrib['id'],
                protein_count= elem.attrib['protein_count'],
                name= elem.attrib['short_name'],
                type= elem.attrib['type'],
                publications= pubs,
                parent_list=parent_ids,
                child_list=child_ids,
                member_list=member_ids,
                )
        )

    return result


def interpro_xrefs(
        db_type: Literal[
                'go',
                'structural',
                'external',
            ],
    ) -> Dict[str, Union[List[str], Dict[str, List[str]]]]:
    """
    Downloads cross-references for each InterPro entry.

    Args
        db_type: Type of the cross-reference databases.

    Returns
        A dictionary; keys are InterPro IDs.
            If 'db_type' is 'go', values are list of GO terms related to the InterPro entry.
            Otherwise values are dictionaries, where keys are database names
            and the values are list of cross references related to the InterPro entry.
    """

    db_type_dict = {
    'go': 'class_list',
    'structural': 'structure_db_links',
    'external':  'external_doc_list',
    }
    db_type_name= db_type_dict[db_type]

    result={}

    url = urls.urls['interpro']['entries']
    path = curl.Curl(
        url,
        silent = False,
        large = False
    ).fileobj.name

    with gzip.open(path, 'rb') as f_in:
        with open(path.split('.gz')[0], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    parser = etree.iterparse(path.split('.gz')[0], events = ('end',), tag='interpro')

    for ev, elem in parser:

        interpro_id=elem.attrib['id']

        if db_type == 'go':

            go_terms = (
                [go.attrib['id'] for go in elem.find(db_type_name)]
                    if elem.find(db_type_name) is not None else
                None
            )

            result[interpro_id] = go_terms

        else:

            other_db_keys={}
            if elem.find(db_type_name) is not None:

                for link in elem.find(db_type_name):

                    if link.attrib['db'] in other_db_keys:

                        other_db_keys[link.attrib['db']].append(link.attrib['dbkey'])

                    else:

                        other_db_keys[link.attrib['db']]=[]
                        other_db_keys[link.attrib['db']].append(link.attrib['dbkey'])

            result[interpro_id] = other_db_keys

    return result



def _interpro_annotations_api(
    
        page_size: int = 200,
        reviewed: bool = True,
        tax_id: str | int = 9606,
    ) ->  dict:
    
    """
    Downloads UniProtKB proteins and the InterPro entries they match.

    Args
        page_size: Number of results returned at a time.
        reviewed: Downloads only reviewed UniprotKB proteins if True,
            Downloads all UniprotKB proteins otherwise.

    Returns
        A dictionary. Keys are Uniprot IDs, values are sets of annotations.
    """

    InterproAnnotation = collections.namedtuple(
        'InterproAnnotation',
        (
            'interpro_id',
            'organism',
            'start',
            'end',
        ),
    )

    annotations = collections.defaultdict(set)
    page = 0
    proteins = (
        'reviewed'
            if reviewed else
        'uniprot'
    )
    base_url = (
        urls.urls['interpro']['annotations'] % (proteins, tax_id, page_size)
    )
    next_page_url = base_url

    while next_page_url:

        c = curl.Curl(
            next_page_url,
            silent = False,
            large = False
        )

        res = inputs_common.json_read(c.result)
        totalrec = int(res['count'])


        _log(
            'Downloading page %u (total: %s).' % (
                page + 1,
                'unknown'
                    if totalrec < 0 else
                str(math.ceil(totalrec / page_size))
            )
        )

        

        for entry in res['results']:

            entry_info = entry['metadata']

            proteins_url = entry.get('proteins_url')


            if not proteins_url:
                raise ValueError(f"Missing proteins_url for {entry_info['accession']}")
        
            c2 = curl.Curl(proteins_url, silent=False, large=False)
            res2 = inputs_common.json_read(c2.result)

            for protein in res2.get('results', []):

                locations = protein.get('entry_protein_locations', [])

                for location in locations:

                    for fragment in location.get('fragments', []):

                        uniprot_id = protein['accession'].upper()

                        annotations[uniprot_id].add(
                            InterproAnnotation(
                                interpro_id=entry_info['accession'],
                                organism=protein.get('organism'),
                                start=int(fragment['start']),
                                end=int(fragment['end']),
                            )
                        )


        next_page_url = res.get('next')
        page += 1

    return annotations

def interpro_annotations(
    page_size: int = 200,
    reviewed: bool = True,
    tax_id: str | int = 9606,
) -> dict:
    """
    Downloads UniProtKB proteins and the InterPro entries they match.

    NOTE:
    - Default implementation uses BULK dataset (protein2ipr.dat.gz) for performance.
    - Falls back to API if bulk file is not available.
    - This behavior can be changed if needed.

    """

    import os
    import gzip
    import collections

    InterproAnnotation = collections.namedtuple(
        'InterproAnnotation',
        ('interpro_id', 'organism', 'start', 'end')
    )

    annotations = collections.defaultdict(set)

    # BULK DATA PATH
    BULK_FILE_PATH = os.getenv(
    "INTERPRO_BULK_PATH",
    r"D:\derssel\CROSSBAR\protein2ipr.dat.gz"
    )

    MAX_LINES = None  # set integer value for debugging

    # =========================
    # BULK MODE
    # =========================
    if os.path.exists(BULK_FILE_PATH):

        filter_ids = set(
            uniprot._all_uniprots(
                organism=tax_id,
                swissprot=reviewed
            )
        )



        with gzip.open(BULK_FILE_PATH, 'rt', encoding='utf-8', errors='ignore') as f:
            for i, line in enumerate(f):



                if MAX_LINES and i >= MAX_LINES:
                    break

                parts = line.strip().split('\t')

                if len(parts) < 6:
                    continue

                try:
                    uniprot_id = parts[0]
                    interpro_id = parts[1]
                    start = int(parts[4])
                    end = int(parts[5])


                except Exception as e:
                    continue

                if uniprot_id not in filter_ids:
                    continue

                annotations[uniprot_id].add(
                    InterproAnnotation(
                        interpro_id=interpro_id,
                        organism=tax_id,  # bulk dataset does not include organism, using input tax_id
                        start=start,
                        end=end
                    )
                )

        return annotations

    # =========================
    # FALLBACK: API MODE
    # =========================
    print("Bulk file not found, falling back to API...")

    return _interpro_annotations_api(
        page_size=page_size,
        reviewed=reviewed,
        tax_id=tax_id
    )

def interpro2go_annotations() -> dict[str, set[tuple]]:
    """
    Downloads GO term annotations for InterPro entries.

    Returns
        Dict of InterPro entries as keys and sets of GO terms as values.
    """

    import requests
    import collections

    url = urls.urls['interpro']['interpro2go']

    response = requests.get(url)
    response.raise_for_status()

    lines = response.text.splitlines()

    if not lines:
        raise ValueError(
            f"Failed to download interpro2go from {url}"
        )

    Interpro2GOAnnotation = collections.namedtuple(
        'Interpro2GOAnnotation',
        (
            'go_term_id',
            'go_term_name',
        ),
    )

    annotations = collections.defaultdict(set)

    for r in lines:

        if not r.startswith('!'):

            r = r.strip()

            try:
                interpro_id = r.split('InterPro:')[1].split(' ')[0]
                go_term_name = r.split('> GO:')[1].split(' ; ')[0]
                go_term_id = r.split('> GO:')[1].split(' ; ')[1]

                annotations[interpro_id].add(
                    Interpro2GOAnnotation(
                        go_term_id=go_term_id,
                        go_term_name=go_term_name
                    )
                )

            except Exception:
                # malformed line → skip
                continue

    return annotations
