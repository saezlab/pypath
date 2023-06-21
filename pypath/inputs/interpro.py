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


def interpro_annotations(
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

            for protein in entry['protein_subset']:

                locations = protein['entry_protein_locations']

                for location in locations:

                    for fragment in location['fragments']:

                        uniprot_id = protein['accession'].upper()
                        annotations[uniprot_id].add(
                            InterproAnnotation(
                                interpro_id = entry_info['accession'],
                                organism = protein['organism'],
                                start = int(fragment['start']),
                                end = int(fragment['end']),
                            )
                        )

        next_page_url = res.get('next')
        page += 1

    return annotations


def interpro2go_annotations() -> dict[str, set[tuple]]:
    """
    Downloads GO term annotations for InterPro entries.

    Returns
        Dict of InterPro entries as keys and sets of GO terms as values.
    """

    url = urls.urls['interpro']['interpro2go']
    c = curl.Curl(url, large = True, silent = False)

    Interpro2GOAnnotation = collections.namedtuple(
        'Interpro2GOAnnotation',
        (
            'go_term_id',
            'go_term_name',
        ),
    )

    annotations = collections.defaultdict(set)

    for r in c.result:

        if not r.startswith('!'):

            r = r.strip()
            interpro_id = r.split('InterPro:')[1].split(' ')[0]
            go_term_name = r.split('> GO:')[1].split(' ; ')[0]
            go_term_id = r.split('> GO:')[1].split(' ; ')[1]

            annotations[interpro_id].add(
                Interpro2GOAnnotation(
                    go_term_id = go_term_id,
                    go_term_name = go_term_name
                )
            )

    return annotations
