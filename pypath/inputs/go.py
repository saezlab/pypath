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

from past.builtins import xrange, range
from future.utils import iteritems

import os
import re
import json
import collections

from lxml import etree

import pypath.resources.urls as urls
import pypath.share.curl as curl
import pypath.share.progress as progress
import pypath.share.session as session
import pypath.share.common as common
import pypath.utils.taxonomy as taxonomy

_logger = session.Logger(name = 'uniprot_input')
_log = _logger._log


HEADER_ACCEPT_TSV = {'Accept': 'text/tsv'}
HEADER_ACCEPT_JSON = {'Accept': 'application/json'}


def go_annotations_uniprot(organism = 9606, swissprot = 'yes'):
    """
    Deprecated, should be removed soon.
    """

    rev = '' if swissprot is None \
        else ' AND reviewed:%s' % swissprot
    query = 'organism:%u%s' % (int(organism), rev)
    url = urls.urls['uniprot_basic']['url']
    get = {'query': query, 'format': 'tab', 'columns': 'id,go-id'}
    c = curl.Curl(url, get = get, silent = False)
    data = c.result

    return dict([(x[0], [go.strip() for go in x[1].split(';')])
                 for x in [x.split('\t') for x in data.split('\n')]
                 if len(x) > 1])


def go_annotations_goa(
        organism = 'human',
        evidence_codes=False,
    ):
    """
    Downloads GO annotation from UniProt GOA.

    Args:
        organism:
            Organism name or NCBI Taxonomy ID.
        evidence_codes:
            Include evidence codes in the output.
    """

    organism = taxonomy.ensure_common_name(organism)

    annot = dict(
        (asp, collections.defaultdict(set))
        for asp in ('C', 'P', 'F')
    )

    url = urls.urls['goa']['ebi_url'] % (organism.upper(), organism.lower())
    c = curl.Curl(url, silent = False, large = True)

    for line in c.result:
        if not line or line[0] == '!':
            continue

        line = line.strip().split('\t')
        if evidence_codes:
            annot[line[8]][line[1]].add((line[4], line[6]))
        else:
            annot[line[8]][line[1]].add(line[4])

    return dict((k, dict(v)) for k, v in iteritems(annot))


# synonym for the default method
go_annotations = go_annotations_goa


def go_annotations_all(
        organism: int | str = 'human',
        fields: str | list[str] | None = None
    ) -> dict[str, set[tuple[str]]]:

    if organism != '*':
        organism = taxonomy.ensure_common_name(organism)

    all_fields = (
        'db',
        'db_object_id',
        'db_object_symbol',
        'qualifier',
        'go_id',
        'reference',
        'evidence_code',
        'with_or_from',
        'aspect',
        'db_object_name',
        'db_object_synonym',
        'db_object_type',
        'taxon_and_interacting_taxon',
        'date',
        'assigned_By',
        'annotation_extension',
        'gene_product_form_id'
    )

    fields = fields or all_fields
    fields = common.to_list(fields)

    if organism in ('*', None):
        url = urls.urls['goa']['ebi_url'] % ('UNIPROT', 'uniprot_gcrp')
    else:
        url = urls.urls['goa']['ebi_url'] % (organism.upper(), organism.lower())

    c = curl.Curl(url, silent = False, large = True)

    result = collections.defaultdict(set)
    record = collections.namedtuple('GoAnnotation', fields)

    for line in c.result:

        if not line.strip() or line[0] == '!':
            continue

        line = dict(zip(all_fields, line.strip().split('\t')))
        result[line['db_object_id']].add(
            record(**dict(zip(fields, (line.get(f, None) for f in fields))))
        )

    return dict(result)


def go_ancestors_goose(aspects = ('C','F','P')):
    """
    Queries the ancestors of GO terms by AmiGO goose.

    Returns dict of sets where keys are GO accessions and values are sets
    of their ancestors.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    aspects_part = ''
    respaces = re.compile(r'[\s\n]+')

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }

    if set(aspects) != {'C', 'F', 'P'}:

        aspects_part = 'WHERE (%s)' % (
            ' OR '.join(
                'term.term_type = "%s"' % ontologies[asp]
                for asp in aspects
            )
        )

    sql_path = os.path.join(common.DATA, 'goose_ancestors.sql')

    with open(sql_path, 'r') as fp:

        query = fp.read()

    query = query % aspects_part
    query = respaces.sub(r' ', query).strip()

    url = urls.urls['goose']['url'] % query

    c = curl.Curl(url, silent = False, large = True)

    ancestors = collections.defaultdict(set)

    for l in c.result:

        l = l.strip().split('\t')
        ancestors[l[0]].add(l[1])

    return ancestors


def go_ancestors_quickgo(aspects = ('C', 'F', 'P')):
    """
    Queries the ancestors of GO terms by QuickGO REST API.

    Returns dict of sets where keys are GO accessions and values are sets
    of their ancestors.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    desc = go_descendants_quickgo(aspects = aspects)

    return go_descendants_to_ancestors(desc)


# synonym for the default method
go_ancestors = go_ancestors_quickgo


def go_descendants_to_ancestors(desc):
    """
    Turns a dict of descendants to dict of ancestors by swapping the
    relationships. This way descendants will be the keys and their ancestors
    will be the values.
    """

    ancestors = {}

    for asp, dct in iteritems(desc):

        ancestors[asp] = collections.defaultdict(set)

        for anc_term, des in iteritems(dct):

            for des_term, rel in des:

                ancestors[asp][des_term].add((anc_term, rel))

        ancestors[asp] = dict(ancestors[asp])

    return ancestors


def go_descendants_goose(aspects = ('C','F','P')):
    """
    Queries descendants of GO terms by AmiGO goose.

    IMPORTANT:
    This is not the preferred method any more to get descendants.
    Recently the preferred method to access GO annotations is
    ``pypath.inputs.go.go_descendants_quickgo``.
    The data in GO MySQL instances has not been updated since Dec 2016.
    Unfortunately the providers ceased to support MySQL, the most flexible
    and highest performance access to GO data. The replacement is Solr
    which is far from providing the same features as MySQL, for example
    it is unable to provide GO graph relationships. Other service is QuickGO
    which is up to date and has nice ways to query the ontology.

    Returns dict of sets where keys are GO accessions and values are sets
    of their descendants.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    desc = collections.defaultdict(set)

    anc = go_ancestors_goose(aspects = aspects)

    for term, ancs in iteritems(anc):
        for terma in ancs:
            desc[terma].add(term)

    return desc


def go_descendants_quickgo(
        aspects = ('C', 'F', 'P'),
        terms = None,
        relations = None,
        quickgo_download_size = 500,
    ):
    """
    Queries descendants of GO terms by QuickGO REST API.

    Returns dict of sets where keys are GO accessions and values are sets
    of their descendants.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    :param dict terms:
        Result from ``go_terms_solr``. If ``None`` the method will be called.
    """


    def download_in_chunks(terms, chunk_size, target = None):

        target = target or collections.defaultdict(set)

        paginator = common.paginate(terms, chunk_size)

        for p, terms_part in enumerate(paginator):

            url = urls.urls['quickgo_rest']['desc'] % (
                ','.join(terms_part),
                '?relations = %s' % relations_part,
            )

            c = curl.Curl(
                url,
                req_headers = HEADER_ACCEPT_JSON,
                silent = True,
                large = True,
            )

            try:
                result = json.load(c.fileobj)

            except json.decoder.JSONDecodeError:

                done = chunk_size * p
                remaining = terms[done:]
                new_chunk_size = chunk_size // 2

                if new_chunk_size < 10:

                    _log(
                        'Failed to download QuickGO, tried to decrease the '
                        'number of terms in each query, went below 10 terms '
                        'per query but still getting erroneous JSON. '
                        'This might be due to very slow network connection. '
                        'You might increase the timeout of CURL. '
                        'But then it will take forever.'
                    )

                    return target

                return download_in_chunks(
                    terms = remaining,
                    chunk_size = new_chunk_size,
                    target = taret,
                )

            for res in result['results']:
                if 'children' not in res:
                    continue

                target[res['id']].update(
                    set(
                        (child['id'], child['relation'])
                        for child in res['children']
                    )
                )

        return target

    desc = {}

    terms = terms or go_terms_quickgo(aspects = aspects)
    relations = relations or ('is_a', 'part_of', 'occurs_in', 'regulates',)

    relations_part = ','.join(relations)

    for asp in aspects:
        desc[asp] = download_in_chunks(
            terms = list(terms[asp].keys()),
            chunk_size = quickgo_download_size,
        )

    return desc


# synonym for the default method
go_descendants = go_descendants_quickgo


def go_terms_solr(aspects = ('C', 'F', 'P')):
    """
    Queries GO terms by AmiGO Solr.

    Returns dict of dicts where upper level keys are one letter codes of the
    aspects `C`, `F` and `P` for cellular_component, molecular_function and
    biological_process, respectively. Lower level keys are GO accessions
    and values are names of the terms.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    reamp = re.compile(r'[\s\n\r]+([&\?])')
    relin = re.compile(r'[\s\n\r]+')

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    terms = dict((a, {}) for a in aspects)

    query = '''
        ?q = document_category:"ontology_class" AND
            idspace:GO AND
            is_obsolete:0
        &rows = 9999999
        &start = 0
        &fl = annotation_class,annotation_class_label,source
    '''

    query = relin.sub(' ', reamp.sub(r'\1', query.strip()))

    # downloading data
    url = urls.urls['golr']['url'] % query

    c = curl.Curl(url, silent = False, large = True)

    # parsing XML by lxml.etree.iterparse
    parser = etree.iterparse(c.fileobj, events = ('start', 'end'))
    root = next(parser)
    used_elements = []

    for ev, elem in parser:
        if ev == 'end' and elem.tag == 'doc':
            asp  = elem.find('.//str[@name="source"]').text
            asp  = ontol_short[asp]

            if asp not in aspects:
                continue

            term = elem.find('.//str[@name="annotation_class"]').text
            name = elem.find('.//str[@name="annotation_class_label"]').text

            terms[asp][term] = name

        used_elements.append(elem)

        # removing used elements to keep memory low
        if len(used_elements) > 1000:
            for _ in xrange(500):
                e = used_elements.pop(0)
                e.clear()

    # closing the XML
    c.fileobj.close()
    del c

    return terms


def go_terms_quickgo(aspects = ('C','F','P')):
    """
    Queries GO terms by the QuickGO REST API.

    Return dict of dicts where upper level keys are one letter codes of the
    aspects `C`, `F` and `P` for cellular_component, molecular_function and
    biological_process, respectively. Lower level keys are GO accessions
    and values are names of the terms.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    result = dict((a, {}) for a in aspects)
    url = urls.urls['quickgo_rest']['terms']
    last_page = 9999999
    this_page = 1
    prg = progress.Progress(
        name = 'Downloading data from QuickGO',
        interval = 1,
    )

    while this_page <= last_page:
        page_url = url % this_page

        c = curl.Curl(page_url, silent = True)

        this_result = json.loads(c.result)
        last_page = this_result['pageInfo']['total']


        for res in this_result['results']:
            if 'aspect' not in res:
                continue

            asp = ontol_short[res['aspect']]

            if res['isObsolete'] or asp not in aspects:
                continue

            result[asp][res['id']] = res['name']

        if prg.total is None:
            prg.set_total(last_page)

        prg.step()

        this_page += 1

    return result


# synonym for the default method
go_terms = go_terms_quickgo


def go_terms_goose(aspects = ('C','F','P')):
    """
    Queries GO terms by AmiGO goose.

    Return dict of dicts where upper level keys are one letter codes of the
    aspects `C`, `F` and `P` for cellular_component, molecular_function and
    biological_process, respectively. Lower level keys are GO accessions
    and values are names of the terms.

    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    """

    aspects_part = ''
    respaces = re.compile(r'[\s\n]+')

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    if set(aspects) != {'C', 'F', 'P'}:
        aspects_part = 'WHERE (%s)' % (
            ' OR '.join(
                'term.term_type = "%s"' % ontologies[asp]
                for asp in aspects
            )
        )

    sql_path = os.path.join(common.DATA, 'goose_terms.sql')

    with open(sql_path, 'r') as fp:
        query = fp.read()

    query = query % aspects_part
    query = respaces.sub(r' ', query).strip()

    url = urls.urls['goose']['url'] % query

    c = curl.Curl(url, silent = False, large = True)

    terms = {'P': {}, 'C': {}, 'F': {}}

    for l in c.result:
        l = l.strip().split('\t')

        if l[1] not in ontol_short:
            continue

        aspect = ontol_short[l[1]]
        terms[aspect][l[2]] = l[0]

    return terms


def go_annotations_quickgo(
        organism = 9606,
        aspects = ('C','F','P'),
        relations = ('is_a', 'part_of'),
    ):
    """
    Queries GO annotations by QuickGO REST API.

    IMPORTANT:
    Recently the preferred method to access GO annotations is
    ``pypath.inputs.go.go_annotations_goa``.
    Contrary to its name QuickGO is super slow, otherwise it should yield
    up to date data, identical to the GOA file.

    Returns terms in dict of dicts and annotations in dict of dicts of sets.
    In both dicts the keys are aspects by their one letter codes.
    In the term dicts keys are GO accessions and values are their names.
    In the annotation dicts keys are UniProt IDs and values are sets
    of GO accessions.

    :param int organism:
        NCBI Taxonomy ID of one organism. Default is human (9606).
    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    :param list uniprots:
        Optionally a list of UniProt IDs. If `None`, results for all proteins
        returned.
    """

    annot = dict((a, collections.defaultdict(set)) for a in aspects)

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    url = urls.urls['quickgo_rest']['annot']

    aspects_part = ','.join(ontologies[a] for a in aspects)
    relations_part = ','.join(relations)

    page = 1

    while True:

        this_url = url % (
            aspects_part, # aspect
            relations_part, # goUsageRelationships
            organism, # taxonId
            page,
        )

        c = curl.Curl(
            url = this_url,
            req_headers = HEADER_ACCEPT_TSV,
            silent = False,
            large = True
        )

        _ = next(c.result) # the header row

        for line in c.result:
            line = line.strip().split('\t')

            if line[3] not in relations:
                continue

            annot[line[5]][line[1]].add(line[4])

        page += 1

    return annot


def go_annotations_solr(
        organism = 9606,
        aspects = ('C', 'F', 'P'),
        references = False,
    ):
    """
    Queries GO annotations by AmiGO Solr.

    Before other methods have been provided to access GO.
    Now this is the preferred method to get annotations.
    Returns terms in dict of dicts and annotations in dict of dicts of sets.
    In both dicts the keys are aspects by their one letter codes.
    In the term dicts keys are GO accessions and values are their names.
    In the annotation dicts keys are UniProt IDs and values are sets
    of GO accessions.

    :param int organism:
        NCBI Taxonomy ID of one organism. Default is human (9606).
    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    :param bool references:
        Retrieve the references (PubMed IDs) for the annotations.
        Currently not implemented.
    """

    reamp = re.compile(r'[\s\n\r]+([&\?])')
    relin = re.compile(r'[\s\n\r]+')

    annot = dict((a, collections.defaultdict(set)) for a in aspects)

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    # assembling the query

    if len(aspects) < 3:
        aspects_part = ' AND (%s)' % (
            ' OR '.join('aspect:%s' % a for a in aspects)
        )

    else:
        aspects_part = ''

    refs_part = ',reference' if references else ''

    query = '''
        ?q = taxon:"NCBITaxon:%u" AND
            type:protein AND
            document_category:annotation AND
            source:UniProtKB%s
        &rows = 9999999
        &start = 0
        &fl = bioentity,annotation_class,aspect%s
    ''' % (
        organism,
        aspects_part,
        refs_part
    )

    query = relin.sub(' ', reamp.sub(r'\1', query.strip()))

    # downloading data
    url = urls.urls['golr']['url'] % query
    c = curl.Curl(url, silent = False, large = True)

    # parsing XML by lxml.etree.iterparse
    parser = etree.iterparse(c.fileobj, events = ('start', 'end'))
    root = next(parser)
    used_elements = []

    for ev, elem in parser:

        if ev == 'end' and elem.tag == 'doc':

            id_ = elem.find('.//str[@name="bioentity"]').text

            if not id_.startswith('UniProtKB:'):
                continue

            asp  = elem.find('.//str[@name="aspect"]').text

            if asp not in aspects:
                continue

            term = elem.find('.//str[@name="annotation_class"]').text
            id_  = id_[10:] # removing the `UniProtKB:` prefix

            # adding the term to the annotation dict
            annot[asp][id_].add(term)

        used_elements.append(elem)

        # removing used elements to keep memory low
        if len(used_elements) > 1000:

            for _ in xrange(500):

                e = used_elements.pop(0)
                e.clear()

    # closing the XML
    c.fileobj.close()
    del c

    return terms, annot


def go_annotations_goose(organism = 9606, aspects = ('C', 'F', 'P'),
                         uniprots = None):
    """
    Queries GO annotations by AmiGO goose.

    IMPORTANT:
    This is not the preferred method any more to get terms and annotations.
    Recently the preferred method to access GO annotations is
    ``pypath.inputs.go.go_annotations_solr``.
    The data in GO MySQL instances has not been updated since Dec 2016.
    Unfortunately the providers ceased to support MySQL, the most flexible
    and highest performance access to GO data. The replacement is Solr
    which is far from providing the same features as MySQL.

    Returns terms in dict of dicts and annotations in dict of dicts of sets.
    In both dicts the keys are aspects by their one letter codes.
    In the term dicts keys are GO accessions and values are their names.
    In the annotation dicts keys are UniProt IDs and values are sets
    of GO accessions.

    :param int organism:
        NCBI Taxonomy ID of one organism. Default is human (9606).
    :param tuple aspects:
        GO aspects: `C`, `F` and `P` for cellular_component,
        molecular_function and biological_process, respectively.
    :param list uniprots:
        Optionally a list of UniProt IDs. If `None`, results for all proteins
        returned.
    """

    aspects_part = ''
    uniprot_part = ''
    respaces = re.compile(r'[\s\n]+')

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }
    ontol_short = dict(reversed(i) for i in ontologies.items())

    if set(aspects) != {'C', 'F', 'P'}:
        aspects_part = '(%s) AND' % (
            ' OR '.join(
                'term.term_type="%s"' % ontologies[asp]
                for asp in aspects
            )
        )

    if uniprots is not None:
        uniprot_part = 'dbxref.xref_key IN (%s) AND' % (
            ','.join('"%s"' % uniprot for uniprot in uniprots)
        )

    sql_path = os.path.join(common.DATA, 'goose_annotations.sql')

    with open(sql_path, 'r') as fp:
        query = fp.read()

    query = query % (organism, aspects_part, uniprot_part)
    query = respaces.sub(r' ', query).strip()

    url = urls.urls['goose']['url'] % query

    c = curl.Curl(url, silent = False, large = True)

    terms = {'P': {}, 'C': {}, 'F': {}}
    annot = {
        'C': collections.defaultdict(set),
        'F': collections.defaultdict(set),
        'P': collections.defaultdict(set),
    }

    for l in c.result:
        l = l.strip().split('\t')

        aspect = ontol_short[l[1]]

        terms[aspect][l[2]] = l[0]
        annot[aspect][l[5]].add(l[2])

    return terms, annot


def get_go_desc(go_ids, organism = 9606):
    """
    Deprecated, should be removed soon.
    """

    go_ids = (
        ','.join(sorted(go_ids))
        if type(go_ids) in {list, tuple, set} else
        go_ids
    )

    url = urls.urls['quickgo_desc']['url'] % (organism, go_ids)
    c = curl.Curl(
        url,
        silent = False,
        large = True,
        req_headers = HEADER_ACCEPT_TSV,
    )
    _ = c.result.readline()

    return set(l.split('\t')[1] for l in c.result)


def get_go_quick(
        organism = 9606,
        slim = False,
        names_only = False,
        aspects = ('C', 'F', 'P'),
    ):
    """
    Deprecated, should be removed soon.

    Loads GO terms and annotations from QuickGO.
    Returns 2 dicts: `names` are GO terms by their IDs,
    `terms` are proteins GO IDs by UniProt IDs.
    """

    ontologies = {
        'C': 'cellular_component',
        'F': 'molecular_function',
        'P': 'biological_process',
    }

    terms = {
        'C': collections.defaultdict(set),
        'F': collections.defaultdict(set),
        'P': collections.defaultdict(set),
    }

    names = {}
    aspects_param = ','.join(sorted(ontologies[a] for a in aspects))

    url = urls.urls['quickgo']['url'] % (
        organism,
        aspects_param,
        '&goUsage = slim' if slim else '',
    )

    c = curl.Curl(
        url,
        silent = False,
        large = True,
        req_headers = HEADER_ACCEPT_TSV,
        keep_failed = True,
    )

    _ = next(c.result)

    for l in c.result:

        l = l.split('\t')

        if not names_only:

            terms[l[5]][l[1]].add(l[4])

    return {'terms': terms, 'names': names}


def get_goslim(url = None):

    rego = re.compile(r'GO:[0-9]{7}')
    url = (
        url
            if isinstance(url, str) else
        urls.urls['goslim_gen']['url']
    )
    c = curl.Curl(url, silent = False)
    data = c.result
    result = []

    for l in data.split('\n'):

        if l.startswith('id:'):

            result += rego.findall(l)

    return result
