#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#  Helps to translate from the mouse data to human data
#
#  Copyright
#  2014-2020
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#                  Olga Ivanova
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

import re
import os
import collections

import pypath.share.session as session_mod
import pypath.share.common as common
import pypath.share.curl as curl
import pypath.resources.urls as urls

_logger = session_mod.Logger(name = 'biomart_input')


# for mouse homologues: Filter name = "with_mmusculus_homolog"
_filter_xml_template = '<Filter name="%s" excluded="0"/>'
_attr_xml_template = '<Attribute name="%s" />'


def biomart_query(
        attrs,
        filters = None,
        transcript = False,
        gene = False,
        dataset = 'hsapiens_gene_ensembl',
    ):

    _attrs = []

    if gene:

        _attrs.append('ensembl_gene_id')

    if transcript:

        _attrs.append('ensembl_transcript_id')

    _attrs.extend(common.to_list(attrs))
    filters = common.to_list(filters)

    record = collections.namedtuple(
        'EnsemblRecord',
        _attrs
    )

    _logger._log(
        'Downloading data from Ensembl Biomart: '
        'dataset=`%s`, '
        '%s'
        'attributes=`%s`.' % (
            dataset,
            (
                'filters=`%s`, ' % ', '.join(filters)
                    if filters else
                ''
            ),
            ', '.join(_attrs),
        )
    )

    rewsp = re.compile(r'\n\s+')

    xml_template_path = os.path.join(common.DATA, 'ensembl_biomart_query.xml')

    with open(xml_template_path, 'r') as fp:

        xml_template = fp.read()

    filter_part = ''.join(
        _filter_xml_template % _filter
        for _filter in filters
    )
    attr_part = ''.join(
        _attr_xml_template % _attr
        for _attr in _attrs
    )

    xml_query = xml_template % (
        dataset,
        filter_part,
        attr_part,
    )
    xml_query = rewsp.sub('', xml_query)

    biomart_url = urls.urls['ensembl']['biomart_url'] % xml_query

    c = curl.Curl(biomart_url, large = True, silent = False)

    success = False

    for line in c.result:

        line = line.strip('\n\r').split('\t')

        success = success or line[0] == '[success]'

        if len(line) == len(record._fields):

            yield record(*line)

    if not success:

        _logger._log(
            'Error: Interrupted transfer while downlading data '
            'from Ensembl Biomart (missing `success` tag).'
        )