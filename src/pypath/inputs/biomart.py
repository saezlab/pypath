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
import pypath.utils.taxonomy as taxonomy

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
    """
    Query the Ensembl Biomart web service.
    Use https://www.ensembl.org/biomart/martview/ to check for attribute
    and dataset names.

    Args:
        attrs (str,list): One or more Ensembl attribute names.
        filters (str,list): One or more Ensembl filter names.
        transcript (bool): Include Ensembl transcript IDs in the result.
        gene (bool): Include Ensembl gene IDs in the result.
        dataset (str): An Ensembl dataset name.

    Yields:
        Named tuples with the requested attributes for each record returned
        by Ensembl Biomart.
    """

    _attrs = []

    if gene:

        _attrs.append('ensembl_gene_id')

    if transcript:

        _attrs.append('ensembl_transcript_id')

    _attrs.extend(common.to_list(attrs))
    filters = common.to_list(filters)

    record = collections.namedtuple(
        'EnsemblRecord',
        _attrs,
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


def biomart_homology(
        source_organism = 9606,
        target_organism = 10090,
    ):
    """
    Retrieves homology data from Ensembl Biomart.

    Returns:
        List of named tuples filtered to genes of the source organism having
        orthologues in the target organism, with homology related fields.
    """

    def ensure_organism(organism):

        organism_ensembl = taxonomy.ensure_ensembl_name(organism)

        if not organism_ensembl:

            msg = 'Could not find Ensembl taxon ID for `%s`.' % str(organism)
            _log(msg)
            raise ValueError(msg)

        return organism_ensembl


    source_organism = ensure_organism(source_organism)
    target_organism = ensure_organism(target_organism)

    attrs = [
        'ensembl_peptide_id'
    ]

    homolog_attrs = [
        'homolog_ensembl_peptide',
        'homolog_ensembl_gene',
        'homolog_orthology_type',
        'homolog_orthology_confidence',
        'homolog_canonical_transcript_protein',
    ]

    homolog_attrs = [
        '%s_%s' % (target_organism, attr)
        for attr in homolog_attrs
    ]

    filters = [
        'with_%s_homolog' % target_organism
    ]

    return list(
        biomart_query(
            attrs = attrs + homolog_attrs,
            filters = filters,
            transcript = True,
            gene = True,
            dataset = '%s_gene_ensembl' % source_organism,
        )
    )
