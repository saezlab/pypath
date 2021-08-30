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

from future.utils import iteritems

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


microarrays = {
    'AFFY HG U95E': 'affy_hg_u95e',
    'AFFY HG G110': 'affy_hc_g110',
    'AFFY HG Focus': 'affy_hg_focus',
    'AFFY HG U133A 2': 'affy_hg_u133a_2',
    'AFFY HG U133A': 'affy_hg_u133a',
    'AFFY HG U133B': 'affy_hg_u133b',
    'AFFY HG U133 Plus 2': 'affy_hg_u133_plus_2',
    'AFFY HG U95A': 'affy_hg_u95a',
    'AFFY HG U95Av2': 'affy_hg_u95av2',
    'AFFY HG U95B': 'affy_hg_u95b',
    'AFFY HG U95D': 'affy_hg_u95d',
    'AFFY HG U95C': 'affy_hg_u95c',
    'AFFY HTA 2 0': 'affy_hta_2_0',
    'AFFY HG HuEx 1 0 st v2': 'affy_huex_1_0_st_v2',
    'AFFY HG HuGeneFL': 'affy_hugenefl',
    'AFFY HG HuGene 1 0 st v1': 'affy_hugene_1_0_st_v1',
    'AFFY HG HuGene 2 0 st v1': 'affy_hugene_2_0_st_v1',
    'AFFY HG PrimeView': 'affy_primeview',
    'PHALANX OneArray': 'phalanx_onearray',
    'ILLUMINA HumanWG 6 V3': 'illumina_humanwg_6_v3',
    'ILLUMINA HumanWG 6 V2': 'illumina_humanwg_6_v2',
    'ILLUMINA HumanWG 6 V1': 'illumina_humanwg_6_v1',
    'ILLUMINA HumanRef 8 V3': 'illumina_humanref_8_v3',
    'ILLUMINA HumanHT 12 V4': 'illumina_humanht_12_v4',
    'ILLUMINA HumanHT 12 V3': 'illumina_humanht_12_v3',
    'AGILENT WholeGenome 4x44k v2': 'agilent_wholegenome_4x44k_v2',
    'CODELINK CODELINK': 'codelink_codelink',
    'AGILENT WholeGenome 4x44k v1': 'agilent_wholegenome_4x44k_v1',
    'AGILENT WholeGenome': 'agilent_wholegenome',
    'AGILENT SurePrint G3 GE 8x60k v2': 'agilent_sureprint_g3_ge_8x60k_v2',
    'AGILENT SurePrint G3 GE 8x60k': 'agilent_sureprint_g3_ge_8x60k',
    'AGILENT GPL6848': 'agilent_gpl6848',
    'AGILENT GPL26966': 'agilent_gpl26966',
    'AGILENT CGH 44b': 'agilent_cgh_44b',
    'AFFY HG U133 X3P': 'affy_u133_x3p',
}


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


def biomart_microarray(array_type):
    """
    Microarray probe identifier mappings.

    Args:
        array_type (str): The microarray model, as shown on the BioMart
            webpage, or the corresponding code. For a full list of available
            identifiers see the ``microarrays`` attribute of this module.

    Returns:
        A dictionary with Ensembl gene, transcript and peptide IDs as keys
        and sets of microarray probe IDs as values.
    """

    _microarrays = dict(
        (k.lower(), v)
        for k, v in iteritems(microarrays)
    )
    all_array_types = set(microarrays.values())

    _array_type = array_type.lower().replace('probe', '').strip()
    _array_type = common.maybe_in_dict(_microarrays, _array_type)

    if _array_type not in all_array_types:

        msg = 'No such array type in Ensembl BioMart: `%s` (%s).' % (
            array_type,
            _array_type
        )
        _logger._log(msg)
        raise ValueError(msg)

    attrs = ['ensembl_peptide_id', _array_type]

    result = collections.defaultdict(set)

    for r in biomart_query(attrs = attrs, transcript = True, gene = True):

        array_probe_id = getattr(r, _array_type)

        if array_probe_id:

            for gene_attr in ('gene', 'transcript', 'peptide'):

                gene_id = getattr(r, 'ensembl_%s_id' % gene_attr)

                if gene_id:

                    result[gene_id].add(array_probe_id)

    return result
