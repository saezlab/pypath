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

from future.utils import iteritems

import re
import os
import json
import collections

import pypath.share.session as session_mod
import pypath.share.common as common
import pypath_common.data as _data
import pypath.share.curl as curl
import pypath.share.settings as settings
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy

_logger = session_mod.Logger(name = 'biomart_input')


# for mouse homologues: Filter name = "with_mmusculus_homolog"
_filter_xml_template = '<Filter name="%s" excluded="0"/>'
_attr_xml_template = '<Attribute name="%s" />'


def biomart_query(
        attrs: str | list[str],
        filters: str | list = None,
        transcript: bool = False,
        peptide: bool = False,
        gene: bool = False,
        dataset: str = 'hsapiens_gene_ensembl',
    ):
    """
    Query the Ensembl Biomart web service.
    Use https://www.ensembl.org/biomart/martview/ to check for attribute
    and dataset names.

    Args
        attrs:
            One or more Ensembl attribute names.
        filters:
            One or more Ensembl filter names.
        transcript:
            Include Ensembl transcript IDs in the result.
        peptide:
            Include Ensembl peptide IDs in the result.
        gene:
            Include Ensembl gene IDs in the result.
        dataset:
            An Ensembl dataset name.

    Yields:
        Named tuples with the requested attributes for each record returned
        by Ensembl Biomart.
    """

    _attrs = []

    if gene:

        _attrs.append('ensembl_gene_id')

    if transcript:

        _attrs.append('ensembl_transcript_id')

    if peptide:

        _attrs.append('ensembl_peptide_id')

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

    xml_template_path = _data.path('ensembl_biomart_query.xml')

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
    c = curl.Curl(
        biomart_url,
        req_headers = [settings.get('user_agent')],
        large = True,
        silent = False,
    )
    success = False

    for line in c.result:

        _line = line.strip('\n\r').split('\t')

        if _line[0] == '[success]':

            success = True
            continue

        if line.strip() and len(_line) == len(record._fields):

            yield record(*_line)

    if not success:

        _logger._log(
            'Error: Interrupted transfer while downlading data '
            'from Ensembl Biomart (missing `success` tag).'
        )


def biomart_homology(
        source_organism: int | str = 9606,
        target_organism: int | str = 10090,
        extra_fields: str | Iterable[str] = 'external_gene_name',
    ):
    """
    Retrieves orthology data from Ensembl Biomart.

    Returns
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

    homolog_attrs = [
        'homolog_ensembl_peptide',
        'homolog_ensembl_gene',
        'homolog_orthology_type',
        'homolog_orthology_confidence',
        'homolog_canonical_transcript_protein',
        'homolog_associated_gene_name',
    ]

    homolog_attrs = [
        '%s_%s' % (target_organism, attr)
        for attr in homolog_attrs
    ] + common.to_list(extra_fields)

    filters = [
        'with_%s_homolog' % target_organism
    ]

    return list(
        biomart_query(
            attrs = homolog_attrs,
            filters = filters,
            transcript = True,
            gene = True,
            peptide = True,
            dataset = '%s_gene_ensembl' % source_organism,
        )
    )



def biomart_microarray_types(organism: int | str = 9606):
    """
    Retrieves a list of available microarray types for an organism.

    Args
        organism:
            Name or ID of an organism.
    """

    organism = taxonomy.ensure_ensembl_name(organism)

    url = urls.urls['ensembl']['arraytypes'] % organism
    c = curl.Curl(url, req_headers = [settings.get('user_agent')])
    result = json.loads(c.result)

    _ = [
        r.update(
            label = '%s %s' % (
                r['vendor'],
                re.sub('[-_]', ' ', r['array']),
            )
        )
        for r in result
    ]

    return result


def biomart_microarray(
        array_type: str,
        gene: bool = True,
        transcript: bool = False,
        peptide: bool = False,
        organism: int | str = 9606
    ):
    """
    Microarray probe identifier mappings.

    Args
        array_type:
            The microarray model, as shown on the BioMart
            webpage, or the corresponding code. For a full list of available
            identifiers see the ``biomart_microarray_types``.
        gene:
            Include the mapping to Ensembl gene IDs.
        transcript:
            Include the mapping to Ensembl transcript IDs.
        peptide:
            Include the mapping to Ensembl peptide IDs.
        organism:
            Name or ID of an organism.

    Returns
        A dictionary with Ensembl gene, transcript and peptide IDs as keys
        and sets of microarray probe IDs as values.
    """

    organism = taxonomy.ensure_ensembl_name(organism)
    array_types = biomart_microarray_types(organism = organism)
    array_types = {
        at['label'].lower().replace(' ', '_')
        for at in array_types
    }

    _array_type = (
        array_type.
        lower().
        replace('probe', '').
        strip().
        replace(' ', '_')
    )

    if _array_type not in array_types:

        msg = 'No such array type in Ensembl BioMart: `%s` (%s).' % (
            array_type,
            _array_type
        )
        _logger._log(msg)
        raise ValueError(msg)

    attrs = [_array_type]
    dataset = '%s_gene_ensembl' % organism

    biomart_result = biomart_query(
        attrs = attrs,
        transcript = transcript,
        gene = gene,
        peptide = peptide,
        dataset = dataset,
    )
    result = collections.defaultdict(set)
    _locals = locals()
    ensembl_attrs = tuple(
        attr
        for attr in ('gene', 'transcript', 'peptide')
        if _locals[attr]
    )

    for r in biomart_result:

        array_probe_id = getattr(r, _array_type)

        if array_probe_id:

            for ensembl_attr in ensembl_attrs:

                ensembl_id = getattr(r, 'ensembl_%s_id' % ensembl_attr)

                if ensembl_id:

                    result[ensembl_id].add(array_probe_id)

    return dict(result)


def biomart_microarrays(
        organism: int | str = 9606,
        vendor: str | set[str] = None,
        gene: bool = True,
        transcript: bool = False,
        peptide: bool = False
    ):
    """
    Microarray probe identifier mappings for multiple microarrays.
    Retrieves probe mappings for all array types for one organism,
    optionally limited to one or more array vendors. Note: depending
    on the number of array models, it can take minutes to download
    the data.

    Args
        organism:
            Name or ID of an organism.
        vendor:
            One or more vendors. None means all vendors. For
            human, possible values are AFFY, ILLUMINA, AGILENT, CODELINK
            and PHALANX.
        gene:
            Include the mapping to Ensembl gene IDs.
        transcript:
            Include the mapping to Ensembl transcript IDs.
        peptide:
            Include the mapping to Ensembl peptide IDs.

    Returns
        A dictionary with Ensembl gene, transcript and peptide IDs as keys
        and sets of tuples with microarray types and probe IDs as values.
    """

    record = collections.namedtuple(
        'Probe',
        (
            'array',
            'probe',
        )
    )

    array_types = biomart_microarray_types(organism = organism)
    vendor = {v.upper() for v in common.to_set(vendor)}

    result = collections.defaultdict(set)

    for at in array_types:

        if not vendor or at['vendor'] in vendor:

            probe_map = biomart_microarray(
                array_type = at['label'],
                organism = organism,
                gene = gene,
                transcript = transcript,
                peptide = peptide,
            )

            for gene_id, probe_ids in iteritems(probe_map):

                for probe_id in probe_ids:

                    result[gene_id].add(
                        record(
                            array = at['label'].lower().replace(' ', '_'),
                            probe = probe_id,
                        )
                    )

    return dict(result)
