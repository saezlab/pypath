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

from typing import Collection, Dict, List, Optional, Union

import glom

import pypath.share.session as session
import pypath.share.common as common
import pypath.utils.taxonomy as taxonomy
import pypath.resources.urls as urls
import pypath.inputs.ebi as ebi
import pypath.inputs.common as inputs_common

_logger = session.Logger(name = 'proteins_input')

VARIATION_SOURCE_TYPES = (
    'large scale study',
    'mixed',
    'uniprot',
)

VARIATION_DB_TYPES = (
    '1000Genomes',
    'ClinVar',
    'cosmic curated',
    'dbSNP',
    'ESP',
    'ExAC',
    'gnomAD',
    'NCI-TCGA',
    'NCI-TCGA Cosmic',
    'TOPMed',
    'UniProt',
)

VARIATION_CONSEQUENCE_TYPES = (
    'frameshift',
    'inframe deletion',
    'insertion',
    'missense',
    'stop gained',
    'stop lost',
)


def variants(
        organism: Union[str, int] = 9606,
        qs: Optional[Dict] = None,
        sourcetype: Optional[Union[str, List[str]]] = ('uniprot', 'mixed'),
        disease: Optional[Union[str, List[str]]] = None,
        consequencetype: Optional[Union[str, List[str]]] = None,
        fields: Optional[inputs_common.GlomFields] = None,
        feature_fields: Optional[inputs_common.GlomFields] = None,
    ) -> List[tuple]:

    _logger._log('Retrieving variant data from from EBI Proteins.')

    ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism)

    if not ncbi_tax_id:

        msg = 'Failed to recognize organism `{}`.'.format(organism)
        _logger.log(msg)
        raise ValueError(msg)

    qs = qs or {}
    qs['taxid'] = ncbi_tax_id

    qs = dict(
        (
            k,
            ','.join(v) if isinstance(v, Collection) else v
        )
        for k, v in qs.items()
    )

    for param in ('sourcetype', 'consequencetype', 'disease'):

        value = locals()[param]

        if value:

            value = common.to_list(value)
            qs[param] = ','.join(value)

    url = urls.urls['proteins']['url'] % 'variation'

    _fields = {
        'uniprot': 'accession',
        'features': (
            'features',
            [{
                'type': 'type',
                'begin': ('begin', int),
                'end': ('end', int),
                'consequence': 'consequenceType',
                'wild_residue': 'wildType',
                'mutated_residue': (
                    glom.Coalesce('mutatedType', default = None)
                ),
                'somatic': ('somaticStatus', bool),
                'evidence': 'sourceType',
            }]
        ),
    }

    fields = inputs_common.glom_fields(fields)
    _fields.update(fields)

    feature_fields = inputs_common.glom_fields(feature_fields)
    _fields['features'][1][0].update(feature_fields)

    result = ebi.ebi_rest(
        url = url,
        qs = qs,
        fields = _fields,
        page_param = 'offset',
        size_param = 'size',
        by_page = False,
        paginate = True,
    )

    _logger._log('Finished retrieving variant data from from EBI Proteins.')

    return result
