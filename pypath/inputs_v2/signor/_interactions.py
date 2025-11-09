#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
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

"""
Parse SIGNOR interaction data and emit Entity records.
"""

from __future__ import annotations

import csv
from collections.abc import Generator

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity as SilverEntity
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    ReferenceTypeCv,
)
from ..tabular_builder import (
    Annotations,
    Column,
    Entities,
    Entity as EntitySchema,
    Identifiers,
    Members,
)

__all__ = ['signor_interactions']


def signor_interactions() -> Generator[SilverEntity, None, None]:
    """
    Download SIGNOR causalTab interactions and yield `SilverEntity` objects.
    """

    opener = download_and_open(
        url='https://signor.uniroma2.it/download_entity.php',
        filename='signor_all_causalTab.txt',
        subfolder='signor',
        query={
            'format': 'causalTab',
            'submit': 'Download',
        },
        post=True,
        large=True,
    )

    mi_processing = {'extract_value': r'(MI:\d+)'}
    signor_processing = {
        'extract_prefix': r'^([^:]+):',
        'extract_value': r'^[^:]+:(.*)',
    }
    uniprot_processing = {'extract_value': r'uniprotkb:([^|"]+)'}
    tax_processing = {'extract_value': r'taxid:([-\d]+)'}
    pubmed_processing = {'extract_value': r'(?i)pubmed:(\d+)'}

    schema = EntitySchema(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=Identifiers(
            Column(
                13,
                delimiter='|',
                processing=signor_processing,
                cv={
                    'signor': IdentifierNamespaceCv.SIGNOR,
                    'signor-interaction': IdentifierNamespaceCv.SIGNOR,
                },
            ),
        ),
        annotations=Annotations(
            Column(11, delimiter='|', processing=mi_processing, cv=lambda value: value.value),
            Column(6, delimiter='|', processing=mi_processing, cv=lambda value: value.value),
            Column(45, delimiter='|', processing=mi_processing, cv=lambda value: value.value),
            Column(44, delimiter='|', processing=mi_processing, cv=lambda value: value.value),
            Column(8, delimiter='|', processing=pubmed_processing, cv=ReferenceTypeCv.PUBMED),
        ),
        members=Members(
            Entities(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=Identifiers(
                    Column(0, delimiter='|', processing=uniprot_processing, cv=IdentifierNamespaceCv.UNIPROT),
                    Column(2, delimiter='|', processing=uniprot_processing, cv=IdentifierNamespaceCv.UNIPROT),
                ),
                annotations=Annotations(
                    Column(16, delimiter='|', processing=mi_processing, cv=lambda value: value.value),
                    Column(18, delimiter='|', processing=mi_processing, cv=lambda value: value.value),
                    Column(9, delimiter='|', processing=tax_processing, cv=IdentifierNamespaceCv.NCBI_TAX_ID),
                ),
            ),
            Entities(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=Identifiers(
                    Column(1, delimiter='|', processing=uniprot_processing, cv=IdentifierNamespaceCv.UNIPROT),
                    Column(3, delimiter='|', processing=uniprot_processing, cv=IdentifierNamespaceCv.UNIPROT),
                ),
                annotations=Annotations(
                    Column(17, delimiter='|', processing=mi_processing, cv=lambda value: value.value),
                    Column(19, delimiter='|', processing=mi_processing, cv=lambda value: value.value),
                    Column(10, delimiter='|', processing=tax_processing, cv=IdentifierNamespaceCv.NCBI_TAX_ID),
                ),
            ),
        ),
    )

    reader = csv.DictReader(opener.result, delimiter='\t')
    for entity in reader:
        yield schema(entity)
