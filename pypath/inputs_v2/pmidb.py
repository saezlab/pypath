"""
Parse PMI-DB data and emit Entity records.

This module converts protein-metabolite interaction information from PMI-DB into
Entity records using the declarative schema pattern.
"""

from __future__ import annotations

import pandas as pd

from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers.base import iter_tsv, _first_handle
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    CV,
    Member,
    MembershipBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    OntologyCv,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

BASE_URL = 'http://easybioai.com/PMIDB/static/%s.txt'
files = ['interaction', 'protein_infor', 'metabolites_infor']

config = ResourceConfig(
    id=ResourceCv.PMIDB,
    name='PMI-DB',
    url='http://easybioai.com/PMIDB/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='33554247',
    primary_category='interactions',
    description=(
        'PMIDB provides both interactions and non-interactions between protein '
        'and metabolite, which not only reduces the experimental cost for '
        'biological experimenters, but also facilitates the construction of '
        'more accurate algorithms for researchers using machine learning.'
    ),
)

download = {
    f: Download(
        url=BASE_URL % f,
        filename=f'{f}.txt',
        subfolder='PMIDB',
        large=True,
        ext='.txt',
        default_mode='r',
    )
    for f in files
}

def parser(opener, header=None, sep='\t'):

    if header is not None and iter(header):

        entries = [line.strip().split(sep) for line in _first_handle(opener)]

        if x := len(entries[0]) - len(header):

            raise KeyError(
                'Length of header is %s than the number of entries'
                % ('smaller' if x > 0 else 'larger')
            )

        yield from [
            {k: v for k, v in zip(header, line)}
            for line in entries
        ]

    else:

        return iter_tsv(opener)

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={},
    map={},
    transform={},
)

schema_interaction = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.UNIPROT,
                        value=f('Uniprot_KB_id')
                    ),
                ),
            )
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.CHEMICAL,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.KEGG_COMPOUND,
                        value=f('KEGG')
                    ),
                ),
            )
        ),
    )
)

# ================================= RESOURCE ===================================

#resource = Resource(
#    config=config,
#    data=Dataset(
#        download=download,
#        mapper=schema,
#        raw_parser=parser,
#    ),
#)

# ================================= REFERENCE ==================================
