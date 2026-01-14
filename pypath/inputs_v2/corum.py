"""
Parse CORUM data and emit Entity records.

This module converts CORUM protein complex data into Entity records using the
declarative schema pattern defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembershipBuilder,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.corum import _raw


config = ResourceConfig(
    id=ResourceCv.CORUM,
    name='CORUM',
    url='https://mips.helmholtz-muenchen.de/corum/',
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='30357367',
    description=(
        'CORUM is a manually curated repository of experimentally characterized '
        'protein complexes from mammalian organisms. Each complex includes '
        'information about subunit composition, functional annotations, and '
        'literature references.'
    ),
)

f = FieldConfig(
    delimiter=';',
)

complexes_schema = EntityBuilder(
    entity_type=EntityTypeCv.COMPLEX,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CORUM, value=f('ComplexID')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('ComplexName')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('PubMed ID', delimiter=';')),
        CV(term=MoleculeAnnotationsCv.FUNCAT, value=f('FunCat description', delimiter=';')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('GO ID', delimiter=';')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(
                    term=IdentifierNamespaceCv.UNIPROT,
                    value=f('subunits(UniProt IDs)', delimiter=';'),
                ),
            ),
        )
    ),
)

resource = Resource(
    config,
    complexes=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/CORUM_allComplexes.txt.zip',
            filename='CORUM_allComplexes.txt.zip',
            subfolder='corum',
            large=True,
            ext='zip',
        ),
        mapper=complexes_schema,
        raw_parser=_raw,
    ),
)
