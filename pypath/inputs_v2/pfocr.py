"""PFOCR figure gene and chemical associations from WikiPathways GMT exports."""

from __future__ import annotations

import re

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Membership
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.pfocr import current_pfocr_filename, current_pfocr_url, iter_gmt


config = ResourceConfig(
    id=ResourceCv.PFOCR,
    name='PFOCR',
    url='https://pfocr.wikipathways.org/',
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    primary_category='pathways',
    pubmed='38007419',
    description=(
        'PFOCR (Pathway Figure OCR) extracts genes and chemicals from pathway '
        'figures in the biomedical literature. This module parses the current '
        'WikiPathways PFOCR GMT exports.'
    ),
)

_CHEBI_RE = re.compile(r'^(?:CHEBI:|chebi:)?(\d+)$')
_PMC_RE = re.compile(r'^(PMC\d+)__')


def association_to_entity(row: dict) -> Entity:
    """Map one PFOCR GMT row to a figure-level association entity."""
    data_type = row.get('data_type') or 'gene'
    member_type = EntityTypeCv.SMALL_MOLECULE if data_type == 'chemical' else EntityTypeCv.PROTEIN
    identifier_type = IdentifierNamespaceCv.CHEBI if data_type == 'chemical' else IdentifierNamespaceCv.ENTREZ

    annotations = [
        Annotation(term=IdentifierNamespaceCv.PFOCR, value=row.get('pfocr_id')),
        Annotation(term=IdentifierNamespaceCv.NAME, value=row.get('caption')),
    ]
    if row.get('taxonomy_id'):
        annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=row['taxonomy_id']))
    pmc_match = _PMC_RE.match(row.get('pfocr_id') or '')
    if pmc_match:
        annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED_CENTRAL, value=pmc_match.group(1)))

    return Entity(
        type=EntityTypeCv.ASSOCIATION,
        identifiers=[],
        annotations=annotations,
        membership=[
            Membership(
                member=Entity(
                    type=EntityTypeCv.CV_TERM,
                    identifiers=[
                        Identifier(type=IdentifierNamespaceCv.PFOCR, value=row['pfocr_id']),
                        Identifier(type=IdentifierNamespaceCv.NAME, value=row.get('caption') or row['pfocr_id']),
                    ],
                ),
            ),
            Membership(
                member=Entity(
                    type=member_type,
                    identifiers=_member_identifiers(row['identifier'], identifier_type),
                    annotations=_member_annotations(identifier_type, row.get('taxonomy_id')),
                ),
            ),
        ],
    )


def _member_identifiers(identifier: str, identifier_type) -> list[Identifier]:
    identifier = identifier.strip()
    if identifier_type is IdentifierNamespaceCv.CHEBI:
        match = _CHEBI_RE.match(identifier)
        value = match.group(1) if match else identifier
        return [Identifier(type=IdentifierNamespaceCv.CHEBI, value=value)]

    return [Identifier(type=IdentifierNamespaceCv.ENTREZ, value=identifier)]


def _member_annotations(identifier_type, taxonomy_id: str | None) -> list[Annotation] | None:
    if identifier_type is IdentifierNamespaceCv.ENTREZ and taxonomy_id:
        return [Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxonomy_id)]
    return None


def _download(data_type: str, species: str = 'Homo_sapiens') -> Download:
    return Download(
        url=lambda **_kwargs: current_pfocr_url(data_type=data_type, species=species),
        filename=lambda **_kwargs: current_pfocr_filename(data_type=data_type, species=species),
        subfolder='pfocr',
        large=True,
        encoding='utf-8',
        default_mode='r',
    )


resource = Resource(
    config,
    human_gene_associations=Dataset(
        download=_download('gene'),
        mapper=association_to_entity,
        raw_parser=lambda opener, **kwargs: iter_gmt(opener, data_type='gene', species='Homo_sapiens', **kwargs),
    ),
    human_chemical_associations=Dataset(
        download=_download('chemical'),
        mapper=association_to_entity,
        raw_parser=lambda opener, **kwargs: iter_gmt(opener, data_type='chemical', species='Homo_sapiens', **kwargs),
    ),
)
