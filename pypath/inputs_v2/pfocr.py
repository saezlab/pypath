"""PFOCR figure mentions of genes and chemicals.

A figure is an InformationContentEntity; mentions never assert a molecular
interaction between co-mentioned entities.
"""

from __future__ import annotations

import re

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.pfocr import (
    current_pfocr_filename,
    current_pfocr_url,
    iter_gmt,
)


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


def association_to_entity(row: dict) -> Relation | None:
    """A figure mentions a gene or chemical; co-mention is not biological interaction."""
    data_type = row.get('data_type') or 'gene'
    if data_type not in {'gene', 'chemical'}:
        raise ValueError(f'Unsupported PFOCR data type: {data_type}')
    identifier = str(row.get('identifier') or '').strip()
    if data_type == 'chemical':
        match = _CHEBI_RE.fullmatch(identifier)
        if not match:
            return None
        identifier = match.group(1)
    elif not identifier.isdigit():
        return None
    if not row.get('pfocr_id'):
        return None
    annotations = []
    taxon = str(row.get('taxonomy_id') or '')
    if data_type == 'gene' and taxon.isdigit() and int(taxon) > 0:
        annotations.append(
            Annotation(term=slots.in_taxon, value='NCBITaxon:' + taxon)
        )
    figure_annotations = []
    pmc_match = _PMC_RE.match(row['pfocr_id'])
    if pmc_match:
        figure_annotations.append(
            Annotation(
                term=slots.publications,
                value='PMC:' + pmc_match.group(1).removeprefix('PMC'),
            )
        )
    return Relation(
        subject=Entity(
            type=model.InformationContentEntity,
            identifiers=[
                Identifier(type=Namespace.PFOCR, value=row['pfocr_id']),
                Identifier(
                    type=Namespace.NAME,
                    value=row.get('caption') or row['pfocr_id'],
                ),
            ],
            annotations=figure_annotations,
        ),
        predicate=slots.mentions,
        object=Entity(
            type=model.ChemicalEntity
            if data_type == 'chemical'
            else model.Gene,
            identifiers=[
                Identifier(
                    type=Namespace.CHEBI
                    if data_type == 'chemical'
                    else Namespace.ENTREZ,
                    value=identifier,
                )
            ],
            annotations=annotations,
        ),
    )


def _download(data_type: str, species: str = 'Homo_sapiens') -> Download:
    return Download(
        url=lambda **_kwargs: current_pfocr_url(
            data_type=data_type, species=species
        ),
        filename=lambda **_kwargs: current_pfocr_filename(
            data_type=data_type, species=species
        ),
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
        raw_parser=lambda opener, **kwargs: iter_gmt(
            opener, data_type='gene', species='Homo_sapiens', **kwargs
        ),
    ),
    human_chemical_associations=Dataset(
        download=_download('chemical'),
        mapper=association_to_entity,
        raw_parser=lambda opener, **kwargs: iter_gmt(
            opener, data_type='chemical', species='Homo_sapiens', **kwargs
        ),
    ),
)
