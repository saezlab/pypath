"""OncoTree (MSK): the cancer type classification, with its hierarchy and cross-references."""

from __future__ import annotations

from typing import Any

from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, Topic

from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers.oncotree import iter_terms, oncotree_record_to_term
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity


config = ResourceConfig(
    id=ResourceCv.ONCOTREE,
    name='OncoTree',
    full='OncoTree cancer type classification',
    url='https://oncotree.mskcc.org/',
    # The web application (checked 2026-10-12) states no license for the data: it asks to cite
    # Kundra et al. 2021 and nothing more. The source repository https://github.com/cBioPortal/oncotree
    # does: its README says "This work is licensed under a Creative Commons Attribution 4.0
    # International License" (http://creativecommons.org/licenses/by/4.0/) and its LICENSE.md is
    # the CC BY 4.0 legal code (GitHub: CC-BY-4.0). This is the repository's license, taken
    # for the OncoTree data the repository holds.
    license=LicenseCV.CC_BY_4_0,
    # Stable releases twice a year (2025-04-08, 2025-10-03); the release read is named in the URL.
    update_category=UpdateCategoryCV.REGULAR,
    # Kundra et al., "OncoTree: A Cancer Classification System for Precision Oncology",
    # JCO Clinical Cancer Informatics 2021 (doi 10.1200/CCI.20.00108).
    pubmed='33625877',
    topics=(Topic.DISEASE_PHENOTYPE,),
    content=(Content.ONTOLOGIES,),
    summary='OncoTree: the cancer types of Memorial Sloan Kettering by code, with their hierarchy, main type, tissue and cross-references.',
    description=(
        'OncoTree is the cancer classification of Memorial Sloan Kettering Cancer Center, '
        'maintained with the cBioPortal community. It arranges 897 tumor types of the release '
        'of 2025-10-03 in a tree of tissues, main types and subtypes, each with a short code '
        '(BRCA, LUAD, ANGS), a name, the main type and tissue it belongs to, cross-references '
        'to the UMLS and the NCI Thesaurus, and the codes it replaced in earlier releases. '
        'Resources cite its codes for the cancer type of a cohort or a sample (IntOGen); this '
        'module provides their names and hierarchy.'
    ),
)

ONTOLOGY_ID = 'oncotree'
MAIN_TYPE = 'oncotree:main_type'
TISSUE = 'oncotree:tissue'

# The codes are bare (BRCA), in no prefix, so every identifier is in this namespace.
_term_entity = ontology_entity_mapper(
    oncotree_record_to_term,
    ontology_id=ONTOLOGY_ID,
    identifier_type=Namespace.ONCOTREE,
)


def map_term(row: dict[str, Any]) -> Entity | None:
    """A tumor type: its code, name, previous codes, cross-references, parent, main type and tissue."""
    entity = _term_entity(row)
    if entity is None:
        return None
    attributes = [
        Annotation(term=term, value=row[key])
        for term, key in ((MAIN_TYPE, 'main_type'), (TISSUE, 'tissue'))
        if row.get(key)
    ]
    return entity._replace(annotations=[*(entity.annotations or []), *attributes] or None)


# The release of 2025-10-03, the latest stable one (the API lists it as
# ``oncotree_latest_stable`` and, under its own name, ``oncotree_2025_10_03``: the same
# content). The tree has the root, TISSUE, which the flat list of tumor types lacks.
RELEASE = 'oncotree_2025_10_03'

resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url=f'https://oncotree.mskcc.org/api/tumorTypes/tree?version={RELEASE}',
            filename=f'{RELEASE}_tree.json',
            subfolder='oncotree',
            large=False,
            ext='json',
        ),
        mapper=map_term,
        raw_parser=iter_terms,
    ),
)
