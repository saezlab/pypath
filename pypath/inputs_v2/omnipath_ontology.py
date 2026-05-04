"""Export the OmniPath ontology as an OBO artifact."""

from __future__ import annotations

import inspect
from collections.abc import Generator
from typing import Any

from pypath.internals import cv_terms
from pypath.internals.cv_terms import CvEnum, LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.ontology_schema import OntologyDocument, OntologyTerm
from pypath.inputs_v2.base import OntologyDataset, Resource, ResourceConfig


config = ResourceConfig(
    id=ResourceCv.OMNIPATH_ONTOLOGY,
    name='OmniPath Ontology',
    url='https://omnipathdb.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    primary_category='ontologies',
    resource_kind='ontology',
    description='OmniPath controlled vocabulary ontology export.',
)


def _format_name(name: str) -> str:
    return name.lower().replace('_', ' ')



def _extract_om_terms() -> list[dict]:
    terms = {}
    parent_terms = {}

    for _, obj in inspect.getmembers(cv_terms):
        if not (inspect.isclass(obj) and issubclass(obj, CvEnum) and obj is not CvEnum):
            continue

        parent_accession = None
        if hasattr(obj, 'parent_cv_term') and obj.parent_cv_term:
            parent_term = obj.parent_cv_term
            if isinstance(parent_term, tuple):
                parent_accession = parent_term[0]
                if parent_accession.startswith('OM:') and parent_accession not in parent_terms:
                    parent_name = parent_term[1] if len(parent_term) > 1 else ''
                    parent_def = parent_term[2] if len(parent_term) > 2 else ''
                    parent_terms[parent_accession] = {
                        'accession': parent_accession,
                        'name': _format_name(parent_name) if parent_name else parent_accession,
                        'definition': parent_def,
                        'is_a': 'MI:0000',
                    }
            elif isinstance(parent_term, str):
                parent_accession = parent_term

        for member in obj:
            accession = member.value
            if not accession.startswith('OM:'):
                continue
            terms[accession] = {
                'accession': accession,
                'name': _format_name(member.name),
                'definition': getattr(member, 'definition', None) or '',
                'is_a': parent_accession if parent_accession else 'MI:0000',
            }

    for parent_acc, parent_data in parent_terms.items():
        if parent_acc not in terms:
            terms[parent_acc] = parent_data

    return list(terms.values())



def _iter_om_terms(_opener=None, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    yield from _extract_om_terms()


def _map_om_term(row: dict[str, Any]) -> OntologyTerm:
    return OntologyTerm(
        id=row['accession'],
        name=row['name'],
        definition=row.get('definition') or None,
        is_a=[row['is_a']] if row.get('is_a') else None,
    )


resource = Resource(
    config,
    ontology=OntologyDataset(
        download=None,
        mapper=_map_om_term,
        raw_parser=_iter_om_terms,
        document=OntologyDocument(ontology='omnipath'),
        extension='obo',
        file_stem='omnipath_mi',
    ),
)
