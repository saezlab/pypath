"""Export the OmniPath ontology as an OBO artifact."""

from __future__ import annotations

import inspect
import re
import urllib.request

from pypath.internals import cv_terms
from pypath.internals.cv_terms import CvEnum, LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.inputs_v2.base import ArtifactDataset, Resource, ResourceConfig


PSI_MI_URL = 'https://raw.githubusercontent.com/HUPO-PSI/psi-mi-CV/master/psi-mi.obo'


config = ResourceConfig(
    id=ResourceCv.OMNIPATH_ONTOLOGY,
    name='OmniPath Ontology',
    url='https://omnipathdb.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    description='Combined PSI-MI and OmniPath controlled vocabulary ontology export.',
)


def _fetch_psi_mi() -> str:
    with urllib.request.urlopen(PSI_MI_URL) as response:
        return response.read().decode('utf-8')



def _fix_malformed_dates(content: str) -> str:
    content = re.sub(r'^date: \d{2}:\d{2}:\d{4}.*\n', '', content, flags=re.MULTILINE)
    content = re.sub(r'^creation_date:.*\n', '', content, flags=re.MULTILINE)
    return content



def _format_name(name: str) -> str:
    return name.lower().replace('_', ' ')



def _escape_obo_string(value: str) -> str:
    if not value:
        return ''
    return value.replace('\\', '\\\\').replace('"', '\\"')



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



def _format_om_terms(terms: list[dict]) -> str:
    lines = []
    for term in sorted(terms, key=lambda t: t['accession']):
        lines.append('[Term]')
        lines.append(f"id: {term['accession']}")
        lines.append(f"name: {term['name']}")
        if term['definition']:
            lines.append(f'def: "{_escape_obo_string(term["definition"])}" []')
        if term['is_a']:
            lines.append(f"is_a: {term['is_a']}")
        lines.append('')
    return '\n'.join(lines)



def render_omnipath_obo(_opener=None, **_kwargs) -> str:
    psi_mi = _fix_malformed_dates(_fetch_psi_mi())
    om_obo = _format_om_terms(_extract_om_terms())
    return psi_mi + '\n' + om_obo


resource = Resource(
    config,
    ontology=ArtifactDataset(
        renderer=render_omnipath_obo,
        extension='obo',
        file_stem='omnipath_mi',
    ),
)
