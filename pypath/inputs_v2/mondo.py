"""MONDO ontology and gene-disease associations from ``mondo.obo``.

The ontology dataset exports MONDO terms and hierarchy as OBO. The annotations
 dataset emits association entities linking MONDO disease terms to HGNC gene
entities for MONDO relationships whose target is an HGNC gene identifier.
"""

from __future__ import annotations

import re
from collections.abc import Generator, Iterable
from typing import Any

from pypath.internals.cv_terms import (
    DiseaseAnnotationCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Membership
from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers.obo import iter_obo, obo_record_to_term


_MONDO_DOWNLOAD = Download(
    url='https://purl.obolibrary.org/obo/mondo.obo',
    filename='mondo.obo',
    subfolder='mondo',
    large=True,
)

# Gene-disease relation predicates observed in MONDO OBO records with HGNC
# targets. ``RO:0004001`` is gain-of-function germline mutation; ``RO:0004004``
# is somatic mutation. The symbolic predicates are used verbatim in MONDO.
GENE_DISEASE_RELATIONS = frozenset({
    'has_material_basis_in_germline_mutation_in',
    'disease_has_basis_in_dysfunction_of',
    'disease_has_basis_in_disruption_of',
    'disease_causes_dysfunction_of',
    'RO:0004001',
    'RO:0004004',
})

_HGNC_TARGET_RE = re.compile(r'^(?:https?://identifiers\.org/hgnc/|HGNC:)(\d+)$', re.I)
_SOURCE_RE = re.compile(r'\bsource="([^"]+)"')


config = ResourceConfig(
    id=ResourceCv.MONDO,
    name='Mondo Disease Ontology',
    url='https://mondo.monarchinitiative.org/',
    pubmed='41052288',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    primary_category='diseases',
    resource_kind='ontology',
    annotation_ontologies=(OntologyCv.MONDO,),
    description=(
        'The Mondo Disease Ontology harmonizes disease concepts across sources. '
        'This module provides the MONDO OBO ontology and gene-disease '
        'associations encoded as MONDO relationships to HGNC genes.'
    ),
)


def iter_gene_disease_associations(
    opener,
    *,
    relations: Iterable[str] | None = None,
    include_logical_definitions: bool = True,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield MONDO disease-to-HGNC-gene relationship records.

    By default this includes the explicit ``relationship`` lines and also
    ``intersection_of`` logical-definition lines. Duplicate disease/gene/relation
    tuples are removed, so logical definitions do not duplicate explicit
    relationships.
    """
    allowed_relations = set(relations or GENE_DISEASE_RELATIONS)
    text = _read_text(opener)
    current: dict[str, Any] | None = None
    seen: set[tuple[str, str, str]] = set()

    for raw_line in text.splitlines():
        line = raw_line.strip()

        if line == '[Term]':
            current = {'id': '', 'name': ''}
            continue
        if line.startswith('['):
            current = None
            continue
        if current is None or not line:
            continue

        tag, sep, value = line.partition(':')
        if not sep:
            continue
        tag = tag.strip()
        value = value.strip()

        if tag == 'id' and not current['id']:
            current['id'] = value
        elif tag == 'name' and not current['name']:
            current['name'] = value
        elif tag == 'relationship' or (include_logical_definitions and tag == 'intersection_of'):
            association = _parse_gene_relation(value, allowed_relations)
            if not association or not current.get('id'):
                continue

            key = (current['id'], association['hgnc_id'], association['relation'])
            if key in seen:
                continue
            seen.add(key)

            yield {
                'mondo_id': current['id'],
                'mondo_name': current.get('name') or current['id'],
                **association,
            }


def gene_disease_association_to_entity(row: dict[str, Any]) -> Entity:
    """Map a MONDO-HGNC association to an association entity."""
    gene_identifiers = [
        Identifier(type=IdentifierNamespaceCv.HGNC, value=row['hgnc_id']),
    ]
    if row.get('gene_symbol'):
        gene_identifiers.append(
            Identifier(type=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=row['gene_symbol'])
        )

    annotations = [Annotation(term=DiseaseAnnotationCv.RELATIONSHIP, value=row['relation'])]
    for source in row.get('sources') or []:
        annotations.append(Annotation(term=DiseaseAnnotationCv.SOURCE, value=source))

    return Entity(
        type=EntityTypeCv.ASSOCIATION,
        identifiers=[],
        annotations=annotations,
        membership=[
            Membership(
                member=Entity(
                    type=EntityTypeCv.CV_TERM,
                    identifiers=[
                        Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=row['mondo_id']),
                    ],
                ),
            ),
            Membership(
                member=Entity(
                    type=EntityTypeCv.PROTEIN,
                    identifiers=gene_identifiers,
                    annotations=[Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value='9606')],
                ),
            ),
        ],
    )


def _parse_gene_relation(value: str, allowed_relations: set[str]) -> dict[str, Any] | None:
    base, _, target_name = value.partition(' ! ')
    parts = base.strip().split(None, 2)
    if len(parts) < 2:
        return None

    relation, target = parts[0], parts[1]
    if relation not in allowed_relations:
        return None

    hgnc_match = _HGNC_TARGET_RE.match(target)
    if not hgnc_match:
        return None

    rest = parts[2] if len(parts) > 2 else ''
    return {
        'relation': relation,
        'hgnc_id': hgnc_match.group(1),
        'gene_symbol': target_name.strip() or None,
        'sources': _dedupe(_SOURCE_RE.findall(rest)),
    }


def _read_text(opener) -> str:
    handle = None
    if opener and getattr(opener, 'result', None):
        handle = next(iter(opener.result.values()), None) if isinstance(opener.result, dict) else opener.result
    if handle is None:
        return ''
    if hasattr(handle, 'seek'):
        handle.seek(0)
    content = handle.read() if hasattr(handle, 'read') else ''.join(handle)
    return content.decode('utf-8', 'ignore') if isinstance(content, bytes) else str(content)


def _dedupe(values: Iterable[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


terms_schema = ontology_entity_mapper(obo_record_to_term, ontology_id='mondo')


resource = Resource(
    config,
    terms=Dataset(
        download=_MONDO_DOWNLOAD,
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
    annotations=Dataset(
        download=_MONDO_DOWNLOAD,
        mapper=gene_disease_association_to_entity,
        raw_parser=iter_gene_disease_associations,
    ),
)
