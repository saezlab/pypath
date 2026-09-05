"""MONDO ontology and gene-disease associations from ``mondo.obo``.

The ontology dataset exports MONDO terms and hierarchy as OBO. The annotations
 dataset emits explicit broad relations linking MONDO disease terms to HGNC gene
entities for MONDO relationships whose target is an HGNC gene identifier.
"""

from __future__ import annotations

from biolink_model.datamodel.model import Gene, OntologyClass, slots
from omnipath_core.naming import Namespace

import re
from collections.abc import Generator, Iterable
from typing import Any

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
)
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
GENE_DISEASE_RELATIONS = frozenset(
    {
        'has_material_basis_in_germline_mutation_in',
        'disease_has_basis_in_dysfunction_of',
        'disease_has_basis_in_disruption_of',
        'disease_causes_dysfunction_of',
        'RO:0004001',
        'RO:0004004',
    }
)

_HGNC_TARGET_RE = re.compile(
    r'^(?:https?://identifiers\.org/hgnc/|HGNC:)(\d+)$', re.I
)
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

    stanza: list[str] = []
    for raw_line in [*text.splitlines(), '[End]']:
        line = raw_line.strip()
        if line.startswith('['):
            if current is not None and not any(
                value == 'is_obsolete: true' for value in stanza
            ):
                for value in stanza:
                    tag, _, relation_text = value.partition(':')
                    if tag != 'relationship' and not (
                        include_logical_definitions and tag == 'intersection_of'
                    ):
                        continue
                    association = _parse_gene_relation(
                        relation_text.strip(), allowed_relations
                    )
                    if not association or not current.get('id'):
                        continue
                    key = (
                        current['id'],
                        association['hgnc_id'],
                        association['relation'],
                    )
                    if key not in seen:
                        seen.add(key)
                        yield {
                            'mondo_id': current['id'],
                            'mondo_name': current.get('name') or current['id'],
                            **association,
                        }
            current = {'id': '', 'name': ''} if line == '[Term]' else None
            stanza = []
        elif current is not None and line:
            stanza.append(line)
            tag, _, value = line.partition(':')
            if tag in ('id', 'name') and not current[tag]:
                current[tag] = value.strip()


def gene_disease_association_to_entity(row: dict[str, Any]) -> Relation | None:
    """Preserve explicit MONDO gene links without inventing a causal crosswalk."""
    if row.get('relation') not in GENE_DISEASE_RELATIONS or not row.get(
        'hgnc_id'
    ):
        return None
    gene_identifiers = [Identifier(type=Namespace.HGNC, value=row['hgnc_id'])]
    if row.get('gene_symbol'):
        gene_identifiers.append(
            Identifier(type=Namespace.GENESYMBOL, value=row['gene_symbol'])
        )
    # Preserve published source predicates as external attributes; symbolic
    # predicates and provenance strings remain in the parsed evidence payload.
    annotations = (
        [Annotation(term=row['relation'], value=f'HGNC:{row["hgnc_id"]}')]
        if row['relation'].startswith('RO:')
        else []
    )
    return Relation(
        subject=Entity(
            type=OntologyClass,
            identifiers=[
                Identifier(type=Namespace.MONDO, value=row['mondo_id']),
                *(
                    [Identifier(type=Namespace.NAME, value=row['mondo_name'])]
                    if row.get('mondo_name')
                    else []
                ),
            ],
        ),
        predicate=slots.associated_with,
        object=Entity(
            type=Gene,
            identifiers=gene_identifiers,
            annotations=[
                Annotation(term=slots.in_taxon, value='NCBITaxon:9606')
            ],
        ),
        annotations=annotations,
    )


def _parse_gene_relation(
    value: str, allowed_relations: set[str]
) -> dict[str, Any] | None:
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
        handle = (
            next(iter(opener.result.values()), None)
            if isinstance(opener.result, dict)
            else opener.result
        )
    if handle is None:
        return ''
    if hasattr(handle, 'seek'):
        handle.seek(0)
    content = handle.read() if hasattr(handle, 'read') else ''.join(handle)
    return (
        content.decode('utf-8', 'ignore')
        if isinstance(content, bytes)
        else str(content)
    )


def _dedupe(values: Iterable[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


terms_schema = ontology_entity_mapper(
    obo_record_to_term, ontology_id='mondo', identifier_type=Namespace.MONDO
)


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
