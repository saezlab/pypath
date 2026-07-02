"""
Parse CellChatDB data and emit Entity records.

This module converts CellChatDB ligand-receptor interactions, complexes and
cofactor groups into Entity records.
"""

from __future__ import annotations

import functools
import re
from typing import Any

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    InterCellAnnotations,
    InteractionMetadataCv,
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Membership
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cellchat import (
    iter_cellchat_cofactors,
    iter_cellchat_cofactor_interactions,
    iter_cellchat_complexes,
    iter_cellchat_interactions,
)


config = ResourceConfig(
    id=ResourceCv.CELLCHAT,
    name='CellChatDB',
    url='https://github.com/jinworks/CellChat',
    license=LicenseCV.GPL_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33597522',
    primary_category='interactions',
    description=(
        'CellChatDB is a manually curated database of literature-supported '
        'cell-cell communication interactions for CellChat, including '
        'ligand-receptor pairs, multi-subunit complexes and signaling '
        'cofactors across multiple species.'
    ),
)


BASE_URL = 'https://raw.githubusercontent.com/jinworks/CellChat/main/data/'

SPECIES = {
    'human': '9606',
    'mouse': '10090',
    'zebrafish': '7955',
}


def _download(species: str) -> Download:
    return Download(
        url=f'{BASE_URL}CellChatDB.{species}.rda',
        filename=f'cellchatdb_{species}.rda',
        subfolder='cellchat',
        default_mode='rb',
        encoding=None,
    )


def _parser(parser: Any, species: str) -> Any:
    return functools.partial(parser, taxon_id=SPECIES[species])


_pmid_re = re.compile(r'PMID:?\s*(\d+)', re.IGNORECASE)
_pmc_re = re.compile(r'\bPMC(?:ID)?\s*:?\s*(\d+)', re.IGNORECASE)


def _values(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value).strip()
    return [text] if text else []


def _first(value: Any) -> str:
    vals = _values(value)
    return vals[0] if vals else ''


def _annotations(*items: Annotation | None) -> list[Annotation] | None:
    out: list[Annotation] = []
    seen: set[tuple[Any, Any, Any]] = set()
    for item in items:
        if item is None:
            continue
        key = (item.term, item.value, item.units)
        if key not in seen:
            out.append(item)
            seen.add(key)
    return out or None


def _annotation(term: Any, value: Any = None) -> Annotation | None:
    if value is None:
        return None
    value = str(value).strip()
    return Annotation(term=term, value=value) if value else None


def _neurotransmitter_annotation(value: Any) -> Annotation | None:
    if str(value or '').strip().upper() != 'TRUE':
        return None
    return Annotation(term=InteractionMetadataCv.NEUROTRANSMITTER_INTERACTION)


def _location_terms(location: str | None) -> list[Annotation]:
    location_lower = (location or '').lower()
    terms: list[Annotation] = []
    if 'secreted' in location_lower:
        terms.append(Annotation(term=InterCellAnnotations.SECRETED))
    if 'membrane' in location_lower:
        terms.append(Annotation(term=InterCellAnnotations.MEMBRANE))
    if 'cytoplasm' in location_lower:
        terms.append(Annotation(term=InterCellAnnotations.CYTOPLASM))
    return terms


def _pubmed_annotations(evidence: str) -> list[Annotation]:
    return [
        Annotation(term=IdentifierNamespaceCv.PUBMED, value=pmid)
        for pmid in sorted(set(_pmid_re.findall(evidence or '')))
    ]


def _pmc_annotations(evidence: str) -> list[Annotation]:
    return [
        Annotation(term=IdentifierNamespaceCv.PUBMED_CENTRAL, value=f'PMC{pmc}')
        for pmc in sorted(set(_pmc_re.findall(evidence or '')))
    ]


def _protein_entity(gene: str, taxon_id: str) -> Entity:
    return Entity(
        type=EntityTypeCv.PROTEIN,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=gene),
        ],
        annotations=[
            Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
        ],
    )


def _normalize_fallback_gene(gene: str, taxon_id: str) -> str:
    gene = gene.strip()
    return gene.lower() if taxon_id == SPECIES['zebrafish'] else gene


def _complex_genes_from_name(name: str, taxon_id: str) -> list[str]:
    parts = [part.strip() for part in name.split('_') if part.strip()]
    if len(parts) < 2:
        return []
    if parts[-1].isdigit():
        parts = parts[:-1]
    if len(parts) < 2:
        return []
    return [_normalize_fallback_gene(part, taxon_id) for part in parts]


def _is_named_group(name: str) -> bool:
    return bool(re.search(r'\s', name.strip()))


def _protein_group_entity(
    *,
    name: str,
    genes: list[str],
    taxon_id: str,
    annotations: list[Annotation] | None = None,
) -> Entity:
    genes = [gene for gene in genes if gene]
    annotations = _annotations(
        Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
        *(annotations or []),
    )

    if not genes:
        genes = _complex_genes_from_name(name, taxon_id)

    if len(genes) == 1:
        identifiers = [
            Identifier(type=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=genes[0]),
        ]
        if name and name != genes[0]:
            identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=name))
        return Entity(
            type=EntityTypeCv.PROTEIN,
            identifiers=identifiers,
            annotations=annotations,
        )

    if genes:
        return Entity(
            type=EntityTypeCv.COMPLEX,
            identifiers=[Identifier(type=IdentifierNamespaceCv.NAME, value=name)],
            annotations=annotations,
            membership=[
                Membership(member=_protein_entity(gene, taxon_id))
                for gene in genes
            ],
        )

    if name and _is_named_group(name):
        return Entity(
            type=EntityTypeCv.COMPLEX,
            identifiers=[Identifier(type=IdentifierNamespaceCv.NAME, value=name)],
            annotations=annotations,
        )

    return Entity(
        type=EntityTypeCv.PROTEIN,
        identifiers=[
            Identifier(
                type=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                value=_normalize_fallback_gene(name, taxon_id),
            ),
            Identifier(type=IdentifierNamespaceCv.NAME, value=name),
        ],
        annotations=annotations,
    )


def _participant_annotations(
    row: dict[str, Any],
    role: str,
    role_term: InterCellAnnotations | None = None,
) -> list[Annotation]:
    annotations = [
        Annotation(term=role_term) if role_term is not None else None,
    ]
    annotations.extend(_location_terms(row.get(f'{role}_location')))
    return [item for item in annotations if item is not None]


def _target_role_term(row: dict[str, Any]) -> InterCellAnnotations | None:
    target_role = _first(row.get('target_role')).lower()
    if target_role == 'ligand':
        return InterCellAnnotations.LIGAND
    if target_role == 'receptor':
        return InterCellAnnotations.RECEPTOR
    return None


def map_cellchat_interaction(row: dict[str, Any]) -> Entity:
    taxon_id = _first(row.get('taxon_id'))
    evidence = _first(row.get('evidence'))

    annotations = _annotations(
        Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
        _neurotransmitter_annotation(row.get('is_neurotransmitter')),
        *_pubmed_annotations(evidence),
        *_pmc_annotations(evidence),
    )

    ligand = _protein_group_entity(
        name=_first(row.get('ligand_name')),
        genes=_values(row.get('ligand_genes')),
        taxon_id=taxon_id,
        annotations=_participant_annotations(
            row,
            'ligand',
            InterCellAnnotations.LIGAND,
        ),
    )
    receptor = _protein_group_entity(
        name=_first(row.get('receptor_name')),
        genes=_values(row.get('receptor_genes')),
        taxon_id=taxon_id,
        annotations=_participant_annotations(
            row,
            'receptor',
            InterCellAnnotations.RECEPTOR,
        ),
    )

    return Entity(
        type=EntityTypeCv.INTERACTION,
        identifiers=[
            Identifier(
                type=IdentifierNamespaceCv.CELLCHAT,
                value=_first(row.get('interaction_name')),
            ),
            Identifier(
                type=IdentifierNamespaceCv.NAME,
                value=_first(row.get('interaction_name_2')) or _first(row.get('interaction_name')),
            ),
        ],
        annotations=annotations,
        membership=[
            Membership(
                member=ligand,
                annotations=[
                    Annotation(term=InterCellAnnotations.LIGAND),
                ],
            ),
            Membership(
                member=receptor,
                annotations=[
                    Annotation(term=InterCellAnnotations.RECEPTOR),
                ],
            ),
        ],
    )


def map_cellchat_cofactor_interaction(row: dict[str, Any]) -> Entity:
    taxon_id = _first(row.get('taxon_id'))
    evidence = _first(row.get('evidence'))
    effect = _first(row.get('effect'))
    cofactor_gene = _first(row.get('cofactor_gene'))

    annotations = _annotations(
        Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
        *_pubmed_annotations(evidence),
        *_pmc_annotations(evidence),
    )

    cofactor = _protein_group_entity(
        name=cofactor_gene,
        genes=[cofactor_gene],
        taxon_id=taxon_id,
    )
    target = _protein_group_entity(
        name=_first(row.get('target_name')),
        genes=_values(row.get('target_genes')),
        taxon_id=taxon_id,
        annotations=_participant_annotations(
            row,
            'target',
            _target_role_term(row),
        ),
    )

    return Entity(
        type=EntityTypeCv.INTERACTION,
        identifiers=[
            Identifier(
                type=IdentifierNamespaceCv.CELLCHAT,
                value=_first(row.get('cofactor_interaction_name')),
            ),
            Identifier(
                type=IdentifierNamespaceCv.NAME,
                value=(
                    f'{cofactor_gene} '
                    f'{effect} '
                    f'{_first(row.get("target_name"))}'
                ),
            ),
        ],
        annotations=annotations,
        membership=[
            Membership(
                member=cofactor,
            ),
            Membership(
                member=target,
                annotations=(
                    [Annotation(term=role_term)]
                    if (role_term := _target_role_term(row)) is not None
                    else None
                ),
            ),
        ],
    )


def map_cellchat_group(row: dict[str, Any]) -> Entity:
    return _protein_group_entity(
        name=_first(row.get('name')),
        genes=_values(row.get('genes')),
        taxon_id=_first(row.get('taxon_id')),
    )


resource = Resource(
    config,
    human_interactions=Dataset(
        download=_download('human'),
        mapper=map_cellchat_interaction,
        raw_parser=_parser(iter_cellchat_interactions, 'human'),
    ),
    mouse_interactions=Dataset(
        download=_download('mouse'),
        mapper=map_cellchat_interaction,
        raw_parser=_parser(iter_cellchat_interactions, 'mouse'),
    ),
    zebrafish_interactions=Dataset(
        download=_download('zebrafish'),
        mapper=map_cellchat_interaction,
        raw_parser=_parser(iter_cellchat_interactions, 'zebrafish'),
    ),
    human_cofactor_interactions=Dataset(
        download=_download('human'),
        mapper=map_cellchat_cofactor_interaction,
        raw_parser=_parser(iter_cellchat_cofactor_interactions, 'human'),
    ),
    mouse_cofactor_interactions=Dataset(
        download=_download('mouse'),
        mapper=map_cellchat_cofactor_interaction,
        raw_parser=_parser(iter_cellchat_cofactor_interactions, 'mouse'),
    ),
    zebrafish_cofactor_interactions=Dataset(
        download=_download('zebrafish'),
        mapper=map_cellchat_cofactor_interaction,
        raw_parser=_parser(iter_cellchat_cofactor_interactions, 'zebrafish'),
    ),
    human_complexes=Dataset(
        download=_download('human'),
        mapper=map_cellchat_group,
        raw_parser=_parser(iter_cellchat_complexes, 'human'),
    ),
    mouse_complexes=Dataset(
        download=_download('mouse'),
        mapper=map_cellchat_group,
        raw_parser=_parser(iter_cellchat_complexes, 'mouse'),
    ),
    zebrafish_complexes=Dataset(
        download=_download('zebrafish'),
        mapper=map_cellchat_group,
        raw_parser=_parser(iter_cellchat_complexes, 'zebrafish'),
    ),
    human_cofactors=Dataset(
        download=_download('human'),
        mapper=map_cellchat_group,
        raw_parser=_parser(iter_cellchat_cofactors, 'human'),
    ),
    mouse_cofactors=Dataset(
        download=_download('mouse'),
        mapper=map_cellchat_group,
        raw_parser=_parser(iter_cellchat_cofactors, 'mouse'),
    ),
    zebrafish_cofactors=Dataset(
        download=_download('zebrafish'),
        mapper=map_cellchat_group,
        raw_parser=_parser(iter_cellchat_cofactors, 'zebrafish'),
    ),
)
