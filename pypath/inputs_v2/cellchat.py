"""
Parse CellChatDB data and emit Entity records.

This module converts CellChatDB ligand-receptor interactions, complexes and
cofactor groups into Entity records.
"""

from __future__ import annotations
import functools
import re
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    NamedThing,
    MacromolecularComplex,
    slots,
    DirectionQualifierEnum,
)
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
    Relation,
)
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
    description='CellChatDB is a manually curated database of literature-supported cell-cell communication interactions for CellChat, including ligand-receptor pairs, multi-subunit complexes and signaling cofactors across multiple species.',
)
BASE_URL = 'https://raw.githubusercontent.com/jinworks/CellChat/main/data/'
SPECIES = {'human': '9606', 'mouse': '10090', 'zebrafish': '7955'}


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


_pmid_re = re.compile('PMID:?\\s*(\\d+)', re.IGNORECASE)
_pmc_re = re.compile('\\bPMC(?:ID)?\\s*:?\\s*(\\d+)', re.IGNORECASE)


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
        key = (format_term(item.term), format_term(item.value), item.units)
        if key not in seen:
            out.append(item)
            seen.add(key)
    return out or None


def _annotation(term: Any, value: Any = None) -> Annotation | None:
    if value is None:
        return None
    value = str(value).strip()
    return Annotation(term=term, value=value) if value else None


def _pubmed_annotations(evidence: str) -> list[Annotation]:
    return [
        Annotation(term=slots.publications, value=f'PMID:{pmid}')
        for pmid in sorted(set(_pmid_re.findall(evidence or '')))
    ]


def _pmc_annotations(evidence: str) -> list[Annotation]:
    return [
        Annotation(term=slots.publications, value=f'PMC:PMC{pmc}')
        for pmc in sorted(set(_pmc_re.findall(evidence or '')))
    ]


def _protein_entity(gene: str, taxon_id: str) -> Entity:
    return Entity(
        type=Protein,
        identifiers=[Identifier(type=Namespace.GENESYMBOL, value=gene)],
        annotations=[
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}')
        ],
    )


def _protein_group_entity(
    *,
    name: str,
    genes: list[str],
    taxon_id: str,
    annotations: list[Annotation] | None = None,
    group_type: type = MacromolecularComplex,
) -> Entity:
    genes = [gene for gene in genes if gene]
    annotations = _annotations(
        Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
        *(annotations or []),
    )
    if len(genes) == 1:
        return Entity(
            type=Protein,
            identifiers=[
                Identifier(type=Namespace.GENESYMBOL, value=genes[0]),
                Identifier(type=Namespace.NAME, value=genes[0]),
            ],
            annotations=annotations,
        )
    return Entity(
        type=group_type if genes else NamedThing,
        identifiers=[Identifier(type=Namespace.NAME, value=name)],
        annotations=annotations,
        membership=[
            Membership(member=_protein_entity(gene, taxon_id)) for gene in genes
        ]
        or None,
    )


def map_cellchat_interaction(row: dict[str, Any]) -> Relation:
    taxon_id = _first(row.get('taxon_id'))
    evidence = _first(row.get('evidence'))
    annotations = _annotations(
        Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
        *_pubmed_annotations(evidence),
        *_pmc_annotations(evidence),
    )
    ligand = _protein_group_entity(
        name=_first(row.get('ligand_name')),
        genes=_values(row.get('ligand_genes')),
        taxon_id=taxon_id,
    )
    receptor = _protein_group_entity(
        name=_first(row.get('receptor_name')),
        genes=_values(row.get('receptor_genes')),
        taxon_id=taxon_id,
    )
    return Relation(
        subject=ligand,
        predicate=slots.interacts_with,
        object=receptor,
        identifiers=[
            Identifier(
                type=Namespace.CELLCHAT,
                value=_first(row.get('interaction_name')),
            ),
            Identifier(
                type=Namespace.NAME,
                value=_first(row.get('interaction_name_2'))
                or _first(row.get('interaction_name')),
            ),
        ],
        annotations=[*(annotations or [])],
    )


def map_cellchat_cofactor_interaction(row: dict[str, Any]) -> Relation:
    taxon_id = _first(row.get('taxon_id'))
    evidence = _first(row.get('evidence'))
    effect = _first(row.get('effect'))
    direction = {
        'stimulation': DirectionQualifierEnum.increased,
        'inhibition': DirectionQualifierEnum.decreased,
    }.get(effect)
    if direction is None:
        raise ValueError(f'Unknown CellChat cofactor effect: {effect!r}')
    cofactor_gene = _first(row.get('cofactor_gene'))
    annotations = _annotations(
        Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
        *_pubmed_annotations(evidence),
        *_pmc_annotations(evidence),
    )
    cofactor = _protein_group_entity(
        name=cofactor_gene, genes=[cofactor_gene], taxon_id=taxon_id
    )
    target = _protein_group_entity(
        name=_first(row.get('target_name')),
        genes=_values(row.get('target_genes')),
        taxon_id=taxon_id,
    )
    return Relation(
        subject=cofactor,
        predicate=slots.affects,
        object=target,
        identifiers=[
            Identifier(
                type=Namespace.CELLCHAT,
                value=_first(row.get('cofactor_interaction_name')),
            ),
            Identifier(
                type=Namespace.NAME,
                value=f'{cofactor_gene} {effect} {_first(row.get("target_name"))}',
            ),
        ],
        annotations=[
            *(annotations or []),
            Annotation(term=slots.object_direction_qualifier, value=direction),
            Annotation(term=slots.original_subject, value=_first(row.get("cofactor_group")) or cofactor_gene),
            Annotation(term=slots.original_predicate, value=_first(row.get("cofactor_role")) or effect),
        ],
    )


def map_cellchat_group(row: dict[str, Any]) -> Entity:
    return _protein_group_entity(
        name=_first(row.get('name')),
        genes=_values(row.get('genes')),
        taxon_id=_first(row.get('taxon_id')),
    )


def map_cellchat_cofactor_group(row: dict[str, Any]) -> Entity:
    return _protein_group_entity(
        name=_first(row.get('name')),
        genes=_values(row.get('genes')),
        taxon_id=_first(row.get('taxon_id')),
        group_type=NamedThing,
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
        mapper=map_cellchat_cofactor_group,
        raw_parser=_parser(iter_cellchat_cofactors, 'human'),
    ),
    mouse_cofactors=Dataset(
        download=_download('mouse'),
        mapper=map_cellchat_cofactor_group,
        raw_parser=_parser(iter_cellchat_cofactors, 'mouse'),
    ),
    zebrafish_cofactors=Dataset(
        download=_download('zebrafish'),
        mapper=map_cellchat_cofactor_group,
        raw_parser=_parser(iter_cellchat_cofactors, 'zebrafish'),
    ),
)
