"""
Parse ICELLNET v2 data and emit Entity records.
"""

from __future__ import annotations
import re
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    MacromolecularComplex,
    slots,
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
from pypath.inputs_v2.parsers.base import iter_csv

config = ResourceConfig(
    id=ResourceCv.ICELLNET,
    name='ICELLNET v2',
    url='https://github.com/soumelis-lab/ICELLNET',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='38490248',
    primary_category='interactions',
    description='ICELLNET is a transcriptomic framework with an expert-curated human ligand-receptor interaction database accounting for multi-subunit ligands and receptors.',
)
download = Download(
    url='https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.csv',
    filename='icellnetdb.csv',
    subfolder='icellnet',
    encoding='utf-8-sig',
)
_pmid_re = re.compile('\\d+')
_spreadsheet_currency_gene_re = re.compile('^(\\d+),00[\\s\\xa0]+(DKK)$')


def _clean(value: Any) -> str:
    return str(value or '').strip()


def _clean_gene(value: Any) -> str:
    value = _clean(value)
    if match := _spreadsheet_currency_gene_re.fullmatch(value):
        return f'{match.group(2)}{match.group(1)}'
    return value


def _components(row: dict[str, Any], prefix: str, count: int) -> list[str]:
    return [
        value
        for idx in range(1, count + 1)
        if (value := _clean_gene(row.get(f'{prefix} {idx}')))
    ]


def _pubmeds(value: Any) -> list[str]:
    return _pmid_re.findall(_clean(value))


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


def _protein(gene: str, annotations: list[Annotation] | None = None) -> Entity:
    return Entity(
        type=Protein,
        identifiers=[Identifier(type=Namespace.GENESYMBOL, value=gene)],
        annotations=[
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{"9606"}'),
            *(annotations or []),
        ],
    )


def _participant(
    name: str, genes: list[str], annotations: list[Annotation] | None = None
) -> Entity:
    if len(genes) == 1:
        return Entity(
            type=Protein,
            identifiers=[Identifier(type=Namespace.GENESYMBOL, value=genes[0])],
            annotations=[
                Annotation(term=slots.in_taxon, value=f'NCBITaxon:{"9606"}'),
                *(annotations or []),
            ],
        )
    return Entity(
        type=MacromolecularComplex,
        identifiers=[Identifier(type=Namespace.NAME, value=name)],
        annotations=[
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{"9606"}'),
            *(annotations or []),
        ],
        membership=[Membership(member=_protein(gene)) for gene in genes],
    )


def _interaction_name(ligands: list[str], receptors: list[str]) -> str:
    return f'{"+".join(ligands)} - {"+".join(receptors)}'


def map_icellnet_interaction(row: dict[str, Any]) -> Relation | None:
    ligands = _components(row, 'Ligand', 4)
    receptors = _components(row, 'Receptor', 5)
    if not ligands or not receptors:
        return None
    name = _interaction_name(ligands, receptors)
    annotations = _annotations(
        Annotation(term=slots.in_taxon, value=f'NCBITaxon:{"9606"}'),
        *[
            Annotation(term=slots.publications, value=f'PMID:{pmid}')
            for pmid in _pubmeds(row.get('Reference'))
        ],
    )
    return Relation(
        subject=_participant('+'.join(ligands), ligands, []),
        predicate=slots.interacts_with,
        object=_participant('+'.join(receptors), receptors, []),
        identifiers=[Identifier(type=Namespace.NAME, value=name)],
        annotations=[*(annotations or [])],
    )


def iter_icellnet_interaction_rows(opener, **kwargs: Any):
    for row in iter_csv(opener, delimiter=';', **kwargs):
        if _components(row, 'Ligand', 4) and _components(row, 'Receptor', 5):
            yield row


def iter_icellnet_complex_rows(opener, **kwargs: Any):
    seen: set[tuple[str, ...]] = set()
    for row in iter_csv(opener, delimiter=';', **kwargs):
        for prefix, count in (('Ligand', 4), ('Receptor', 5)):
            genes = _components(row, prefix, count)
            if len(genes) <= 1:
                continue
            key = tuple(genes)
            if key in seen:
                continue
            seen.add(key)
            yield {'name': '+'.join(genes), 'genes': genes, 'taxon_id': '9606'}


def map_icellnet_complex(row: dict[str, Any]) -> Entity:
    return _participant(_clean(row.get('name')), list(row.get('genes') or []))


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=map_icellnet_interaction,
        raw_parser=iter_icellnet_interaction_rows,
    ),
    complexes=Dataset(
        download=download,
        mapper=map_icellnet_complex,
        raw_parser=iter_icellnet_complex_rows,
    ),
)
