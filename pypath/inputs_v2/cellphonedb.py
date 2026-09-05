"""
Parse CellPhoneDB data and emit Entity records.

This module converts CellPhoneDB interactions and complexes into Entity
records using the declarative schema pattern.
"""

from __future__ import annotations
import re
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    MacromolecularComplex,
    ChemicalEntity,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.silver_schema import Annotation, Entity, Identifier
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembershipBuilder,
    MembersFromList,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_csv

config = ResourceConfig(
    id=ResourceCv.CELLPHONEDB,
    name='CellPhoneDB',
    url='https://www.cellphonedb.org/',
    license=LicenseCV.MIT,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='40133495',
    primary_category='interactions',
    description='CellPhoneDB is a publicly available repository of curated receptors, ligands and their interactions, designed to enable the analysis of cell-cell communication from single-cell transcriptomics data.',
)
BASE_URL = (
    'https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/'
)
download_interactions = Download(
    url=BASE_URL + 'interaction_input.csv',
    filename='cellphonedb_interactions.csv',
    subfolder='cellphonedb',
    ext='csv',
)
download_complexes = Download(
    url=BASE_URL + 'complex_input.csv',
    filename='cellphonedb_complexes.csv',
    subfolder='cellphonedb',
    ext='csv',
)
download_proteins = Download(
    url=BASE_URL + 'protein_input.csv',
    filename='cellphonedb_proteins.csv',
    subfolder='cellphonedb',
    ext='csv',
)
UNIPROT_ACC_RE = re.compile(
    '^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$'
)
HUMAN_TAXON_ID = '9606'
SYNTHETIC_METABOLITE_SYSTEM_RE = re.compile('^(.+?)_by[A-Za-z0-9].*')


def _extract_pmid(token: str) -> str | None:
    """Extract PubMed ID from a token."""
    m = re.search('PMID:?\\s*(\\d+)', token, re.IGNORECASE)
    return m.group(1) if m else None


def _extract_pmc(token: str) -> str | None:
    """Extract PubMed Central ID from a token."""
    m = re.search('PMC\\s*(\\d+)', token, re.IGNORECASE)
    return f'PMC{m.group(1)}' if m else None


def _extract_comment(token: str) -> str | None:
    """Return the token if it's not a PMID or PMC."""
    if re.search('PMID|PMC', token, re.IGNORECASE):
        return None
    return token.strip()


def _source_split(row: dict[str, Any]) -> list[str]:
    """Split the source column into individual tokens using regex."""
    return re.split('[;,]\\s*', row.get('source') or '')


def _extract_uniprot_acc(val: str) -> str | None:
    """Return the value if it matches the UniProt accession pattern."""
    return val if UNIPROT_ACC_RE.match(val) else None


def _extract_non_uniprot(val: str) -> str | None:
    """Return the value if it does NOT match the UniProt accession pattern."""
    return val if not UNIPROT_ACC_RE.match(val) else None


def _synthetic_metabolite_name(val: str) -> str | None:
    """Extract the molecule name from CellPhoneDB synthetic metabolite systems."""
    m = SYNTHETIC_METABOLITE_SYSTEM_RE.match(val or '')
    return m.group(1) if m else None


def _is_synthetic_metabolite_system(val: str) -> bool:
    return _synthetic_metabolite_name(val) is not None


def _extract_partner_name(val: str) -> str | None:
    if UNIPROT_ACC_RE.match(val):
        return None
    return _synthetic_metabolite_name(val) or val


def _get_partner_type(col: str) -> Any:
    """
    Determine entity type for a partner.

    CellPhoneDB encodes some metabolite systems as names such as
    ``Glutamate_byGLS2_and_SLC1A1``. These are not protein complexes: the
    biological entity is the metabolite, while the suffix is a resource-specific
    label about supporting proteins. We preserve that label as annotation
    instead of turning the proteins into complex members.
    """

    def _type_selector(row: dict[str, Any]) -> type:
        val = row.get(col, '')
        if _is_synthetic_metabolite_system(val):
            return ChemicalEntity
        return Protein if UNIPROT_ACC_RE.match(val) else MacromolecularComplex

    return _type_selector


f = FieldConfig(
    extract={
        'pmid': _extract_pmid,
        'pmc': _extract_pmc,
        'comment': _extract_comment,
        'uniprot_acc': _extract_uniprot_acc,
        'non_uniprot': _extract_non_uniprot,
        'partner_name': _extract_partner_name,
    }
)
partner_a_builder = EntityBuilder(
    entity_type=_get_partner_type('partner_a'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f('partner_a', extract='uniprot_acc')),
        CV(term=Namespace.NAME, value=f('partner_a', extract='partner_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=lambda row: None
            if _is_synthetic_metabolite_system(row.get('partner_a', ''))
            else f'NCBITaxon:{HUMAN_TAXON_ID}',
        )
    ),
)
partner_b_builder = EntityBuilder(
    entity_type=_get_partner_type('partner_b'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f('partner_b', extract='uniprot_acc')),
        CV(term=Namespace.NAME, value=f('partner_b', extract='partner_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=lambda row: None
            if _is_synthetic_metabolite_system(row.get('partner_b', ''))
            else f'NCBITaxon:{HUMAN_TAXON_ID}',
        )
    ),
)
interactions_schema = RelationBuilder(
    subject=partner_a_builder,
    predicate=slots.interacts_with,
    object=partner_b_builder,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('interactors'))
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.publications,
            value=f(
                _source_split, extract='pmid', transform=lambda v: f'PMID:{v}'
            ),
        ),
        CV(
            term=slots.publications,
            value=f(
                _source_split, extract='pmc', transform=lambda v: f'PMC:{v}'
            ),
        ),
    ),
)
protein_complexes_schema = EntityBuilder(
    entity_type=MacromolecularComplex,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('complex_name'))
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=Protein,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.UNIPROT,
                    value=f(
                        lambda row: [
                            row.get(f'uniprot_{i}')
                            for i in range(1, 6)
                            if row.get(f'uniprot_{i}')
                        ]
                    ),
                )
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=slots.in_taxon, value=f'NCBITaxon:{HUMAN_TAXON_ID}')
            ),
        )
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=f'NCBITaxon:{HUMAN_TAXON_ID}'),
    ),
)


def _identifier(type_: object, value: object) -> Identifier | None:
    value = str(value or '').strip()
    return Identifier(type=type_, value=value) if value else None


def _annotation(term: object, value: object = None) -> Annotation | None:
    value = str(value or '').strip() if value is not None else None
    return Annotation(term=term, value=value) if value else None


def _identifiers(*items: Identifier | None) -> list[Identifier]:
    out: list[Identifier] = []
    seen: set[tuple[object, str]] = set()
    for item in items:
        if item is None:
            continue
        key = (item.type, item.value)
        if key in seen:
            continue
        out.append(item)
        seen.add(key)
    return out


def _annotations(*items: Annotation | None) -> list[Annotation] | None:
    out: list[Annotation] = []
    seen: set[tuple[object, object, object]] = set()
    for item in items:
        if item is None:
            continue
        key = (format_term(item.term), format_term(item.value), item.units)
        if key in seen:
            continue
        out.append(item)
        seen.add(key)
    return out or None


def complexes_schema(row: dict[str, Any]) -> Entity:
    name = str(row.get('complex_name') or '').strip()
    metabolite_name = _synthetic_metabolite_name(name)
    if not metabolite_name:
        return protein_complexes_schema(row)
    return Entity(
        type=ChemicalEntity,
        identifiers=_identifiers(_identifier(Namespace.NAME, metabolite_name)),
        annotations=_annotations(),
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=iter_csv,
    ),
    complexes=Dataset(
        download=download_complexes,
        mapper=complexes_schema,
        raw_parser=iter_csv,
    ),
)
