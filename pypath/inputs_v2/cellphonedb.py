"""
Parse CellPhoneDB data and emit Entity records.

This module converts CellPhoneDB interactions and complexes into Entity
records using the declarative schema pattern.
"""

from __future__ import annotations
from collections.abc import Generator
import re
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    MacromolecularComplex,
    ChemicalEntity,
    OntologyClass,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from omnipath_core.silver_schema import format_term
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    EntityRef,
    Identifier,
    OntologyRelation,
    Relation,
)
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
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(
        Content.LIGAND_RECEPTOR_INTERACTIONS,
        Content.PROTEIN_COMPLEXES,
        Content.MOLECULAR_ANNOTATIONS,
    ),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Ligand–receptor interactions, multisubunit protein complexes, and the '
        'roles of the proteins and complexes (receptor, secreted, transmembrane).'
    ),
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
    r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?(?:-PRO_\d+)?$'
)
HUMAN_TAXON_ID = '9606'
SYNTHETIC_METABOLITE_SYSTEM_RE = re.compile('^(.+?)_by[A-Za-z0-9].*')
CHEBI_RE = re.compile(r'^CHEBI:\d+$')
# A metabolite system ("Glutamate_byGLS_and_SLC1A1") has the ChEBI term of its
# metabolite in secreted_desc, or nothing. The ChEBI terms of the metabolites
# written without one, found by name or synonym in ChEBI.
METABOLITE_CHEBI = {
    '5-HT': 'CHEBI:28790',
    'Acetylcholine': 'CHEBI:15355',
    'Adenosine': 'CHEBI:16335',
    'Adrenaline': 'CHEBI:28918',
    'b-PEA': 'CHEBI:18397',
    'DHEAsulfate': 'CHEBI:16814',
    'Dehydroepiandrosterone': 'CHEBI:28689',
    'Desmosterol': 'CHEBI:17737',
    'Desoxycorticosterone': 'CHEBI:16973',
    'Dihydrotestosterone': 'CHEBI:16330',
    'Dopamine': 'CHEBI:18243',
    'GABA': 'CHEBI:16865',
    'Glutamate': 'CHEBI:16015',
    'Glycine': 'CHEBI:15428',
    'LTC4': 'CHEBI:16978',
    'LTD4': 'CHEBI:28666',
    'Noradrenaline': 'CHEBI:18357',
    'PGD2': 'CHEBI:15555',
    # Testosterone_byHSD17B12 has a receptor name in secreted_desc; the other
    # testosterone systems have this term.
    'Testosterone': 'CHEBI:17347',
}
# Systems of two metabolites. "SerotoninDopamin_byDDC" has the term of
# serotonin only; DDC makes both.
MIXED_METABOLITES = {'AldosteroneCorticosterone', 'SerotoninDopamin'}


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


def metabolite_chebi(name: str, secreted_desc: str = '') -> str | None:
    """The ChEBI term of the metabolite of a metabolite system, if it has one.

    Several terms in secreted_desc (progesterone with the other products of
    HSD3B1) are no single metabolite.
    """
    metabolite = _synthetic_metabolite_name(name)
    if not metabolite or metabolite in MIXED_METABOLITES:
        return None
    terms = [t for t in re.split(r'[;,\s]+', secreted_desc or '') if CHEBI_RE.match(t)]
    if terms:
        return terms[0] if len(terms) == 1 else None
    return METABOLITE_CHEBI.get(metabolite)


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
        CV(term=Namespace.CHEBI, value=f('partner_a_chebi')),
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
        CV(term=Namespace.CHEBI, value=f('partner_b_chebi')),
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
        identifiers=_identifiers(
            _identifier(Namespace.CHEBI, metabolite_chebi(name, row.get('secreted_desc'))),
            _identifier(Namespace.NAME, metabolite_name),
        ),
        annotations=_annotations(),
    )


def iter_interactions(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The interactions, with the ChEBI term of a metabolite partner from the
    complex table."""
    chebi = {
        row['complex_name']: metabolite_chebi(row['complex_name'], row.get('secreted_desc'))
        for row in iter_csv(download_complexes.open())
    }
    for row in iter_csv(opener, **kwargs):
        yield {
            **row,
            **{
                f'{side}_chebi': chebi.get(row.get(side, '')) or metabolite_chebi(row.get(side, ''))
                for side in ('partner_a', 'partner_b')
            },
        }


# Roles of the proteins and complexes: the flags of protein_input.csv and
# complex_input.csv, which say which side of an interaction a partner can be
# on, and the classes of the receptors and secreted proteins ("Cytokine",
# "TNF_receptor_add"). A class is below its role in the role_terms dataset; a
# partner is associated with its classes, or with the role where it has none.
# The two roles with classes have a namespace each: a class name can be a
# receptor class and a secreted class.
ROLE_NAMESPACE = 'cellphonedb_role'
ROLES = ('receptor', 'secreted', 'transmembrane', 'peripheral', 'integrin')
CLASSIFIED_ROLES = {
    'receptor': ('receptor_desc', 'cellphonedb_receptor_class'),
    'secreted': ('secreted_desc', 'cellphonedb_secreted_class'),
}


def _role_classes(value: str) -> list[str]:
    """The classes in a receptor_desc or secreted_desc cell, as one spelling:
    lower case, words joined by "_", without the curators' "add" mark. Classes
    are separated by "|"; text after ";" is a note ("nuclear_receptor;needs
    transporter")."""
    classes = []
    for item in (value or '').split('|'):
        item = item.split(';')[0].strip()
        if not item or item.upper() == 'FALSE' or item.startswith('http') or CHEBI_RE.match(item):
            continue
        words = [w for w in re.split(r'[\s,_]+', item.lower()) if w and w != 'add']
        item = '_'.join(words)
        if item and item not in classes:
            classes.append(item)
    return classes


def _roles(row: dict[str, Any]) -> Generator[tuple[str, str | None], None, None]:
    """(role, class) of each role of a row; the class is None where the role
    has none."""
    metabolite = _is_synthetic_metabolite_system(row.get('complex_name') or '')
    for role in ROLES:
        if str(row.get(role) or '').upper() != 'TRUE':
            continue
        column = CLASSIFIED_ROLES.get(role, ('', ''))[0]
        # The description of a metabolite system is its ChEBI term, if any.
        classes = [] if metabolite else _role_classes(row.get(column))
        if classes:
            yield from ((role, c) for c in classes)
        else:
            yield role, None


def _role_ref(role: str, role_class: str | None) -> EntityRef:
    if role_class:
        return EntityRef(OntologyClass, CLASSIFIED_ROLES[role][1], role_class)
    return EntityRef(OntologyClass, ROLE_NAMESPACE, role)


def _protein_entity(row: dict[str, Any]) -> Entity:
    value = str(row.get('uniprot') or '').strip()
    if UNIPROT_ACC_RE.match(value):
        identifier = Identifier(type=Namespace.UNIPROT, value=value)
    else:
        # HLA genes are written without the dash ("HLAA"); IFNA1 by its symbol.
        if value.startswith('HLA') and '-' not in value:
            value = f'HLA-{value[3:]}'
        identifier = Identifier(type=Namespace.GENESYMBOL, value=value)
    return Entity(
        type=Protein,
        identifiers=[identifier],
        annotations=[Annotation(term=slots.in_taxon, value=f'NCBITaxon:{HUMAN_TAXON_ID}')],
    )


def iter_roles(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One row for each role of each protein or complex."""
    for row in iter_csv(opener, **kwargs):
        for role, role_class in _roles(row):
            yield {**row, 'role': role, 'role_class': role_class}


def iter_role_terms(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The roles, and the classes found in the proteins and complexes tables."""
    seen = set()
    for role in ROLES:
        seen.add((role, None))
        yield {'role': role, 'role_class': None}
    for rows in (iter_csv(opener, **kwargs), iter_csv(download_complexes.open())):
        for row in rows:
            for term in _roles(row):
                if term not in seen:
                    seen.add(term)
                    yield {'role': term[0], 'role_class': term[1]}


def map_role_term(row: dict[str, Any]) -> Entity:
    """A role, or a class below its role."""
    ref = _role_ref(row['role'], row['role_class'])
    relations = None
    if row['role_class']:
        relations = [
            OntologyRelation(
                predicate=slots.subclass_of,
                object=_role_ref(row['role'], None),
                ontology_id='cellphonedb',
            )
        ]
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=ref.identifier_type, value=ref.identifier),
            Identifier(type=Namespace.NAME, value=ref.identifier),
        ],
        ontology_relations=relations,
    )


def _role_relation(subject: Entity, row: dict[str, Any]) -> Relation:
    return Relation(
        subject=subject,
        predicate=slots.associated_with,
        object=_role_ref(row['role'], row['role_class']),
    )


def map_protein_role(row: dict[str, Any]) -> Relation:
    """A protein associated with one of its roles."""
    return _role_relation(_protein_entity(row), row)


def map_complex_role(row: dict[str, Any]) -> Relation:
    """A complex or metabolite system associated with one of its roles."""
    return _role_relation(complexes_schema(row), row)


resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=iter_interactions,
    ),
    complexes=Dataset(
        download=download_complexes,
        mapper=complexes_schema,
        raw_parser=iter_csv,
    ),
    role_terms=Dataset(
        download=download_proteins,
        mapper=map_role_term,
        raw_parser=iter_role_terms,
    ),
    protein_roles=Dataset(
        download=download_proteins,
        mapper=map_protein_role,
        raw_parser=iter_roles,
    ),
    complex_roles=Dataset(
        download=download_complexes,
        mapper=map_complex_role,
        raw_parser=iter_roles,
    ),
)
