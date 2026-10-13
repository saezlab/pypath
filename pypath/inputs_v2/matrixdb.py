"""
Parse MatrixDB and emit its ECM protein list and its curated interactions.

MatrixDB is the database of interactions of the extracellular matrix (ECM). Three
datasets are read from its download page:

- ``terms``: the list of ECM and secreted proteins, as one term of MatrixDB;
- ``ecm_proteins``: each human, reviewed protein of the list
  (``ecm_proteins.csv``: UniProtKB proteins of all species that MatrixDB selects
  by Gene Ontology cellular component terms and UniProt keywords), associated
  with the term. The GO terms and keywords that select the protein and the
  number of MatrixDB interactions of the protein are attributes of the
  association;
- ``interactions``: the manually curated interactions of MatrixDB (release 4.0,
  ``matrixdb_CORE``) in PSI-MI TAB 2.7: proteins, protein fragments, chemicals
  (glycosaminoglycans, cations, lipids), complexes and multimers, with detection
  methods, interaction types, publications, participant roles and affinity
  parameters.

The legacy module read the lists of ECM, membrane and secreted proteins as three
lists; MatrixDB now offers one (``ecm_proteins.csv``) and the membrane and secreted
lists are gone. The legacy module kept the human Swiss-Prot proteins of a list;
this module keeps the same proteins, by the species and the UniProt dataset that
the file states for each (2,835 of its 49,978 proteins, file of 2026-10-12: 41,145
are unreviewed TrEMBL entries and 44,455 are of the 15 other species). The release that adds the interactions imported from IMEx databases
(``matrixdb_all``, 208 MB) is not read: those interactions are the ones the IMEx
sources already provide.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.molecular_forms import molecular_form_from_identifiers
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2 import _mitab
from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.matrixdb import (
    iter_ecm_proteins,
    iter_interactions,
    iter_list_terms,
)
from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    EntityRef,
    Identifier,
    Relation,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.MATRIXDB,
    name='MatrixDB',
    url='https://matrixdb.expasy.org/',
    # https://matrixdb.expasy.org/download: "All data provided on this page is
    # licensed under the Creative Commons Attribution License (CC BY 4.0)".
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='39558161',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.PHYSICAL_INTERACTIONS, Content.MOLECULAR_ANNOTATIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED, KnowledgeOrigin.EXPERIMENTAL),
    summary=(
        'Manually curated interactions of extracellular matrix proteins, '
        'glycosaminoglycans and other biomolecules, and the list of human ECM '
        'and secreted proteins.'
    ),
    description=(
        'MatrixDB is a database of the interactions of the extracellular matrix '
        '(ECM). It collects the interactions of ECM proteins, proteoglycans, '
        'glycosaminoglycans and polysaccharides with each other and with cell '
        'surface proteins, from the literature and in house experiments, and '
        'is a member of the IMEx consortium. It selects ECM and secreted '
        'proteins from UniProtKB by Gene Ontology terms and UniProt keywords and '
        'by the matrisome divisions and categories.'
    ),
)

# ================================ ECM PROTEINS ================================

LIST_NAMESPACE = 'matrixdb_protein_list'
GO_COMPONENT = 'matrixdb:go_cellular_component'
KEYWORD = 'matrixdb:uniprot_keyword'
INTERACTION_COUNT = 'matrixdb:interaction_count'

# The criteria table of https://matrixdb.expasy.org/download ("MatrixDB selects
# extracellular matrix (ECM) and secreted proteins according to the following
# criteria"): the names that ``ecm_proteins.csv`` writes for the GO cellular
# component terms and the UniProtKB keywords, and their accessions there.
GO_TERMS = {
    'Extracellular matrix': 'GO:0031012',
    'Basement membrane': 'GO:0005604',
    'Extracellular space': 'GO:0005615',
    'Extracellular region': 'GO:0005576',
}
KEYWORDS = {
    'Basement membrane': 'KW-0084',
    'Extracellular matrix': 'KW-0272',
    'Extracellular space/secreted': 'KW-0964',
}
UNIPROT_DATASETS = ('Swiss-Prot', 'UniProtKB unreviewed (TrEMBL)')
# The proteins of the list that are kept: reviewed entries of human, as in the
# legacy module (see the module docstring).
KEPT_DATASET = 'Swiss-Prot'
KEPT_TAXON = '9606'
LIST_NAME = 'ECM and secreted proteins'


def map_list_term(row: dict[str, Any]) -> Entity:
    """The list of ECM and secreted proteins as a term."""
    return Entity(
        type=model.OntologyClass,
        identifiers=[
            Identifier(type=LIST_NAMESPACE, value=row['list']),
            Identifier(type=Namespace.NAME, value=LIST_NAME),
        ],
        annotations=[
            Annotation(
                term=slots.description,
                value=(
                    'The extracellular matrix (ECM) and secreted proteins that '
                    'MatrixDB selects from UniProtKB by Gene Ontology cellular '
                    'component terms (extracellular matrix, basement membrane, '
                    'extracellular space, extracellular region) and UniProt '
                    'keywords (basement membrane, extracellular matrix, '
                    'extracellular space/secreted).'
                ),
            )
        ],
    )


def _crosswalk(names: list[str], table: dict[str, str], what: str) -> list[str]:
    unknown = [name for name in names if name not in table]
    if unknown:
        raise ValueError(f'Unknown MatrixDB {what}: {unknown}')
    return [table[name] for name in names]


def _interaction_count(row: dict[str, Any]):
    value = measurement(row['interaction_count'], source_field='interaction_count')
    if value is None:
        raise ValueError(
            f'MatrixDB {row["accession"]}: interaction count {row["interaction_count"]!r}'
        )
    return value


def map_ecm_protein(row: dict[str, Any]) -> Relation | None:
    """A protein associated with the list; its selection criteria are attributes.

    None for a protein that is not kept: an unreviewed entry or one of another
    species than human.
    """
    if row['dataset'] not in UNIPROT_DATASETS:
        raise ValueError(f'Unknown MatrixDB UniProt dataset: {row["dataset"]!r}')
    if row['dataset'] != KEPT_DATASET or row['taxon'] != KEPT_TAXON:
        return None
    annotations = [
        Annotation(term=GO_COMPONENT, value=accession)
        for accession in _crosswalk(row['go_terms'], GO_TERMS, 'GO term')
    ]
    annotations += [
        Annotation(term=KEYWORD, value=accession)
        for accession in _crosswalk(row['keywords'], KEYWORDS, 'UniProt keyword')
    ]
    annotations += [
        # How many interactions MatrixDB has for the protein: a Measurement, a
        # count without a unit.
        Annotation(term=INTERACTION_COUNT, value=_interaction_count(row)),
    ]
    return Relation(
        subject=Entity(
            type=model.Protein,
            identifiers=[Identifier(type=Namespace.UNIPROT, value=row['accession'])],
            annotations=[
                Annotation(term=slots.in_taxon, value=f'NCBITaxon:{row["taxon"]}')
            ],
        ),
        predicate=slots.associated_with,
        object=EntityRef(model.OntologyClass, LIST_NAMESPACE, 'ecm_proteins'),
        annotations=annotations,
    )


# ================================ INTERACTIONS ================================

MATRIXDB = 'matrixdb'
MATRIXDB_INTERACTION = 'matrixdb_interaction'
EXPANSION_METHOD = 'matrixdb:expansion_method'
# MatrixDB identifies the biomolecules that are not proteins itself; the page
# defines the types by the prefix of the identifier.
_MATRIXDB_TYPES = {
    'PFRAG': model.Protein,  # protein fragment
    'MULT': model.MacromolecularComplex,  # multimer
    'GAG': model.ChemicalEntity,  # glycosaminoglycan
    'CAT': model.ChemicalEntity,  # cation
    'LIP': model.ChemicalEntity,  # lipid
    'SMALLMOL': model.ChemicalEntity,  # small molecule
    'SPEP': model.ChemicalEntity,  # synthetic peptide
}
_DATABASE_TYPES = {
    'uniprotkb': model.Protein,
    'chebi': model.ChemicalEntity,
    'complex portal': model.MacromolecularComplex,
}
_NAMESPACES = {**_mitab.MITAB_NAMESPACES, MATRIXDB: MATRIXDB}

# Interaction types of the file that are not covered by the shared rules for direct
# interactions and associations: colocalization is a different assertion.
_COLOCALIZATION = 'MI:0403'


def _primary(row: dict[str, Any], side: str):
    items = _mitab.parse_mitab_cell(row[f'id_{side}'])
    if len(items) != 1:
        raise ValueError(f'MatrixDB interactor without one identifier: {row[f"id_{side}"]!r}')
    return items[0]


def _entity_type(side: str):
    def _type(row: dict[str, Any]):
        item = _primary(row, side)
        if item.database == MATRIXDB:
            prefix = item.accession.split('_')[0]
            if prefix not in _MATRIXDB_TYPES:
                raise ValueError(f'Unknown MatrixDB biomolecule: {item.accession!r}')
            return _MATRIXDB_TYPES[prefix]
        if item.database not in _DATABASE_TYPES:
            raise ValueError(f'Unknown MatrixDB interactor database: {item.database!r}')
        return _DATABASE_TYPES[item.database]

    return _type


def _identifier_pairs(side: str):
    """The primary identifier, the MatrixDB identifier, then names and gene names."""

    def _pairs(row: dict[str, Any]) -> list[tuple[Any, str]]:
        pairs = []
        for column in (f'id_{side}', f'alt_ids_{side}'):
            pairs += _mitab.identifier_pairs(row[column], _NAMESPACES)
        for item in _mitab.parse_mitab_cell(row[f'aliases_{side}']):
            if item.description == 'short label':
                pairs.append((Namespace.NAME, item.accession))
            elif item.description == 'gene name':
                pairs.append((Namespace.GENESYMBOL, item.accession))
            else:
                raise ValueError(f'Unknown MatrixDB alias: {item!r}')
        return list(dict.fromkeys(pairs))

    return _pairs


def _participant_form(side: str):
    pairs = _identifier_pairs(side)

    def _form(row: dict[str, Any]) -> dict | None:
        return molecular_form_from_identifiers(
            [
                {'ns': str(ns), 'id': accession}
                for ns, accession in pairs(row)
                if ns == Namespace.UNIPROT
            ]
        )

    return _form


def _participant(side: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=_entity_type(side),
        molecular_form=_participant_form(side),
        identifiers=IdentifiersBuilder(CV.from_pairs(_identifier_pairs(side))),
        annotations=AnnotationsBuilder(
            *_mitab.participant_annotations(
                taxon=f'taxid_{side}',
                biological_role=f'biological_role_{side}',
                experimental_role=f'experimental_role_{side}',
                identification_method=f'identification_method_{side}',
            ),
        ),
    )


def _predicate(row: dict[str, Any]):
    """
    Predicate for the PSI-MI interaction types of the row.

    A direct interaction or one of its kinds (cleavage, enzymatic reactions,
    self interaction) is a direct physical interaction, an association or
    physical association an interaction, and colocalization alone colocalization.
    Any other type gives None.
    """
    types = _mitab.mi_terms(row['interaction_types'])
    if set(types) == {_COLOCALIZATION}:
        return slots.colocalizes_with
    return _mitab.physical_predicate(types)


def _interaction_ids(row: dict[str, Any]) -> list[tuple[Any, str]]:
    pairs = []
    for item in _mitab.parse_mitab_cell(row['interaction_ids']):
        if item.database == MATRIXDB:
            pairs.append((MATRIXDB_INTERACTION, item.accession))
        elif item.database in ('intact', 'imex'):
            pairs.append((_mitab.MITAB_NAMESPACES[item.database], item.accession))
        else:
            raise ValueError(f'Unknown MatrixDB interaction identifier: {item!r}')
    return list(dict.fromkeys(pairs))


_PARAMETER_FIELDS = {'kd': 'kd', 'kon': 'kon', 'koff': 'koff'}


def _parameters(row: dict[str, Any]) -> list:
    """
    The kinetic parameters of the interaction as measurements.

    ``kd:5e-08(molar)``; ``kon`` (per mole per second) and ``koff`` (per second)
    are the others. A range (``kd:1e-08~5e-08(molar)``) is not a single value and
    is left out.
    """
    result = []
    for item in _mitab.parse_mitab_cell(row['parameters']):
        if item.database not in _PARAMETER_FIELDS:
            raise ValueError(f'Unknown MatrixDB interaction parameter: {item!r}')
        value = measurement(
            item.accession, unit=item.description, source_field=item.database
        )
        if value is not None:
            result.append(value)
    return result


def _expansion(row: dict[str, Any]) -> list[str]:
    return _mitab.mi_terms(row['expansion_methods'])


interactions_schema = RelationBuilder(
    subject=_participant('a'),
    predicate=_predicate,
    object=_participant('b'),
    identifiers=IdentifiersBuilder(CV.from_pairs(_interaction_ids)),
    annotations=AnnotationsBuilder(
        *_mitab.interaction_annotations(
            detection_method='detection_methods',
            interaction_type='interaction_types',
            publications='publications',
            host_organism='host_organisms',
        ),
        CV(term=slots.has_quantitative_value, value=_parameters),
        # "spoke expansion": the pair is one of the pairs that a bait and each of
        # its preys of a multi-participant experiment make.
        CV(term=EXPANSION_METHOD, value=_expansion),
    ),
)


def exclusion_reason(row: dict[str, Any]) -> str | None:
    """Why a row gives no relation; None for a row that is mapped."""
    if _mitab.is_negative(row['negative']):
        return 'negative result'
    if _predicate(row) is None:
        raise ValueError(
            f'Unknown MatrixDB interaction type: {row["interaction_types"]!r}'
        )
    return None


def map_interaction(row: dict[str, Any]):
    return None if exclusion_reason(row) else interactions_schema(row)


# =================================== RESOURCE =================================

# matrixdb.univ-lyon1.fr, the address in the legacy module, now points to the same
# server under another name whose certificate does not match; the site is at
# matrixdb.expasy.org.
BASE_URL = 'https://matrixdb.expasy.org/downloads/'
ecm_proteins_download = Download(
    url=BASE_URL + 'ecm_proteins.csv',
    filename='matrixdb_ecm_proteins.csv',
    subfolder='matrixdb',
    large=True,
)

resource = Resource(
    config,
    terms=Dataset(
        download=ecm_proteins_download,
        mapper=map_list_term,
        raw_parser=iter_list_terms,
    ),
    ecm_proteins=Dataset(
        download=ecm_proteins_download,
        mapper=map_ecm_protein,
        raw_parser=iter_ecm_proteins,
    ),
    interactions=Dataset(
        download=Download(
            url=BASE_URL + 'matrixdb_CORE_4_0.tab.zip',
            filename='matrixdb_CORE_4_0.tab.zip',
            subfolder='matrixdb',
            large=True,
            ext='zip',
        ),
        mapper=map_interaction,
        raw_parser=iter_interactions,
    ),
)
