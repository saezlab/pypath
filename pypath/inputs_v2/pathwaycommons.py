"""
Parse the Pathway Commons binary relations and emit Relation records.

Pathway Commons integrates pathway and interaction databases into one BioPAX
model, and derives binary relations from it with a set of rules (its "binary
SIF" types, documented at https://www.pathwaycommons.org/pc2/formats). This
module reads the extended SIF of all data sources for human genes
(``pc-hgnc.txt.gz`` of release 14): one line is one relation of one type
between two participants, with the data sources it is derived from, its
PubMed IDs and its pathways. Participants are identified by HGNC symbol (genes
and gene products) or ChEBI identifier (small molecules).

The relation types are mapped one by one (``RELATION_TYPES``), each from what
Pathway Commons says it means. The SIF types state no sign, so no relation has
a direction of effect. A type that relates reactions rather than molecules, or
that has no faithful predicate, is not mapped.
"""

from __future__ import annotations

import re
from typing import Any, NamedTuple

from biolink_model.datamodel.model import (
    Gene,
    GeneOrGeneProductOrChemicalEntityAspectEnum as Aspect,
    ChemicalEntity,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.pathwaycommons import iter_interactions
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.PATHWAYCOMMONS,
    name='PathwayCommons',
    full='Pathway Commons',
    url='https://www.pathwaycommons.org/',
    # pathwaycommons.org, FAQ "Is it free to use?": "Pathway Commons
    # distributes pathway information with the intellectual property
    # restrictions of the source database; Only databases that are freely
    # available for academics are included." The data sources page and
    # metadata.json give each source an availability, "free" or "academic"
    # (academic: HumanCyc, HPRD, CORUM, MSigDB, miRTarBase, DrugBank, CTD,
    # KEGG).
    license=LicenseCV.COMPOSITE,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='31647099',
    topics=(Topic.SIGNALING, Topic.GENE_REGULATION, Topic.PHARMACOLOGY),
    content=(
        Content.PHYSICAL_INTERACTIONS,
        Content.CAUSAL_REGULATION,
        Content.COMPOUND_TARGET_INTERACTIONS,
    ),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Binary relations between human genes, gene products and small '
        'molecules that Pathway Commons derives from the pathway and '
        'interaction databases it integrates, with their data sources and '
        'publications.'
    ),
    description=(
        'Pathway Commons is a collection of publicly available pathway and '
        'molecular interaction data from many databases (Reactome, NCI-PID, '
        'PANTHER, HumanCyc, KEGG, BioGRID, IntAct, HPRD, CTD, DrugBank and '
        'others), converted to the BioPAX standard, normalized and merged. '
        'Its binary SIF files reduce the BioPAX model to pairwise relations '
        'of defined types, such as "controls-state-change-of" and '
        '"in-complex-with"; the extended format adds the data sources, '
        'PubMed IDs and pathways of each relation.'
    ),
)

# Release 14 (21 May 2024) is the latest release in the download directory.
# The directory is download.baderlab.org/PathwayCommons/PC2/; releases before
# 12 also had one SIF file for each data source, which later releases do not.
VERSION = 14
_URL = (
    'https://download.baderlab.org/PathwayCommons/PC2/'
    f'v{VERSION}/pc-hgnc.txt.gz'
)
# The file holds human genes, identified by HGNC symbol.
_TAXON = 'NCBITaxon:9606'


class _Type(NamedTuple):
    """What a relation type of Pathway Commons asserts."""

    predicate: Any
    # ``object_aspect_qualifier`` of the effect, where the type names one.
    aspect: Any = None
    # The class of the participants: 'molecule' for a gene or gene product,
    # 'chemical' for a small molecule.
    subject: str = 'molecule'
    object: str = 'molecule'


# Relation type -> predicate. The descriptions are those of the Pathway
# Commons documentation ("Types of Binary Relations"); None leaves the type
# out. A type that is not here raises.
RELATION_TYPES: dict[str, _Type | None] = {
    # "Proteins are participants of the same MolecularInteraction."
    'interacts-with': _Type(slots.interacts_with),
    # "Proteins are members of the same complex." The source expands a
    # complex to the pairs of its members; the pair says nothing about
    # contact between the two.
    'in-complex-with': _Type(slots.in_complex_with),
    # "First protein controls a reaction that changes the state of the
    # second protein." Which state is not given.
    'controls-state-change-of': _Type(slots.affects),
    # "... changes the phosphorylation status of the second protein": the
    # type covers phosphorylation and dephosphorylation, and the sign of the
    # control is not given.
    'controls-phosphorylation-of': _Type(slots.affects, Aspect.phosphorylation),
    # "... changes the cellular location of the second protein."
    'controls-transport-of': _Type(slots.affects, Aspect.localization),
    # "First protein controls a conversion or a template reaction that
    # changes expression of the second protein." No sign.
    'controls-expression-of': _Type(slots.affects, Aspect.expression),
    # "A small molecule has an effect on the protein state."
    'chemical-affects': _Type(slots.affects, subject='chemical'),
    # "The protein controls a reaction of which the small molecule is an
    # output."
    'controls-production-of': _Type(slots.affects, object='chemical'),
    # "The protein controls a reaction that changes cellular location of the
    # small molecule."
    'controls-transport-of-chemical': _Type(
        slots.affects, Aspect.localization, object='chemical'
    ),
    # Types left out:
    # "First protein controls a reaction whose output molecule is input to
    # another reaction controled by the second protein": a relation between
    # two reactions, with proteins as their controllers.
    'catalysis-precedes': None,
    # "The small molecule is consumed by a reaction that is controled by a
    # protein": no predicate for the consumption of a molecule by a reaction
    # that a protein controls.
    'consumption-controlled-by': None,
    # "Small molecules are input to a biochemical reaction": small molecules
    # of one reaction, not a statement between the two.
    'reacts-with': None,
    # "A reaction consumes a small molecule to produce another small
    # molecule": a relation between the participants of a reaction.
    'used-to-produce': None,
}

# ``INTERACTION_DATA_SOURCE``: every data source the file names, in the file's
# spelling ("pid" is NCI-PID). A name that is not here raises.
DATA_SOURCES = {
    name: name
    for name in (
        'BIND', 'Biofactoid', 'BioGRID', 'CORUM', 'CTD', 'DIP', 'DrugBank',
        'HPRD', 'HumanCyc', 'INOH', 'IntAct', 'KEGG', 'MSigDB', 'NetPath',
        'PANTHER', 'Pathbank', 'PhosphoSite', 'Reactome', 'Recon X', 'pid',
    )
}

# BioPAX reference types of a participant. A participant with a protein
# reference is a protein (its symbol also names the gene and transcript the
# reference types of the same participant can show); one with only a small
# molecule reference is a small molecule; one with only DNA and RNA references
# has no protein in the file and is a gene.
_CHEMICAL_TYPE = 'SmallMoleculeReference'
_PROTEIN_TYPE = 'ProteinReference'
_NUCLEIC_ACID_TYPES = frozenset(
    {'DnaReference', 'DnaRegionReference', 'RnaReference', 'RnaRegionReference'}
)

_CHEBI = re.compile(r'(?i)chebi:(\d+)')
_UNIPROT_ACCESSION = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
)


def participant_kind(participant: dict[str, Any]) -> str:
    """'chemical', 'protein' or 'gene' by the reference types of the participant."""
    types = set(participant['types'])
    if types == {_CHEMICAL_TYPE}:
        return 'chemical'
    if _CHEMICAL_TYPE in types or not types:
        raise ValueError(
            f'Unknown Pathway Commons participant types: {participant["types"]!r}'
        )
    if _PROTEIN_TYPE in types:
        if not types <= _NUCLEIC_ACID_TYPES | {_PROTEIN_TYPE}:
            raise ValueError(
                f'Unknown Pathway Commons participant types: {participant["types"]!r}'
            )
        return 'protein'
    if not types <= _NUCLEIC_ACID_TYPES:
        raise ValueError(
            f'Unknown Pathway Commons participant types: {participant["types"]!r}'
        )
    return 'gene'


def _identifiers(participant: dict[str, Any], kind: str) -> list[Identifier]:
    if kind == 'chemical':
        # ``chebi:15377`` and ``CHEBI:15377`` both name ChEBI 15377. One
        # participant has a PubChem identifier in the ChEBI form.
        match = _CHEBI.fullmatch(participant['id'])
        identifiers = (
            [Identifier(type=Namespace.CHEBI, value=f'CHEBI:{match.group(1)}')]
            if match
            else []
        )
        identifiers.extend(
            Identifier(type=Namespace.NAME, value=name)
            for name in dict.fromkeys(participant['names'])
        )
        return identifiers
    identifiers = [Identifier(type=Namespace.GENESYMBOL, value=participant['id'])]
    if kind == 'protein':
        # One UniProt accession, or none where the participant has several.
        accessions = list(
            dict.fromkeys(
                accession
                for database, accession in participant['unification_xrefs']
                if database == 'uniprot'
            )
        )
        if len(accessions) == 1 and _UNIPROT_ACCESSION.fullmatch(accessions[0]):
            identifiers.insert(
                0, Identifier(type=Namespace.UNIPROT, value=accessions[0])
            )
    return identifiers


_ENTITY_TYPES = {'chemical': ChemicalEntity, 'protein': Protein, 'gene': Gene}


def _participant(participant: dict[str, Any]) -> Entity:
    kind = participant_kind(participant)
    return Entity(
        type=_ENTITY_TYPES[kind],
        identifiers=_identifiers(participant, kind),
        annotations=(
            None
            if kind == 'chemical'
            else [Annotation(term=slots.in_taxon, value=_TAXON)]
        ),
    )


def _relation_type(row: dict[str, Any]) -> _Type | None:
    name = row['interaction_type']
    if name not in RELATION_TYPES:
        raise ValueError(f'Unknown Pathway Commons relation type: {name!r}')
    return RELATION_TYPES[name]


def _data_sources(row: dict[str, Any]) -> list[Annotation]:
    annotations = []
    for name in row['data_sources']:
        if name not in DATA_SOURCES:
            raise ValueError(f'Unknown Pathway Commons data source: {name!r}')
        annotations.append(
            Annotation(term=slots.supporting_data_source, value=DATA_SOURCES[name])
        )
    return annotations


def _publications(row: dict[str, Any]) -> list[Annotation]:
    # A value that is not a PubMed ID is not a publication.
    return [
        Annotation(term=slots.publications, value=curie)
        for curie in dict.fromkeys(map(pmid_curie, row['pubmed']))
        if curie
    ]


def exclusion_reason(row: dict[str, Any]) -> str | None:
    """Why a line gives no relation; None for a line that is mapped.

    A participant of another class than the type expects raises.
    """
    relation_type = _relation_type(row)
    if relation_type is None:
        return 'relation type without a faithful predicate'
    for side, expected in (
        ('participant_a', relation_type.subject),
        ('participant_b', relation_type.object),
    ):
        kind = participant_kind(row[side])
        if (kind == 'chemical') != (expected == 'chemical'):
            raise ValueError(
                f'Pathway Commons {row["interaction_type"]!r} relation with a '
                f'{kind} as {side}: {row[side]["id"]!r}'
            )
        if not _identifiers(row[side], kind):
            return 'participant without an identifier'
    return None


# The SIF relation type as the file writes it. Several types become the same
# predicate (``controls-state-change-of``, ``chemical-affects`` and
# ``controls-production-of`` are all a bare ``affects``); the source's own
# type keeps them apart.
_INTERACTION_TYPE = 'pathwaycommons:interaction_type'


def map_interaction(row: dict[str, Any]) -> Relation | None:
    """One line of the extended SIF as a relation; None for a type left out."""
    if exclusion_reason(row):
        return None
    relation_type = _relation_type(row)
    return Relation(
        subject=_participant(row['participant_a']),
        predicate=relation_type.predicate,
        object=_participant(row['participant_b']),
        annotations=[
            *(
                [
                    Annotation(
                        term=slots.object_aspect_qualifier, value=relation_type.aspect
                    )
                ]
                if relation_type.aspect is not None
                else []
            ),
            Annotation(term=_INTERACTION_TYPE, value=row['interaction_type']),
            *_data_sources(row),
            *_publications(row),
        ],
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url=_URL,
            filename=f'pc{VERSION}-hgnc.txt.gz',
            subfolder='pathwaycommons',
            large=True,
            ext='gz',
        ),
        mapper=map_interaction,
        raw_parser=iter_interactions,
    ),
)
