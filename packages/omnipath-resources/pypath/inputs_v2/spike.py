"""
Parse the SPIKE database export and emit complexes, families and relations.

SPIKE (Signaling Pathway Integrated Knowledge Engine) stores human signaling
pathways in an entity-relationship model. Its entities are genes (one object
for a gene, its RNA and its protein, identified by the Entrez Gene ID), small
molecules (ChEBI), complexes, families and "events" (processes, stimuli and
other general entities). Its relations are

- regulations: a source, an effect (promote/activate, inhibit/repress or
  unknown) and a target, with the biochemical mechanism;
- interactions: undirected physical links found by an experimental method;
- containments: the members of complexes and families.

Complexes and families are emitted with their stated members, and take part
in relations as themselves: nothing is expanded to members or to pairs.

A regulation gives one or two relations:

- its effect on the target, with the direction the source states. The target
  can be a complex or a family while the mechanism acts on one member (the
  "physical target"); the relation is with the target.
- for a mechanism that is a modification (phosphorylation, ubiquitination,
  cleavage ...), the modification of the physical target. A regulation with
  an unknown effect and such a mechanism ("ATM phosphorylates MCM2") gives
  this relation only.

The target of a regulation can be another regulation or an interaction: "CHEK2
inhibits the inhibition of TP53 by MDM2, by phosphorylating TP53". The effect
is then on that relation, not on the physical target, and only what the
source does to the physical target is emitted.

Every relation has a quality level from 1 (focused biochemical studies) to 4
(high-throughput screens and imports from other databases). The level is kept
on each relation and no level is left out: there is one dataset per level, so
that restricting the build to some levels is a choice of datasets. Elements
flagged as deleted are not emitted.

The file is the snapshot of 25 January 2012; SPIKE has not been updated since.
"""

from __future__ import annotations

import functools
import re
from typing import Any

from biolink_model.datamodel.model import (
    ChemicalEntity,
    DirectionQualifierEnum,
    GeneFamily,
    GeneOrGeneProductOrChemicalEntityAspectEnum,
    MacromolecularComplex,
    MicroRNA,
    NamedThing,
    Protein,
    ProteinFamily,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.spike import (
    QUALITY_LEVELS,
    iter_groups,
    iter_interactions,
    iter_regulations,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.SPIKE,
    name='SPIKE',
    full='Signaling Pathway Integrated Knowledge Engine',
    url='https://www.cs.tau.ac.il/~spike/',
    # The database is "freely available for download"; the license of the
    # SPIKE software, which covers the data in it, is for academic use.
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='21097778',
    topics=(Topic.SIGNALING,),
    content=(
        Content.CAUSAL_REGULATION,
        Content.PHYSICAL_INTERACTIONS,
        Content.PROTEIN_COMPLEXES,
        Content.ENZYME_PTM_RELATIONSHIPS,
    ),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'Signed regulations with their mechanisms, physical interactions, '
        'complexes and protein families of human signaling pathways.'
    ),
    description=(
        'SPIKE is a database of highly curated human signaling pathways from '
        'Tel Aviv University, focused on the DNA damage response, cell cycle, '
        'programmed cell death, stress-activated transcription factors, MAPK '
        'and immune signaling and hearing. Curators record regulations with '
        'their effect, biochemical mechanism, publications and a quality '
        'level. The database also holds relations imported from KEGG, '
        'Reactome and other pathway databases, and protein-protein '
        'interactions from IntAct and large two-hybrid and '
        'immunoprecipitation screens.'
    ),
)

download_database = functools.partial(
    Download,
    url='https://www.cs.tau.ac.il/~spike/download/LatestSpikeDB.xml.zip',
    filename='LatestSpikeDB.xml.zip',
    subfolder='spike',
    ext='zip',
    needed=['LatestSpikeDB.xml'],
    encoding=None,
    default_mode='rb',
)

# SPIKE holds human genes only.
TAXON = 'NCBITaxon:9606'
# SPIKE accessions (``TAU:C-000033``) have no Namespace member yet.
SPIKE_NAMESPACE = 'spike'

# ``effect`` of a regulation -> direction of the effect on the target.
EFFECTS = {
    '1': DirectionQualifierEnum.increased,  # promote/activate
    '2': DirectionQualifierEnum.decreased,  # inhibit/repress
    '3': None,  # unknown
}

_ASPECT = GeneOrGeneProductOrChemicalEntityAspectEnum
_UP = DirectionQualifierEnum.increased
_DOWN = DirectionQualifierEnum.decreased
# ``mechanism`` of a regulation that is a modification of the physical
# target -> the aspect modified and its direction. The spelling is the
# source's.
MODIFICATIONS = {
    'Phosphorylation': (_ASPECT.phosphorylation, _UP),
    'Dephosphorylation': (_ASPECT.phosphorylation, _DOWN),
    'Ubiquitination': (_ASPECT.ubiquitination, _UP),
    'Polyubiquitination': (_ASPECT.ubiquitination, _UP),
    'Deubiqitination': (_ASPECT.ubiquitination, _DOWN),
    'Methylation': (_ASPECT.methylation, _UP),
    'Demethylation': (_ASPECT.methylation, _DOWN),
    'Acetylation': (_ASPECT.acetylation, _UP),
    'Deacetylation': (_ASPECT.acetylation, _DOWN),
    'Sumoylation': (_ASPECT.sumoylation, _UP),
    'Desumoylation': (_ASPECT.sumoylation, _DOWN),
    'Hydroxylation': (_ASPECT.hydroxylation, _UP),
    'Dehydroxylation': (_ASPECT.hydroxylation, _DOWN),
    'Glycosylation': (_ASPECT.glycosylation, _UP),
    'S-glutathionylation': (_ASPECT.glutathionylation, _UP),
    'Deglutathionylation': (_ASPECT.glutathionylation, _DOWN),
    'PARylation': (getattr(_ASPECT, 'ADP-ribosylation'), _UP),
    'Proteolytic Processing': (_ASPECT.cleavage, _UP),
    'Protein Degradation': (_ASPECT.degradation, _UP),
}
# ``mechanism`` that names the level at which the effect is exerted -> the
# aspect of the target that the effect increases or decreases.
EFFECT_ASPECTS = {
    'Transcription Regulation': _ASPECT.expression,
    'mRNA stability': _ASPECT.stability,
}
# The other values of ``mechanism``. Neddylation has no aspect in Biolink.
OTHER_MECHANISMS = frozenset(
    {
        'Physical Interaction',
        'Translational Regulation',
        'Indirect',
        'Neddylation',
        'Deneddylation',
        'Other',
        'N/A',
    }
)
MECHANISMS = frozenset({*MODIFICATIONS, *EFFECT_ASPECTS, *OTHER_MECHANISMS})

# ``biologicalAssay`` of an interaction -> PSI-MI interaction detection method.
ASSAYS = {
    'Y2H': 'MI:0018',  # two hybrid
    'IP': 'MI:0019',  # coimmunoprecipitation
    'CHIP': 'MI:0402',  # chromatin immunoprecipitation assay
    'Other': None,
}

# ``dataSource``: the database a record was imported from, in the source's
# spelling. The other values are not databases.
DATA_SOURCES = {
    'IntAct DB': 'IntAct DB',
    'Kegg': 'Kegg',
    'Reactome': 'Reactome',
    'PID': 'PID',
    'Curie': 'Curie',
    'EntrezGene': 'EntrezGene',
    'Published research': None,
    'PromoterAnalysis': None,
    "Shiloh's lab": None,
    "Avraham's lab": None,
    'Naor lab': None,
}

# The statements of a regulation.
EFFECT = 'effect'
MODIFICATION = 'modification'
ACTION = 'action'

_GROUP_TYPES = {'Complex': MacromolecularComplex, 'Family': ProteinFamily}
# ``XRef/@db`` codes: 1 is Entrez Gene, 4 is ChEBI.
_ENTREZ_DB = '1'
_CHEBI_DB = '4'
_number = re.compile(r'[1-9]\d*')
# HGNC reserves the root ``MIR`` for microRNA genes; the XML does not mark
# the genes that do not code for a protein.
_mirna_symbol = re.compile(r'MIR(?:LET)?\d')


def _names(entity: dict[str, Any]) -> list[Identifier]:
    """Short name, and the long name as a synonym when it differs."""
    name, long_name = entity.get('name'), entity.get('long_name')
    identifiers = []
    if name:
        identifiers.append(Identifier(type=Namespace.NAME, value=name))
    if long_name and long_name != name:
        identifiers.append(Identifier(type=Namespace.SYNONYM, value=long_name))
    return identifiers


def _xref(entity: dict[str, Any], db: str) -> str | None:
    xref_id = entity.get('xref_id') or ''
    if entity.get('xref_db') == db and _number.fullmatch(xref_id):
        return xref_id
    return None


def _microrna_genes(group: dict[str, Any]) -> bool:
    """Whether every member of a group is a microRNA gene."""
    members = group.get('member_kinds') or ()
    return bool(members) and all(
        member.get('kind') == 'Gene'
        and _mirna_symbol.match(member.get('name') or '')
        for member in members
    )


def participant(entity: dict[str, Any] | None) -> Entity | None:
    """The entity record of a SPIKE entity; None if deleted or not an entity."""
    if not entity or entity.get('is_deleted'):
        return None
    kind = entity.get('kind')
    accession = entity.get('accession')
    if kind == 'Gene':
        # A gene stands for the gene, its RNA and its protein.
        symbol = entity.get('name')
        entrez = _xref(entity, _ENTREZ_DB)
        if not (entrez or symbol):
            return None
        identifiers = []
        if entrez:
            identifiers.append(Identifier(type=Namespace.ENTREZ, value=entrez))
        if symbol:
            identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=symbol))
        long_name = entity.get('long_name')
        if long_name and long_name != symbol:
            identifiers.append(Identifier(type=Namespace.NAME, value=long_name))
        return Entity(
            type=MicroRNA if _mirna_symbol.match(symbol or '') else Protein,
            identifiers=identifiers,
            annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        )
    if kind == 'SignallingMolecule':
        chebi = _xref(entity, _CHEBI_DB)
        identifiers = [Identifier(type=Namespace.CHEBI, value=chebi)] if chebi else []
        identifiers.extend(_names(entity))
        return Entity(type=ChemicalEntity, identifiers=identifiers) if identifiers else None
    if kind == 'Event':
        entity_type = NamedThing
    elif kind == 'Group':
        group_type = entity.get('group_type')
        if group_type not in _GROUP_TYPES:
            raise ValueError(f'Unknown SPIKE group type: {group_type!r}')
        entity_type = _GROUP_TYPES[group_type]
        if entity_type is ProteinFamily and _microrna_genes(entity):
            entity_type = GeneFamily
    else:
        # A regulation or an interaction, or a reference to nothing.
        return None
    if not accession:
        return None
    return Entity(
        type=entity_type,
        identifiers=[
            Identifier(type=SPIKE_NAMESPACE, value=accession),
            *_names(entity),
        ],
    )


def _data_source(row: dict[str, Any]) -> list[Annotation]:
    data_source = row.get('data_source')
    if data_source not in DATA_SOURCES:
        raise ValueError(f'Unknown SPIKE data source: {data_source!r}')
    database = DATA_SOURCES[data_source]
    if database is None:
        return []
    return [Annotation(term=slots.supporting_data_source, value=database)]


def _description(row: dict[str, Any]) -> list[Annotation]:
    description = row.get('description')
    return [Annotation(term=slots.description, value=description)] if description else []


def _evidence(row: dict[str, Any]) -> list[Annotation]:
    """Quality level, publications, data source and curator's note of a relation."""
    level = row.get('integrity')
    if level not in QUALITY_LEVELS:
        raise ValueError(f'Unknown SPIKE quality level: {level!r}')
    return [
        Annotation(term=slots.has_confidence_level, value=level),
        *(
            Annotation(term=slots.publications, value=f'PMID:{pmid}')
            for pmid in row.get('pubmed') or ()
            if is_pmid(pmid)
        ),
        *_data_source(row),
        *_description(row),
    ]


def regulation_statements(row: dict[str, Any]) -> list[str]:
    """Which relations a regulation gives; none for a deleted regulation."""
    mechanism, effect = row.get('mechanism'), row.get('effect')
    if mechanism not in MECHANISMS:
        raise ValueError(f'Unknown SPIKE mechanism: {mechanism!r}')
    if effect not in EFFECTS:
        raise ValueError(f'Unknown SPIKE effect: {effect!r}')
    if row.get('is_deleted'):
        return []
    modification = mechanism in MODIFICATIONS
    on_entity = (row.get('target') or {}).get('kind') not in ('Regulation', 'Interaction')
    statements = []
    if on_entity and (EFFECTS[effect] is not None or not modification):
        statements.append(EFFECT)
    if modification:
        statements.append(MODIFICATION)
    elif not on_entity:
        statements.append(ACTION)
    return statements


def map_regulation(row: dict[str, Any]) -> Relation | None:
    """One statement of a regulation: the source affects the target."""
    statement = row.get('STATEMENT')
    if statement is None or row.get('is_deleted'):
        return None
    mechanism = row.get('mechanism')
    if statement == EFFECT:
        target = row.get('target')
        aspect, direction = EFFECT_ASPECTS.get(mechanism), EFFECTS[row.get('effect')]
    elif statement == MODIFICATION:
        target = row.get('physical_target')
        aspect, direction = MODIFICATIONS[mechanism]
    elif statement == ACTION:
        target = row.get('physical_target')
        aspect = direction = None
    else:
        raise ValueError(f'Unknown statement of a SPIKE regulation: {statement!r}')
    subject, object_ = participant(row.get('source')), participant(target)
    if subject is None or object_ is None:
        return None
    annotations = []
    if direction is not None:
        annotations.append(
            Annotation(term=slots.object_direction_qualifier, value=direction)
        )
    if aspect is not None:
        annotations.append(
            Annotation(term=slots.object_aspect_qualifier, value=aspect)
        )
    return Relation(
        subject=subject,
        predicate=slots.affects,
        object=object_,
        identifiers=[Identifier(type=SPIKE_NAMESPACE, value=row['accession'])],
        annotations=[*annotations, *_evidence(row)],
    )


def map_interaction(row: dict[str, Any]) -> Relation | None:
    """An undirected physical interaction of two SPIKE entities."""
    if row.get('is_deleted'):
        return None
    assay = row.get('biological_assay')
    if assay not in ASSAYS:
        raise ValueError(f'Unknown SPIKE biological assay: {assay!r}')
    first, second = participant(row.get('protein_a')), participant(row.get('protein_b'))
    if first is None or second is None:
        return None
    method = ASSAYS[assay]
    return Relation(
        subject=first,
        predicate=slots.interacts_with,
        object=second,
        identifiers=[Identifier(type=SPIKE_NAMESPACE, value=row['accession'])],
        annotations=[
            *(
                [Annotation(term=slots.has_evidence_of_type, value=method)]
                if method
                else []
            ),
            *_evidence(row),
        ],
    )


def map_group(row: dict[str, Any]) -> Entity | None:
    """A complex or a family with its members; None if deleted or without members."""
    group = participant(row)
    if group is None:
        return None
    members = [
        member
        for member in map(participant, row.get('members') or ())
        if member is not None
    ]
    if not members:
        return None
    return group._replace(
        annotations=[*_data_source(row), *_description(row)] or None,
        membership=[Membership(member=member) for member in members],
    )


def _relations(raw_parser, mapper, level: str) -> Dataset:
    return Dataset(
        download=download_database(),
        mapper=mapper,
        raw_parser=functools.partial(raw_parser, level=level),
    )


def _regulations(level: str) -> Dataset:
    return _relations(
        functools.partial(iter_regulations, statements_of=regulation_statements),
        map_regulation,
        level,
    )


def _groups(group_type: str) -> Dataset:
    return Dataset(
        download=download_database(),
        mapper=map_group,
        raw_parser=functools.partial(iter_groups, group_type=group_type),
    )


# One dataset per quality level, so that leaving out the low levels is a
# choice of datasets. Together they are all the relations of the file.
resource = Resource(
    config,
    complexes=_groups('Complex'),
    protein_families=_groups('Family'),
    **{
        f'regulations_quality_{level}': _regulations(level)
        for level in QUALITY_LEVELS
    },
    **{
        f'interactions_quality_{level}': _relations(iter_interactions, map_interaction, level)
        for level in QUALITY_LEVELS
    },
)
