"""
Parse the NetPath PSI-MI export and emit physical interaction records.

NetPath distributed one PSI-MI XML 2.5 file per signaling pathway. Each
``<entry>`` holds one interaction of type "physical interaction" (MI:0218)
with its experiments (publications, host organism, detection method) and its
participants. The enzyme catalysis and transport reactions of NetPath are not
in the PSI-MI export.

One parsed row is one interaction in one experiment. An interaction of one
or two distinct interactors becomes an ``interacts_with`` relation; an
interaction of three or more is what NetPath calls a "complex interaction"
(the proteins were found in one complex, PMID 20067622) and becomes a
``MacromolecularComplex`` with its members. Nothing is expanded to pairs.

netpath.org has answered 503 since 2023, so the files come from the OmniPath
rescued-data server: the batch export of 3 September 2011 (26 pathways) and a
copy of the "browse" page that names the pathways.
"""

from __future__ import annotations

import functools
import re
from typing import Any

from biolink_model.datamodel.model import (
    MacromolecularComplex,
    Pathway,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.netpath import (
    UNIPROT_RE,
    iter_netpath_complexes,
    iter_netpath_interactions,
    iter_netpath_pathways,
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
    Identifier,
    Membership,
    Relation,
)

# Identifier systems without a Namespace member yet.
_HPRD = 'hprd'
_OMIM = 'omim'
_NETPATH = 'netpath'

# PSI-MI participant attribute without a Biolink slot; values are MI terms.
_EXPERIMENTAL_ROLE = 'psi_mi:experimental_role'

# interactionType -> predicate. The export has a single type.
INTERACTION_TYPE_PREDICATES = {
    'MI:0218': slots.interacts_with,  # physical interaction
}
# interactorType -> entity type.
INTERACTOR_TYPES = {
    'MI:0326': Protein,  # protein
}
# NetPath writes a missing OMIM number as "0".
_ACCESSION_PATTERNS = {
    'hprd': re.compile(r'\d{5}'),
    'entrez': re.compile(r'[1-9]\d*'),
    'refseq': re.compile(r'[NXY]P_\d+(?:\.\d+)?'),
    # The parser moves the values of the "GenBank" field that have the form
    # of a UniProt accession to ``uniprot``; the rest are GenBank protein
    # accessions (three letters, five to seven digits) or are rejected.
    'genbank': re.compile(r'[A-Z]{3}\d{5,7}(?:\.\d+)?'),
    'uniprot': UNIPROT_RE,
    'omim': re.compile(r'[1-9]\d{5}'),
}
_ACCESSION_NAMESPACES = (
    ('entrez', Namespace.ENTREZ),
    ('refseq', Namespace.REFSEQ_PROTEIN),
    ('hprd', _HPRD),
    ('genbank', Namespace.GENBANK),
    ('uniprot', Namespace.UNIPROT),
    ('omim', _OMIM),
)
# Interactors carry the identifiers of the human HPRD entry whatever organism
# NetPath gives them: the same HPRD and Entrez Gene records appear under
# "Human", "Mouse", "Rat" and "Mammalia". The organism is the species of the
# protein used in the experiment, so a taxon is stated only where it agrees
# with the identifiers.
_HUMAN_TAXON = '9606'


config = ResourceConfig(
    id=ResourceCv.NETPATH,
    name='NetPath',
    url='http://www.netpath.org/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_2_5,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='20067622',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    topics=(Topic.SIGNALING,),
    content=(
        Content.PHYSICAL_INTERACTIONS,
        Content.PROTEIN_COMPLEXES,
        Content.PATHWAYS,
    ),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Literature-curated physical interactions of signaling pathways, '
        'with their pathway, publications and detection methods.'
    ),
    description=(
        'NetPath is a manually curated resource of human signal transduction '
        'pathways, developed by the Pandey lab at Johns Hopkins University '
        'and the Institute of Bioinformatics, Bangalore. For each pathway, '
        'curators collected physical interactions, enzyme catalysis and '
        'transport reactions and transcriptionally regulated genes from the '
        'literature. The molecular interaction subset was distributed in '
        'PSI-MI format; the website is no longer online.'
    ),
)


def _accessions(participant: dict[str, Any], field: str) -> list[str]:
    """Well-formed accessions of one cross-reference field."""
    pattern = _ACCESSION_PATTERNS[field]
    return list(
        dict.fromkeys(
            value.strip()
            for value in participant.get(field) or []
            if pattern.fullmatch(value.strip())
        )
    )


def _participant_entity(participant: dict[str, Any]) -> Entity:
    if participant.get('type_mi') not in INTERACTOR_TYPES:
        raise ValueError(
            f'Unknown NetPath interactor type: {participant.get("type_mi")!r} '
            f'({participant.get("type")!r})'
        )
    identifiers = [
        Identifier(type=namespace, value=value)
        for field, namespace in _ACCESSION_NAMESPACES
        for value in _accessions(participant, field)
    ]
    # The short label of a NetPath molecule is its gene symbol.
    symbol = (participant.get('symbol') or '').strip()
    name = (participant.get('name') or '').strip()
    if symbol:
        identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=symbol))
    if name and name != symbol:
        identifiers.append(Identifier(type=Namespace.NAME, value=name))
    identifiers.extend(
        Identifier(type=Namespace.SYNONYM, value=synonym)
        for synonym in participant.get('synonyms') or []
        if synonym not in (symbol, name)
    )
    annotations = []
    if participant.get('taxon') == _HUMAN_TAXON:
        annotations.append(
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{_HUMAN_TAXON}')
        )
    return Entity(
        type=INTERACTOR_TYPES[participant['type_mi']],
        identifiers=identifiers,
        annotations=annotations,
    )


def _role_annotations(participant: dict[str, Any]) -> list[Annotation]:
    return [
        Annotation(term=_EXPERIMENTAL_ROLE, value=role)
        for role in participant.get('experimental_role_mi') or []
        if role
    ]


def _evidence_annotations(row: dict[str, Any]) -> list[Annotation]:
    """Method, publications and pathway of one observation."""
    annotations = []
    if row.get('detection_method_mi'):
        annotations.append(
            Annotation(
                term=slots.has_evidence_of_type,
                value=row['detection_method_mi'],
            )
        )
    annotations.extend(
        Annotation(term=slots.publications, value=curie)
        for pmid in row.get('pubmed') or []
        if (curie := pmid_curie(pmid))
    )
    if row.get('pathway_id'):
        annotations.append(
            Annotation(term=slots.has_topic, value=row['pathway_id'])
        )
    return annotations


def _host_annotations(row: dict[str, Any]) -> list[Annotation]:
    """The organism in which the experiment was done."""
    host = str(row.get('host_taxon') or '')
    if host.isdigit() and int(host) > 0:
        return [Annotation(term=slots.in_taxon, value=f'NCBITaxon:{host}')]
    return []


def _predicate(row: dict[str, Any]):
    interaction_type = row.get('interaction_type_mi')
    if interaction_type not in INTERACTION_TYPE_PREDICATES:
        raise ValueError(
            f'Unknown NetPath interaction type: {interaction_type!r} '
            f'({row.get("interaction_type")!r})'
        )
    return INTERACTION_TYPE_PREDICATES[interaction_type]


def map_netpath_interaction(row: dict[str, Any]) -> Relation | None:
    """A physical interaction of two proteins, or of a protein with itself.

    A source interaction with a single interactor is a homo-oligomer
    ("Interaction between IKBKG and IKBKG"): subject and object are the same.
    """
    participants = row.get('participants') or []
    if not 1 <= len(participants) <= 2:
        return None
    predicate = _predicate(row)
    first, second = participants[0], participants[-1]
    annotations = [*_evidence_annotations(row), *_host_annotations(row)]

    def entity(participant: dict[str, Any]) -> Entity:
        built = _participant_entity(participant)
        return built._replace(
            annotations=[*built.annotations, *_role_annotations(participant)]
        )

    return Relation(
        subject=entity(first),
        predicate=predicate,
        object=entity(second),
        annotations=annotations,
    )


def map_netpath_complex(row: dict[str, Any]) -> Entity | None:
    """The complex of a NetPath interaction with three or more interactors.

    The source gives the complex no accession and no name of its own; the
    members define it. As for an interaction, the organism in which the
    experiment was done is the complex's taxon; the members carry the taxon
    of their own identifiers (see ``_participant_entity``).
    """
    participants = row.get('participants') or []
    if len(participants) < 3:
        return None
    _predicate(row)
    return Entity(
        type=MacromolecularComplex,
        identifiers=[],
        annotations=[*_evidence_annotations(row), *_host_annotations(row)],
        membership=[
            Membership(
                member=_participant_entity(participant),
                annotations=_role_annotations(participant) or None,
            )
            for participant in participants
        ],
    )


def map_netpath_pathway(row: dict[str, Any]) -> Entity | None:
    pathway_id = (row.get('pathway_id') or '').strip()
    if not pathway_id:
        return None
    identifiers = [Identifier(type=_NETPATH, value=pathway_id)]
    name = (row.get('name') or '').strip()
    if name:
        identifiers.append(Identifier(type=Namespace.NAME, value=name))
    return Entity(type=Pathway, identifiers=identifiers)


_psimi_download = functools.partial(
    Download,
    # The 3 September 2011 batch export, 26 pathways. Was at
    # http://www.netpath.org/data/batch/NetPath_PSI-MI.zip
    url='https://rescued.omnipathdb.org/NetPath_PSI-MI.zip',
    filename='NetPath_PSI-MI.zip',
    subfolder='netpath',
    large=False,
    ext='zip',
    default_mode='rb',
    encoding=None,
)

resource = Resource(
    config,
    interactions=Dataset(
        download=_psimi_download(),
        mapper=map_netpath_interaction,
        raw_parser=iter_netpath_interactions,
    ),
    complexes=Dataset(
        download=_psimi_download(),
        mapper=map_netpath_complex,
        raw_parser=iter_netpath_complexes,
    ),
    pathways=Dataset(
        download=Download(
            # A copy of http://www.netpath.org/browse
            url='https://rescued.omnipathdb.org/netpath-browse.html',
            filename='netpath-browse.html',
            subfolder='netpath',
            large=False,
        ),
        mapper=map_netpath_pathway,
        raw_parser=iter_netpath_pathways,
    ),
)
