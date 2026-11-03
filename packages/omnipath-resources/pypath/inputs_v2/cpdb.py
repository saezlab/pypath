"""
Parse the ConsensusPathDB human protein interactions and emit records.

ConsensusPathDB integrates the interaction data of other databases and
distributes the protein interaction part of its network as one table
(``ConsensusPathDB_human_PPI.gz``). A line is one interaction with the
databases that report it, its publications, its participants and, for a
pair, a confidence score. The data access page of ConsensusPathDB says how
to read the participants:

- a line with one protein is a self-interaction, and becomes an
  ``interacts_with`` relation of the protein with itself;
- a line with two proteins is an interaction, and becomes an
  ``interacts_with`` relation;
- a line with more than two proteins is, in the words of the page, a
  complex interaction: the members are stated, and the line becomes a
  ``MacromolecularComplex`` with its members. It is never expanded to pairs.

The file does not say that a pair interacts directly, so no relation is
``directly_physically_interacts_with``. The databases that contribute the
interaction are repeated as ``supporting_data_source``.

A participant known only by genomic identifiers can stand for several UniProt
entries (the entry names are concatenated with dots, see the parser). Such a
participant is not one protein, and a line that has one is not mapped: for a
complex the others alone would misstate the composition.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import MacromolecularComplex, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2._measurements import measurement as _measurement
from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cpdb import iter_interactions
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.CPDB,
    name='CPDB',
    full='ConsensusPathDB',
    url='http://cpdb.molgen.mpg.de/',
    # The site (cpdb.molgen.mpg.de/CPDB/rlFrame, "Licensing information"):
    # "The use of ConsensusPathDB is free for academic users. Commercial users
    # should contact Dr. Ralf Herwig (herwig [at] molgen.mpg.de). Interaction
    # data from ConsensusPathDB is available under the license terms of each
    # of the contributing databases listed above."
    # The data access page adds that the complete integrated network cannot be
    # released "due to several licensing issues".
    license=LicenseCV.COMPOSITE,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='23143270',
    topics=(),
    content=(Content.PHYSICAL_INTERACTIONS, Content.PROTEIN_COMPLEXES),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'Human protein interactions, self-interactions and complex '
        'interactions that ConsensusPathDB integrates from other interaction '
        'databases, with the reporting databases, publications and, for '
        'pairs, a confidence score.'
    ),
    description=(
        'ConsensusPathDB (CPDB), developed at the Max Planck Institute for '
        'Molecular Genetics in Berlin, integrates interaction networks of '
        'Homo sapiens (protein interactions, genetic, metabolic, signaling, '
        'gene regulatory and drug-target interactions) and biochemical '
        'pathways from currently 31 public resources, merging identical '
        'interactions and recording the source databases of each. Only the '
        'protein interaction part of the network is available for download. '
        'Binary protein interactions are scored with an integrative '
        'confidence method (IntScore).'
    ),
)

# ConsensusPathDB release 35 (5 June 2021) is the file at this address; the
# file carries no other release marker.
_URL = 'http://cpdb.molgen.mpg.de/download/ConsensusPathDB_human_PPI.gz'
_TAXON = 'NCBITaxon:9606'

# ``source_databases``: every database the file names, in the file's spelling
# ("Manual upload" is what ConsensusPathDB writes for interactions its curators
# uploaded). A name that is not here raises.
SOURCE_DATABASES = {
    name: name
    for name in (
        'BIND', 'BioCarta', 'Biogrid', 'CORUM', 'DIP', 'HPRD', 'HumanCyc',
        'INOH', 'InnateDB', 'IntAct', 'MIPS-MPPI', 'MINT', 'Manual upload',
        'MatrixDB', 'NetPath', 'PDB', 'PDZBase', 'PID', 'PINdb',
        'PhosphoPOINT', 'Reactome', 'Spike',
    )
}

_UNIPROT_ACCESSION = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+)?'
)
_UNIPROT_ENTRY = re.compile(r'[A-Z0-9]+_[A-Z0-9]+')
_HGNC = re.compile(r'HGNC:(\d+)')
_NUMBER = re.compile(r'[1-9]\d*')
_ENSEMBL_GENE = re.compile(r'ENSG\d{11}')


def _is_ambiguous(participant: dict[str, list[str]]) -> bool:
    """Whether the participant stands for several UniProt entries."""
    return len(participant['uniprot_entry']) > 1 or len(participant['uniprot_id']) > 1


def _identifiers(participant: dict[str, list[str]]) -> list[Identifier]:
    """UniProt accession and entry name, then the gene identifiers.

    The gene identifiers (name, HGNC, Entrez, Ensembl) are of the gene of the
    protein; where the file gives several genes for one protein, none of them
    is used, since any of them could be another protein's.
    """
    identifiers = [
        Identifier(type=Namespace.UNIPROT, value=value)
        for value in participant['uniprot_id']
        if _UNIPROT_ACCESSION.fullmatch(value)
    ]
    identifiers.extend(
        Identifier(type=Namespace.UNIPROT_ENTRY, value=value)
        for value in participant['uniprot_entry']
        if _UNIPROT_ENTRY.fullmatch(value)
    )
    for key, namespace, pattern in (
        ('gene_name', Namespace.GENESYMBOL, None),
        ('hgnc_id', Namespace.HGNC, _HGNC),
        ('entrez_gene', Namespace.ENTREZ, _NUMBER),
        ('ensembl_gene', Namespace.ENSG, _ENSEMBL_GENE),
    ):
        values = participant[key]
        if len(values) != 1 or not values[0]:
            continue
        if pattern is None:
            identifiers.append(Identifier(type=namespace, value=values[0]))
            continue
        match = pattern.fullmatch(values[0])
        if match:
            identifiers.append(
                Identifier(
                    type=namespace,
                    value=match.group(1) if match.groups() else values[0],
                )
            )
    return identifiers


def _protein(participant: dict[str, list[str]]) -> Entity:
    return Entity(
        type=Protein,
        identifiers=_identifiers(participant),
        annotations=[Annotation(term=slots.in_taxon, value=_TAXON)],
    )


def _data_sources(row: dict[str, Any]) -> list[Annotation]:
    annotations = []
    for name in row['source_databases']:
        if name not in SOURCE_DATABASES:
            raise ValueError(f'Unknown ConsensusPathDB source database: {name!r}')
        annotations.append(
            Annotation(term=slots.supporting_data_source, value=SOURCE_DATABASES[name])
        )
    return annotations


def _publications(row: dict[str, Any]) -> list[Annotation]:
    # A token that is not a PubMed ID (``unassigned1304``, a long number, a
    # ``rtd:`` or other label) is no publication.
    return [
        Annotation(term=slots.publications, value=curie)
        for curie in dict.fromkeys(map(pmid_curie, row['publications']))
        if curie
    ]


def _confidence(row: dict[str, Any]) -> list[Annotation]:
    measurement = _measurement(
        row['confidence'], source_field='interaction_confidence'
    )
    return (
        [Annotation(term=slots.has_confidence_score, value=measurement)]
        if measurement
        else []
    )


def exclusion_reason(row: dict[str, Any]) -> str | None:
    """Why a line gives no record; None for a line that is mapped."""
    participants = row['participants']
    if not participants:
        return 'no participant'
    if any(_is_ambiguous(participant) for participant in participants):
        return 'participant with several UniProt entries'
    if any(not _identifiers(participant) for participant in participants):
        return 'participant without a well-formed identifier'
    return None


def map_interaction(row: dict[str, Any]) -> Relation | None:
    """A self-interaction (one protein) or an interaction (two) as a relation."""
    participants = row['participants']
    if len(participants) > 2 or exclusion_reason(row):
        return None
    return Relation(
        subject=_protein(participants[0]),
        predicate=slots.interacts_with,
        object=_protein(participants[-1]),
        annotations=[
            *_data_sources(row),
            *_publications(row),
            *_confidence(row),
        ],
    )


def map_complex(row: dict[str, Any]) -> Entity | None:
    """A complex interaction (more than two proteins) as a complex.

    ConsensusPathDB gives the complex neither an accession nor a name; its
    members define it. The confidence column is ``NA`` for these lines.
    """
    participants = row['participants']
    if len(participants) < 3 or exclusion_reason(row):
        return None
    return Entity(
        type=MacromolecularComplex,
        identifiers=[],
        annotations=[
            *_data_sources(row),
            *_publications(row),
            *_confidence(row),
        ],
        membership=[
            Membership(member=_protein(participant)) for participant in participants
        ],
    )


def _download() -> Download:
    return Download(
        url=_URL,
        filename='ConsensusPathDB_human_PPI.gz',
        subfolder='cpdb',
        large=True,
        ext='gz',
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=_download(),
        mapper=map_interaction,
        raw_parser=iter_interactions,
    ),
    complexes=Dataset(
        download=_download(),
        mapper=map_complex,
        raw_parser=iter_interactions,
    ),
)
