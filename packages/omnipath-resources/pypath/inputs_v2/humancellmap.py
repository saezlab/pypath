"""
Parse the Human Cell Map and emit the predicted localizations of its preys.

The Human Cell Map (Go et al. 2021) is a BioID proximity-labelling map of the
HEK293 cell: 192 baits that mark compartments of the cell were used to label
the proteins near them, and the proteins that were significant for a bait (the
preys) were placed in the map by two methods, a non-negative matrix
factorization (NMF) and a spatial analysis of functional enrichment (SAFE).
Each prey has the compartments of the map that each method predicts for it. The
file read here is the list of the 4,424 significant preys.

A prey is related to each compartment label of each method with ``located_in``,
with the method, and the rank (NMF) or domain (SAFE) that the prey is assigned
to. The labels are the map's own, and are emitted as terms. The bait-prey
proximity interactions of the map (the SAINTexpress report, ``saint-latest.txt``,
81 MB) are not imported: the legacy reader did not read them.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.humancellmap import (
    NMF,
    SAFE,
    iter_localizations,
    iter_terms,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Relation

config = ResourceConfig(
    id=ResourceCv.HUMANCELLMAP,
    name='HumanCellMap',
    full='Human Cell Map',
    url='https://humancellmap.org/',
    # "Data from the human cell map is licensed under a Creative Commons
    # Attribution 4.0 International License" https://humancellmap.org/help/download
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='34079125',
    topics=(Topic.SUBCELLULAR_LOCALIZATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL, KnowledgeOrigin.COMPUTATIONALLY_PREDICTED),
    summary=(
        'Subcellular localizations of human proteins predicted from a BioID '
        'proximity-labelling map of the HEK293 cell, by two methods (NMF and SAFE).'
    ),
    description=(
        'The Human Cell Map is a proximity-dependent biotinylation (BioID) map '
        'of a human cell. 192 baits that mark compartments, organelles and '
        'structures of the HEK293 cell were used to label neighbouring '
        'proteins, and the significant preys were assigned a localization from '
        'the map by non-negative matrix factorization and by spatial analysis '
        'of functional enrichment.'
    ),
)

# The map is of human cells (HEK293), and its preys are human proteins.
HUMAN = 'NCBITaxon:9606'
# The labels of the map are its own.
LOCALIZATION_NAMESPACE = 'humancellmap_localization'
METHODS = (NMF, SAFE)

# Source attributes of the statement: the method, and the rank of the NMF or
# the domain of the SAFE network the prey is assigned to (numbers of the map).
_METHOD = 'humancellmap:method'
_GROUPS = {NMF: 'humancellmap:nmf_rank', SAFE: 'humancellmap:safe_domain'}


def _prey(row: dict[str, Any]) -> Entity:
    """A protein for a prey with a UniProt accession, a gene for the rest."""
    if row['uniprot']:
        entity_type = Protein
        identifiers = [Identifier(type=Namespace.UNIPROT, value=row['uniprot'])]
    else:
        entity_type = Gene
        identifiers = []
    if row['entrez']:
        identifiers.append(Identifier(type=Namespace.ENTREZ, value=row['entrez']))
    if row['hgnc']:
        identifiers.append(Identifier(type=Namespace.HGNC, value=row['hgnc']))
    identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=row['symbol']))
    return Entity(
        type=entity_type,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=HUMAN)],
    )


def map_localization(row: dict[str, Any]) -> Relation:
    """A prey located in a compartment of the map, by one method."""
    if row['method'] not in METHODS:
        raise ValueError(f'Unknown Human Cell Map method: {row["method"]!r}')
    return Relation(
        subject=_prey(row),
        predicate=slots.located_in,
        object=Entity(
            type=OntologyClass,
            identifiers=[Identifier(type=LOCALIZATION_NAMESPACE, value=row['label'])],
        ),
        annotations=[
            Annotation(term=_METHOD, value=row['method']),
            Annotation(term=_GROUPS[row['method']], value=row['group']),
        ],
    )


def map_term(row: dict[str, Any]) -> Entity:
    """A compartment label of the map."""
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=LOCALIZATION_NAMESPACE, value=row['label']),
            Identifier(type=Namespace.NAME, value=row['label']),
        ],
    )


# Version 1 of the files, which is also the latest (the two files are identical).
preys_download = Download(
    url='https://humancellmap.org/resources/downloads/preys-v1.txt',
    filename='humancellmap_preys-v1.txt',
    subfolder='humancellmap',
    large=True,
)

resource = Resource(
    config,
    terms=Dataset(
        download=preys_download,
        mapper=map_term,
        raw_parser=iter_terms,
    ),
    localizations=Dataset(
        download=preys_download,
        mapper=map_localization,
        raw_parser=iter_localizations,
    ),
)
