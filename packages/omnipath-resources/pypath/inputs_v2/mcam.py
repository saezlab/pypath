"""
Parse the Mammalian Cell Adhesion Molecule (MCAM) database list of human cell adhesion molecules.

MCAM is a database of cell adhesion molecules (CAMs) of mouse, human and rat,
with their sequences, functions and expression. The human list, 148 proteins,
is a plain list: one term, the human cell adhesion molecules of MCAM, and each
listed protein is associated with it.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.mcam import ORGANISM, iter_list, iter_proteins
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Association,
    Entity,
    EntityRef,
    Identifier,
)

config = ResourceConfig(
    id=ResourceCv.MCAM,
    name='MCAM',
    full='Mammalian Cell Adhesion Molecule Database',
    url='https://app1.unmc.edu/mcam/',
    # The web site states no license (the article, Cancer Informatics 2008, is CC BY 3.0).
    # The legacy registry (resources.json) records "Unspecified NC-SA", a placeholder.
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='19259402',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    # The article (Cancer Informatics 2008, "Data collection", "Evaluation of data
    # and classification of CAMs"): "The entries which did not have CAM related
    # annotations in UniProtKB/Swiss-Prot were validated manually for CAMs using
    # PubMed literature searches. Entries not validated as CAMs were removed".
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='The list of human cell adhesion molecules of MCAM.',
    description=(
        'MCAM is a Web-based database of mammalian cell adhesion molecules '
        '(CAMs), built to accelerate the identification of functional CAMs. It '
        'gives mouse, human and rat CAMs with their sequence information and '
        'characteristics such as gene functions and expression patterns in '
        'normal and tumor tissues and cell lines. The list of human CAMs is '
        'imported.'
    ),
)

download = Download(
    url='https://app1.unmc.edu/mcam/exel.cfm?seq=HumanUniprot',
    filename='mcam_human_uniprot.tsv',
    subfolder='mcam',
    default_mode='r',
)

LIST_NAMESPACE = 'mcam_list'
LIST_ID = 'human_cell_adhesion_molecule'
LIST_NAME = 'MCAM human cell adhesion molecule'
TAXON = 'NCBITaxon:9606'


def map_protein(row: dict[str, Any]) -> Entity:
    """A listed protein, associated with the list. Its accessions are those of
    its UniProt entry: the first is the primary, the others are secondary."""
    if row['organism'] != ORGANISM:
        raise ValueError(f'MCAM protein of organism {row["organism"]!r}')
    primary, *secondary = row['accessions']
    return Entity(
        type=Protein,
        identifiers=[
            Identifier(type=Namespace.UNIPROT, value=primary),
            *(Identifier(type=Namespace.UNIPROT_SEC, value=a) for a in secondary),
            Identifier(type=Namespace.UNIPROT_ENTRY, value=row['entry_name']),
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        associations=[
            Association(
                object=EntityRef(
                    type=OntologyClass,
                    identifier_type=LIST_NAMESPACE,
                    identifier=LIST_ID,
                ),
                predicate=slots.associated_with,
            )
        ],
    )


def map_term(row: dict[str, Any]) -> Entity:
    """The list as a term. Every row of the file belongs to the one list, so the
    term is the same for each row; the dataset keeps it to the first."""
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=LIST_NAMESPACE, value=LIST_ID),
            Identifier(type=Namespace.NAME, value=LIST_NAME),
        ],
        annotations=[
            Annotation(
                term=slots.description,
                value='The human cell adhesion molecules of the MCAM database.',
            )
        ],
    )


resource = Resource(
    config,
    terms=Dataset(
        download=download,
        mapper=map_term,
        raw_parser=iter_list,
    ),
    proteins=Dataset(
        download=download,
        mapper=map_protein,
        raw_parser=iter_proteins,
    ),
)
