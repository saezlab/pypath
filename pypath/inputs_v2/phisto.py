from __future__ import annotations

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.parsers.phisto import iter_interactions # my parser

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2._publications import pmid_curie

def _text(row: dict[str, Any], column: str) -> str | None:
    """Return a stripped cell value, or None for a blank/missing cell."""
    value = str(row.get(column) or '').strip()
    return value or None

# we need to add topic Host-microbe in resource_tags.yaml.
config = ResourceConfig(
    id=ResourceCv.PHISTO,
    name='PHISTO',
    url='https://phisto.org/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='24103452', # change this for phisto
    topics=(Topic.DISEASE_PHENOTYPE,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL, KnowledgeOrigin.LITERATURE_CURATED),
    summary='Pathogen-host interaction data with identifiers and experimental evidence.',
    description='PHISTO is a comprehensive database on pathogen-host interactions. It contains data on interactions between pathogens and their host organisms, with detailed classification of both pathogens and host species.',
)

download_interactions = Download(
    url='https://raw.githubusercontent.com/AlicjaWojciechowska/phi_base/main/phi_data.csv', # I ll need help to change it, they don't have a normal download link, sorry for this
    filename='phi_data.csv',
    subfolder='phisto', # where to cache the file
    large=True,
    encoding='iso-8859-1',
)

def _pathogen_protein() -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.UNIPROT,
                value=lambda row: _text(row, 'Uniprot ID [pathogen]'),
            ),
            CV(
                term=Namespace.UNIPROT_ENTRY,
                value=lambda row: _text(row, 'Uniprot ID [human]'),
            ),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.in_taxon, # this is a slot for ncbi it
                value=lambda row: (
                    f"NCBITaxon:{_text(row, 'Taxonomy ID')}"
                    if _text(row, 'Taxonomy ID')
                    else None
                ),
            ),
            CV(
                term=slots.in_taxon_label, # this is a slot for human-readable name
                value=lambda row: _text(row, 'Pathogen'),
            ),
        ),
    )


def _human_protein() -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.UNIPROT,
                value=lambda row: _text(row, 'Human UniProt ID'),
            ),
            CV(
                term=Namespace.UNIPROT_ENTRY,
                value=lambda row: _text(row, 'Human Protein'),
            ),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value='NCBITaxon:9606'),
        ),
    )


interactions_schema = RelationBuilder(
    subject=_pathogen_protein(),
    predicate=slots.interacts_with,
    object=_human_protein(),
    annotations=AnnotationsBuilder(
        CV(
            term='phisto:experimental_method',
            value=lambda row: _text(row, 'Experimental Method'),
        ),
        CV(
            term=slots.publications,
            value=lambda row: pmid_curie(_text(row, 'Pubmed ID')),
        ),
    ),
)


resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=iter_interactions,
    ),
)
