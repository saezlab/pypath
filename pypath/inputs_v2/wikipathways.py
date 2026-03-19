"""
Parse WikiPathways RDF data and emit Entity records.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    ParticipantMetadataCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.wikipathways import _raw, current_rdf_url


config = ResourceConfig(
    id=ResourceCv.WIKIPATHWAYS,
    name='WikiPathways',
    url='https://www.wikipathways.org/',
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='22073070',
    description=(
        'WikiPathways is a community-curated pathway database. '
        'This inputs_v2 module parses the current RDF pathway export and '
        'emits pathway and directed interaction entities.'
    ),
)


entity_type_map = {value.value: value for value in EntityTypeCv}

f = FieldConfig(
    map={
        'entity_type': lambda value: entity_type_map.get(value, EntityTypeCv.PHYSICAL_ENTITY),
    },
)


pathways_schema = EntityBuilder(
    entity_type=EntityTypeCv.PATHWAY,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.WIKIPATHWAYS, value=f('pathway_id')),
        CV(term=IdentifierNamespaceCv.WIKIPATHWAYS_VERSION, value=f('pathway_version_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('title')),
        CV(term=IdentifierNamespaceCv.SYSTEMATIC_NAME, value=f('pathway_uri')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('description')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('taxon_id')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed_ids', delimiter=';')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('ontology_terms', delimiter=';')),
    ),
)


def _member(prefix: str, role) -> Member:
    return Member(
        entity=EntityBuilder(
            entity_type=f(f'{prefix}_entity_type', map='entity_type'),
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.NAME, value=f(f'{prefix}_label')),
                CV(term=IdentifierNamespaceCv.SYSTEMATIC_NAME, value=f(f'{prefix}_uri')),
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f(f'{prefix}_uniprot', delimiter=';')),
                CV(term=IdentifierNamespaceCv.ENTREZ, value=f(f'{prefix}_entrez', delimiter=';')),
                CV(term=IdentifierNamespaceCv.ENSEMBL, value=f(f'{prefix}_ensembl', delimiter=';')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f(f'{prefix}_chebi', delimiter=';')),
                CV(term=IdentifierNamespaceCv.HMDB, value=f(f'{prefix}_hmdb', delimiter=';')),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f(f'{prefix}_kegg_compound', delimiter=';')),
                CV(
                    term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                    value=f(f'{prefix}_pubchem_compound', delimiter=';'),
                ),
                CV(term=IdentifierNamespaceCv.HGNC, value=f(f'{prefix}_hgnc', delimiter=';')),
            ),
        ),
        annotations=AnnotationsBuilder(
            CV(term=role),
        ),
    )


interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('interaction_local_id')),
        CV(term=IdentifierNamespaceCv.SYSTEMATIC_NAME, value=f('interaction_uri')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('interaction_types', delimiter=';')),
        CV(term=IdentifierNamespaceCv.WIKIPATHWAYS, value=f('pathway_id')),
        CV(term=IdentifierNamespaceCv.WIKIPATHWAYS_VERSION, value=f('pathway_version_id')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('taxon_id')),
    ),
    membership=MembershipBuilder(
        _member('source', ParticipantMetadataCv.SOURCE),
        _member('target', ParticipantMetadataCv.TARGET),
    ),
)


resource = Resource(
    config,
    pathways=Dataset(
        download=Download(
            url=current_rdf_url,
            filename='wikipathways_rdf_wp.zip',
            subfolder='wikipathways',
            large=True,
            ext='zip',
            default_mode='rb',
        ),
        mapper=pathways_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='pathways', **kwargs),
    ),
    interactions=Dataset(
        download=Download(
            url=current_rdf_url,
            filename='wikipathways_rdf_wp.zip',
            subfolder='wikipathways',
            large=True,
            ext='zip',
            default_mode='rb',
        ),
        mapper=interactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='interactions', **kwargs),
    ),
)
