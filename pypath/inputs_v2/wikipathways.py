"""
Parse WikiPathways RDF data and emit interaction-centric Entity records.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    OntologyAnnotationCv,
    OntologyCv,
    ParticipantMetadataCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
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
    pubmed='37941138',
    primary_category='pathways',
    annotation_ontologies=(OntologyCv.WIKIPATHWAYS,),
    description=(
        'WikiPathways is a community-curated pathway database. '
        'This inputs_v2 module parses the current RDF pathway export and '
        'emits pathway and directed interaction entities.'
    ),
)


entity_type_map = {value.value: value for value in EntityTypeCv}
_TAXON_SCOPED_ENTITY_TYPES = {
    EntityTypeCv.PROTEIN,
    EntityTypeCv.GENE,
    EntityTypeCv.RNA,
    EntityTypeCv.DNA,
}

f = FieldConfig(
    extract={
        'ensembl_id': r'^(ENS[A-Z0-9]*\d+(?:\.\d+)?)$',
        'uniprot_id': r'^((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?)$',
        'entrez_id': r'^(\d+)$',
        'chebi': r'^(?:CHEBI:)?(\d+)$',
    },
    map={
        'entity_type': lambda value: entity_type_map.get(value, EntityTypeCv.PHYSICAL_ENTITY),
    },
)


def _member_taxon_id(prefix: str):
    def _value(row):
        entity_type = entity_type_map.get(
            row.get(f'{prefix}_entity_type'),
            EntityTypeCv.PHYSICAL_ENTITY,
        )
        return row.get('taxon_id') if entity_type in _TAXON_SCOPED_ENTITY_TYPES else None

    return _value


pathways_schema = EntityBuilder(
    entity_type=EntityTypeCv.PATHWAY,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.WIKIPATHWAYS, value=f('pathway_id')),
        CV(term=IdentifierNamespaceCv.WIKIPATHWAYS_VERSION, value=f('pathway_version_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('title')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('organism_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('taxon_id')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed_ids', delimiter=';')),
        CV(term=OntologyAnnotationCv.DEFINITION, value=f('description')),
    ),
)


def _member(prefix: str, role) -> Member:
    return Member(
        entity=EntityBuilder(
            entity_type=f(f'{prefix}_entity_type', map='entity_type'),
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.NAME, value=f(f'{prefix}_label')),
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f(f'{prefix}_uniprot', delimiter=';', extract='uniprot_id')),
                CV(term=IdentifierNamespaceCv.ENTREZ, value=f(f'{prefix}_entrez', delimiter=';', extract='entrez_id')),
                CV(term=IdentifierNamespaceCv.ENSEMBL, value=f(f'{prefix}_ensembl', delimiter=';', extract='ensembl_id')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f(f'{prefix}_chebi', delimiter=';', extract='chebi')),
                CV(term=IdentifierNamespaceCv.HMDB, value=f(f'{prefix}_hmdb', delimiter=';')),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f(f'{prefix}_kegg_compound', delimiter=';')),
                CV(
                    term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                    value=f(f'{prefix}_pubchem_compound', delimiter=';'),
                ),
                CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f(f'{prefix}_hgnc', delimiter=';')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=_member_taxon_id(prefix)),
            ),
            associations=AssociationsBuilder(
                AssociationBuilder(
                    object_entity_type=EntityTypeCv.PATHWAY,
                    object_identifier_type=IdentifierNamespaceCv.WIKIPATHWAYS,
                    object_identifier=f('pathway_term_accession'),
                ),
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
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('interaction_types', delimiter=';')),
        CV(term=IdentifierNamespaceCv.WIKIPATHWAYS_VERSION, value=f('pathway_version_id')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('taxon_id')),
    ),
    associations=AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=EntityTypeCv.PATHWAY,
            object_identifier_type=IdentifierNamespaceCv.WIKIPATHWAYS,
            object_identifier=f('pathway_id'),
        ),
    ),
    membership=MembershipBuilder(
        _member('source', ParticipantMetadataCv.SOURCE),
        _member('target', ParticipantMetadataCv.TARGET),
    ),
)


download = Download(
    url=current_rdf_url,
    filename='wikipathways_rdf_wp.zip',
    subfolder='wikipathways',
    large=True,
    ext='zip',
    default_mode='rb',
)


resource = Resource(
    config,
    pathways=Dataset(
        download=download,
        mapper=pathways_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='pathways', **kwargs),
    ),
    interactions=Dataset(
        download=download,
        mapper=interactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='interactions', **kwargs),
    ),
)
