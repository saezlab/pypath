"""
Parse WikiPathways RDF data and emit interaction-centric Entity records.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    OntologyCv,
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
from pypath.internals.ontology_builder import OntologyBuilder
from pypath.internals.ontology_schema import OntologyDocument
from pypath.inputs_v2.base import Dataset, Download, OntologyDataset, Resource, ResourceConfig
from pypath.inputs_v2.parsers.wikipathways import _raw, current_rdf_url


config = ResourceConfig(
    id=ResourceCv.WIKIPATHWAYS,
    name='WikiPathways',
    url='https://www.wikipathways.org/',
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='22073070',
    primary_category='pathways',
    annotation_ontologies=(OntologyCv.WIKIPATHWAYS,),
    description=(
        'WikiPathways is a community-curated pathway database. '
        'This inputs_v2 module parses the current RDF pathway export and '
        'emits directed interaction entities plus a standalone pathway ontology export.'
    ),
)


entity_type_map = {value.value: value for value in EntityTypeCv}

f = FieldConfig(
    extract={
        'ensembl_id': r'^(ENS[A-Z0-9]*\d+(?:\.\d+)?)$',
        'uniprot_id': r'^((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?)$',
        'entrez_id': r'^(\d+)$',
    },
    map={
        'entity_type': lambda value: entity_type_map.get(value, EntityTypeCv.PHYSICAL_ENTITY),
    },
)


pathway_ontology_schema = OntologyBuilder(
    id='id',
    name='name',
    definition='definition',
    synonyms=f('synonyms', delimiter=';'),
    comments=f('comments', delimiter=';'),
    xrefs=f('xrefs', delimiter=';'),
)

pathway_ontology_document = OntologyDocument(
    ontology='wikipathways',
    default_namespace='wikipathways',
    remark='WikiPathways pathway ontology exported from the current RDF pathway archive via pypath.',
)


def _member(prefix: str, role) -> Member:
    return Member(
        entity=EntityBuilder(
            entity_type=f(f'{prefix}_entity_type', map='entity_type'),
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.NAME, value=f(f'{prefix}_label')),
                CV(term=IdentifierNamespaceCv.SYSTEMATIC_NAME, value=f(f'{prefix}_uri')),
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f(f'{prefix}_uniprot', delimiter=';', extract='uniprot_id')),
                CV(term=IdentifierNamespaceCv.ENTREZ, value=f(f'{prefix}_entrez', delimiter=';', extract='entrez_id')),
                CV(term=IdentifierNamespaceCv.ENSEMBL, value=f(f'{prefix}_ensembl', delimiter=';', extract='ensembl_id')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f(f'{prefix}_chebi', delimiter=';')),
                CV(term=IdentifierNamespaceCv.HMDB, value=f(f'{prefix}_hmdb', delimiter=';')),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f(f'{prefix}_kegg_compound', delimiter=';')),
                CV(
                    term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                    value=f(f'{prefix}_pubchem_compound', delimiter=';'),
                ),
                CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f(f'{prefix}_hgnc', delimiter=';')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('taxon_id')),
                CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('pathway_term_accession')),
                CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('pathway_ontology_terms', delimiter=';')),
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
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('pathway_term_accession')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('pathway_ontology_terms', delimiter=';')),
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
    pathway_ontology=OntologyDataset(
        download=download,
        mapper=pathway_ontology_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='pathway_terms', **kwargs),
        document=pathway_ontology_document,
        extension='obo',
        file_stem='wikipathways',
    ),
    interactions=Dataset(
        download=download,
        mapper=interactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='interactions', **kwargs),
    ),
)
