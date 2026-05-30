"""
Parse HPO data and emit Entity records.

This module converts Human Phenotype Ontology (HPO) gene-phenotype annotation
data into Entity records using the declarative schema pattern.

The HPO provides standardized vocabulary for describing phenotypic abnormalities
and their associations with genes.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    OntologyCv,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.inputs_v2.parsers.obo import iter_obo, obo_record_to_term


config = ResourceConfig(
    id=ResourceCv.HPO,
    name='Human Phenotype Ontology',
    url='https://hpo.jax.org/',
    license=LicenseCV.HPO,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33264411',
    primary_category='phenotypes',
    resource_kind='ontology',
    annotation_ontologies=(OntologyCv.HUMAN_PHENOTYPE_ONTOLOGY,),
    description=(
        'The Human Phenotype Ontology (HPO) provides a standardized vocabulary '
        'for describing phenotypic abnormalities encountered in human disease. '
        'This module provides gene-phenotype associations linking genes to HPO terms.'
    ),
)

f = FieldConfig()

annotations_schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.ENTREZ, value=f('ncbi_gene_id')),
        CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('gene_symbol')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value='9606'),
    ),
    associations=AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=EntityTypeCv.CV_TERM,
            object_identifier_type=IdentifierNamespaceCv.CV_TERM_ACCESSION,
            object_identifier=f('hpo_id'),
        ),
    ),
)


terms_schema = ontology_entity_mapper(obo_record_to_term, ontology_id='hpo')


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/hp.obo',
            filename='hp.obo',
            subfolder='hpo',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
    annotations=Dataset(
        download=Download(
            url='https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt',
            filename='genes_to_phenotype.txt',
            subfolder='hpo',
            large=True,
        ),
        mapper=annotations_schema,
        raw_parser=iter_tsv,
    ),
)
