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
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv


config = ResourceConfig(
    id=ResourceCv.HPO,
    name='Human Phenotype Ontology',
    url='https://hpo.jax.org/',
    license=LicenseCV.HPO,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33264411',
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
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('hpo_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('hpo_name')),
    ),
)

resource = Resource(
    config,
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
