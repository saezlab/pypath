"""PHI-base pathogen-host interaction phenotype data. Generally, the authors really like plants. 
But there is some info on homo sapiens too. We will take a subset of PHI-base which covers associations of 
taxon - its protein - phenotype in the host [human readable text with links to PHI-base identifiers]. """

from __future__ import annotations

import json

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phibase import iter_phi_phenotypes
from pypath.share.downloads import download_and_open


config = ResourceConfig(
    id=ResourceCv.PHIBASE,
    name='PHI-base',
    full='Pathogen-Host Interactions Database',
    url='https://phi5.phi-base.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='39588765',
    topics=(Topic.DISEASE_PHENOTYPE,),
    content=(Content.TRAIT_ASSOCIATIONS,),
    knowledge_origin=(
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.LITERATURE_CURATED,
    ),
    summary=(
        'Experimentally supported pathogen-host phenotype annotations for '
        'pathogen genes and genotypes.'
    ),
    description=(
        'PHI-base is a manually curated database of experimentally verified '
        'pathogenicity, virulence and effector genes from fungal, oomycete and '
        'bacterial pathogens infecting animal, plant, fungal and insect hosts.'
    ),
)


def _latest_release(force_refresh: bool = False) -> tuple[str, str]:
    """Return the URL and filename of the latest PHI-base 5 ZIP archive."""

    opener = download_and_open(
        'https://zenodo.org/api/records?q=conceptrecid:10722192&sort=mostrecent&size=1',
        'latest_release.json',
        'phibase',
        large=False,
        force_download=force_refresh,
    )
    record = json.loads(opener.result)['hits']['hits'][0]
    archive = next(file for file in record['files'] if file['key'].endswith('.zip'))

    return archive['links']['self'], archive['key']


f = FieldConfig()


def _pathogen_protein() -> EntityBuilder:
    """Build the pathogen protein and its source and taxonomy annotations."""

    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=f('pathogen_uniprot')),
            CV(term=Namespace.GENESYMBOL, value=f('pathogen_gene_name')),
        ),
        annotations=AnnotationsBuilder(
            CV(
                term=slots.in_taxon,
                value=lambda row: (
                    f"NCBITaxon:{row['pathogen_taxon_id']}"
                    if row.get('pathogen_taxon_id')
                    else None
                ),
            ),
            CV(
                term=slots.in_taxon_label,
                value=f('pathogen_organism'),
            ),
            CV(term='phibase:phig_id', value=f('pathogen_phig_id')),
            CV(term=slots.description, value=f('pathogen_product')),
            CV(term='phibase:strain', value=f('pathogen_strain')),
        ),
    )


phi_phenotype_schema = RelationBuilder(
    subject=_pathogen_protein(),
    predicate=slots.associated_with,
    object=EntityBuilder(
        entity_type=OntologyClass,
        identifiers=IdentifiersBuilder(
            CV(term='phibase', value=f('high_level_term')),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=f('publication')),
        CV(
            term=slots.in_taxon,
            value=lambda row: (
                f"NCBITaxon:{row['host_taxon_id']}"
                if row.get('host_taxon_id')
                else None
            ),
        ),
        CV(term=slots.in_taxon_label, value=f('host_organism')),
        CV(term='phibase:host_strain', value=f('host_strain')),
        CV(term='phibase:host_genotype_id', value=f('host_genotype_id')),
        CV(term='phibase:metagenotype_id', value=f('metagenotype_id')),
        CV(term='phibase:phenotype_term_id', value=f('phenotype_term_id')),
        CV(term='phibase:annotation_extension', value=f('annotation_extensions')),
        CV(term='phibase:condition', value=f('conditions')),
        CV(term='phibase:evidence_code', value=f('evidence_code')),
        CV(term='phibase:figure', value=f('figure')),
        CV(term='phibase:phi4_id', value=f('phi4_ids')),
        CV(term='phibase:curation_session_id', value=f('curation_session_id')),
        CV(term='phibase:curation_status', value=f('curation_status')),
        CV(term='phibase:checked', value=f('checked')),
    ),
)


resource = Resource(
    config,
    phi_phenotypes=Dataset(
        download=Download(
            url=lambda force_refresh=False, **_: _latest_release(force_refresh)[0],
            filename=lambda force_refresh=False, **_: _latest_release(force_refresh)[1],
            subfolder='phibase',
            large=True,
            ext='zip',
        ),
        mapper=phi_phenotype_schema,
        raw_parser=iter_phi_phenotypes,
    ),
)
