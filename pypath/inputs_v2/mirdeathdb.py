"""
Map miRDeathDB microRNA-target observations to Relation records.

miRDeathDB records the microRNAs that act in programmed cell death (apoptosis,
autophagy) and the genes they target, one literature-curated observation per row
with the organism, tissue and PubMed ID. Each row becomes a relation from the
microRNA to the target gene.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, MicroRNA, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.mirdeathdb import iter_interactions
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.MIRDEATHDB,
    name='miRDeathDB',
    url='http://www.rna-world.org/mirdeathdb/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='22743998',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Literature-curated microRNA-target gene pairs from studies of apoptosis and autophagy.',
    description=(
        'miRDeathDB is a database of microRNAs involved in programmed cell death '
        'and of their target genes. Each entry gives the microRNA, the target '
        'gene, the death pathway (apoptosis or autophagy), the organism, the '
        'tissue and the publication. It was the predecessor of ncRDeathDB, and '
        'its last update was in 2012.'
    ),
)

# The site is gone; the file is the copy that OmniPath rescued.
download_interactions = Download(
    url='https://rescued.omnipathdb.org/miRDeathDB_all_data.txt',
    filename='miRDeathDB_all_data.txt',
    subfolder='mirdeathdb',
)

# Every value of the ``Organism`` column with its ``tax_id`` column.
ORGANISMS = {
    'human': '9606',
    'mouse': '10090',
    'rat': '10116',
    'Drosophila': '7227',
    'Xenopus': '262014',
}
# The miRBase organism code of each taxon above.
MIRBASE_CODES = {
    '9606': 'hsa',
    '10090': 'mmu',
    '10116': 'rno',
    '7227': 'dme',
    '262014': 'xtr',
}

f = FieldConfig(
    transform={'pmid': lambda v: f'PMID:{v}'},
)


def taxon(row: dict[str, Any]) -> str:
    """The taxon of the observation: the source's ``tax_id``, which agrees with ``Organism``."""
    if ORGANISMS.get(row['organism']) != row['tax_id']:
        raise ValueError(
            f'Unknown organism and taxon: {row["organism"]!r}, {row["tax_id"]!r}'
        )
    return f'NCBITaxon:{row["tax_id"]}'


def mirna_taxon(row: dict[str, Any]) -> str | None:
    """The taxon of the microRNA: that of the observation.

    A microRNA identified by a miRBase product of another organism (a viral or
    a human microRNA studied in a cell line of another species) has no taxon here.
    """
    result = taxon(row)
    if row['mirbase_mature'] and row['mature_organism_codes'] != [MIRBASE_CODES[row['tax_id']]]:
        return None
    return result


# A symbol with exactly one miRBase product is that product (its accession and
# name, the symbol as a synonym). A symbol with several products stands for the
# microRNA the source names; the products are not asserted one by one.
microrna = EntityBuilder(
    entity_type=MicroRNA,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.MIRBASE_MATURE, value=f('mirbase_mature')),
        CV(
            term=Namespace.NAME,
            value=lambda row: row['mirbase_mature_name'] or row['rna_symbol'],
        ),
        CV(
            term=Namespace.SYNONYM,
            value=lambda row: row['rna_symbol'] if row['mirbase_mature_name'] else None,
        ),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=mirna_taxon)),
)

# The Entrez gene ID identifies the gene; the symbol and the names come after.
gene = EntityBuilder(
    entity_type=Gene,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.ENTREZ, value=f('entrez_id')),
        CV(term=Namespace.HGNC, value=f('hgnc')),
        CV(term=Namespace.ENSG, value=f('ensembl_gene')),
        CV(term=Namespace.GENESYMBOL, value=f('gene_symbol')),
        CV(term=Namespace.NAME, value=f('gene_name')),
        CV(term=Namespace.SYNONYM, value=f('gene_synonyms')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=taxon)),
)


def _participant(builder: EntityBuilder, key: str):
    def build(row: dict[str, Any]) -> Any:
        # A row without the participant's name (or, for the gene, symbol and Entrez ID) has no relation.
        return builder.build(row) if row.get(key) else None

    return build


# The relation says that the microRNA targets the gene. The pathway (the cell
# death process the study is about), the tissue and the expression direction of
# the microRNA are not properties of that relation and are not carried.
interactions_relation = RelationBuilder(
    subject=_participant(microrna, 'rna_symbol'),
    predicate=slots.interacts_with,
    object=lambda row: (
        gene.build(row) if row['entrez_id'] or row['gene_symbol'] else None
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=f('pmids', transform='pmid')),
    ),
)



def interactions_mapper(row: dict[str, Any]):
    """The relation of a row, or None if the row names no target gene. An unknown organism raises."""
    # Errors raised inside the builders' callbacks are not propagated, so the
    # crosswalk is checked here.
    taxon(row)
    return interactions_relation(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_mapper,
        raw_parser=iter_interactions,
    ),
)
