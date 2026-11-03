"""
Map ncRDeathDB non-coding RNA-gene observations to Relation records.

ncRDeathDB, the successor of miRDeathDB, records the non-coding RNAs (microRNAs,
lncRNAs, snoRNAs) that act in programmed cell death and the genes they interact
with in 12 species. Each row with a gene becomes a relation from the RNA to the
gene. The literature-curated observations and the homology predictions that the
database adds for other species are separate datasets.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, MicroRNA, NoncodingRNAProduct, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.ncrdeathdb import (
    iter_homology_predictions,
    iter_interactions,
)
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
    id=ResourceCv.NCRDEATHDB,
    name='ncRDeathDB',
    url='http://www.rna-society.org/ncrdeathdb/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='26431463',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Literature-curated interactions of microRNAs, lncRNAs and snoRNAs with genes '
        'in apoptosis, autophagy and necrosis, and homology predictions for other species.'
    ),
    description=(
        'ncRDeathDB, the successor of miRDeathDB, archives non-coding RNAs that act '
        'in programmed cell death (apoptosis, autophagy, necrosis) and the genes they '
        'interact with, from the literature, in 12 species. For microRNAs it adds '
        'predictions for further species, made from the curated human entries by '
        'homology.'
    ),
)

# The original, http://www.rna-society.org/ncrdeathdb/data/allNcRNACelldeathData.txt,
# answers "Invalid Request" (HTTP 416) to the download client, which does not send
# a browser's User-Agent. OmniPath's rescued copy is byte-identical to it
# (SHA-256 5614ee98...feb5, 2486291 bytes).
download_table = Download(
    url='https://rescued.omnipathdb.org/ncRDeathDB_allNcRNACelldeathData.txt',
    filename='ncRDeathDB_allNcRNACelldeathData.txt',
    subfolder='ncrdeathdb',
    encoding='iso-8859-1',
)

# Every value of the ``Organism`` column with its ``tax_id`` column.
ORGANISMS = {
    'Homo sapiens': '9606',
    'Mus musculus': '10090',
    'Rattus norvegicus': '10116',
    'Macaca mulatta': '9544',
    'Pan troglodytes': '9598',
    'Bos taurus': '9913',
    'Canis lupus familiaris': '9615',
    'Gallus gallus': '9031',
    'Danio rerio': '7955',
    'Drosophila melanogaster': '7227',
    'Caenorhabditis elegans': '6239',
    'Xenopus': '262014',
}
# The miRBase organism code of each taxon above.
MIRBASE_CODES = {
    '9606': 'hsa',
    '10090': 'mmu',
    '10116': 'rno',
    '9544': 'mml',
    '9598': 'ptr',
    '9913': 'bta',
    '9615': 'cfa',
    '9031': 'gga',
    '7955': 'dre',
    '7227': 'dme',
    '6239': 'cel',
    '262014': 'xtr',
}
# Every value of the ``RNA Category`` column. Biolink has no class for lncRNAs
# or snoRNAs; both are non-coding RNA products.
RNA_TYPES = {
    'miRNA': MicroRNA,
    'lncRNA': NoncodingRNAProduct,
    'snoRNA': NoncodingRNAProduct,
}

f = FieldConfig()


def taxon(row: dict[str, Any]) -> str:
    """The taxon of the observation: the source's ``tax_id``, which agrees with ``Organism``."""
    if ORGANISMS.get(row['organism']) != row['tax_id']:
        raise ValueError(
            f'Unknown organism and taxon: {row["organism"]!r}, {row["tax_id"]!r}'
        )
    return f'NCBITaxon:{row["tax_id"]}'


def rna_taxon(row: dict[str, Any]) -> str | None:
    """The taxon of the RNA: that of the observation.

    A microRNA identified by a miRBase product of another organism (a viral or
    a human microRNA studied in a cell line of another species) has no taxon here.
    """
    result = taxon(row)
    if row.get('mirbase_mature') and row['mature_organism_codes'] != [MIRBASE_CODES[row['tax_id']]]:
        return None
    return result


def rna_type(row: dict[str, Any]):
    if row['rna_category'] not in RNA_TYPES:
        raise ValueError(f'Unknown RNA category: {row["rna_category"]!r}')
    return RNA_TYPES[row['rna_category']]


# A microRNA symbol with exactly one miRBase product is that product (its
# accession and name, the symbol as a synonym). Any other symbol, and every
# lncRNA and snoRNA, is identified by its name alone.
ncrna = EntityBuilder(
    entity_type=rna_type,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.MIRBASE_MATURE, value=lambda row: row.get('mirbase_mature')),
        CV(
            term=Namespace.NAME,
            value=lambda row: row.get('mirbase_mature_name') or row['rna_symbol'],
        ),
        CV(
            term=Namespace.SYNONYM,
            value=lambda row: row['rna_symbol'] if row.get('mirbase_mature_name') else None,
        ),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=rna_taxon)),
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

# The relation says that the RNA interacts with the gene. The death pathway, the
# tissue and the expression direction of the RNA are not properties of that
# relation and are not carried.
relation = RelationBuilder(
    subject=lambda row: ncrna.build(row) if row['rna_symbol'] else None,
    predicate=slots.interacts_with,
    object=lambda row: gene.build(row) if row['entrez_id'] or row['gene_symbol'] else None,
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=lambda row: [f'PMID:{p}' for p in row['pmids']]),
        CV(term=slots.description, value=lambda row: row['description']),
    ),
)



def interactions_mapper(row: dict[str, Any]):
    """The relation of a row, or None if the row names no RNA or gene. An unknown organism or RNA category raises."""
    # Errors raised inside the builders' callbacks are not propagated, so the
    # crosswalks are checked here.
    taxon(row)
    rna_type(row)
    return relation(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=download_table,
        mapper=interactions_mapper,
        raw_parser=iter_interactions,
    ),
    homology_predictions=Dataset(
        download=download_table,
        mapper=interactions_mapper,
        raw_parser=iter_homology_predictions,
    ),
)
