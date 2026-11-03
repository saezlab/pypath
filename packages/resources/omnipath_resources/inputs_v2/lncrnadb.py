"""
Map lncRNAdb associations to Relation records.

lncRNAdb is a reference database of functionally characterised lncRNAs. Each entry
lists the molecules the lncRNA is associated with (proteins, genes, transcripts), with
a description and the publication. Each association becomes a relation from the
lncRNA to the partner.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import (
    Gene,
    NoncodingRNAProduct,
    Protein,
    RNAProduct,
    Transcript,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.lncrnadb import iter_associations
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
    id=ResourceCv.LNCRNADB,
    name='lncRNAdb',
    url='https://web.archive.org/web/20190505115633/http://lncrnadb.org/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.NAR_OPEN_ACCESS,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='25332394',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Literature-curated associations of functional lncRNAs with proteins, genes and transcripts.',
    description=(
        'lncRNAdb is a reference database of long non-coding RNAs with a '
        'characterised function. Its entries give the nomenclature, the species and '
        'the associations of the lncRNA with proteins, genes and other transcripts, '
        'each with a description and the publication. The site is gone; its domain '
        'now belongs to an unrelated website.'
    ),
)

# The REST service (http://lncrna.com/rest/all/nomenclature/species/association) and
# the site (lncrnadb.org) are gone: both domains now serve unrelated pages. The
# file is the copy that OmniPath rescued.
download_associations = Download(
    url='https://rescued.omnipathdb.org/lncrnadb_association.xml',
    filename='lncrnadb_association.xml',
    subfolder='lncrnadb',
    default_mode='rb',
    encoding=None,
)

# Every ``ComponentType`` value of an association that names its partner (the
# source writes some in capitals).
PARTNER_TYPES = {
    'Protein': Protein,
    'PROTEIN': Protein,
    'protein': Protein,
    'Gene': Gene,
    'Transcript': Transcript,
    'RNA': RNAProduct,
    'ncRNA': NoncodingRNAProduct,
}
# The ``Species`` value of every entry that has a single species and
# associations, with its taxon. ``Human herpesvirus`` is no taxon but a group of
# viruses, so it has none. An entry with several species ties no association to
# one of them, and has no taxon.
SPECIES = {
    'Homo sapiens (Human)': '9606',
    'Mus musculus (Mouse)': '10090',
    'Rattus norvegicus (Rat)': '10116',
    "Saccharomyces cerevisiae (Baker's yeast)": '4932',
    'Schizosaccharomyces pombe (Fission yeast)': '4896',
    'Arabidopsis thaliana (Mouse-ear cress)': '3702',
    'Xenopus laevis (African clawed frog)': '8355',
    'Gallus Gallus (Chicken)': '9031',
    'Dictyostelium discoideum (Slime mold)': '44689',
    'Paramecium tetraurelia': '5888',
    'Oxytricha trifallax (Sterkiella histriomuscorum)': '1172189',
    'Neurospora crassa': '5141',
    'Human herpesvirus': None,
}

f = FieldConfig()


def taxon(row: dict[str, Any]) -> str | None:
    """The taxon of the lncRNA of an entry with exactly one species."""
    if len(row['species']) != 1:
        return None
    species = row['species'][0]
    if species not in SPECIES:
        raise ValueError(f'Unknown species: {species!r}')
    return f'NCBITaxon:{SPECIES[species]}' if SPECIES[species] else None


def partner_type(row: dict[str, Any]):
    if row['partner_type'] not in PARTNER_TYPES:
        raise ValueError(f'Unknown partner type: {row["partner_type"]!r}')
    return PARTNER_TYPES[row['partner_type']]


# The source gives the lncRNA a name and aliases (some are Ensembl gene IDs).
lncrna = EntityBuilder(
    entity_type=NoncodingRNAProduct,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('lncrna')),
        CV(term=Namespace.SYNONYM, value=f('aliases')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=taxon)),
)

# The partner is known by its name only.
partner = EntityBuilder(
    entity_type=partner_type,
    identifiers=IdentifiersBuilder(CV(term=Namespace.NAME, value=f('partner'))),
)

# The source says only that the partner is associated with the lncRNA; the
# description tells how, in prose.
associations_relation = RelationBuilder(
    subject=lncrna,
    predicate=slots.interacts_with,
    object=lambda row: partner.build(row) if row['partner_type'] and row['partner'] else None,
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=lambda row: [f'PMID:{p}' for p in row['pmids']]),
        CV(term=slots.description, value=f('description')),
    ),
)



def associations_mapper(row: dict[str, Any]):
    """The relation of an association, or None if it names no partner. An unknown species or partner type raises."""
    # Errors raised inside the builders' callbacks are not propagated, so the
    # crosswalks are checked here.
    taxon(row)
    if row['partner_type'] or row['partner']:
        partner_type(row)
    return associations_relation(row)


resource = Resource(
    config,
    associations=Dataset(
        download=download_associations,
        mapper=associations_mapper,
        raw_parser=iter_associations,
    ),
)
