"""
Parse miRTarBase microRNA - target gene interactions and emit Relation records.

miRTarBase collects experimentally validated interactions between mature
microRNAs and their target genes, curated from the literature. Its tables list
one experimental support block per line: the miRNA name, the target gene symbol
and Entrez Gene ID, the organism of each, the experiments, the PubMed ID and a
support type that says whether the block supports a functional interaction.

A functional block becomes a relation, the miRNA ``affects`` the target gene.
The tables name no direction and no aspect of that effect, so there is none in
the record. A non-functional block is a negative result and produces no
relation.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, MicroRNA, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.mirtarbase import iter_interactions
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.MIRTARBASE,
    name='miRTarBase',
    url='https://awi.cuhk.edu.cn/miRTarBase/',
    # The LICENSE file of the download directory says the data are "provided
    # at no cost in the public domain" (checked when the module was written,
    # October 2026); the legacy registry recorded an academic-only license.
    license=LicenseCV.PUBLIC,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='39578692',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Experimentally supported microRNA - target gene interactions of all '
        'species with their publications and support strength; '
        'non-functional results are left out.'
    ),
    description=(
        'miRTarBase is a database of experimentally validated interactions '
        'between microRNAs and their target genes, curated from the '
        'literature by the Warshel Institute for Computational Biology of the '
        'Chinese University of Hong Kong, Shenzhen. Each interaction lists '
        'the experiments that support it, such as reporter assays, Western '
        'blots, qPCR, microarrays or sequencing, with the publication, and '
        'whether the block supports a functional interaction, strongly or '
        'weakly, or a non-functional one. This module imports the tables of '
        'strong evidence (reporter assay or Western blot) and of weak '
        'evidence other than CLIP-Seq. The table of weak CLIP-Seq evidence is '
        'not imported.'
    ),
)

# The release the download page calls 11.0 still serves its files from the
# directory of release 10.0.
_BASE_URL = 'https://awi.cuhk.edu.cn/miRTarBase/downloads/files/10.0/'

# Support Type -> whether the block supports an interaction. The "(Weak)"
# variants are the source's own label for weaker support.
FUNCTIONAL = frozenset({'Functional MTI', 'Functional MTI (Weak)'})
# A negative result: it must never become an interaction.
NON_FUNCTIONAL = frozenset({'Non-Functional MTI', 'Non-Functional MTI (Weak)'})

# Organism codes of the "Species" columns -> NCBI Taxonomy. The target codes
# agree with the organism of the Entrez Gene record of each target. The codes
# of the four miRNA-only organisms (hcmv, mdv1, sja, vun) are miRBase codes.
TAXA = {
    'ath': '3702',  # Arabidopsis thaliana
    'bmo': '7091',  # Bombyx mori
    'bta': '9913',  # Bos taurus
    'cel': '6239',  # Caenorhabditis elegans
    'cfa': '9615',  # Canis lupus familiaris
    'cgr': '10029',  # Cricetulus griseus
    'chi': '9925',  # Capra hircus
    'dme': '7227',  # Drosophila melanogaster
    'dre': '7955',  # Danio rerio
    'ebv': '10376',  # human gammaherpesvirus 4 (Epstein-Barr virus)
    'eca': '9796',  # Equus caballus
    'gga': '9031',  # Gallus gallus
    'ggo': '9595',  # Gorilla gorilla gorilla
    'gma': '3847',  # Glycine max
    'hcmv': '10359',  # Human betaherpesvirus 5
    'hsa': '9606',  # Homo sapiens
    'kshv': '37296',  # Human gammaherpesvirus 8
    'mdv1': '10390',  # Gallid alphaherpesvirus 2
    'mml': '9544',  # Macaca mulatta
    'mmu': '10090',  # Mus musculus
    'oar': '9940',  # Ovis aries
    'ocu': '9986',  # Oryctolagus cuniculus
    'ppa': '9597',  # Pan paniscus
    'ptr': '9598',  # Pan troglodytes
    'rno': '10116',  # Rattus norvegicus
    'sja': '6182',  # Schistosoma japonicum
    'sly': '4081',  # Solanum lycopersicum
    'ssc': '9823',  # Sus scrofa
    'tgu': '59729',  # Taeniopygia guttata
    'vsv': '11277',  # Vesicular stomatitis Indiana virus
    'vun': '3917',  # Vigna unguiculata
    'xtr': '8364',  # Xenopus tropicalis
}


def _taxon(column: str):
    def curie(row: dict[str, Any]) -> str:
        code = row.get(column)
        if code not in TAXA:
            raise ValueError(f'Unknown miRTarBase organism code: {code!r}')
        return f'NCBITaxon:{TAXA[code]}'

    return curie


f = FieldConfig(
    extract={
        # Stricter than the shared PubMed identifier rule (1 to 8 digits): the
        # file writes 0 and 1 where a block has no publication (MIRT733936,
        # hsa-miR-665 and SDC1, has 1, which is PubMed's 1975 paper on a formate
        # assay in methanol poisoning), and every other number of the tables has
        # seven or eight digits.
        'pmid': r'^([1-9]\d{6,7})$',
    },
    transform={
        'pmid_curie': lambda v: f'PMID:{v}',
    },
)

# The source names the miRNA by its mature name only (no miRBase accession) and
# does not translate it between miRBase versions; the -3p and -5p forms stay
# separate names. The organism code is the column's, whatever the name's prefix.
mirna = EntityBuilder(
    entity_type=MicroRNA,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('mirna')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=_taxon('mirna_species')),
    ),
)

# The Entrez Gene ID and the symbol identify the gene, not its transcript.
target_gene = EntityBuilder(
    entity_type=Gene,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.GENESYMBOL, value=f('target_gene')),
        CV(term=Namespace.ENTREZ, value=f('entrez_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=_taxon('target_species')),
    ),
)

interactions_schema = RelationBuilder(
    subject=mirna,
    predicate=slots.affects,
    object=target_gene,
    annotations=AnnotationsBuilder(
        # The source's own label of the support strength, as written.
        CV(term=slots.has_confidence_level, value=f('support_type')),
        CV(term=slots.publications, value=f('pmid', extract='pmid', transform='pmid_curie')),
    ),
)


def map_interaction(row: dict[str, Any]):
    """The relation of a functional block; None for a non-functional one."""
    support_type = row.get('support_type')
    if support_type in NON_FUNCTIONAL:
        return None
    if support_type not in FUNCTIONAL:
        raise ValueError(f'Unknown miRTarBase support type: {support_type!r}')
    return interactions_schema(row)


def _dataset(filename: str) -> Dataset:
    return Dataset(
        download=Download(
            url=_BASE_URL + filename,
            filename=filename,
            subfolder='mirtarbase',
            large=True,
        ),
        mapper=map_interaction,
        raw_parser=iter_interactions,
    )


resource = Resource(
    config,
    strong_evidence=_dataset('miRTarBase_SE_WR.csv'),
    weak_evidence_other=_dataset('miRTarBase_WE_Other.csv'),
)
