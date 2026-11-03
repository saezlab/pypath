"""
Map the LncRNADisease interaction table to Relation records.

LncRNADisease is a database of lncRNAs associated with diseases. Its first version
also published a table of literature-curated interactions of lncRNAs with proteins,
RNAs and DNA (and of lncRNA genes with transcription factors); this module reads
that table. Each row becomes a relation whose predicate follows the interaction
class of the source: binding, regulation or co-expression.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import (
    Gene,
    NoncodingRNAProduct,
    Protein,
    RNAProduct,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.lncdisease import iter_interactions
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
    id=ResourceCv.LNCDISEASE,
    name='LncRNADisease',
    url='http://www.rnanut.net/lncrnadisease/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='23175614',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Literature-curated interactions of lncRNAs with proteins, RNAs and DNA, and of lncRNA genes with transcription factors.',
    description=(
        'LncRNADisease collects lncRNAs associated with diseases from the '
        'literature. Version 1 also published a table of the interactions of '
        'lncRNAs with other molecules: binding, regulation and co-expression at the '
        'RNA level (proteins, RNAs) and at the DNA level (genes, and transcription '
        'factors at the promoters of lncRNA genes), with the species and the '
        'publication. This module reads that table; the associations with diseases '
        'that the current version (v3.0) offers for download are not part of it.'
    ),
)

# The first server (210.73.221.6) is gone and the current site, LncRNADisease v3.0,
# no longer publishes this table; the file is the copy that OmniPath rescued.
download_interactions = Download(
    url='https://rescued.omnipathdb.org/LncRNADisease_data2.txt',
    filename='LncRNADisease_data2.txt',
    subfolder='lncdisease',
)

# Every value of the ``species`` column with its taxon. The strain names are
# those of the NCBI taxonomy records of the strains.
SPECIES = {
    'Human': '9606',
    'Mouse': '10090',
    'Saccharomyces cerevisiae S288c': '559292',
    'Drosophila melanogaster': '7227',
    'Dictyostelium discoideum AX4': '352472',
    'Candida albicans SC5314': '237561',
    'Caenorhabditis elegans': '6239',
    'Arabidopsis thaliana': '3702',
}
# Biolink has no class for a lncRNA; the first partner of an ``RNA-`` level is
# the lncRNA, of a ``DNA-`` level the gene of the lncRNA.
FIRST_TYPES = {'RNA': NoncodingRNAProduct, 'DNA': Gene}
# The second partner. ``TF`` is a transcription factor, a protein.
SECOND_TYPES = {'Protein': Protein, 'TF': Protein, 'RNA': RNAProduct, 'DNA': Gene}

# Every value of the ``interaction_class`` column with the kind of relation.
CLASSES = {
    'binding': 'binding',
    'Binding': 'binding',
    # Both are stated; only the binding is carried.
    'binding/Regulatory': 'binding',
    'Binding/Regulatory': 'binding',
    'Co-expression': 'coexpression',
    'Regulatory': 'regulation',
}
# At this level the source (and the paper) say that the transcription factor
# regulates the lncRNA gene.
TF_LEVEL = 'DNA-TF'

# Accession forms of the ``accession`` column. A RefSeq transcript, a GenBank
# nucleotide record or a GenBank GI number identifies an RNA; a RefSeq genomic
# region (NG_, NC_) does not.
_transcript = re.compile(r'((?:NR|NM|XR|XM)_\d+(?:\.\d+)?)')
_genomic = re.compile(r'(?:NG|NC)_\d+(?:\.\d+)?')
_genbank = re.compile(r'[A-Z]{1,2}\d{5,6}(?:\.\d+)?')
_gi = re.compile(r'\d+')


def accession(row: dict[str, Any], namespace: Namespace) -> str | None:
    """The sequence accession of the lncRNA, if it is one of ``namespace``.

    The accession cell sometimes continues with a query string (``NG_001085.2?from=1``).
    """
    value = (row['accession'] or '').split('?')[0]
    if not value or _genomic.fullmatch(value):
        return None
    for pattern, ns in (
        (_transcript, Namespace.REFSEQ),
        (_genbank, Namespace.GENBANK),
        (_gi, Namespace.GENBANK_GI),
    ):
        if pattern.fullmatch(value):
            return value if ns == namespace else None
    raise ValueError(f'Unknown sequence accession: {row["accession"]!r}')


def taxon(row: dict[str, Any]) -> str | None:
    if row['species'] is None:
        return None
    if row['species'] not in SPECIES:
        raise ValueError(f'Unknown species: {row["species"]!r}')
    return f'NCBITaxon:{SPECIES[row["species"]]}'


def first_type(row: dict[str, Any]):
    if row['first_type'] not in FIRST_TYPES:
        raise ValueError(f'Unknown first partner type: {row["first_type"]!r}')
    return FIRST_TYPES[row['first_type']]


def second_type(row: dict[str, Any]):
    if row['second_type'] not in SECOND_TYPES:
        raise ValueError(f'Unknown second partner type: {row["second_type"]!r}')
    return SECOND_TYPES[row['second_type']]


f = FieldConfig()

# The lncRNA (or, at the DNA level, its gene): the accession of the transcript
# identifies only an RNA.
lncrna = EntityBuilder(
    entity_type=first_type,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.REFSEQ, value=lambda row: accession(row, Namespace.REFSEQ) if row['first_type'] == 'RNA' else None),
        CV(term=Namespace.GENBANK, value=lambda row: accession(row, Namespace.GENBANK) if row['first_type'] == 'RNA' else None),
        CV(term=Namespace.GENBANK_GI, value=lambda row: accession(row, Namespace.GENBANK_GI) if row['first_type'] == 'RNA' else None),
        CV(term=Namespace.NAME, value=f('lncrna')),
        CV(term=Namespace.SYNONYM, value=f('aliases')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=taxon)),
)

# The interacting molecule is known by its name only.
partner = EntityBuilder(
    entity_type=second_type,
    identifiers=IdentifiersBuilder(CV(term=Namespace.NAME, value=f('partner'))),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=taxon)),
)

annotations = AnnotationsBuilder(
    CV(term=slots.publications, value=lambda row: [f'PMID:{p}' for p in row['pmids']]),
    CV(term=slots.description, value=f('description')),
)


def _relation(predicate, subject, obj) -> RelationBuilder:
    return RelationBuilder(
        subject=subject,
        predicate=predicate,
        object=obj,
        annotations=annotations,
    )


def _present(builder, key: str):
    return lambda row: builder.build(row) if row[key] else None


lncrna_first = _present(lncrna, 'lncrna')
partner_second = _present(partner, 'partner')

# Binding is physical and has no direction. Co-expression is symmetric. A
# regulatory interaction has no stated direction, except at the DNA-TF level.
binding = _relation(slots.interacts_with, lncrna_first, partner_second)
coexpression = _relation(slots.coexpressed_with, lncrna_first, partner_second)
association = _relation(slots.associated_with, lncrna_first, partner_second)
regulation_of_gene = _relation(slots.affects, partner_second, lncrna_first)


def interactions_mapper(row: dict[str, Any]):
    """The relation of a row; None when the source gives no level (``N/A``), which leaves the partners' types unknown."""
    taxon(row)
    if row['level'] is None:
        return None
    # Errors raised inside the builders' callbacks are not propagated, so the
    # crosswalks are checked here.
    first_type(row)
    second_type(row)
    if row['first_type'] == 'RNA':
        for namespace in (Namespace.REFSEQ, Namespace.GENBANK, Namespace.GENBANK_GI):
            accession(row, namespace)
    if row['interaction_class'] not in CLASSES:
        raise ValueError(f'Unknown interaction class: {row["interaction_class"]!r}')
    kind = CLASSES[row['interaction_class']]
    if kind == 'binding':
        builder = binding
    elif kind == 'coexpression':
        builder = coexpression
    else:
        builder = regulation_of_gene if row['level'] == TF_LEVEL else association
    return builder(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_mapper,
        raw_parser=iter_interactions,
    ),
)
