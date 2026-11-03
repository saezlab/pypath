"""
Parse the census of human transcription factors and emit the list and its genes.

Vaquerizas et al. 2009 assembled the human transcription factors (TFs) from
the DNA-binding domains of InterPro, the literature and other TF databases,
and classified each locus by the evidence that it encodes a TF. Supplementary
Table 2 lists the loci with their class, Ensembl gene, HGNC symbol,
DNA-binding InterPro domains and families, and tissue specificity. The module
emits

- ``terms``: the census as one term, the list of human TFs, which describes
  the classes;
- ``transcription_factors``: each gene of classes a, b, c and other, a member
  of the list, with its class as the confidence level and its tissues;
- ``interpro_entries``: each such gene associated with the InterPro entries
  of its DNA-binding domains and DNA-binding families.

The genes of class x are the loci that the census judged not to be TFs: they
are statements that a gene is not a TF and are not emitted.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import Gene, OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.tfcensus import (
    iter_census,
    iter_census_interpro,
    iter_census_list,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    EntityRef,
    Identifier,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.TFCENSUS,
    name='TFcensus',
    full='Census of Human Transcription Factors',
    url='https://doi.org/10.1038/nrg2538',
    # The supplementary file states no license. License as recorded in the
    # legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='19274049',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Human transcription factor genes with the census class of evidence, '
        'DNA-binding InterPro entries and tissue specificity.'
    ),
    description=(
        'The census of Vaquerizas, Kummerfeld, Teichmann and Luscombe '
        'catalogues the sequence-specific DNA-binding transcription factors '
        'of the human genome, assembled from the DNA-binding domains of '
        'InterPro, the literature and TF databases and curated by hand, with '
        'their expression across tissues. The repertoire of the paper (1,391 '
        'TFs) is the classes a, b and other; the module also carries the '
        'class c of possible TFs.'
    ),
)

# The census is of the human genome.
TAXON = 'NCBITaxon:9606'
# The list has no Namespace member.
LIST_NAMESPACE = 'tfcensus_list'
LIST_ID = 'human_transcription_factors'
# The classes of the census and what they mean, as the Supplementary Material
# of the paper defines them (section 1, "Identification of TF repertoire"):
# "Genes classified as 'a' or 'b' are probable TFs for which there is
# experimental evidence of regulatory function in any mammalian organism ('a')
# or have an equivalent protein domain arrangement ('b'); those classified as
# 'c' are possible TFs that contain non-promiscuous InterPro DNA-binding
# domains [...] but for which we do not have further functional evidence."
# Class 'other' holds "27 curated probable TFs from other sources [...]
# containing undefined DNA-binding domains". The repertoire of the paper
# (1,391 TFs) is the classes a, b and other; c is not in it.
CLASS_MEANINGS = {
    'a': 'probable TFs with experimental evidence of regulatory function in a mammal',
    'b': 'probable TFs with a protein domain arrangement equivalent to that of a class a TF',
    'c': (
        'possible TFs with non-promiscuous InterPro DNA-binding domains '
        '(found only in TFs) and no further functional evidence'
    ),
    'other': (
        'probable TFs curated from other sources (Gene Ontology, TRANSFAC) '
        'with undefined DNA-binding domains'
    ),
}
NOT_TF_CLASS = 'x'
# Which InterPro column of the table an association comes from, with the
# source's name of the column.
INTERPRO_COLUMN = 'tfcensus:interpro_column'
# The tissue-specificity labels of the table, as written.
TISSUE = 'tfcensus:tissue'
_ENSEMBL_GENE = re.compile(r'^ENSG\d{11}$')
_INTERPRO_ENTRY = re.compile(r'^IPR\d{6}$')


def _is_member(row: dict[str, Any]) -> bool:
    """Whether the gene is in the list: any class but x. Unknown classes raise."""
    census_class = row['census_class']
    if census_class not in CLASS_MEANINGS and census_class != NOT_TF_CLASS:
        raise ValueError(f'Unknown TF census class in {row!r}')
    return census_class != NOT_TF_CLASS


def map_list(row: dict[str, Any]) -> Entity:
    """The census as the term for the list of human TFs, with its classes."""
    classes = '; '.join(
        f'class {name}: {meaning}' for name, meaning in CLASS_MEANINGS.items()
    )
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=LIST_NAMESPACE, value=LIST_ID),
            Identifier(
                type=Namespace.NAME,
                value='Human transcription factors (TF census)',
            ),
        ],
        annotations=[
            Annotation(
                term=slots.description,
                value=(
                    'Loci that encode human transcription factors in the '
                    'census of Vaquerizas et al. 2009, in classes by the '
                    f'evidence for each: {classes}. The 1,391 TFs of the '
                    'paper are the classes a, b and other; class c is a '
                    'possible TF. The census also has the class x for '
                    'unlikely TFs (predicted genes, promiscuous DNA-binding '
                    'domains, an established function outside '
                    'transcription); those loci are not members of the list.'
                ),
            )
        ],
    )


def _gene(row: dict[str, Any]) -> Entity:
    if not _ENSEMBL_GENE.match(row['ensembl']):
        raise ValueError(f'Malformed TF census Ensembl gene ID in {row!r}')
    identifiers = [Identifier(type=Namespace.ENSG, value=row['ensembl'])]
    if row['symbol']:
        identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=row['symbol']))
    return Entity(
        type=Gene,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def map_transcription_factor(row: dict[str, Any]) -> Relation | None:
    """A gene of the list, with its class as the confidence level."""
    if not _is_member(row):
        return None
    return Relation(
        subject=_gene(row),
        predicate=slots.associated_with,
        object=EntityRef(OntologyClass, LIST_NAMESPACE, LIST_ID),
        annotations=[
            Annotation(term=slots.has_confidence_level, value=row['census_class']),
            *(Annotation(term=TISSUE, value=tissue) for tissue in row['tissues']),
        ],
    )


def map_interpro_entry(row: dict[str, Any]) -> Relation | None:
    """A gene of the list and an InterPro entry of its DNA-binding domains."""
    if not _is_member(row):
        return None
    if not _INTERPRO_ENTRY.match(row['interpro']):
        raise ValueError(f'Malformed InterPro ID in {row!r}')
    return Relation(
        subject=_gene(row),
        predicate=slots.associated_with,
        object=EntityRef(OntologyClass, Namespace.INTERPRO, row['interpro']),
        annotations=[Annotation(term=INTERPRO_COLUMN, value=row['interpro_column'])],
    )


# The Springer supplement is Supplementary Table 2 of the paper (the legacy
# registry calls it S3).
download = Download(
    url=(
        'https://static-content.springer.com/esm/art%3A10.1038%2Fnrg2538/'
        'MediaObjects/41576_2009_BFnrg2538_MOESM6_ESM.txt'
    ),
    filename='41576_2009_BFnrg2538_MOESM6_ESM.txt',
    subfolder='tfcensus',
    large=True,
)

resource = Resource(
    config,
    terms=Dataset(download=download, mapper=map_list, raw_parser=iter_census_list),
    transcription_factors=Dataset(
        download=download, mapper=map_transcription_factor, raw_parser=iter_census
    ),
    interpro_entries=Dataset(
        download=download, mapper=map_interpro_entry, raw_parser=iter_census_interpro
    ),
)
