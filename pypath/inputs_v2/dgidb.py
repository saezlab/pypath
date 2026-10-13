"""
DGIdb gene categories and drug-gene interactions.

DGIdb aggregates the drug-gene interaction claims and the druggable gene
category claims of many source databases and maps their gene and drug names
to concepts (HGNC, ChEMBL, RxNorm...). This module reads three of its
downloads:

- ``categories.tsv``: a gene category claim names a gene as its source does and
  a category (``KINASE``, ``DRUGGABLE GENOME``) of DGIdb's vocabulary. The gene
  is joined from ``genes.tsv``, which maps each claim to an HGNC gene. The
  statement is the gene ``associated_with`` the category, which is its own
  entity (the ``terms`` dataset); the database of the claim is the supporting
  data source.
- ``interactions.tsv``: a drug-gene interaction claim of a source database with
  its interaction types. The type decides the predicate and the qualifiers
  (``TYPES``).
"""

from __future__ import annotations

from functools import partial
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import (
    CausalMechanismQualifierEnum as Mechanism,
    DirectionQualifierEnum as Direction,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.dgidb import (
    iter_categories,
    iter_category_terms,
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

# DGIdb's own site stays on the data of June 2024 under /data/latest/ (a
# different file layout), the releases are on GitHub. The module is pinned to
# the release its audit counts are those of (33,889 category claims, 122,980
# interaction type rows); ``releases/latest/download/`` is a moving target, and
# served the same three files (byte for byte) when this was pinned.
RELEASE = '2026-09'
RELEASE_URL = f'https://github.com/dgidb/dgidb-data/releases/download/{RELEASE}/'

config = ResourceConfig(
    id=ResourceCv.DGIDB,
    name='DGIdb',
    full='Drug-Gene Interaction database',
    url='https://dgidb.org/',
    # The site (https://dgidb.org/about/overview/data-accessibility): "The source
    # code for DGIdb is open-source and made available under the MIT license. [...]
    # The data used in DGIdb is all open access"; it names no license for the
    # data, and the downloads page notes that some sources restrict
    # redistribution. MIT is the license of the legacy registry (resources.json).
    license=LicenseCV.MIT,
    update_category=UpdateCategoryCV.REGULAR,
    # DGIdb 5.0, Nucleic Acids Research 2024.
    pubmed='37953380',
    topics=(Topic.PHARMACOLOGY,),
    content=(Content.MOLECULAR_ANNOTATIONS, Content.COMPOUND_TARGET_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Druggable gene categories and drug-gene interactions with their types, '
        f'aggregated from many source databases (data release {RELEASE}).'
    ),
    description=(
        'The Drug-Gene Interaction database (DGIdb) aggregates, normalizes and '
        'curates drug-gene interaction claims and druggable gene category claims '
        'from publications and expert resources, with a focus on cancer. Genes and '
        'drugs are mapped to HGNC, ChEMBL, RxNorm and other concepts. Each claim '
        'keeps the source database it comes from; DGIdb adds an interaction score '
        f'that ranks the claims. This module reads the data release {RELEASE} of '
        'the dgidb-data repository on GitHub (DGIdb v5.0.13).'
    ),
)


def _download(name: str) -> Download:
    return Download(
        url=f'{RELEASE_URL}{name}.tsv',
        filename=f'dgidb_{name}.tsv',
        subfolder=f'dgidb/release-{RELEASE}',
        large=True,
    )


download_categories = _download('categories')
download_genes = _download('genes')
download_interactions = _download('interactions')

# Namespace of the gene categories, a vocabulary of DGIdb.
CATEGORY_NAMESPACE = 'dgidb_category'

# Every source database of ``categories.tsv`` and ``interactions.tsv`` as the
# files write it (release ``RELEASE``), with the label carried as the supporting
# data source. A source that is not here raises.
CATEGORY_SOURCES = {
    name: name
    for name in (
        'BaderLab', 'CIViC', 'COSMIC', 'CarisMolecularIntelligence', 'FoundationOneGenes',
        'GO', 'GuideToPharmacology', 'HingoraniCasas', 'HopkinsGroom', 'HumanProteinAtlas',
        'IDG', 'MskImpact', 'Oncomine', 'Pharos', 'RussLampel', 'Tempus', 'dGene',
    )
}
INTERACTION_SOURCES = {
    name: name
    for name in (
        'COSMIC', 'CGI', 'CIViC', 'CKB-CORE', 'CancerCommons', 'ChEMBL', 'ClearityFoundationBiomarkers',
        'ClearityFoundationClinicalTrial', 'ClinPGx', 'DTC', 'DoCM', 'Drug Repurposing Hub', 'FDA',
        'GuideToPharmacology', 'MOAlmanac', 'MyCancerGenome', 'MyCancerGenomeClinicalTrial', 'NCI',
        'OncoKB', 'PRISM', 'TALC', 'TEND', 'TTD', 'TdgClinicalTrial',
    )
}

# Interaction type -> (predicate, direction of the effect on the gene, mechanism).
# The types and their meaning are those of DGIdb
# (https://dgidb.org/about/overview/types-and-directionality; the definitions
# are also in its GraphQL API, ``interactionClaimTypes``). A type is an effect
# on the target (``affects``) only where its definition says the drug changes
# the target's activity or expression; one that says the drug binds the target,
# or that is the label of an immune therapy, is an interaction. DGIdb's own
# "directionality" column of that table is a loose grouping and is not used
# where the definition does not say the same. The mechanisms are those of
# GtoPdb's actions (``guidetopharma.py``) where DGIdb's type is the same
# notion, and the Biolink value of the same name otherwise. These are the 16
# types of release ``RELEASE``; a type that is not here raises.
TYPES = {
    # "decreases its expression or activity"
    'inhibitor': (slots.affects, Direction.decreased, Mechanism.inhibition),
    # "binds to a target receptor and activates the receptor"
    'agonist': (slots.affects, Direction.increased, Mechanism.agonism),
    # "activates a biological response from a target"
    'activator': (slots.affects, Direction.increased, Mechanism.activation),
    # "Antagonist interactions are sometimes referred to as blocker interactions"
    'blocker': (slots.affects, Direction.decreased, Mechanism.antagonism),
    # "induces a pharmacological response opposite to that of the agonist"
    'inverse agonist': (slots.affects, Direction.decreased, Mechanism.inverse_agonism),
    # "regulates or changes the activity of its target"
    'modulator': (slots.affects, None, Mechanism.modulation),
    # "increases activity of the target enzyme"
    'positive modulator': (slots.affects, Direction.increased, Mechanism.positive_modulation),
    # "negatively regulates the amount or activity of its target"
    'negative modulator': (slots.affects, Direction.decreased, Mechanism.negative_modulation),
    # "enhances the sensitivity of the target to the target's ligands"
    'potentiator': (slots.affects, None, Mechanism.potentiation),
    # "promotes degeneration of the target protein through cleaving of the peptide bonds"
    'cleavage': (slots.affects, Direction.decreased, Mechanism.cleavage),
    # "binds to an mRNA target to inhibit translation"
    'antisense oligonucleotide': (
        slots.affects, Direction.decreased, Mechanism.antisense_oligonucleotide_inhibition,
    ),
    # "physically binding to their target"
    'binder': (slots.interacts_with, None, Mechanism.binding),
    # "specifically binds the target molecule"; the functions that follow
    # (neutralization, cytotoxicity) are not stated for the claim.
    'antibody': (slots.interacts_with, None, None),
    # "stimulate or restore an immune response to their target"
    'vaccine': (slots.interacts_with, None, None),
    # "used to upregulate or downregulate the immune system": no effect on the gene is stated.
    'immunotherapy': (slots.interacts_with, None, None),
    # "a label given by the reporting source to an interaction that doesn't belong to other types"
    'other/unknown': (slots.interacts_with, None, None),
}
# A claim without a type is an interaction of unknown kind.
UNTYPED = (slots.interacts_with, None, None)

# DGIdb concept prefixes -> namespace of the accession written after the colon
# (``hgnc:10452``, ``chembl:CHEMBL592316``, ``rxcui:12574``).
GENE_NAMESPACES = {
    'hgnc': Namespace.HGNC,
    'ncbigene': Namespace.ENTREZ,
    'ensembl': Namespace.ENSG,
}
DRUG_NAMESPACES = {
    'chembl': Namespace.CHEMBL,
    'drugbank': Namespace.DRUGBANK,
    'iuphar.ligand': Namespace.GUIDETOPHARMA,
    'rxcui': Namespace.RXCUI,
    'ncit': Namespace.NCIT,
    'wikidata': Namespace.WIKIDATA,
    'hemonc': Namespace.HEMONC,
    'drugsatfda.nda': Namespace.DRUGSATFDA_NDA,
    'drugsatfda.anda': Namespace.DRUGSATFDA_ANDA,
    'chemidplus': Namespace.CHEMIDPLUS,
}
HUMAN = 'NCBITaxon:9606'
# DGIdb gene names are HGNC symbols.
HGNC_PREFIX = 'hgnc'

# The sources are checked in the mappers (the builders return None for an
# unknown key); the crosswalk here gives the label carried.
f = FieldConfig(
    map={
        'category_source': CATEGORY_SOURCES,
        'interaction_source': INTERACTION_SOURCES,
    },
)


def concept(value: str | None, namespaces: dict[str, Any]) -> tuple[str, str] | None:
    """(prefix, accession) of a concept id as written; None for none, an unknown prefix raises."""
    if value is None:
        return None
    prefix, colon, accession = value.partition(':')
    if not colon or not accession or prefix not in namespaces:
        raise ValueError(f'Unknown DGIdb concept: {value!r}')
    return prefix, accession


def _accession(column: str, namespaces: dict[str, Any], namespace):
    """The accession of ``column`` when its concept is of ``namespace``."""

    def value(row: dict[str, Any]) -> str | None:
        found = concept(row[column], namespaces)
        return found[1] if found and namespaces[found[0]] == namespace else None

    return value


def _gene_taxon(row: dict[str, Any]) -> str | None:
    # An HGNC gene is human; the file says nothing of the species of the others.
    found = concept(row['gene_concept_id'], GENE_NAMESPACES)
    return HUMAN if found and found[0] == HGNC_PREFIX else None


gene = EntityBuilder(
    entity_type=model.Gene,
    identifiers=IdentifiersBuilder(
        *(
            CV(term=namespace, value=_accession('gene_concept_id', GENE_NAMESPACES, namespace))
            for namespace in GENE_NAMESPACES.values()
        ),
        CV(term=Namespace.GENESYMBOL, value=f('gene_name')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_gene_taxon)),
)

# The drug has the concepts the file gives and its name: ``drug_name`` is the
# name of the concept (``CHEMBL:CHEMBL592316`` where the concept has no name,
# which is no name), ``drug_claim_name`` the name the source gave.
_PLACEHOLDER_NAME = 'CHEMBL:'


def _drug_name(row: dict[str, Any]) -> str | None:
    name = row['drug_name']
    return None if name is None or name.startswith(_PLACEHOLDER_NAME) else name


drug = EntityBuilder(
    entity_type=model.ChemicalEntity,
    identifiers=IdentifiersBuilder(
        *(
            CV(term=namespace, value=_accession('drug_concept_id', DRUG_NAMESPACES, namespace))
            for namespace in dict.fromkeys(DRUG_NAMESPACES.values())
        ),
        CV(term=Namespace.NAME, value=_drug_name),
        CV(term=Namespace.NAME, value=lambda row: None if _drug_name(row) else row['drug_claim_name']),
        CV(term=Namespace.SYNONYM, value=lambda row: row['drug_claim_name'] if _drug_name(row) else None),
    ),
)

category = EntityBuilder(
    entity_type=model.OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=CATEGORY_NAMESPACE, value=f('category')),
        CV(term=Namespace.NAME, value=f('category')),
    ),
)

# The statement's object is the same entity as the one of the ``terms`` dataset.
category_object = EntityBuilder(
    entity_type=model.OntologyClass,
    identifiers=IdentifiersBuilder(CV(term=CATEGORY_NAMESPACE, value=f('category'))),
)

category_schema = RelationBuilder(
    subject=gene,
    predicate=slots.associated_with,
    object=category_object,
    annotations=AnnotationsBuilder(
        CV(term=slots.supporting_data_source, value=f('source', map='category_source')),
    ),
)

# DGIdb's own scores. The interaction score ranks the claims of DGIdb (the
# evidence of the claim, and how specific the drug and the gene are); it is
# not a probability, and the scale is open.
_SCORES = {
    'dgidb:interaction_score': 'interaction_score',
    'dgidb:drug_specificity_score': 'drug_specificity_score',
    'dgidb:gene_specificity_score': 'gene_specificity_score',
    'dgidb:evidence_score': 'evidence_score',
}
# Properties of the drug, repeated on each of its claims: the file writes
# ``true`` or ``false`` (or leaves the cell empty), carried as written through
# this crosswalk; any other value raises (``map_interaction``).
YES_NO = {'true': 'true', 'false': 'false'}
_FLAGS = {
    'dgidb:drug_is_approved': 'drug_is_approved',
    'dgidb:drug_is_immunotherapy': 'drug_is_immunotherapy',
    'dgidb:drug_is_antineoplastic': 'drug_is_antineoplastic',
}


def _type(row: dict[str, Any]):
    return TYPES[row['interaction_type']] if row['interaction_type'] else UNTYPED


def _score(column: str):
    return lambda row: measurement(row[column], source_field=column)


def _flag(column: str):
    return lambda row: None if row[column] is None else YES_NO[row[column]]


interaction_schema = RelationBuilder(
    subject=drug,
    predicate=lambda row: _type(row)[0],
    object=gene,
    annotations=AnnotationsBuilder(
        CV(term=slots.object_direction_qualifier, value=lambda row: _type(row)[1]),
        CV(term=slots.causal_mechanism_qualifier, value=lambda row: _type(row)[2]),
        CV(term=slots.supporting_data_source, value=f('source', map='interaction_source')),
        *(CV(term=term, value=_score(column)) for term, column in _SCORES.items()),
        *(CV(term=term, value=_flag(column)) for term, column in _FLAGS.items()),
    ),
)

def _check_source(row: dict[str, Any], sources: dict[str, str]) -> None:
    if row['source'] not in sources:
        raise ValueError(f'Unknown DGIdb source database: {row["source"]!r}')


def map_category(row: dict[str, Any]):
    """The gene's statement for a category claim; None for a claim that maps to no gene."""
    if row['category'] is None:
        raise ValueError(f'DGIdb category claim without a category: {row!r}')
    _check_source(row, CATEGORY_SOURCES)
    if concept(row['gene_concept_id'], GENE_NAMESPACES) is None:
        return None
    return category_schema(row)


def map_interaction(row: dict[str, Any]):
    """The statement of an interaction claim; None for a claim that maps to no gene."""
    _check_source(row, INTERACTION_SOURCES)
    if row['interaction_type'] and row['interaction_type'] not in TYPES:
        raise ValueError(f'Unknown DGIdb interaction type: {row["interaction_type"]!r}')
    for flag in _FLAGS.values():
        if row[flag] is not None and row[flag] not in YES_NO:
            raise ValueError(f'Unknown DGIdb value of {flag}: {row[flag]!r}')
    concept(row['drug_concept_id'], DRUG_NAMESPACES)
    if concept(row['gene_concept_id'], GENE_NAMESPACES) is None:
        return None
    return interaction_schema(row)


resource = Resource(
    config,
    terms=Dataset(
        download=download_categories,
        mapper=category,
        raw_parser=iter_category_terms,
    ),
    categories=Dataset(
        download=download_categories,
        mapper=map_category,
        raw_parser=partial(iter_categories, genes_download=download_genes),
    ),
    interactions=Dataset(
        download=download_interactions,
        mapper=map_interaction,
        raw_parser=iter_interactions,
    ),
)
