"""
Gene sets of the Molecular Signatures Database (MSigDB), human collections.

MSigDB is a collection of annotated gene sets for gene set enrichment analysis.
The module imports the human collections of release 2026.1 from the JSON files
(HGNC gene symbols and the fields of each set) that the Broad Institute's data
server offers without a login, as two datasets:

* ``terms``: the collections and sub-collections (MSigDB's own codes, as
  ``C2:CP:REACTOME``) and the gene sets, each gene set below the collection of
  its file;
* ``genes``: a statement for each membership of a gene in a gene set, made on
  the gene (one record for each gene, with all its gene sets).

Each gene set is a source term: its name is its identifier. The term also
carries

* MSigDB's stable accession of the set (``M5890``, the ``systematicName`` field)
  as a second identifier, in the namespace ``msigdb``;
* the page of the gene set at the MSigDB site as ``msigdb:url``;
* the publication the set comes from (the ``pmid`` field, filled for all hallmark,
  PID, miRDB, GTRD, 3CA, VAX and C8 sets, nearly all CGP, C6 and ImmuneSigDB
  sets and part of the legacy TFT sets) as ``publications``;
* what the set was made from (the ``exactSource`` field) as a cross-reference,
  where it is an accession of the Gene Ontology (``GO:0006915``, the GO sets of
  C5), the Human Phenotype Ontology (``HP:0001507``, the HPO sets of C5),
  Reactome (``R-HSA-109581``) or WikiPathways (``WP100``). The other contents
  of the field (a UniProt accession for GTRD, a KEGG MEDICUS network id, the
  table or figure of a publication, a data set) are not accessions of those
  vocabularies and are left out. The gene set is not related to the ontology
  term: which predicate fits is not decided.

The JSON files hold the same gene sets and the same gene lists as the GMT files
of the data server, which are therefore not read.

Collections left out for their license (see ``config``): ``C2:CP:BIOCARTA`` and
``C2:CP:KEGG_LEGACY``. The aggregated files of MSigDB (``c2.all``, ``c2.cp``,
``c3.all``, ``c3.mir``, ``c3.tft``, ``c4.all``, ``c5.all``, ``c5.go``,
``c7.all`` and the whole-database file) repeat the sets of the files read and
are not read. The mouse collections (``Mm`` releases) are not imported.
"""

from __future__ import annotations

import functools
import re

from biolink_model.datamodel.model import Gene, OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.msigdb import iter_genes, iter_terms
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import EntityRef, OntologyRelation
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    AssociationBuilder,
    AssociationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)

config = ResourceConfig(
    id=ResourceCv.MSIGDB,
    name='MSigDB',
    full='Molecular Signatures Database',
    url='https://www.gsea-msigdb.org/gsea/msigdb/',
    # https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp (checked 2026-10-12), for
    # MSigDB v6.0 to v7.5.1 and v2022.1 and above: "The contents of MSigDB are protected by
    # copyright (c) 2004-2025 Broad Institute, Inc., Massachusetts Institute of Technology,
    # and Regents of the University of California, subject to the terms and conditions of
    # the Creative Commons Attribution 4.0 International License. Some MSigDB gene sets are
    # distributed with further restrictions." The page lists them:
    # - KEGG MEDICUS ("KEGG_MEDICUS_" sets of C2:CP:KEGG_MEDICUS): Creative Commons
    #   Attribution-ShareAlike 4.0 International License. Imported.
    # - Other KEGG ("KEGG_" sets of C2:CP:KEGG_LEGACY): "under license to the Broad
    #   Institute, Inc., with qualified permission to include in this release", conditioned
    #   on KEGG's academic and commercial licensing policies. Left out.
    # - BioCarta ("BIOCARTA_" sets of C2:CP:BIOCARTA): "under license to the Broad Institute,
    #   Inc. with qualified permission to include in this release", subject to BioCarta's
    #   disclaimer. Left out.
    # - AAAS/STKE Cell Signaling Database ("ST_"): removed from MSigDB v7.0 on; not in release
    #   2026.1.
    # The imported collections are therefore CC BY 4.0 with one CC BY-SA 4.0 subset, which
    # no single LicenseCV member describes.
    license=LicenseCV.COMPOSITE,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='21546393',
    topics=(Topic.MOLECULAR_CLASSIFICATION,),
    content=(Content.MOLECULAR_ANNOTATIONS, Content.PATHWAYS),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary='Human gene sets (hallmark, curated, regulatory target, cell type and others) of MSigDB, with their accession, publication and ontology or pathway source.',
    description=(
        'The Molecular Signatures Database (MSigDB) is a collection of '
        'annotated gene sets for use with gene set enrichment analysis. The '
        'human collections group about 35,000 gene sets: hallmark sets, '
        'positional sets, curated sets (canonical pathways and chemical and '
        'genetic perturbations), regulatory target sets of transcription '
        'factors and microRNAs, computational sets, ontology sets, oncogenic '
        'signatures, immunologic signatures, cell type signatures and '
        'computational perturbation signatures. Each gene set has a stable '
        'MSigDB accession, a page at the MSigDB site and, where there is one, '
        'the publication it comes from. The gene sets of the Gene Ontology, '
        'the Human Phenotype Ontology, Reactome and WikiPathways name the term '
        'they were made from.'
    ),
)

# ================================== DOWNLOAD ==================================

# The data server of the Broad Institute lists the release directories and serves
# the files without a login. Next to the GMT files (symbols or Entrez ids) it has
# one JSON file for each collection with the same sets and gene lists and the
# fields of the sets (accession, publication, exact source). The legacy URLs
# (www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=...) need a
# registered login. ``Hs`` is the human release; the mouse one is ``2026.1.Mm``.
RELEASE = '2026.1.Hs'
BASE_URL = f'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{RELEASE}/'


def _download(stem: str) -> Download:
    return Download(
        url=f'{BASE_URL}{stem}.v{RELEASE}.json',
        filename=f'{stem}.v{RELEASE}.json',
        subfolder='msigdb',
        # Each file is one JSON object, read whole (at most 10 MB).
        large=False,
    )


# Collection code -> title, description, parent code and file stem, after
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp. A collection without
# a file groups sub-collections whose files are read. Order: parents first.
_COLLECTION_SPECS = (
    ('H', 'hallmark gene sets', None, 'h.all',
     'Gene sets that summarize specific well-defined biological states or processes and display coherent expression.'),
    ('C1', 'positional gene sets', None, 'c1.all',
     'Gene sets corresponding to human chromosome cytogenetic bands.'),
    ('C2', 'curated gene sets', None, None,
     'Gene sets curated from online pathway databases and the biomedical literature.'),
    ('C2:CGP', 'chemical and genetic perturbations', 'C2', 'c2.cgp',
     'Expression signatures of genetic and chemical perturbations.'),
    ('C2:CP', 'canonical pathways', 'C2', None,
     'Gene sets from pathway databases, canonical representations of biological processes.'),
    ('C2:CP:KEGG_MEDICUS', 'KEGG MEDICUS subset of CP', 'C2:CP', 'c2.cp.kegg_medicus',
     'Canonical pathways gene sets derived from the KEGG MEDICUS pathway database.'),
    ('C2:CP:PID', 'PID subset of CP', 'C2:CP', 'c2.cp.pid',
     'Canonical pathways gene sets derived from the PID pathway database.'),
    ('C2:CP:REACTOME', 'Reactome subset of CP', 'C2:CP', 'c2.cp.reactome',
     'Canonical pathways gene sets derived from the Reactome pathway database.'),
    ('C2:CP:WIKIPATHWAYS', 'WikiPathways subset of CP', 'C2:CP', 'c2.cp.wikipathways',
     'Canonical pathways gene sets derived from the WikiPathways pathway database.'),
    ('C3', 'regulatory target gene sets', None, None,
     'Gene sets of potential targets of regulation by transcription factors or microRNAs.'),
    ('C3:MIR', 'microRNA targets', 'C3', None,
     'Gene sets of predicted targets of microRNAs.'),
    ('C3:MIR:MIRDB', 'miRDB subset of MIR', 'C3:MIR', 'c3.mir.mirdb',
     'High-confidence gene-level predictions of human miRNA targets catalogued by miRDB.'),
    ('C3:MIR:MIR_LEGACY', 'MIR_LEGACY subset of MIR', 'C3:MIR', 'c3.mir.mir_legacy',
     'Older gene sets of genes sharing putative target sites of human mature miRNAs in their 3\'-UTRs.'),
    ('C3:TFT', 'transcription factor targets', 'C3', None,
     'Gene sets of predicted targets of transcription factors.'),
    ('C3:TFT:GTRD', 'GTRD subset of TFT', 'C3:TFT', 'c3.tft.gtrd',
     'Genes that share predicted transcription factor binding sites of the GTRD around the TSS for a transcription factor.'),
    ('C3:TFT:TFT_LEGACY', 'TFT_LEGACY subset of TFT', 'C3:TFT', 'c3.tft.tft_legacy',
     'Older gene sets of genes sharing upstream cis-regulatory motifs, potential transcription factor binding sites.'),
    ('C4', 'computational gene sets', None, None,
     'Gene sets defined by mining large collections of cancer-oriented expression data.'),
    ('C4:3CA', 'Curated Cancer Cell Atlas', 'C4', 'c4.3ca',
     'Gene sets mined from the metaprograms of the Curated Cancer Cell Atlas (3CA).'),
    ('C4:CGN', 'cancer gene neighborhoods', 'C4', 'c4.cgn',
     'Gene sets defined by expression neighborhoods centered on 380 cancer-associated genes.'),
    ('C4:CM', 'cancer modules', 'C4', 'c4.cm',
     'Gene sets (modules) significantly changed in a variety of cancer conditions, defined by Segal et al. 2004.'),
    ('C5', 'ontology gene sets', None, None,
     'Gene sets of the genes annotated by the same ontology term.'),
    ('C5:GO', 'Gene Ontology gene sets', 'C5', None,
     'Gene sets derived from the Gene Ontology.'),
    ('C5:GO:BP', 'BP subset of GO', 'C5:GO', 'c5.go.bp',
     'Gene sets derived from the Gene Ontology Biological Process ontology.'),
    ('C5:GO:CC', 'CC subset of GO', 'C5:GO', 'c5.go.cc',
     'Gene sets derived from the Gene Ontology Cellular Component ontology.'),
    ('C5:GO:MF', 'MF subset of GO', 'C5:GO', 'c5.go.mf',
     'Gene sets derived from the Gene Ontology Molecular Function ontology.'),
    ('C5:HPO', 'Human Phenotype Ontology gene sets', 'C5', 'c5.hpo',
     'Gene sets derived from the Human Phenotype Ontology.'),
    ('C6', 'oncogenic signature gene sets', None, 'c6.all',
     'Gene sets that represent signatures of cellular pathways which are often dis-regulated in cancer.'),
    ('C7', 'immunologic signature gene sets', None, None,
     'Gene sets that represent cell states and perturbations within the immune system.'),
    ('C7:IMMUNESIGDB', 'ImmuneSigDB subset of C7', 'C7', 'c7.immunesigdb',
     'Gene sets of chemical and genetic perturbations of the immune system, curated from published studies in human and mouse immunology.'),
    ('C7:VAX', 'vaccine response gene sets', 'C7', 'c7.vax',
     'Gene sets curated by the Human Immunology Project Consortium (HIPC) on human transcriptomic immune responses to vaccinations.'),
    ('C8', 'cell type signature gene sets', None, 'c8.all',
     'Gene sets of curated cluster markers for cell types identified in single-cell sequencing studies of human tissue.'),
    ('C9', 'computational perturbation signature gene sets', None, 'c9.all',
     'Gene sets of the coordinate transcriptional responses to expression of specific genes, from computational analysis of DepMap CRISPR screen data and CCLE expression profiles.'),
)

COLLECTIONS = {
    code: {
        'title': title,
        'description': description,
        'parent': parent,
        'download': _download(stem) if stem else None,
    }
    for code, title, parent, stem, description in _COLLECTION_SPECS
}

# =================================== SCHEMA ===================================

COLLECTION = 'msigdb_collection'
GENE_SET = 'msigdb_gene_set'
# MSigDB's own accession system (``M5890``, the systematic name of a gene set).
ACCESSION = 'msigdb'
# The page of a gene set at the MSigDB site (the ``msigdbURL`` field).
URL = 'msigdb:url'
# The human files are of human gene symbols.
TAXON = 'NCBITaxon:9606'

f = FieldConfig()


def _gene_set_relations(row) -> list[OntologyRelation]:
    return [
        OntologyRelation(
            predicate=slots.subclass_of,
            object=EntityRef(OntologyClass, COLLECTION, row['collection']),
        )
    ]


def _collection_relations(row) -> list[OntologyRelation]:
    if not row.get('parent'):
        return []
    return [
        OntologyRelation(
            predicate=slots.subclass_of,
            object=EntityRef(OntologyClass, COLLECTION, row['parent']),
        )
    ]


collection_schema = EntityBuilder(
    entity_type=OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=COLLECTION, value=f('id')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.description, value=f('description')),
    ),
    ontology_relations=_collection_relations,
)

# What a gene set was made from (``exactSource``), where it is an accession of
# one of these vocabularies: the exact form of the accession, written as the
# module of the vocabulary writes its identifiers. The field also holds a UniProt
# accession (C3:TFT:GTRD), a KEGG MEDICUS network id (N00001), the table or figure
# of a publication, a data set and descriptions of a comparison; none of them is
# an accession of these vocabularies, and the term has no cross-reference for them.
_SOURCE_ACCESSIONS = (
    re.compile(r'GO:\d{7}'),  # go.py: Namespace.GO
    re.compile(r'HP:\d{7}'),  # hpo.py: Namespace.HPO
    re.compile(r'R-HSA-[1-9]\d*'),  # reactome.py: Namespace.REACTOME, the stable id
    re.compile(r'WP[1-9]\d*'),  # wikipathways.py: Namespace.WIKIPATHWAYS, the pathway id
)


def source_xref(value: str | None) -> str | None:
    """The cross-reference for an ``exactSource`` that is an accession, else None."""
    if value and any(pattern.fullmatch(value) for pattern in _SOURCE_ACCESSIONS):
        return value
    return None


# The term's identity is its first identifier, the name of the gene set, which is
# what the memberships of the ``genes`` dataset refer to. The systematic name is
# a second identifier of the same set.
gene_set_schema = EntityBuilder(
    entity_type=OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=GENE_SET, value=f('id')),
        CV(term=Namespace.NAME, value=f('id')),
        CV(term=ACCESSION, value=f('systematic_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=URL, value=f('url')),
        CV(term=slots.publications, value=lambda row: pmid_curie(row.get('pmid'))),
        CV(term=slots.xref, value=lambda row: source_xref(row.get('exact_source'))),
    ),
    ontology_relations=_gene_set_relations,
)


def map_term(row: dict):
    """The entity of a collection or of a gene set."""
    kind = row.get('kind')
    if kind == 'collection':
        return collection_schema(row)
    if kind == 'gene_set':
        if row.get('collection') not in COLLECTIONS:
            raise ValueError(f'Unknown MSigDB collection: {row.get("collection")!r}')
        return gene_set_schema(row)
    raise ValueError(f'Unknown MSigDB term kind: {kind!r}')


# Membership has nothing attached in the files, so it is an association of the
# gene. A gene set is a list of genes chosen by one criterion (membership in a
# pathway, a signature, a target set), so the gene is associated with the set
# and nothing causal, expressed or located is stated.
gene_schema = EntityBuilder(
    entity_type=Gene,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.GENESYMBOL, value=f('gene')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=TAXON)),
    associations=AssociationsBuilder(
        AssociationBuilder(
            predicate=slots.associated_with,
            object_entity_type=OntologyClass,
            object_identifier_type=GENE_SET,
            object_identifier=f('gene_sets'),
        ),
    ),
)

# ================================= RESOURCE ===================================

resource = Resource(
    config,
    terms=Dataset(
        download=None,
        mapper=map_term,
        raw_parser=functools.partial(iter_terms, collections=COLLECTIONS),
    ),
    genes=Dataset(
        download=None,
        mapper=gene_schema,
        raw_parser=functools.partial(iter_genes, collections=COLLECTIONS),
    ),
)

# ================================= REFERENCE ==================================

# "HALLMARK_TNFA_SIGNALING_VIA_NFKB": {"collection": "H", "systematicName": "M5890",
#   "pmid": "26771021", "exactSource": "", "externalDetailsURL": "",
#   "msigdbURL": "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_TNFA_SIGNALING_VIA_NFKB",
#   "geneSymbols": ["ABCA1", "ACKR3", ...], ...}
# (one object for a collection: gene set name -> fields)
