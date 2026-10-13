"""
Parse the COMPLEAT human protein complexes and emit Entity records.

COMPLEAT (Harvard Medical School, Perrimon lab and the Drosophila RNAi
Screening Center) is a protein complex enrichment analysis tool. Its backend
is a resource of complexes of human, fly and yeast genes, assembled in two
ways, which the file keeps apart in its ``Literature`` / ``Predicted`` column:

- ``Literature``: complexes collected from CORUM, PINdb, CYC2008, Gene
  Ontology, SignaLink, flyReactome and a Drosophila AP-MS study, and mapped to
  human genes by orthology (DIOPT) where they come from another species. The
  column ``database`` names one source database for each complex; it does not
  list all sources of a complex that was merged from several.
- ``Predicted``: complexes found by clustering high-quality protein-protein
  interaction networks (CFinder, NetworkBlast).

The two are different datasets here, ``literature_complexes`` and
``predicted_complexes``, so that predicted complexes are a choice of dataset
and not a hidden filter. The members are Entrez Gene IDs of the human
resource, as the file gives them.

The legacy input read the same file (12 columns, no header), turned the
members into UniProt accessions and merged complexes with equal members. This
module keeps each source complex, with its own accession, and the Entrez
Gene IDs.

The file is the older export of the tool, kept on rescued.omnipathdb.org.
The download of the tool today (flyrnai.org/compleat, ``ComplexDownload``)
returns a different file of 6 unlabelled columns, with 9,881 instead of
9,703 complexes (the 178 more are complexes whose only source is a KEGG
module, which are metabolic pathway modules) and the sources of a complex as
unexplained codes.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel.model import MacromolecularComplex, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.compleat import iter_compleat_complexes
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)

config = ResourceConfig(
    id=ResourceCv.COMPLEAT,
    name='COMPLEAT',
    full='COMPLEAT protein COMPLex Enrichment Analysis Tool',
    url='https://www.flyrnai.org/compleat/',
    # The site and the paper (not open access) state no license.
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_SA_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='23443684',
    content=(Content.PROTEIN_COMPLEXES,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Human protein complexes collected from the literature databases and '
        'predicted from interaction networks, with their member genes.'
    ),
    description=(
        'COMPLEAT is an online tool that analyses high-throughput data sets by '
        'protein complex enrichment. Its resource of complexes of human, fly and '
        'yeast genes was assembled from complexes annotated in CORUM, PINdb, '
        'CYC2008, Gene Ontology, SignaLink, flyReactome and a Drosophila AP-MS '
        'study, mapped across species by orthology, and from complexes predicted '
        'from high-quality protein-protein interaction networks with CFinder and '
        'NetworkBlast.'
    ),
)

# ``database`` -> the database it names (``supporting_data_source``). The
# values that are not databases give none.
DATABASES = {
    'CORUM': 'CORUM',
    'CYC2008': 'CYC2008',
    'PINdb': 'PINdb',
    'Gene Ontology': 'Gene Ontology',
    'SignaLink': 'SignaLink',
    'flyReactome': 'flyReactome',
    'Additional publication': None,
    'Predicted from PPI network': None,
}

# COMPLEAT accessions (``HC1``) have no Namespace member yet.
COMPLEAT_NAMESPACE = 'compleat'
# The file is the human resource of the tool.
TAXON = 'NCBITaxon:9606'


def map_complex(row: dict[str, Any]) -> Entity | None:
    """A complex with its Entrez Gene members; None without members."""
    members = row.get('members') or []
    if not members:
        return None
    identifiers = [Identifier(type=COMPLEAT_NAMESPACE, value=row['compleat_id'])]
    if row.get('name'):
        identifiers.append(Identifier(type=Namespace.NAME, value=row['name']))
    if row['database'] not in DATABASES:
        raise ValueError(f'Unknown COMPLEAT database: {row["database"]!r}')
    database = DATABASES[row['database']]
    return Entity(
        type=MacromolecularComplex,
        identifiers=identifiers,
        annotations=[
            *(
                [Annotation(term=slots.supporting_data_source, value=database)]
                if database
                else []
            ),
            *(
                Annotation(term=slots.has_evidence_of_type, value=term)
                for term in row['detection_methods']
            ),
            *(
                Annotation(term=slots.publications, value=f'PMID:{pmid}')
                for pmid in row['publications']
            ),
        ]
        or None,
        membership=[
            Membership(
                member=Entity(
                    type=Protein,
                    identifiers=[Identifier(type=Namespace.ENTREZ, value=member)],
                    annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
                )
            )
            for member in members
        ],
    )


# A copy of the file the tool served as ComplexDownload?requestType=complexDownload&org=Human
# (http://www.flyrnai.org/compleat/ComplexDownload, last modified 2021-08-20);
# the tool serves another format now.
download_complexes = functools.partial(
    Download,
    url='https://rescued.omnipathdb.org/compleat.tsv',
    filename='compleat.tsv',
    subfolder='compleat',
    large=False,
)

resource = Resource(
    config,
    literature_complexes=Dataset(
        download=download_complexes(),
        mapper=map_complex,
        raw_parser=functools.partial(iter_compleat_complexes, evidence='Literature'),
    ),
    predicted_complexes=Dataset(
        download=download_complexes(),
        mapper=map_complex,
        raw_parser=functools.partial(iter_compleat_complexes, evidence='Predicted'),
    ),
)
