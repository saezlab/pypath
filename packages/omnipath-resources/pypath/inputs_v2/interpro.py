"""
Parse the InterPro entries and emit them as terms with their hierarchy, GO terms and proteins.

InterPro classifies protein sequences into families, domains, repeats and
sites by integrating the signatures of its member databases. The module emits

- ``terms``: each InterPro entry (``IPR000008``, C2 domain) as a term with its
  name, its type and its parents in the entry hierarchy;
- ``go_terms``: the GO terms that InterPro assigns to entries (interpro2go);
- ``protein_entries``: each reviewed human protein associated with the
  InterPro entries that match it, from UniProt's cross-references.

The entries and the hierarchy are those of the current InterPro release (the
files of ``current_release`` on the InterPro FTP site), the proteins are the
reviewed human proteins of the current UniProt release. The file with every
protein's matches (``protein2ipr.dat.gz``, 12 GB) is not read, and the
positions of the matches are not represented: UniProt's table has none, and
the InterPro API gives them one protein at a time.
"""

from __future__ import annotations

from functools import partial
from typing import Any

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.interpro import (
    iter_entries,
    iter_go_terms,
    iter_protein_entries,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Association,
    Entity,
    EntityRef,
    Identifier,
    OntologyRelation,
)

config = ResourceConfig(
    id=ResourceCv.INTERPRO,
    name='InterPro',
    url='https://www.ebi.ac.uk/interpro/',
    # "all of the InterPro, Pfam, PRINTS and SFLD downloadable data provided on
    # the InterPro website is freely available under CC0 1.0 Universal (CC0
    # 1.0) Public Domain Dedication"
    # (https://interpro-documentation.readthedocs.io/en/latest/license.html).
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='39565202',
    topics=(Topic.MOLECULAR_CLASSIFICATION,),
    content=(Content.ONTOLOGIES, Content.MOLECULAR_ANNOTATIONS),
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary=(
        'InterPro protein families, domains and sites with their hierarchy and '
        'GO terms, and the entries that match each reviewed human protein.'
    ),
    description=(
        'InterPro provides functional analysis of proteins by classifying '
        'them into families and predicting domains and important sites. It '
        'integrates the predictive models, signatures, of its member '
        'databases (Pfam, PROSITE, CATH-Gene3D, PANTHER, SMART and others) '
        'into entries that describe the same family, domain or site. The '
        'entries, their hierarchy and their GO terms are read from InterPro; '
        'the statements about proteins are not: they are the InterPro '
        'cross-references of UniProt for the reviewed human proteins (UniProt '
        'takes them from the matches of InterPro), without the positions of '
        'the matches.'
    ),
)

NAMESPACE = Namespace.INTERPRO
ONTOLOGY_ID = 'interpro'
# The type of an entry as ``entry.list`` writes it (Domain, Family, Repeat ...).
ENTRY_TYPE = 'interpro:type'
# The proteins are the reviewed proteins of the human proteome.
TAXON = 'NCBITaxon:9606'
RELEASE = 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/'


def _entry_ref(accession: str) -> EntityRef:
    return EntityRef(OntologyClass, NAMESPACE, accession)


def map_entry(row: dict[str, Any]) -> Entity:
    """An InterPro entry with its name, type and parents."""
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=NAMESPACE, value=row['accession']),
            Identifier(type=Namespace.NAME, value=row['name']),
        ],
        annotations=[Annotation(term=ENTRY_TYPE, value=row['type'])],
        ontology_relations=[
            OntologyRelation(
                predicate=slots.subclass_of,
                object=_entry_ref(parent),
                ontology_id=ONTOLOGY_ID,
            )
            for parent in row['parents']
        ],
    )


def map_go_terms(row: dict[str, Any]) -> Entity:
    """An InterPro entry associated with its GO terms."""
    return Entity(
        type=OntologyClass,
        identifiers=[Identifier(type=NAMESPACE, value=row['accession'])],
        associations=[
            Association(
                object=EntityRef(OntologyClass, Namespace.GO, term),
                predicate=slots.associated_with,
            )
            for term in row['go_terms']
        ],
    )


def map_protein(row: dict[str, Any]) -> Entity | None:
    """A reviewed human protein associated with the InterPro entries that match it."""
    if not row['entries']:
        return None
    return Entity(
        type=Protein,
        identifiers=[Identifier(type=Namespace.UNIPROT, value=row['accession'])],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        associations=[
            Association(object=_entry_ref(entry), predicate=slots.associated_with)
            for entry in row['entries']
        ],
    )


_entry_list = Download(
    url=RELEASE + 'entry.list',
    filename='interpro_entry.list',
    subfolder='interpro',
    large=True,
)
_tree = Download(
    url=RELEASE + 'ParentChildTreeFile.txt',
    filename='interpro_ParentChildTreeFile.txt',
    subfolder='interpro',
    large=True,
)
_interpro2go = Download(
    url=RELEASE + 'interpro2go',
    filename='interpro_interpro2go.txt',
    subfolder='interpro',
    large=True,
)
# UniProt's cross-references to InterPro, which UniProt takes from InterPro's
# matches: every reviewed (Swiss-Prot) human protein, one line each.
_proteins = Download(
    url=(
        'https://rest.uniprot.org/uniprotkb/stream'
        '?compressed=true'
        '&format=tsv'
        '&query=(organism_id:9606) AND (reviewed:true)'
        '&fields=accession,xref_interpro'
    ),
    filename='interpro_uniprot_human_reviewed.tsv.gz',
    subfolder='interpro',
    large=True,
    encoding='utf-8',
    default_mode='r',
    ext='gz',
)

resource = Resource(
    config,
    terms=Dataset(
        download=_entry_list,
        mapper=map_entry,
        raw_parser=partial(iter_entries, tree_download=_tree),
    ),
    go_terms=Dataset(
        download=_interpro2go, mapper=map_go_terms, raw_parser=iter_go_terms
    ),
    protein_entries=Dataset(
        download=_proteins, mapper=map_protein, raw_parser=iter_protein_entries
    ),
)
