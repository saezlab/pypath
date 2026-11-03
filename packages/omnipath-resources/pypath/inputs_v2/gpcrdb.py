"""
The GPCR classification of GPCRdb: classes, ligand types and receptor families.

GPCRdb classifies the G protein-coupled receptors in a tree of four levels:
the receptor class (Class A (Rhodopsin), Class B1 (Secretin), ...), the ligand
type (Aminergic receptors, Peptide receptors, ...), the receptor family
(5-Hydroxytryptamine receptors, ...) and the receptor. This module emits the
three levels above the receptors as protein families, each with the members
the outline puts below it: a class has its ligand types, a ligand type its
families and a family its receptors, as the UniProt accessions the outline
gives for the human, rat and mouse receptor. The family tree of GPCRdb's
service gives the identifiers of the families (its slugs), which the outline
lacks.

The legacy input mapped each human or mouse receptor to its class, ligand type
("family") and receptor family ("subfamily") as three strings; the rat receptors
were not read.
"""

from __future__ import annotations

from functools import partial
from typing import Any

from biolink_model.datamodel.model import Protein, ProteinFamily, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.gpcrdb import iter_families
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)

config = ResourceConfig(
    id=ResourceCv.GPCRDB,
    name='GPCRdb',
    url='https://gpcrdb.org/',
    # The site and its data repository (protwis/gpcrdb_data) state no license
    # (checked 2026-10-12). License as recorded in the legacy registry
    # (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33270898',
    topics=(Topic.SIGNALING,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'GPCR classes, ligand types and receptor families, each with the '
        'human, rat and mouse receptors below it.'
    ),
    description=(
        'GPCRdb is the reference database of G protein-coupled receptors. It '
        'classifies the receptors by class, ligand type and receptor family '
        '(following the IUPHAR classification) and holds their sequences, '
        'structures, ligands, mutations and signalling data.'
    ),
)

# ================================== DOWNLOAD ==================================

# The outline of the receptors of GPCRdb's data repository (the legacy URL).
download_outline = Download(
    url=(
        'https://raw.githubusercontent.com/protwis/gpcrdb_data/master/'
        'protein_data/proteins_and_families.txt'
    ),
    filename='gpcrdb_proteins_and_families.txt',
    subfolder='gpcrdb',
    large=True,
)
# GPCRdb's family tree, with the slug of every node.
download_tree = Download(
    url='https://gpcrdb.org/services/proteinfamily/',
    filename='gpcrdb_proteinfamily.json',
    subfolder='gpcrdb',
    large=True,
)

# =================================== SCHEMA ===================================

# GPCRdb's slugs number the nodes of the tree, and a slug starts with the slug
# of its parent (``001_001_001``).
FAMILY = 'gpcrdb_family'


def _family(node: dict[str, Any]) -> Entity:
    return Entity(
        type=ProteinFamily,
        identifiers=[
            Identifier(FAMILY, node['slug']),
            Identifier(Namespace.NAME, node['name']),
        ],
    )


def _receptor(member: dict[str, Any]) -> Entity:
    return Entity(
        type=Protein,
        identifiers=[Identifier(Namespace.UNIPROT, member['accession'])],
        annotations=[Annotation(slots.in_taxon, f'NCBITaxon:{member["taxon"]}')],
    )


def map_family(row: dict[str, Any]) -> Entity:
    """A class, a ligand type or a family, with the members the outline lists."""
    if row['level'] == 'family':
        members = [_receptor(member) for member in row['members']]
    elif row['level'] in ('class', 'ligand_type'):
        members = [_family(child) for child in row['children']]
    else:
        raise ValueError(f'Unknown GPCRdb level: {row["level"]!r}')
    entity = _family(row)
    return entity._replace(
        membership=[
            Membership(member=member, predicate=slots.has_member)
            for member in members
        ]
    )


# ================================= RESOURCE ===================================


resource = Resource(
    config,
    families=Dataset(
        download=download_outline,
        mapper=map_family,
        raw_parser=partial(iter_families, tree=download_tree),
    ),
)

# ================================== REFERENCE =================================

# proteins_and_families.txt, a receptor line (the human gene, then rat, mouse):
#   "gpcr",IUPHAR family number,family name,IUPHAR receptor number,receptor name,
#   "","","",aliases,HGNC id,symbol,gene name,locus,RefSeq mRNA,RefSeq protein,
#   UniProt accession,Entrez Gene id,<id>, rat: symbol,name,locus,mRNA,protein,
#   accession,Entrez,<MGI id>, mouse: symbol,name,locus,mRNA,protein,accession,Entrez
#   e.g. "gpcr","1","5-Hydroxytryptamine receptors","1","5-HT<sub>1A</sub> receptor",
#        "","","","ADRBRL1|5-HT1A|...","5286","HTR1A",...,"P08908","3350","2845",
#        "Htr1a",...,"P19327","24473","MGI:96273","Htr1a",...,"Q64264","15550"
