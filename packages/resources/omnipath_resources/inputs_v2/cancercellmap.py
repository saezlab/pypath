"""
Parse the Cancer Cell Map binary network of Pathway Commons and emit Relation records.

The Cancer Cell Map was curated in BioPAX. Its own site is gone; the last
release of Pathway Commons 1 (October 2011) distributes it as a "tab delimited
network" of two tables: the edges, which name two Pathway Commons records, a
binary interaction rule and PubMed IDs, and the nodes, which give each record's
gene symbol, UniProt accessions, Entrez Gene IDs, type and organism.

    CPATH_RECORD_ID_A  INTERACTION_TYPE  CPATH_RECORD_ID_B  GENE_SYMBOL_A  GENE_SYMBOL_B  INTERACTION_DATA_SOURCE  INTERACTION_PUBMED_ID
    101318             INTERACTS_WITH    140665             MAP3K4         AXIN1          CELL_MAP                 PUBMED:15262978;

The interaction types are the rules by which Pathway Commons reduced BioPAX to
pairs. Only two of the four in this file state something about the pair itself;
see ``INTERACTION_TYPES``.
"""

from __future__ import annotations
from collections.abc import Callable
from functools import partial
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2._publications import PMID_PATTERN
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cancercellmap import iter_interactions

config = ResourceConfig(
    id=ResourceCv.CANCERCELLMAP,
    name='CancerCellMap',
    full='The Cancer Cell Map',
    url='https://download.baderlab.org/PathwayCommons/PC1/last_release-2011/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_2_5,
    update_category=UpdateCategoryCV.STATIC,
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS, Content.CAUSAL_REGULATION),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Protein interactions and state changes of ten cancer signaling pathways, with publications.',
    description='The Cancer Cell Map is a set of ten human cancer related signaling pathways (among them EGFR, TGF-beta receptor, Wnt, NOTCH, Hedgehog and androgen receptor) curated from the literature at the Institute of Bioinformatics, Bangalore, and reviewed at the Memorial Sloan-Kettering Cancer Center, which released it in 2006. Its own website no longer exists; the data is distributed by Pathway Commons, whose first version exported each source as a network of pairwise relationships inferred from the BioPAX pathways.',
)
_BASE_URL = (
    'https://download.baderlab.org/PathwayCommons/PC1/last_release-2011/'
    'tab_delim_network/by_source/'
)
download_edges = Download(
    url=_BASE_URL + 'cell-map-edge-attributes.txt.zip',
    filename='cell-map-edge-attributes.txt.zip',
    subfolder='cancercellmap',
    ext='zip',
    needed=['cell-map-edge-attributes.txt'],
)
download_nodes = Download(
    url=_BASE_URL + 'cell-map-node-attributes.txt.zip',
    filename='cell-map-node-attributes.txt.zip',
    subfolder='cancercellmap',
    ext='zip',
    needed=['cell-map-node-attributes.txt'],
)
# The binary interaction rules of Pathway Commons 1 that occur in the file and
# the predicate each maps to; None for a rule that is not a statement about the
# two proteins themselves.
#
# INTERACTS_WITH: "The entities participate in an interaction. Controllers are
#   not included. This interaction is undirected."
# STATE_CHANGE: "The first entity controls a reaction that changes the state of
#   the second entity, e.g. by phosphorylation or other posttranslational
#   modification, or by a change in subcellular location. This interaction is
#   directed." Neither the kind of change nor its sign is given.
# IN_SAME_COMPONENT: "Two entities belong to the same molecular complex. This
#   does not necessarily mean they interact directly." It is every pair of
#   members of a complex, without the complex and without publications.
# REACTS_WITH: "The entities participate in a conversion as substrates or
#   products." A protein paired with itself is its own conversion; two proteins
#   share a conversion, which need not bring them into contact.
INTERACTION_TYPES = {
    'INTERACTS_WITH': slots.interacts_with,
    'STATE_CHANGE': slots.affects,
    'IN_SAME_COMPONENT': None,
    'REACTS_WITH': None,
}
# Every node of the file is a protein.
NODE_TYPES = {'protein': Protein}


def _predicate(row: dict[str, Any]) -> Any:
    interaction_type = row.get('INTERACTION_TYPE')
    if interaction_type not in INTERACTION_TYPES:
        raise ValueError(
            f'Unknown Cancer Cell Map interaction type: {interaction_type!r}'
        )
    return INTERACTION_TYPES[interaction_type]


def _node_type(side: str) -> Callable[[dict[str, Any]], Any]:
    def node_type(row: dict[str, Any]) -> Any:
        value = row.get(f'NODE_TYPE_{side}')
        if value not in NODE_TYPES:
            raise ValueError(f'Unknown Cancer Cell Map node type: {value!r}')
        return NODE_TYPES[value]

    return node_type


f = FieldConfig(
    extract={
        # One node has a GenBank protein accession (``BAC32181``) in this column.
        'uniprot': r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$',
        'number': r'^(\d+)$',
        # Six rows have ``null:null;PUBMED:Literature;`` instead.
        'pmid': rf'^PUBMED:({PMID_PATTERN})$',
    },
    transform={
        'taxon_curie': lambda v: f'NCBITaxon:{v}',
        'pmid_curie': lambda v: f'PMID:{v}',
    },
)


def _participant(side: str) -> EntityBuilder:
    # A record lists the accessions of one protein in no particular order:
    # the current one, earlier ones merged into it and unreviewed entries.
    # A record with no accession, gene and symbol yields no participant.
    return EntityBuilder(
        entity_type=_node_type(side),
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.UNIPROT,
                value=f(f'UNIPROT_ACCESSION_{side}', delimiter=':', extract='uniprot'),
            ),
            CV(
                term=Namespace.ENTREZ,
                value=f(f'ENTREZ_GENE_ID_{side}', delimiter=':', extract='number'),
            ),
            CV(term=Namespace.GENESYMBOL, value=f(f'GENE_SYMBOL_{side}')),
        ),
        annotations=AnnotationsBuilder(
            # Human, mouse or rat; 34 records are of "Mammalia" (40674).
            CV(
                term=slots.in_taxon,
                value=f(f'NCBI_TAX_ID_{side}', extract='number', transform='taxon_curie'),
            ),
        ),
    )


# For STATE_CHANGE the first record is the controller and the second the
# protein whose state changes. The relation gets no direction or aspect
# qualifier: the rule covers any change of state.
interactions_schema = RelationBuilder(
    subject=_participant('A'),
    predicate=_predicate,
    object=_participant('B'),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.publications,
            value=f(
                'INTERACTION_PUBMED_ID',
                delimiter=';',
                extract='pmid',
                transform='pmid_curie',
            ),
        ),
    ),
)


def interactions_mapper(row: dict[str, Any]) -> Any:
    """Relation of an edge, None for a rule that is not mapped or a mirror row."""
    if _predicate(row) is None or row.get('MIRROR'):
        return None
    return interactions_schema.build(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=download_edges,
        mapper=interactions_mapper,
        raw_parser=partial(iter_interactions, nodes_download=download_nodes),
    ),
)
