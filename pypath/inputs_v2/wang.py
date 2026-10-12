"""
Parse the Human Signaling Network of the Wang Lab and emit Relation records.

Version 6 of the network (January 2014) is one file with a table of nodes
(node number, gene name, Entrez Gene ID, molecular function, cellular location)
and a table of links between node numbers. The notes at its top define the
three link types:

    +,activation
    _,inhibition
    0,physical interaction

The network grew from the human signaling network of Cui et al. 2007
(Supplementary Table 9 of that publication has the same layout); version 6
supersedes it and is the only version read here.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import (
    DirectionQualifierEnum,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.wang import iter_network
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Relation
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.WANG,
    name='Wang',
    full='Human Signaling Network version 6',
    url='http://www.cancer-systemsbiology.org/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_NC_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='18091723',
    topics=(Topic.SIGNALING,),
    content=(Content.CAUSAL_REGULATION, Content.PHYSICAL_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Activating, inhibiting and physical links between human signaling '
        'proteins.'
    ),
    description=(
        'The Human Signaling Network is the manually curated signaling '
        'network of the laboratory of Edwin Wang. It started as the network '
        'of Cui et al. 2007, curated from the BioCarta pathway diagrams and '
        'joined with the CA1 neuron network and the Cancer Cell Map, and was '
        'extended in yearly versions from further pathway databases and '
        'publications. Version 6 of January 2014 has about 6,300 proteins and '
        '63,000 relations, each an activation, an inhibition or a physical '
        'interaction. It gives no reference for the individual relations.'
    ),
)

# The laboratory's pages that offered the file (www.bri.nrc.ca/wang) are gone.
# This is OmniPath's copy of ``HumanSignalingNetwork_v6.csv``; it is identical
# to the one the Internet Archive took from the original address in 2016.
download = Download(
    url='https://rescued.omnipathdb.org/WangLab_HumanSignalingNetwork_v6-format2.csv',
    filename='WangLab_HumanSignalingNetwork_v6-format2.csv',
    subfolder='wang',
)

# The file is titled the human signaling network and its nodes are human genes.
TAXON = 'NCBITaxon:9606'
# Link type -> predicate and direction of the effect on the target. The notes
# of the file define ``_`` as inhibition and its links write ``-``. Activation
# and inhibition are effects of the source node on the target node; the file
# does not say what is changed, so no aspect is set. A physical interaction
# has no direction; the file lists many of them in both orientations, and the
# parser marks the second listing as a mirror, which is not mapped.
EFFECTS = {
    '+': (slots.affects, DirectionQualifierEnum.increased),
    '-': (slots.affects, DirectionQualifierEnum.decreased),
    '_': (slots.affects, DirectionQualifierEnum.decreased),
    '0': (slots.interacts_with, None),
}
# Nodes that are not one gene although the file gives them the Entrez Gene ID
# of a gene or locus with the same symbol or alias. The links of these nodes
# are not emitted: the ID identifies another entity and the file has no
# identifier of the node itself.
#
# Small molecules: their function in the file is ``Molecule`` or ``Lipid``, the
# table of Cui et al. 2007 names them as below, and their links are those of
# the molecule (adenylate cyclases activate and phosphodiesterases inhibit
# ``CAMP``).
#
# Families of proteins: the name is that of a family, not of a gene.
# ``PKC`` (function ``Kinase``) is filed under 50818, the dystonia 10 locus
# (DYT10, alias PKC; a replaced record), ``ACTIN`` (function ``Cytoskeleton``)
# under 280979, the beta-actin gene of the cow, and ``CYP3A`` under 1574, the
# withdrawn NCBI record of the cytochrome P450 subfamily 3A.
NOT_GENES = {
    '8': 'AA',  # arachidonic acid
    '353': 'AMP',  # written APRT in version 6
    '820': 'CAMP',  # cyclic AMP; the gene is cathelicidin
    '51761': 'ATP',  # written ATP8A2 in version 6
    '92170': 'GTP',  # written MTG1 in version 6
    '50818': 'PKC',  # protein kinase C family; DYT10, a dystonia locus
    '280979': 'ACTIN',  # actin family; Bos taurus ACTB
    '1574': 'CYP3A',  # cytochrome P450 subfamily 3A; a withdrawn record
}


# The table of Cui et al. 2007 gives the nodes that are not genes the ID 0.
_entrez_re = re.compile(r'[1-9]\d*')


def _effect(row: dict[str, Any]) -> tuple[Any, Any]:
    link_type = row.get('link_type')
    if link_type not in EFFECTS:
        raise ValueError(f'Unknown Wang network link type: {link_type!r}')
    return EFFECTS[link_type]


f = FieldConfig()


def _node(role: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.ENTREZ, value=f(f'{role}_entrez')),
            CV(term=Namespace.GENESYMBOL, value=f(f'{role}_symbol')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=TAXON),
        ),
    )


# The source node of the link is the subject and its target node the object.
interactions_schema = RelationBuilder(
    subject=_node('source'),
    predicate=lambda row: _effect(row)[0],
    object=_node('target'),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.object_direction_qualifier,
            value=lambda row: _effect(row)[1],
        ),
    ),
)


def map_interaction(row: dict[str, Any]) -> Relation | None:
    """A link of the network; None if one of its nodes is not a gene or the link repeats one."""
    _effect(row)
    if row.get('mirror'):
        return None
    for role in ('source', 'target'):
        entrez = row.get(f'{role}_entrez') or ''
        if entrez in NOT_GENES or not _entrez_re.fullmatch(entrez):
            return None
    return interactions_schema(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=map_interaction,
        raw_parser=iter_network,
    ),
)
