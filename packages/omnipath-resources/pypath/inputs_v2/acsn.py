"""
Parse the ACSN 2.0 binary relations between proteins and emit Relation records.

The Atlas of Cancer Signalling Network is a set of CellDesigner maps: a
reaction network in the process description notation. ACSN also publishes its
own projection of the master map onto pairs of proteins, a table with one row
per pair of map entity names, the kinds of relation found between the two and
the PubMed IDs of the reactions involved. This module ports that table as the
source gives it; it does not read the maps.

Participants. The table names its nodes by the entity names of the maps: HUGO
symbols where the entity is one gene product, and other names (most of them
marked with a final ``*``) for aliases and for generic nodes that stand for
several gene products. ACSN publishes a second file with the HUGO names of
every entity; nothing is read from the spelling of a name. A node is

- one gene product, named by the symbol, when its name is among its HUGO names;
- one gene product, with the entity name kept as a name, when the file gives
  it exactly one HUGO name;
- a family with the listed gene products as members when the file gives it
  several HUGO names, none of them its own name. Nothing is expanded to pairs.

Relations. ``INTERACTION_TYPE`` is the set of CellDesigner reaction and
modification types that connect the pair in the map. The crosswalk in
``INTERACTION_TYPES`` is exact; a value outside it is an error. A row with
several kinds of type gives one relation per kind, each with the PubMed IDs
of the row (the source does not say which publication supports which type).
The table has one row per unordered pair. The first entity is the one that
modifies the reaction of the second in the master map for about 90% of the
modification types; for about 4% the map has the modification in the other
direction only and for 5% in both, and the table still has one row. Relations
are emitted in the order of the table.

acsn.curie.fr no longer resolves (October 2026). The files come from the copy
of the ACSN site at acsn.vincent-noel.fr, the host to which navicell.curie.fr
now redirects. The table is identical to the last Internet Archive capture of
the original (18 March 2025).
"""

from __future__ import annotations

from collections.abc import Iterable
from functools import partial
import re
from typing import Any

from biolink_model.datamodel.model import (
    DirectionQualifierEnum,
    GeneFamily,
    MicroRNA,
    Protein,
    ProteinFamily,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.acsn import iter_relations
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.ACSN,
    name='ACSN',
    full='Atlas of Cancer Signalling Network',
    url='http://acsn.vincent-noel.fr/ACSN2/ACSN2.html',
    # The site asks for a citation and names no license.
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='26192618',
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS, Content.CAUSAL_REGULATION),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Binary relations between proteins derived from the ACSN 2.0 cancer '
        'signaling maps: complex formation and regulation of reactions.'
    ),
    description=(
        'The Atlas of Cancer Signalling Network (ACSN) is a collection of '
        'manually curated maps of the molecular mechanisms of cancer, drawn '
        'with CellDesigner as a reaction network: cell cycle and DNA repair, '
        'cell survival, regulated cell death, telomere maintenance, EMT and '
        'senescence, invasion and motility, and maps of the tumor '
        'microenvironment. ACSN 2.0 has more than 8,000 reactions based on '
        'about 4,500 publications. Besides the maps, ACSN distributes a '
        'projection of the master map onto pairs of proteins.'
    ),
)

_BASE_URL = 'http://acsn.vincent-noel.fr/ACSN2/downloads/'
# Was https://acsn.curie.fr/ACSN2/downloads/ (same file names).
download_relations = Download(
    url=_BASE_URL + 'ACSN2_binary_relations_between_proteins_with_PMID.txt',
    filename='ACSN2_binary_relations_between_proteins_with_PMID.txt',
    subfolder='acsn',
)
download_names = Download(
    url=_BASE_URL + 'ACSN2_HUGO_Correspondence.gmt',
    filename='ACSN2_HUGO_Correspondence.gmt',
    subfolder='acsn',
)

# ACSN is a human resource: its entities are annotated with HUGO names.
TAXON = 'NCBITaxon:9606'
# What a relation of the table says about the pair.
ASSOCIATION = 'association'
MODULATION = 'modulation'
POSITIVE_INFLUENCE = 'positive_influence'
NEGATIVE_INFLUENCE = 'negative_influence'

# Every value of ``INTERACTION_TYPE`` (the cell is a ``;`` separated set) and
# the statement it gives. The values are CellDesigner types.
INTERACTION_TYPES = {
    # The two take part in an association reaction, are the two components
    # of a complex (PROTEIN_INTERACTION) or two of the components of a larger
    # one (COMPLEX_EXPANSION, the source's own expansion of its complexes).
    'HETERODIMER_ASSOCIATION': ASSOCIATION,
    'PROTEIN_INTERACTION': ASSOCIATION,
    'COMPLEX_EXPANSION': ASSOCIATION,
    # The first modifies a reaction in which the second is a reactant or a
    # product. The sign of the modification (catalysis, stimulation,
    # inhibition) is its effect on the reaction: the reaction can be the
    # transcription, the binding, the modification or the degradation of the
    # second, so the table does not say whether the second is increased or
    # decreased. UNKNOWN_ marks a modification drawn as uncertain.
    'CATALYSIS': MODULATION,
    'UNKNOWN_CATALYSIS': MODULATION,
    'PHYSICAL_STIMULATION': MODULATION,
    'TRIGGER': MODULATION,
    'MODULATION': MODULATION,
    'INHIBITION': MODULATION,
    'UNKNOWN_INHIBITION': MODULATION,
    # Reduced notation: an arrow drawn from the first entity to the second.
    'POSITIVE_INFLUENCE': POSITIVE_INFLUENCE,
    'UNKNOWN_POSITIVE_INFLUENCE': POSITIVE_INFLUENCE,
    'NEGATIVE_INFLUENCE': NEGATIVE_INFLUENCE,
    'UNKNOWN_NEGATIVE_INFLUENCE': NEGATIVE_INFLUENCE,
    # The first is a reactant and the second a product of one reaction. Most
    # such pairs are two names of one entity (a gene and its protein, a
    # protein and its cleaved form, a protein and its family): no relation
    # between two proteins is stated.
    'STATE_TRANSITION': None,
    'UNKNOWN_TRANSITION': None,
    'KNOWN_TRANSITION_OMITTED': None,
    'TRANSCRIPTION': None,
    'TRANSLATION': None,
    'TRANSPORT': None,
    'DISSOCIATION': None,
}
# Statement -> predicate and direction of the effect on the second entity.
STATEMENTS = {
    ASSOCIATION: (slots.interacts_with, None),
    MODULATION: (slots.affects, None),
    POSITIVE_INFLUENCE: (slots.affects, DirectionQualifierEnum.increased),
    NEGATIVE_INFLUENCE: (slots.affects, DirectionQualifierEnum.decreased),
}

# An approved HUGO symbol. The correspondence file also has a dozen remains
# of the annotation syntax (``EGFRHUGO:EGFR``, ``:KIR2DL2``, ``6709``).
_symbol = re.compile(r'[A-Z][A-Z0-9]*(?:orf[0-9]+)?(?:-[A-Z0-9]+)*')
# HGNC reserves the root ``MIR`` for microRNA genes.
_mirna_symbol = re.compile(r'MIR(?:LET)?\d')


def _gene_product(symbol: str, name: str | None = None) -> Entity:
    identifiers = [Identifier(type=Namespace.GENESYMBOL, value=symbol)]
    if name and name != symbol:
        identifiers.append(Identifier(type=Namespace.NAME, value=name))
    return Entity(
        type=MicroRNA if _mirna_symbol.match(symbol) else Protein,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def participant(name: str, hugo: Iterable[str]) -> Entity | None:
    """The entity of a node of the table, from its name and its HUGO names."""
    name = (name or '').strip()
    if not name:
        return None
    symbols = [symbol for symbol in hugo if _symbol.fullmatch(symbol)]
    if name in symbols:
        return _gene_product(name)
    if len(symbols) == 1:
        return _gene_product(symbols[0], name)
    taxon = [Annotation(term=slots.in_taxon, value=TAXON)]
    if not symbols:
        # No usable HUGO name: the entity keeps the name the map gives it.
        return Entity(
            type=Protein,
            identifiers=[Identifier(type=Namespace.NAME, value=name)],
            annotations=taxon,
        )
    members = [_gene_product(symbol) for symbol in symbols]
    return Entity(
        type=(
            GeneFamily
            if all(member.type is MicroRNA for member in members)
            else ProteinFamily
        ),
        identifiers=[Identifier(type=Namespace.NAME, value=name)],
        annotations=taxon,
        membership=[Membership(member=member) for member in members],
    )


def map_relation(row: dict[str, Any]) -> Relation | None:
    """One statement about an ordered pair of ACSN entities."""
    statement = row.get('STATEMENT')
    if statement is None:
        return None
    predicate, direction = STATEMENTS[statement]
    first = participant(row.get('INTERACTOR1'), row.get('HUGO1') or ())
    second = participant(row.get('INTERACTOR2'), row.get('HUGO2') or ())
    if first is None or second is None:
        return None
    annotations = []
    if direction is not None:
        annotations.append(
            Annotation(term=slots.object_direction_qualifier, value=direction)
        )
    annotations.extend(
        Annotation(term=slots.publications, value=curie)
        # A pair without a publication has ``N/A``; the table also has a few
        # numbers that are longer than a PubMed ID or start with a zero.
        for pmid in (row.get('PMID') or '').split(';')
        if (curie := pmid_curie(pmid))
    )
    return Relation(
        subject=first,
        predicate=predicate,
        object=second,
        annotations=annotations,
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=download_relations,
        mapper=map_relation,
        raw_parser=partial(
            iter_relations,
            names_download=download_names,
            interaction_types=INTERACTION_TYPES,
            statements=STATEMENTS,
        ),
    ),
)
