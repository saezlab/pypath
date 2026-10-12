"""
Parse the interaction list of the macrophage pathway map and emit Relation records.

Additional file 2 of Raza et al. 2010 lists the interactions drawn on the
authors' map of macrophage activation pathways. Each row names two interacting
partners by official gene symbol and Entrez Gene ID, with the kind of
interactant and its label on the map, the interaction type, its location and
one supporting PubMed ID; an interaction with several publications takes
several rows.

Only rows between two single proteins, or a protein and a gene, are mapped.
The cells of partners drawn as complexes sometimes list all members and
sometimes only the member that takes part, so they cannot be read as complexes.
"""

from __future__ import annotations
from collections.abc import Callable
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    DirectionQualifierEnum,
    Gene,
    GeneOrGeneProductOrChemicalEntityAspectEnum,
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
from pypath.inputs_v2.parsers.macrophage import iter_interactions

config = ResourceConfig(
    id=ResourceCv.MACROPHAGE,
    name='Macrophage',
    full='Integrated map of macrophage pathways',
    url='https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-4-63',
    # The article is under CC BY 2.0, which has no term yet.
    license=LicenseCV.CC_BY_2_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='20470404',
    topics=(Topic.SIGNALING,),
    content=(
        Content.CAUSAL_REGULATION,
        Content.PHYSICAL_INTERACTIONS,
        Content.ENZYME_PTM_RELATIONSHIPS,
    ),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Binding, activation, inhibition, modification and gene expression relationships between proteins of macrophage activation pathways, with publications.',
    description='Raza et al. 2010 constructed an integrated map of the pathways of macrophage activation: Toll-like receptor and other pathogen recognition signaling, NF-kB, interferon and MAP kinase signaling, apoptosis, antigen presentation and the ubiquitin-proteasome system. The interactions on the map were curated manually from the literature and are listed, with their components named by human gene symbols, the type of the interaction, its cellular location and the supporting publications, in a supplementary table of the article.',
)
download = Download(
    url='https://static-content.springer.com/esm/art%3A10.1186%2F1752-0509-4-63/MediaObjects/12918_2010_452_MOESM2_ESM.XLS',
    filename='12918_2010_452_MOESM2_ESM.XLS',
    subfolder='macrophage',
    default_mode='rb',
    encoding=None,
)
# "All components have been depicted using standard human gene nomenclature
# (HGNC)" (PMID 20470404); the Gene IDs are those of the human genes. The
# pathways are "based on data derived from a range of different cellular
# systems and mammalian species (human and mouse)": the ``Species`` column
# names the species or cell line of the experiment. A row without a species
# is human by the nomenclature. Where the experiment's species is not human,
# no taxon is stated for the proteins: the identifiers are those of the human
# genes, and a mouse or hamster taxon on them would contradict the
# identifiers (the rule ``netpath.py`` follows).
TAXON = 'NCBITaxon:9606'
SPECIES_TAXA = {
    None: TAXON,
    'Human': TAXON,
    'Mouse': None,
    'Human epithelial kidney cell line': TAXON,
    'Human monocytic cell line': TAXON,
    'Chinese hamster lung fibroblast cell line': None,
    'Human epithelial kidney cell line, Chinese hamster lung fibroblast cell line': None,
}
# Interactant types that name one gene product or one gene. The other values
# of the column (complexes, DNA elements, drugs...) are checked by the parser.
INTERACTANT_TYPES = {'Protein': Protein, 'Gene': Gene}

_DIRECTION = DirectionQualifierEnum
_ASPECT = GeneOrGeneProductOrChemicalEntityAspectEnum
_INTERACTION = (slots.interacts_with, None, None)
_PHOSPHORYLATION = (slots.affects, _DIRECTION.increased, _ASPECT.phosphorylation)
_DEPHOSPHORYLATION = (slots.affects, _DIRECTION.decreased, _ASPECT.phosphorylation)
_ACTIVATION = (slots.affects, _DIRECTION.increased, _ASPECT.activity)
_INHIBITION = (slots.affects, _DIRECTION.decreased, _ASPECT.activity)
_EXPRESSION = (slots.affects, None, _ASPECT.expression)
_CLEAVAGE = (slots.affects, _DIRECTION.increased, _ASPECT.cleavage)
# Every value of the interaction type column, as the sheet writes it, and the
# predicate, direction and aspect of the effect of partner 1 on partner 2.
# None: the value states no relation between the two partners.
INTERACTION_TYPES = {
    'Binding': _INTERACTION,
    'Binbing': _INTERACTION,
    'Activation': _ACTIVATION,
    # The partner is activated; that it oligomerizes is not represented.
    'Activation:oligomerization': _ACTIVATION,
    'Inhibition': _INHIBITION,
    'Inhibition by Binding': _INHIBITION,
    'Enhancement': (slots.affects, _DIRECTION.increased, None),
    'Phosphorylation': _PHOSPHORYLATION,
    'Phosohorylation': _PHOSPHORYLATION,
    'phosphrylation': _PHOSPHORYLATION,
    'autoPhosphorylation': _PHOSPHORYLATION,
    'auto-Phosphorylation': _PHOSPHORYLATION,
    'dePhosphorylation': _DEPHOSPHORYLATION,
    'Inhibition (of Phosphorylation)': _DEPHOSPHORYLATION,
    'Ubiquitination': (slots.affects, _DIRECTION.increased, _ASPECT.ubiquitination),
    'acetylation': (slots.affects, _DIRECTION.increased, _ASPECT.acetylation),
    'sumo-modification': (slots.affects, _DIRECTION.increased, _ASPECT.sumoylation),
    # Partner 1 is read as the protease. The sheet and the article state no
    # role of the partners; in four mapped rows (LMNA with CASP6, BID with
    # CASP8, two publications each) the substrate is partner 1, and the
    # relation says that the substrate cleaves the protease. Not corrected.
    'Cleavage': _CLEAVAGE,
    'autoCleavage': _CLEAVAGE,
    'Degradation': (slots.affects, _DIRECTION.increased, _ASPECT.degradation),
    # Partner 1 moves partner 2 to another compartment; no sign is stated.
    'Translocation': (slots.affects, None, _ASPECT.transport),
    'Activation of Translocation': (slots.affects, _DIRECTION.increased, _ASPECT.transport),
    # The value does not say that expression goes up; repression has its own.
    'Gene Expression': _EXPRESSION,
    'gene experssion': _EXPRESSION,
    'Gene Expression- Inhibition': (slots.affects, _DIRECTION.decreased, _ASPECT.expression),
    # A complex falls apart; partner 2 is a product or empty.
    'Dissociation': None,
    # Where a protein is, or what it does; both partners are the same protein.
    'localisation': None,
    'localization': None,
    'function': None,
    # A lipid kinase and its substrate.
    'catalysis': None,
    'N/A': None,
    None: None,
}
# Interaction types that a protein can have with itself. In other rows with
# the same gene on both sides, partner 2 is the product of partner 1
# (``proIL-18`` and ``IL18``) or the row describes one protein.
SELF_INTERACTION_TYPES = frozenset(
    {'Binding', 'Binbing', 'autoPhosphorylation', 'auto-Phosphorylation', 'autoCleavage'}
)


def _taxon(row: dict[str, Any]) -> str | None:
    species = row.get('Species')
    if species not in SPECIES_TAXA:
        raise ValueError(f'Unknown macrophage species: {species!r}')
    return SPECIES_TAXA[species]


def _effect(row: dict[str, Any]) -> tuple[Any, Any, Any] | None:
    interaction_type = row.get('Interaction Type')
    if interaction_type not in INTERACTION_TYPES:
        raise ValueError(
            f'Unknown macrophage interaction type: {interaction_type!r}'
        )
    return INTERACTION_TYPES[interaction_type]


f = FieldConfig(
    extract={
        'gene_id': r'^(\d+)$',
        # An official symbol alone. Cells that add an alias or a form
        # (``CHUK(IKKα)``, ``proCASP1``) are identified by the Gene ID only.
        'symbol': r'^([A-Z][A-Z0-9]*(?:-[A-Z0-9]+)*)$',
        # Two cells end in a comma; three hold numbers too long for a PubMed ID.
        'pmid': rf'^({PMID_PATTERN}),?$',
    },
    transform={'pmid_curie': lambda v: f'PMID:{v}'},
)


def _participant(partner: str) -> Callable[[dict[str, Any]], Any]:
    kind = f'{partner} Interactant Type'
    gene_id = f(f'{partner} Gene ID', extract='gene_id')
    builder = EntityBuilder(
        entity_type=lambda row: INTERACTANT_TYPES[row[kind]],
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.ENTREZ, value=gene_id),
            CV(
                term=Namespace.GENESYMBOL,
                value=f(f'{partner} Official Gene Symbol', extract='symbol'),
            ),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=_taxon),
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        # One protein or gene with its Gene ID. Complexes, DNA elements,
        # pathogen components, drugs and cells listing several genes are not
        # participants this mapping can identify.
        if row.get(kind) not in INTERACTANT_TYPES or not gene_id.extract(row):
            return None
        return builder.build(row)

    return build


# Partner 1 is the subject and partner 2 the object.
interactions_schema = RelationBuilder(
    subject=_participant('Partner 1'),
    predicate=lambda row: _effect(row)[0],
    object=_participant('Partner 2'),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.object_direction_qualifier,
            value=lambda row: _effect(row)[1],
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=lambda row: _effect(row)[2],
        ),
        CV(
            term=slots.publications,
            value=f('NCBI-PubMed ID', extract='pmid', transform='pmid_curie'),
        ),
        CV(term=slots.description, value=f('Interaction-description')),
    ),
)


def interactions_mapper(row: dict[str, Any]) -> Any:
    """Relation of a row; None where the type states none or a partner is not one protein or gene."""
    # Checked on every row; a check inside the builders would be swallowed.
    _taxon(row)
    if _effect(row) is None:
        return None
    if (
        row.get('Partner 1 Gene ID') == row.get('Partner 2 Gene ID')
        and row.get('Interaction Type') not in SELF_INTERACTION_TYPES
    ):
        return None
    return interactions_schema.build(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=interactions_mapper,
        raw_parser=iter_interactions,
    ),
)
