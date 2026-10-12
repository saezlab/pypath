"""
Parse the hippocampal CA1 neuron signaling network and emit Relation records.

The network is file S1 of the supporting online material of Ma'ayan et al.
2005. Each line is one link between two components with one supporting
publication; a link with several publications takes several lines:

    Source Human Mouse Function Location Target Human Mouse Function Location Effects Type of interaction PubMedID
    CAMKII Q9UQM7 P11798 Kinase Cytosol SYNAPSIN P17600 O88935 Vesicle Vesicles _ Phosphorylation 7890745

A component has a name, a human and a mouse Swiss-Prot accession, a class and
a cellular location. The effect of the source on the target is activating
(``+``), inhibitory (``_``) or neutral (``0``).
"""

from __future__ import annotations
from collections.abc import Callable
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    DirectionQualifierEnum,
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
from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.ca1 import iter_links

config = ResourceConfig(
    id=ResourceCv.CA1,
    name='CA1',
    url='https://www.science.org/doi/10.1126/science.1108876',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='16099987',
    topics=(Topic.SIGNALING,),
    content=(Content.CAUSAL_REGULATION, Content.PHYSICAL_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Activating, inhibitory and neutral links between signaling proteins of the hippocampal CA1 neuron, with publications.',
    description="Ma'ayan et al. 2005 assembled the signaling network of a mammalian hippocampal CA1 neuron from the experimental literature: 545 components, from ligands and receptors to the transcriptional, translational, secretory and electrical machinery, and 1,259 links between them. The authors read the primary publications of every link to establish that the interaction is direct and functionally relevant, and classified it as activating, inhibitory or neutral, and by its biochemical mechanism.",
)
# The publisher's supplement link and the copy in PubMed Central are behind bot
# protection; this is the publisher's archive, kept by OmniPath.
download = Download(
    url='https://rescued.omnipathdb.org/maayan2005_som_external_files.zip',
    filename='maayan2005_som_external_files.zip',
    subfolder='ca1',
    default_mode='rb',
    encoding=None,
    ext='zip',
    needed=['S1.txt'],
)
# The accession columns used are the ones headed ``Human``.
TAXON = 'NCBITaxon:9606'
# Effects -> predicate and direction of the effect on the target. A neutral
# link "does not specify directionality between components" (scaffolding and
# anchoring), so it is an interaction without a direction. The sign is the
# functional effect; the mechanism in ``Type of interaction`` is independent of
# it, so no aspect is set.
EFFECTS = {
    '+': (slots.affects, DirectionQualifierEnum.increased),
    '_': (slots.affects, DirectionQualifierEnum.decreased),
    '0': (slots.interacts_with, None),
}


def _effect(row: dict[str, Any]) -> tuple[Any, Any]:
    effect = row.get('Effects')
    if effect not in EFFECTS:
        raise ValueError(f'Unknown CA1 effect: {effect!r}')
    return EFFECTS[effect]


f = FieldConfig(
    extract={'pmid': PMID_EXTRACT},
    transform={'pmid_curie': lambda v: f'PMID:{v}'},
)


def _participant(role: str) -> Callable[[dict[str, Any]], Any]:
    accession = f'{role} Human'
    builder = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=f(accession)),
            # The network's own name of the component, e.g. ``CAMKII``.
            CV(term=Namespace.NAME, value=f(role)),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=TAXON),
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        # A component without an accession is known by a name of the network
        # only (``GLYCINE``, ``CALCIUM``, ``CRE``, ``EGF``): the table does not
        # say what kind of molecule it is, so the link is not emitted.
        return builder.build(row) if row.get(accession) else None

    return build


# The source of the link is the subject and its target the object.
interactions_schema = RelationBuilder(
    subject=_participant('Source'),
    predicate=lambda row: _effect(row)[0],
    object=_participant('Target'),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.object_direction_qualifier,
            value=lambda row: _effect(row)[1],
        ),
        CV(
            term=slots.publications,
            value=f('PubMedID', extract='pmid', transform='pmid_curie'),
        ),
    ),
)
resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=interactions_schema,
        raw_parser=iter_links,
    ),
)
