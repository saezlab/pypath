"""Human solute carriers and the substrates and ions they transport.

The substrate-based ontology of Meixner et al. 2020 annotates 446 human SLCs with
their substrates and coupled ions, mapped to ChEBI, with references, transport
mechanism, substrate class and subcellular localization. Legacy OmniPath served
these pairs as its ``SLC`` transporter-substrate network (``slc.slc_interactions``);
this is the same content. The bioparadigms SLC tables are a separate resource
(``slctables``).
"""

from __future__ import annotations

import functools

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.interaction_profiles import TRANSPORT_QUALIFIERS
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from omnipath_resources.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from omnipath_resources.inputs_v2.parsers import slc as _parsers
from omnipath_resources.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from omnipath_resources.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.SLC,
    name='SLC',
    url='https://opendata.cemm.at/gsflab/slcontology/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='32697042',
    content=(Content.MOLECULAR_TRANSPORT,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Human solute carriers and the ChEBI-mapped substrates and ions they transport.',
    description=(
        'A substrate-based ontology for human solute carriers (Meixner et al. 2020): '
        'the substrates and coupled ions of 446 human SLC transporters, mapped to '
        'ChEBI, with literature references, transport mechanism, substrate class '
        'and subcellular localization.'
    ),
)

# The rescued copies of Tables EV1 and EV2 of the article, as the legacy module used.
_annotation_download = Download(
    url='https://rescued.omnipathdb.org/meixner2020-S1.xlsx',
    filename='meixner2020-S1.xlsx',
    subfolder='slc',
    large=True,
    encoding=None,
    default_mode='rb',
)
_mapping_download = Download(
    url='https://rescued.omnipathdb.org/meixner2020-S2.xlsx',
    filename='meixner2020-S2.xlsx',
    subfolder='slc',
    large=True,
    encoding=None,
    default_mode='rb',
)

f = FieldConfig(extract={'chebi': r'CHEBI:(\d+)'})

# Source attributes of a pair, kept as the tables write them.
_SOURCE_ATTRIBUTES = (
    ('slc:role', 'role'),
    ('slc:transport_mechanism', 'transport_mechanism'),
    ('slc:substrate_class', 'substrate_class'),
    ('slc:subcellular_localization', 'localization'),
)

transport_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=model.Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('genesymbol')),
            CV(term=Namespace.ENTREZ, value=f('entrez')),
            CV(term=Namespace.ENSG, value=f('ensg')),
        ),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value='NCBITaxon:9606')),
    ),
    predicate=slots.affects,
    object=EntityBuilder(
        entity_type=model.ChemicalEntity,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.CHEBI, value=f('chebi', extract='chebi')),
            CV(term=Namespace.NAME, value=f('substrate_name')),
        ),
        annotations=AnnotationsBuilder(),
    ),
    annotations=AnnotationsBuilder(
        *(CV(term=term, value=value) for term, value in TRANSPORT_QUALIFIERS),
        CV(term=slots.publications, value=lambda row: ['PMID:' + p for p in row.get('pmids') or ()]),
        *(CV(term=term, value=f(field)) for term, field in _SOURCE_ATTRIBUTES),
    ),
)

resource = Resource(
    config,
    transports=Dataset(
        download=_annotation_download,
        mapper=transport_schema,
        raw_parser=functools.partial(
            _parsers.transport_substrates,
            mapping_download=_mapping_download,
        ),
    ),
)
