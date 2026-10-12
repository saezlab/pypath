"""
Parse the PhosphoNetworks high-resolution network and emit Relation records.

Source file (``highResolutionNetwork.csv``; despite its name it is tab
separated). A line starting with ``>`` names a substrate by gene symbol; the
lines below it, up to the next such line, each assign one kinase to one
phosphorylation site of that substrate:

    >ABI2
    PPTQKPPSPPMSGKG    S177    CDK5    2.0987
    PPTQKPPSPPMSGKG    S177    DYRK4   3.5732

The columns are the sequence window of the site (seven residues on each side
of it, "-" padding a window that runs past an end of the sequence), the site
as residue letter and position, the kinase's gene symbol and a score. Each
line is one site-level observation.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phosphonetworks import iter_network
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

# The site states that the network is of human proteins; the file has no
# organism column.
PHOSPHONETWORKS_TAX_ID = '9606'


def _substrate_form(row: dict[str, Any]) -> dict | None:
    """The phosphorylated residue, a feature of the substrate in this record.

    The position is the source's; the file names no sequence it refers to, so
    the coordinate has no reference sequence. The description is the site and
    its sequence window as written.
    """
    if row.get('position') is None:
        return None
    return {
        'modifications': [
            {
                'term': 'phosphorylated residue',
                'residue': row['residue'],
                'position': row['position'],
                'end_position': row['position'],
                'coordinate_reference': {
                    'identifier': None,
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
                'description': f'{row["site"]} ({row["window"]})',
            }
        ],
    }


config = ResourceConfig(
    id=ResourceCv.PHOSPHONETWORKS,
    name='PhosphoNetworks',
    url='https://phosphonetworks.org/',
    # "Free to academic and non-profit organizations"; commercial users are
    # asked to contact the authors before downloading.
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='24227675',
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS,),
    # Kinase-substrate pairs come from phosphorylation reactions on protein
    # microarrays, refined with a Bayesian model and combined with known pairs;
    # the sites are mass spectrometry sites on those substrates, assigned to a
    # kinase when they match its motif. The file does not mark which kinase-site
    # line rests on which of these.
    knowledge_origin=(
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Human kinase-substrate phosphorylation sites of the PhosphoNetworks '
        'high-resolution network, with the sequence window and a motif score.'
    ),
    description=(
        'PhosphoNetworks is a database of human phosphorylation networks from '
        'the labs of Heng Zhu, Jin Zhang and Jiang Qian. Phosphorylation '
        'reactions of 289 human kinases on protein microarrays of 4,191 human '
        'proteins gave kinase-substrate relationships, which a Bayesian '
        'analysis refined and which were combined with known relationships. '
        'Phosphorylation motifs were derived for the kinases and matched '
        'against phosphorylation sites found by mass spectrometry on the '
        'substrates, so that the network connects kinases to specific '
        'phosphorylation sites. The high-resolution network lists each '
        'kinase-site pair with the sequence window of the site and a score.'
    ),
)

f = FieldConfig()


def _protein(column: str, molecular_form=None) -> EntityBuilder:
    return EntityBuilder(
        entity_type=model.Protein,
        molecular_form=molecular_form,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(column)),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{PHOSPHONETWORKS_TAX_ID}'),
        ),
    )


enzyme_substrate_schema = RelationBuilder(
    subject=_protein('enzyme'),
    predicate=slots.affects,
    # The phosphorylated residue is in the substrate's molecular form.
    object=_protein('substrate', molecular_form=_substrate_form),
    annotations=AnnotationsBuilder(
        # The kinase increases the substrate's phosphorylation; the file says
        # nothing about the substrate's activity.
        CV(
            term=slots.object_direction_qualifier,
            value=model.DirectionQualifierEnum.increased,
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=model.GeneOrGeneProductOrChemicalEntityAspectEnum.phosphorylation,
        ),
        # The score of the last column, 0.44 to 16 in the file; the paper
        # calls it the motif match score. The file has no explicit scale.
        CV(
            term=slots.has_confidence_score,
            value=lambda row: measurement(row.get('score'), source_field='score'),
        ),
    ),
)

resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=Download(
            url='https://phosphonetworks.org/download/highResolutionNetwork.csv',
            filename='highResolutionNetwork.csv',
            subfolder='phosphonetworks',
            large=True,
        ),
        mapper=enzyme_substrate_schema,
        raw_parser=iter_network,
    ),
)
