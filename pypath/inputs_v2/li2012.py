"""
Parse the curated human phosphotyrosine circuits of Li et al. 2012.

Supplemental File S1 lists 583 literature-derived "tyrosine kinase circuits":
a tyrosine kinase family (the writer), a substrate with a phosphotyrosine
position, and a family of SH2 or PTB domain proteins that binds the
phosphorylated site (the reader):

    substrate family  protein/pY position  TK writers  SH2/PTB Readers
    BCR               BCR/177              ABL         GRB2 SH2

The paper assembled the circuits from curated kinase-substrate and
substrate-reader interactions, so a site with two writers and two readers
gives four rows. This module reads the two kinds of interaction back out:
one record per distinct writer and site, and one per distinct site and
reader. Writers and readers are protein families, whose members are in
Supplemental File S3. The predicted circuits (Files S2 and S4) are not read.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.li2012 import (
    FAMILY,
    SITE,
    WRITER,
    iter_families,
    iter_kinase_substrate,
    iter_phosphosite_binding,
)
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

# The data set is the human phosphotyrosine signaling network; the sheets
# have no organism column.
LI2012_TAX_ID = '9606'


def _phosphotyrosine_form(row: dict[str, Any]) -> dict | None:
    """The phosphotyrosine site, a feature of the substrate in this record.

    Every site of the data set is a tyrosine. The paper names no sequence
    the position refers to, so the coordinate has no reference sequence.
    """
    if row.get('position') is None:
        return None
    return {
        'modifications': [
            {
                'term': 'phosphorylated residue',
                'residue': 'Y',
                'position': row['position'],
                'end_position': row['position'],
                'coordinate_reference': {
                    'identifier': None,
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
                'description': row[SITE],
            }
        ],
    }


config = ResourceConfig(
    id=ResourceCv.LI2012,
    name='Li2012',
    url='https://doi.org/10.1101/gr.128819.111',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_NC_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='22194470',
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS, Content.PHYSICAL_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Curated tyrosine kinase-substrate sites and SH2/PTB domain binding '
        'to phosphotyrosine sites.'
    ),
    description=(
        'Li et al. 2012 analysed the evolution of the human phosphotyrosine '
        'signaling network and its use in cancer. The study compiled 583 '
        'literature-derived tyrosine kinase circuits, each a tyrosine kinase '
        'family, a substrate phosphotyrosine site and an SH2 or PTB domain '
        'protein family that binds the phosphorylated site, from the '
        'literature, Phospho.ELM and PhosphoSitePlus, and about 50,000 '
        'predicted circuits.'
    ),
)

f = FieldConfig()


def _family(column: str, molecular_form=None) -> EntityBuilder:
    # A family has only its name in the source; its members are in File S3.
    return EntityBuilder(
        entity_type=model.ProteinFamily,
        molecular_form=molecular_form,
        identifiers=IdentifiersBuilder(CV(term=Namespace.NAME, value=f(column))),
    )


def _substrate() -> EntityBuilder:
    return EntityBuilder(
        entity_type=model.Protein,
        # The phosphotyrosine is in the substrate's molecular form.
        molecular_form=_phosphotyrosine_form,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('substrate')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{LI2012_TAX_ID}'),
        ),
    )


kinase_substrate_schema = RelationBuilder(
    subject=_family(WRITER),
    predicate=slots.affects,
    object=_substrate(),
    annotations=AnnotationsBuilder(
        # The kinase family increases the substrate's phosphorylation; the
        # sheet says nothing about the substrate's activity.
        CV(
            term=slots.object_direction_qualifier,
            value=model.DirectionQualifierEnum.increased,
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=model.GeneOrGeneProductOrChemicalEntityAspectEnum.phosphorylation,
        ),
    ),
)

# The domain of the reader family that binds the phosphotyrosine, "SH2" or
# "PTB" as the circuit sheet writes it after the family name ("GRB2 SH2"). A
# family has no sequence to put the domain on, and a region in the family's
# molecular form would make a different entity of the family than the one
# File S3 defines with its members, so the domain is an attribute of the
# relation. Biolink has no slot for it.
_BINDING_DOMAIN = 'li2012:binding_domain'

# The sheet says that the SH2 or PTB domain of a reader family binds the
# phosphorylated tyrosine. It gives no interaction type or method, so this is
# an interaction, not a stated direct one, and no effect on the substrate is
# asserted.
phosphosite_binding_schema = RelationBuilder(
    subject=_family('reader_family'),
    predicate=slots.interacts_with,
    object=_substrate(),
    annotations=AnnotationsBuilder(
        CV(term=_BINDING_DOMAIN, value=f('reader_domain')),
    ),
)


def map_family(row: dict[str, Any]) -> Entity | None:
    """A kinase or reader family with its member proteins, by Entrez Gene ID."""
    name = row.get(FAMILY)
    if not name:
        return None
    return Entity(
        type=model.ProteinFamily,
        identifiers=[Identifier(type=Namespace.NAME, value=name)],
        membership=[
            Membership(
                member=Entity(
                    type=model.Protein,
                    identifiers=[
                        Identifier(type=Namespace.ENTREZ, value=member['entrez']),
                        Identifier(type=Namespace.GENESYMBOL, value=member['symbol']),
                    ],
                    annotations=[
                        Annotation(
                            term=slots.in_taxon,
                            value=f'NCBITaxon:{LI2012_TAX_ID}',
                        )
                    ],
                ),
                predicate=slots.has_member,
            )
            for member in row.get('members') or ()
        ],
    )


# The supplement is no longer served by the journal (genome.cshlp.org); this
# is OmniPath's copy of Supplementary_files_S1-S5.xls.
download = Download(
    url='https://rescued.omnipathdb.org/li2012-supplementary-S1-S5.xls',
    filename='li2012-supplementary-S1-S5.xls',
    subfolder='li2012',
    large=True,
    encoding=None,
    default_mode='rb',
)

resource = Resource(
    config,
    kinase_substrate=Dataset(
        download=download,
        mapper=kinase_substrate_schema,
        raw_parser=iter_kinase_substrate,
    ),
    phosphosite_binding=Dataset(
        download=download,
        mapper=phosphosite_binding_schema,
        raw_parser=iter_phosphosite_binding,
    ),
    families=Dataset(
        download=download,
        mapper=map_family,
        raw_parser=iter_families,
    ),
)
