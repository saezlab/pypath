"""
Parse the KEA kinase-substrate phosphosite database and emit Relation records.

Source file (tab separated, no header, one kinase-site-reference-database
combination per line; read by ``parsers/kea.py``):

    ABL1    ABL1_Y226    10964922    phosphosite

Each line is one site-level observation: the kinase's gene symbol, the
substrate's gene symbol with the phosphorylated residue and its position, a
reference and the database KEA took the line from.
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.kea import iter_kea
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

# The values of the fourth column. NetworKIN lines are kinase predictions for
# known sites and cite the NetworKIN publication; the others are curated
# databases. The label is kept as the source writes it.
SOURCE_DATABASES = frozenset(
    {'hprd', 'kinexus', 'networkin', 'phosphoelm', 'phosphosite'}
)
# KEA converted mouse and rat records to human gene symbols (PMID 19176546),
# so every participant is a human protein. The file has no organism column.
KEA_TAX_ID = '9606'


def _source_database(row: dict[str, object]) -> str:
    value = str(row.get('source_database') or '')
    if value not in SOURCE_DATABASES:
        raise ValueError(f'Unknown KEA source database: {value!r}')
    return value


def _substrate_form(row: dict[str, object]) -> dict | None:
    """The phosphorylated residue, a feature of the substrate in this record.

    The position is the source's; KEA names no sequence it refers to, so the
    coordinate has no reference sequence.
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
                'description': row['substrate_site'],
            }
        ],
    }


config = ResourceConfig(
    id=ResourceCv.KEA,
    name='KEA',
    url='https://www.maayanlab.net/KEA2/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='19176546',
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary='Kinase-substrate phosphorylation sites with references and source databases.',
    description=(
        'KEA (Kinase Enrichment Analysis) is a kinase enrichment tool of the '
        "Ma'ayan Laboratory. Its background database consolidates "
        'kinase-substrate interactions from PhosphoSitePlus, HPRD, '
        'Phospho.ELM, NetworKIN and Kinexus PhosphoNET, with mouse and rat '
        'records converted to human gene symbols. The phosphosite-level '
        'table of KEA2 lists each kinase with the substrate residue it '
        'phosphorylates, a reference and the database of origin.'
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
            CV(term=slots.in_taxon, value=f'NCBITaxon:{KEA_TAX_ID}'),
        ),
    )


# A line without a kinase or without a parsable site yields no relation: the
# participant has no identifier, so the builder returns None.
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
        CV(
            term=slots.publications,
            value=f('pmid', transform=lambda pmid: f'PMID:{pmid}'),
        ),
        CV(term=slots.supporting_data_source, value=_source_database),
    ),
)

resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=Download(
            url=(
                'https://www.maayanlab.net/KEA2/gsl/'
                'kinase_substrate_phospho-site_level_pmid_resource_database.txt'
            ),
            filename='kinase_substrate_phospho-site_level_pmid_resource_database.txt',
            subfolder='kea',
            large=True,
        ),
        mapper=enzyme_substrate_schema,
        raw_parser=iter_kea,
    ),
)
