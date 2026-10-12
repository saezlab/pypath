"""
Map the kinase-substrate phosphorylation sites of MIMP to Relation records.

MIMP's phosphorylation data file is not a prediction: it holds the
experimentally identified human phosphosites MIMP combined from
PhosphoSitePlus, Phospho.ELM and HPRD, each matched to the longest CCDS
isoform of its gene, and for a part of the sites the kinases those databases
name. ``parsers/mimp.py`` describes the file and reads it.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.mimp import iter_kinase_sites
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

# MIMP collected human phosphosites and matched them to human CCDS proteins
# (PMID 25938373). The file has no organism column.
MIMP_TAX_ID = '9606'


def _substrate_form(row: dict[str, Any]) -> dict | None:
    """The phosphorylated residue, a feature of the substrate in this record.

    The position counts residues of the protein that the RefSeq transcript
    in ``gene_acc`` encodes, the isoform MIMP matched the site to. The
    description is the site with its flanking residues as written.
    """
    if row.get('site_position') is None or not row.get('gene_acc'):
        return None
    return {
        'modifications': [
            {
                'term': 'phosphorylated residue',
                'residue': row['residue'],
                'position': row['site_position'],
                'end_position': row['site_position'],
                'coordinate_reference': {
                    'identifier': {'ns': 'refseq', 'id': row['gene_acc']},
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
                'description': row['flank'],
            }
        ],
    }


config = ResourceConfig(
    id=ResourceCv.MIMP,
    name='MIMP',
    url='https://mimp.baderlab.org/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.LGPL_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='25938373',
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Human kinase-substrate phosphorylation sites that MIMP combined '
        'from PhosphoSitePlus, Phospho.ELM and HPRD.'
    ),
    description=(
        'MIMP (Mutations IMpact on Phosphorylation) is a machine learning '
        'method that predicts whether a missense variant in or next to a '
        'phosphorylation site disrupts the site or creates a new one for a '
        'kinase. Its kinase specificity models are built from '
        'experimentally identified human phosphosites and their kinases, '
        'which MIMP collected from PhosphoSitePlus, Phospho.ELM and HPRD '
        'and matched to the longest CCDS isoform of each protein.'
    ),
)

f = FieldConfig()

enzyme_substrate_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=lambda row: (
            model.ProteinFamily if row.get('enzyme_group') else model.Protein
        ),
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('enzyme_genesymbol')),
            CV(term=Namespace.NAME, value=f('enzyme_name')),
            CV(term=Namespace.NAME, value=f('enzyme_group')),
        ),
        annotations=AnnotationsBuilder(
            # A gene symbol is human; a group or another name has no species.
            CV(
                term=slots.in_taxon,
                value=f(
                    'enzyme_genesymbol',
                    transform=lambda _symbol: f'NCBITaxon:{MIMP_TAX_ID}',
                ),
            ),
        ),
    ),
    predicate=slots.affects,
    object=EntityBuilder(
        entity_type=model.Protein,
        # The phosphorylated residue is in the substrate's molecular form.
        molecular_form=_substrate_form,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('gene')),
            # The transcript of the isoform, a nucleotide accession.
            CV(term=Namespace.REFSEQ, value=f('gene_acc')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{MIMP_TAX_ID}'),
        ),
    ),
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
        # The databases MIMP found the site in. The file does not say which
        # of them names the kinase, and gives a number of publications
        # (npmid) but no PubMed ID.
        CV(term=slots.supporting_data_source, value=lambda row: row.get('sources')),
    ),
)

resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=Download(
            url=(
                'https://mimp.baderlab.org/fetch_data/'
                'phosphorylation_data.tab/phosphorylation_data.tab'
            ),
            filename='phosphorylation_data.tab',
            subfolder='mimp',
            large=True,
        ),
        mapper=enzyme_substrate_schema,
        raw_parser=iter_kinase_sites,
    ),
)
