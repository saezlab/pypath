"""
Map the enzyme-substrate sites of dbPTM to Relation records.

dbPTM integrates modification sites from other databases. Only its
releases up to dbPTM 3.0 named the enzyme of a site, taken from the
"by <enzyme>" notes of Swiss-Prot. The dbPTM 3.0 files are gone from the
dbPTM hosts, so the module reads the table OmniPath rescued from them
(``rescued.omnipathdb.org/old_dbptm.tab``), which is a product of the legacy
pypath reader and not a file of dbPTM; ``parsers/dbptm.py`` describes it and
reads it. A site without an enzyme is not an enzyme-substrate relation and is
not read.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.dbptm import BACTERIAL_NAMES, iter_enzyme_sites
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

_ASPECT = model.GeneOrGeneProductOrChemicalEntityAspectEnum
# The modification types of the rows that name an enzyme, with the aspect of
# the substrate the enzyme increases and the PSI-MOD name of the modified
# residue. The other types of the table occur only without an enzyme.
_MODIFICATIONS = {
    'phosphorylation': (_ASPECT.phosphorylation, 'phosphorylated residue'),
    'methylation': (_ASPECT.methylation, 'methylated residue'),
    'acetylation': (_ASPECT.acetylation, 'acetylated residue'),
    'palmitoylation': (_ASPECT.palmitoylation, 'palmitoylated residue'),
}
# The table holds the human sites of dbPTM: the substrate of every row that
# names an enzyme is a reviewed human UniProt entry (checked in UniProt,
# October 2026). The table has no organism column.
DBPTM_TAX_ID = '9606'


def _substrate_form(row: dict[str, Any]) -> dict | None:
    """The modified residue, a feature of the substrate in this record.

    dbPTM maps its sites to the UniProt entry of the substrate, so the
    position refers to that entry's sequence as it was at the time. The
    description is the sequence window of the site as written.
    """
    if row.get('position') is None or not row.get('substrate'):
        return None
    return {
        'modifications': [
            {
                'term': _MODIFICATIONS[row['typ']][1],
                'residue': row.get('residue'),
                'position': row['position'],
                'end_position': row['position'],
                'coordinate_reference': {
                    'identifier': {'ns': 'uniprot', 'id': row['substrate']},
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
                'description': row.get('instance') or None,
            }
        ],
    }


def _enzyme_taxon(row: dict[str, Any]) -> str | None:
    """Human for a protein of the human table; none for a group or a bacterial name."""
    if row.get('enzyme_group') or row.get('enzyme') in BACTERIAL_NAMES:
        return None
    return f'NCBITaxon:{DBPTM_TAX_ID}'


config = ResourceConfig(
    id=ResourceCv.DBPTM,
    name='dbPTM',
    url='https://biomics.lab.nycu.edu.tw/dbPTM/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.NAR_OPEN_ACCESS,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='23193290',
    topics=(Topic.SIGNALING,),
    content=(Content.ENZYME_PTM_RELATIONSHIPS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Human modification sites with the enzyme Swiss-Prot names for them, '
        'from a rescued table of dbPTM 3.0.'
    ),
    description=(
        'dbPTM is an integrated resource of protein post-translational '
        'modifications. It collects experimentally verified modification '
        'sites from other databases and from the literature and maps them '
        'to UniProtKB entries. Releases up to dbPTM 3.0 also gave the '
        'enzyme of a site where Swiss-Prot names one; the files of later '
        'releases list sites without enzymes.'
    ),
)

f = FieldConfig()

# A line whose enzyme label does not name one enzyme yields no relation: the
# enzyme has no identifier, so the builder returns None.
enzyme_substrate_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=lambda row: (
            model.ProteinFamily if row.get('enzyme_group') else model.Protein
        ),
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=f('enzyme_uniprot')),
            CV(term=Namespace.GENESYMBOL, value=f('enzyme_genesymbol')),
            CV(term=Namespace.NAME, value=f('enzyme_name')),
            CV(term=Namespace.NAME, value=f('enzyme_group')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=_enzyme_taxon),
        ),
    ),
    predicate=slots.affects,
    object=EntityBuilder(
        entity_type=model.Protein,
        # The modified residue is in the substrate's molecular form.
        molecular_form=_substrate_form,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=f('substrate')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{DBPTM_TAX_ID}'),
        ),
    ),
    annotations=AnnotationsBuilder(
        # The enzyme increases this modification of the substrate; the table
        # says nothing about the substrate's activity.
        CV(
            term=slots.object_direction_qualifier,
            value=model.DirectionQualifierEnum.increased,
        ),
        CV(
            term=slots.object_aspect_qualifier,
            value=lambda row: _MODIFICATIONS[row['typ']][0],
        ),
        CV(
            term=slots.publications,
            value=lambda row: [f'PMID:{pmid}' for pmid in row.get('pmids') or ()],
        ),
        CV(term=slots.supporting_data_source, value=f('source')),
    ),
)

resource = Resource(
    config,
    enzyme_substrate=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/old_dbptm.tab',
            filename='old_dbptm.tab',
            subfolder='dbptm',
            large=True,
        ),
        mapper=enzyme_substrate_schema,
        raw_parser=iter_enzyme_sites,
    ),
)
