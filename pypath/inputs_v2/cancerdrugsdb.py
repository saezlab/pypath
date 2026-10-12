"""
Parse the CancerDrugs_DB table of licensed cancer drugs and emit Relation records.

Source file (``cancerdrugsdb.txt``, tab separated, one row per drug):

    Product  EMA  FDA  EN  Other  WHO  Year  Generic  DrugBank ID  ATC  ChEMBL  Indications  Targets  Last Update
    Abemaciclib  Y  Y  N    N  2017  N  <a href=...>DB12001</a>  L01EF03  <a href=...>CHEMBL3301610</a>  Advanced Breast Cancer; ...  KRAS; ERBB2; CDK4; ...  2025-04-02

A drug row lists its target genes in one cell, separated by semicolons. The
module maps each drug-target pair: the drug as a chemical, the target as a
protein identified by its gene symbol. The table gives no mechanism, so the
pair is an association. Approvals and indications are not mapped.
"""

from __future__ import annotations

import re

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cancerdrugsdb import iter_drug_targets
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

# Namespace of the ATC codes; the registry has none yet.
ATC_NAMESPACE = 'atc'
# An ATC code of the fifth level names a substance. The shorter codes of the
# column ("L02BA") name a class of substances and identify no chemical.
_ATC_SUBSTANCE_REGEX = re.compile(r'^[A-Z]\d{2}[A-Z]{2}\d{2}$')

config = ResourceConfig(
    id=ResourceCv.CANCERDRUGSDB,
    name='CancerDrugsDB',
    url='https://data.tp53.org.uk/cancerdrugs.php',
    # Neither the site nor the publication states a license for the data; the
    # paper calls the database "open access".
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='33776770',
    topics=(Topic.PHARMACOLOGY,),
    content=(Content.COMPOUND_TARGET_INTERACTIONS,),
    # The targets come from DGIdb, which collects them from DrugBank, Drug
    # Target Commons and PharmGKB.
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Licensed cancer drugs with the target genes listed for them.',
    description=(
        'CancerDrugs_DB is a curated list of drugs licensed for the treatment '
        'of cancer, made by the Anticancer Fund from the data of the NCI, the '
        'FDA, the EMA and other agencies, and now maintained by the George '
        'Pantziarka TP53 Trust. Supportive-care and diagnostic agents and '
        'investigational drugs are not included. For each drug the table gives '
        'its approvals by agency and year, DrugBank, ATC and ChEMBL '
        'identifiers, indications from DrugBank and target genes from the '
        'Drug-Gene Interaction database (DGIdb).'
    ),
)

f = FieldConfig()


def _atc_substance_codes(row: dict) -> list[str]:
    return [code for code in row['atc_codes'] if _ATC_SUBSTANCE_REGEX.match(code)]


drug_builder = EntityBuilder(
    entity_type=model.ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.DRUGBANK, value=lambda row: row['drugbank_ids']),
        CV(term=Namespace.CHEMBL, value=lambda row: row['chembl_ids']),
        CV(term=ATC_NAMESPACE, value=_atc_substance_codes),
        CV(term=Namespace.NAME, value=f('product')),
    ),
)

# The table lists targets without saying how the drug acts on them, so the
# pair is an association, with no direction or mechanism. The file does not
# name the target's species, so the protein has no taxon.
drug_target_schema = RelationBuilder(
    subject=drug_builder,
    predicate=slots.associated_with,
    object=EntityBuilder(
        entity_type=model.Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('target')),
        ),
    ),
)

resource = Resource(
    config,
    drug_targets=Dataset(
        download=Download(
            url='https://data.tp53.org.uk/cancerdrugsdb.txt',
            filename='cancerdrugsdb.txt',
            subfolder='cancerdrugsdb',
            large=True,
        ),
        mapper=drug_target_schema,
        raw_parser=iter_drug_targets,
    ),
)
