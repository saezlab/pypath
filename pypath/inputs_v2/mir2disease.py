"""
Parse the miR2Disease validated miRNA - target pairs and emit Relation records.

miR2Disease is a database of microRNA deregulation in human diseases. Its
download of experimentally verified pairs says that they were "validated by
luciferase reporter experiments" and "extracted manually from the corresponding
reference". Each line names a human miRNA, the target and the paper.

The miRNA ``affects`` the target. The table states no direction and no aspect
of that effect, so there is none in the record. The target is named as the
paper names it (``p53``, ``caspase-3``, ``KRAS``), which is a gene or protein
name and not always a gene symbol, so it is a name identifier. The table has no
organism column: the organism code of the miRNA name gives the miRNA's taxon and
the target has none. The miRNA - disease table of the database is not imported.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, MicroRNA, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.mir2disease import iter_targets
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

config = ResourceConfig(
    id=ResourceCv.MIR2DISEASE,
    name='miR2Disease',
    url='http://watson.compbio.iupui.edu:8080/miR2Disease/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_NC_2_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='18927107',
    topics=(Topic.GENE_REGULATION, Topic.DISEASE_PHENOTYPE),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Human microRNA - target pairs validated by luciferase reporter '
        'experiments, curated from the literature.'
    ),
    description=(
        'miR2Disease is a manually curated database of microRNA deregulation '
        'in human diseases, built at the Harbin Institute of Technology and '
        'Indiana University. Each entry links a miRNA to a disease with its '
        'expression pattern, the detection method, experimentally verified '
        'target genes and the reference. This module imports the download '
        'of the verified miRNA - target pairs, which the database says were '
        'validated by luciferase reporter experiments and extracted manually '
        'from the references. The miRNA - disease entries are not imported.'
    ),
)

# Organism code of the miRNA name -> NCBI Taxonomy. All 805 names of the table
# start with ``hsa-``.
TAXA = {
    'hsa': '9606',  # Homo sapiens
}


def _mirna_taxon(row: dict[str, Any]) -> str:
    code = row.get('mirna_species')
    if code not in TAXA:
        raise ValueError(f'Unknown miR2Disease organism code: {code!r}')
    return f'NCBITaxon:{TAXA[code]}'


f = FieldConfig()

# The name is kept as the source writes it; it is not translated between miRBase
# versions and the -3p and -5p forms stay separate.
mirna = EntityBuilder(
    entity_type=MicroRNA,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('mirna')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=_mirna_taxon),
    ),
)

# The source names the target by a gene or protein name of the paper, with no
# accession and no organism.
target = EntityBuilder(
    entity_type=Gene,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('target')),
    ),
)

# A row without a target name yields no relation: the target has no identifier.
interactions_schema = RelationBuilder(
    subject=mirna,
    predicate=slots.affects,
    object=target,
)

resource = Resource(
    config,
    mirna_targets=Dataset(
        download=Download(
            url='http://watson.compbio.iupui.edu:8080/miR2Disease/download/miRtar.txt',
            filename='miRtar.txt',
            subfolder='mir2disease',
            large=True,
            encoding='iso-8859-1',
        ),
        mapper=interactions_schema,
        raw_parser=iter_targets,
    ),
)
