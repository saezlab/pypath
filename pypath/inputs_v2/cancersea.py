"""
CancerSEA gene sets of cancer functional states.

CancerSEA defines 14 functional states of cancer cells (stemness, invasion,
hypoxia...). For each state the authors collected a signature, a gene set,
from databases (HCMDB, Cyclebase, StemMapper, GO, MSigDB) and the literature
(the help page: "We manually collected 14 function state signatures from
databases [...] and published literatures"), and scored it in single-cell
datasets. The download page lists the states and links one text file of genes
per state: the signature. A gene of the file is ``associated_with`` the
state; the state is its own entity (the ``terms`` dataset), with the
description of the page.
"""

from __future__ import annotations

from functools import partial

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cancersea import iter_gene_sets, iter_states
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.CANCERSEA,
    name='CancerSEA',
    full='Cancer Single-cell State Atlas',
    url='http://biocc.hrbmu.edu.cn/CancerSEA/',
    # The footer of the site: "CancerSEA is covered by Creative Commons
    # Attribution-NonCommercial CC BY-NC License" (no version is named; the
    # legacy registry, resources.json, records CC BY-NC 4.0).
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.STATIC,
    # Yuan et al., "CancerSEA: a cancer single-cell state atlas", Nucleic Acids
    # Research 2019.
    pubmed='30329142',
    topics=(Topic.DISEASE_PHENOTYPE,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Gene sets (signatures) of 14 functional states of cancer cells, collected from databases and the literature.',
    description=(
        'CancerSEA is a database of the functional states of cancer cells at the '
        'single-cell level. It profiles 14 functional states (angiogenesis, '
        'apoptosis, cell cycle, differentiation, DNA damage, DNA repair, EMT, '
        'hypoxia, inflammation, invasion, metastasis, proliferation, quiescence, '
        'stemness) across single-cell RNA-seq datasets of cancers. The signature '
        'of each state, the gene set this module imports, was collected by the '
        'authors from databases and the literature; the state activities in the '
        'single-cell datasets are not imported.'
    ),
)

# The module reads OmniPath's rescued copy of the page and files, as the legacy
# module did. The live site (http://biocc.hrbmu.edu.cn/CancerSEA/) answers on
# plain http only (its https certificate has expired). Its 14 state files are
# byte-identical to the copy, and the table of states of its page is the same;
# the live page lists two more experiments (EXP0073, EXP0074) than the copy,
# which this module does not read (checked 2026-10-12).
BASE_URL = 'https://rescued.omnipathdb.org/'
# The file names of the states of the page. A state that is not here, or one
# that has gone, raises.
STATE_FILES = (
    'Angiogenesis', 'Apoptosis', 'Cell_Cycle', 'DNA_damage', 'DNA_repair', 'Differentiation',
    'EMT', 'Hypoxia', 'Inflammation', 'Invasion', 'Metastasis', 'Proliferation', 'Quiescence',
    'Stemness',
)
download_page = Download(
    url=f'{BASE_URL}CancerSEA.html',
    filename='CancerSEA-goDownload.html',
    subfolder='cancersea',
    large=False,
)
download_states = {
    name: Download(
        url=f'{BASE_URL}cancersea/{name}.txt',
        filename=f'CancerSEA-{name}.txt',
        subfolder='cancersea',
    )
    for name in STATE_FILES
}

# The states are CancerSEA's own vocabulary.
STATE_NAMESPACE = 'cancersea_state'

f = FieldConfig()

state_terms = EntityBuilder(
    entity_type=model.OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=STATE_NAMESPACE, value=f('state')),
        CV(term=Namespace.NAME, value=f('state')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.description, value=f('description'))),
)

# The gene: its Ensembl gene identifier and symbol, as the file gives them. The
# cells are of human cancers.
state_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=model.Gene,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.ENSG, value=f('ensembl_id')),
            CV(term=Namespace.GENESYMBOL, value=f('symbol')),
        ),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value='NCBITaxon:9606')),
    ),
    predicate=slots.associated_with,
    object=EntityBuilder(
        entity_type=model.OntologyClass,
        identifiers=IdentifiersBuilder(CV(term=STATE_NAMESPACE, value=f('state'))),
    ),
)

resource = Resource(
    config,
    terms=Dataset(
        download=download_page,
        mapper=state_terms,
        raw_parser=iter_states,
    ),
    genes=Dataset(
        download=download_page,
        mapper=state_schema,
        raw_parser=partial(iter_gene_sets, downloads=download_states),
    ),
)
