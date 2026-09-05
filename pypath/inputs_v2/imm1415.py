"""
Parse iMM1415 data and emit Entity records.

This module converts information from the BiGG Model iMM1415 into Entity
records using the declarative schema pattern.
"""

import json
from pypath.inputs_v2.base import ResourceConfig, Download, Resource, Dataset
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    ChemicalEntity,
    slots,
)
from omnipath_core.naming import Namespace

URL = 'http://bigg.ucsd.edu/static/models/iMM1415.json'
config = ResourceConfig(
    id=ResourceCv.IMM1415,
    name='iMM1415',
    url='http://bigg.ucsd.edu/models/iMM1415',
    license=LicenseCV.BIGG,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='20959003',
    primary_category='metabolism',
    description='A high-quality mouse genome-scale metabolic reconstruction, iMM1415(Mus Musculus, 1415 genes)',
)


def parser(opener, **kwargs):
    def process_metabolite(metabolite):
        proc = {k: metabolite.get(k, '') for k in ('id', 'name')}
        proc['compartment'] = comps.get(metabolite.get('compartment'), '')
        proc.update(
            {
                k: (metabolite.get('annotation') or {}).get(k, [])
                for k in annot_keys
                if k in (metabolite.get('annotation') or {}).keys()
            }
        )
        return proc

    result = json.loads(''.join(opener.result))
    comps = result['compartments']
    metabolites = result['metabolites']
    annot_keys = {
        k for m in metabolites for k in (m.get('annotation') or {}).keys()
    }
    yield from [process_metabolite(met) for met in metabolites]


download = Download(
    url=URL,
    filename='iMM1415.json',
    subfolder='bigg',
    large=True,
    ext='.json',
    default_mode='r',
)
f = FieldConfig(
    extract={'metacyc': 'META:(.*)', 'chebi': 'CHEBI:(.*)'},
    map={},
    transform={},
)
schema = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.BIGG_METABOLITE, value=f('bigg.metabolite')),
        CV(term=Namespace.METACYC, value=f('biocyc', extract='metacyc')),
        CV(term=Namespace.CHEBI, value=f('chebi', extract='chebi')),
        CV(term=Namespace.ENVIPATH, value=f('envipath')),
        CV(term=Namespace.HMDB, value=f('hmdb')),
        CV(term=Namespace.INCHIKEY, value=f('inchi_key')),
        CV(term=Namespace.KEGG, value=f('kegg.compound')),
        CV(term=Namespace.KEGG_DRUG, value=f('kegg.drug')),
        CV(term=Namespace.KEGG_GLYCAN, value=f('kegg.glycan')),
        CV(term=Namespace.LIPIDMAPS, value=f('lipidmaps')),
        CV(term=Namespace.METANETX, value=f('metanetx.chemical')),
        CV(term=Namespace.REACTOME, value=f('reactome.compound')),
        CV(term=Namespace.SABIORK_COMPOUND, value=f('sabiork')),
        CV(term=Namespace.SEED_COMPOUND, value=f('seed.compound')),
        CV(term=Namespace.SWISSLIPIDS, value=f('slm')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.has_topic, value=f('sbo'))),
)
resource = Resource(
    config=config,
    metabolites=Dataset(download=download, mapper=schema, raw_parser=parser),
)
