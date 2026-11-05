"""iMM1415 mouse metabolic model: its own metabolites, plus the shared BiGG reaction content.

iMM1415 (Sigurdsson et al. 2010) is a mouse genome-scale metabolic
reconstruction, distributed as BiGG JSON exactly like Recon3D's. Its
``reactions``, ``enzyme_complexes`` and ``genes`` datasets come from the
shared factory in :mod:`omnipath_resources.inputs_v2._bigg`, with the mouse taxon and a
GPR identifier prefix of their own.

The ``metabolites`` dataset stays hand-written and unfiltered: it lists
every metabolite the model carries, not only the ones no reaction mentions
(``orphan_filter=False`` is therefore never passed to the shared factory —
this dataset is not built from it at all). The 16-namespace cross-reference
set and the per-metabolite ``located_in`` compartment relation are wider
than what a reaction participant row carries, and iMM1415 had no
`reactions` dataset at all until this module existed to filter against. A
metabolite that is also a reaction participant is therefore described
twice after this module's reaction content lands: once here, standalone,
and once per reaction, as a participant. That duplication is accepted
rather than resolved by a filter that would cost the wider identifier
coverage and the compartment relation to buy it back.
"""

import json

from biolink_model.datamodel.model import ChemicalEntity, OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.reactions import go_compartment
from omnipath_core.resource_metadata import Content, Topic

from omnipath_resources.inputs_v2._bigg import bigg_resource
from omnipath_resources.inputs_v2._source_context import compartment_name
from omnipath_resources.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from omnipath_resources.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from omnipath_resources.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    AssociationBuilder,
    AssociationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from omnipath_resources.resources import urls

MOUSE_TAXON_ID = 'NCBITaxon:10090'
GPR_PREFIX = 'imm1415.gpr_and'

config = ResourceConfig(
    id=ResourceCv.IMM1415,
    name='iMM1415',
    url='http://bigg.ucsd.edu/models/iMM1415',
    license=LicenseCV.BIGG,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='20959003',
    topics=(Topic.METABOLISM,),
    content=(Content.BIOCHEMICAL_REACTIONS,),
    summary='Mouse metabolic model entities and metabolic network information.',
    description='A high-quality mouse genome-scale metabolic reconstruction, iMM1415(Mus Musculus, 1415 genes)',
)


def parser(opener, **kwargs):
    def process_metabolite(metabolite):
        proc = {k: metabolite.get(k, '') for k in ('id', 'name')}
        proc['compartment'] = comps.get(metabolite.get('compartment'), '')
        proc['compartment_code'] = metabolite.get('compartment') or ''
        proc['formula'] = metabolite.get('formula') or ''
        proc['charge'] = metabolite.get('charge')
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
    url=urls.urls['imm1415']['bigg_json'],
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
    annotations=AnnotationsBuilder(
        CV(term='chemrof:charge', value=f('charge')),
        CV(term=slots.has_chemical_formula, value=f('formula')),
    ),
    # A metabolite of the model is located in its compartment. The model names
    # its compartments in its own way ("mitochondria", "peroxisome/glyoxysome");
    # the codes are those of Recon3D and Human-GEM, so they take the same GO terms.
    associations=AssociationsBuilder(
        AssociationBuilder(
            predicate=slots.located_in,
            object_entity_type=OntologyClass,
            object_identifier_type=Namespace.GO,
            object_identifier=lambda row: go_compartment(
                compartment_name(row.get('compartment_code'))
            ),
        ),
    ),
)

_bigg_datasets = bigg_resource(download=download, taxon=MOUSE_TAXON_ID, gpr_prefix=GPR_PREFIX)
reactions_schema = _bigg_datasets['reactions'].mapper
enzyme_complexes_schema = _bigg_datasets['enzyme_complexes'].mapper

resource = Resource(
    config=config,
    metabolites=Dataset(download=download, mapper=schema, raw_parser=parser),
    reactions=_bigg_datasets['reactions'],
    enzyme_complexes=_bigg_datasets['enzyme_complexes'],
    genes=_bigg_datasets['genes'],
)
