"""
The BiGG model iMM1415 of mouse metabolism: metabolites, reactions and genes.

Metabolites keep their full BiGG cross-references and compartments. Reactions,
their participants and gene rules, AND groups of genes, and the model genes are
read as Recon3D's are, since both are BiGG models: BiGG reaction and metabolite
identifiers are shared by every BiGG model, so a reaction iMM1415 shares with
Recon3D is one entity, and what is iMM1415's own stays apart: its genes (mouse
Entrez identifiers, taxon 10090), its flux bounds and its datasets.
"""

from functools import partial
import json
from omnipath_resources.inputs_v2.base import ResourceConfig, Download, Resource, Dataset
from omnipath_resources.inputs_v2._source_context import compartment_name
from omnipath_resources.inputs_v2.parsers.recon3d import _raw as _bigg_raw
from omnipath_resources.inputs_v2.recon3d import (
    BiggModel,
    _enzyme_complexes_schema,
    _genes_schema,
    _reactions_schema,
)
from omnipath_resources.internals.tabular_builder import (
    CV,
    AssociationBuilder,
    AssociationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from omnipath_resources.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import ChemicalEntity, OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.reactions import go_compartment
from omnipath_core.resource_metadata import Content, Topic

URL = 'http://bigg.ucsd.edu/static/models/iMM1415.json'
config = ResourceConfig(
    id=ResourceCv.IMM1415,
    name='iMM1415',
    url='http://bigg.ucsd.edu/models/iMM1415',
    license=LicenseCV.BIGG,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='20959003',
    topics=(Topic.METABOLISM,),
    content=(Content.BIOCHEMICAL_REACTIONS,),
    summary='Mouse metabolic model metabolites, reactions with stoichiometry, and gene associations.',
    description='A high-quality mouse genome-scale metabolic reconstruction, iMM1415(Mus Musculus, 1415 genes)',
)


def parser(opener, **kwargs):
    def process_metabolite(metabolite):
        proc = {k: metabolite.get(k, '') for k in ('id', 'name')}
        proc['compartment'] = comps.get(metabolite.get('compartment'), '')
        proc['compartment_code'] = metabolite.get('compartment') or ''
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
    # No annotations: every metabolite has the SBO term of a simple chemical.
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
IMM1415 = BiggModel('iMM1415', '10090', 'imm1415.gpr_and')
reactions_schema = _reactions_schema(IMM1415)
genes_schema = _genes_schema(IMM1415)


def enzyme_complexes_schema(row: dict):
    """iMM1415's AND groups of genes."""
    return _enzyme_complexes_schema(row, IMM1415)


resource = Resource(
    config=config,
    metabolites=Dataset(download=download, mapper=schema, raw_parser=parser),
    reactions=Dataset(
        download=download,
        mapper=reactions_schema,
        raw_parser=partial(_bigg_raw, data_type='reactions'),
    ),
    enzyme_complexes=Dataset(
        download=download,
        mapper=enzyme_complexes_schema,
        raw_parser=partial(_bigg_raw, data_type='enzyme_complexes'),
    ),
    genes=Dataset(
        download=download,
        mapper=genes_schema,
        raw_parser=partial(_bigg_raw, data_type='genes'),
    ),
)
