"""
Parse iMM1415 data and emit Entity records.

This module converts information from the BiGG Model iMM1415 into Entity
records using the declarative schema pattern.
"""

import json

from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    MoleculeAnnotationsCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

URL = 'http://bigg.ucsd.edu/static/models/iMM1415.json'

config = ResourceConfig(
    id=ResourceCv.IMM1415,
    name='iMM1415',
    url='http://bigg.ucsd.edu/models/iMM1415',
    license=LicenseCV.BIGG,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='20959003',
    primary_category='metabolism',
    description=(
        'A high-quality mouse genome-scale metabolic reconstruction, iMM1415'
        '(Mus Musculus, 1415 genes)'
    ),
)

# ================================== DOWNLOAD ==================================

def parser(opener, **kwargs):

    def process_metabolite(metabolite):

        proc = {k: metabolite.get(k, '') for k in ('id', 'name')}

        # Processing compartment into human-readable
        proc['compartment'] = comps.get(metabolite.get('compartment'), '')

        # Processing annotations
        proc.update({
            k: metabolite.get('annotation').get(k, [])
            for k in annot_keys
            if k in metabolite.get('annotation', {}).keys()
        })

        return proc


    result = json.loads(''.join(opener.result))

    comps = result['compartments']
    metabolites = result['metabolites']
    annot_keys = {k for m in metabolites for k in m['annotation'].keys()}

    yield from [process_metabolite(met) for met in metabolites]


download = Download(
    url=URL,
    filename='iMM1415.json',
    subfolder='bigg',
    large=True,
    ext='.json',
    default_mode='r',
)

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={
        'metacyc': r"META:(.*)",
        'chebi': r"CHEBI:(.*)",
    },
    map={},
    transform={},
)

schema = EntityBuilder(
    entity_type=EntityTypeCv.CHEMICAL,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.BIGG_METABOLITE,
            value=f('bigg.metabolite')
        ),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(
            term=IdentifierNamespaceCv.METACYC,
            value=f('biocyc', extract='metacyc')
        ),
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('chebi', extract='chebi')),
        CV(term=IdentifierNamespaceCv.ENVIPATH, value=f('envipath')),
        CV(term=IdentifierNamespaceCv.HMDB, value=f('hmdb')),
        CV(
            term=IdentifierNamespaceCv.STANDARD_INCHI_KEY,
            value=f('inchi_key')
        ),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('kegg.compound')),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('kegg.drug')),
        CV(term=IdentifierNamespaceCv.KEGG_GLYCAN, value=f('kegg.glycan')),
        CV(term=IdentifierNamespaceCv.LIPIDMAPS, value=f('lipidmaps')),
        CV(term=IdentifierNamespaceCv.METANETX, value=f('metanetx.chemical')),
        CV(term=IdentifierNamespaceCv.REACTOME_ID, value=f('reactome.compound')),
        CV(term=IdentifierNamespaceCv.SABIORK_COMPOUND, value=f('sabiork')),
        CV(term=IdentifierNamespaceCv.SEED_COMPOUND, value=f('seed.compound')),
        CV(term=IdentifierNamespaceCv.SWISSLIPIDS, value=f('slm')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=MoleculeAnnotationsCv.SUBCELLULAR_LOCATION,
            value=f('compartment')
        ),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('sbo')),
    ),
)

# ================================= RESOURCE ===================================

resource = Resource(
    config=config,
    metabolites=Dataset(
        download=download,
        mapper=schema,
        raw_parser=parser,
    ),
)

# ================================= REFERENCE ==================================
#{
#    'id': '5forthf_c',
#    'name': '5-Formiminotetrahydrofolate',
#    'compartment': 'cytosol',
#    'chebi': [
#        'CHEBI:18603',
#        'CHEBI:15639',
#        'CHEBI:12126',
#        'CHEBI:57456',
#        'CHEBI:2056'
#    ],
#    'inchi_key': ['YCWUVLPMLLBDCU-STQMWFEESA-L'],
#    'seed.compound': ['cpd00502'],
#    'bigg.metabolite': ['5forthf'],
#    'sabiork': ['1913'],
#    'sbo': 'SBO:0000247',
#    'reactome.compound': ['70919'],
#    'metanetx.chemical': ['MNXM915'],
#    'kegg.compound': ['C00664'],
#    'hmdb': ['HMDB01534'],
#    'biocyc': ['META:CPD-671']
#}
#
# All possible annotations/identifiers:
#
# 'id', # Basically [bigg.metabolite]_[compartment]
# X 'name',
# X 'compartment',
# X 'bigg.metabolite',
# X 'biocyc',
# X 'chebi',
# X 'envipath',
# X 'hmdb',
# X 'inchi_key',
# X 'kegg.compound',
# X 'kegg.drug',
# X 'kegg.glycan',
# X 'lipidmaps',
# X 'metanetx.chemical',
# X 'reactome.compound',
# X 'sabiork',
# X 'sbo',
# X 'seed.compound',
# X 'slm'
