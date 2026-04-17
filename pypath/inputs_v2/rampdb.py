"""
Parse RaMP-DB data and emit Entity records.

This module converts annotations of lipids and metabolites into Entity records
using the declarative schema pattern.
"""

import os
import re
import requests
from functools import partial
from pathlib import Path

from bs4 import BeautifulSoup

from pypath.inputs_v2.parsers.base import iter_sqlite
from pypath.inputs_v2.base import ResourceConfig, Download, Resource, Dataset
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeSubtypeCv,
    MoleculeAnnotationsCv,
    AssayAnnotationsCv,
)

# =================================== SET-UP ===================================

BASE_URL = 'https://github.com/ncats/RaMP-DB/raw/refs/heads/%s/db/'

def get_ramp_latest_ver(branch='main'):
    '''
    Retrieves the URL for the latest version of the database SQL file
    '''

    res = requests.get(BASE_URL % branch)
    soup = BeautifulSoup(res.text, 'html.parser')
    files = sorted({
        f.text for f in soup.find_all(title=re.compile("\\.sqlite.gz$"))
    })

    ver = re.search(r'(\d+\.\d+\.\d+)', files[-1]).group(1)

    return ver, BASE_URL % branch + files[-1]

VERSION, URL = get_ramp_latest_ver()

config = ResourceConfig(
    id=ResourceCv.RAMPDB,
    name='RaMP-DB',
    url='https://rampdb.nh.gov/',
    license=LicenseCV.GPL_2_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='36373969',
    primary_category='pathways',
    description=(
        'RaMP-DB (Relational database of Metabolomic Pathways) is a '
        'multi-sourced integrated database with comprehensive annotations on '
        'biological pathways, structure/chemistry, disease and ontology '
        'annotations for genes, proteins, and metabolites.'
    )
)

# ================================== DOWNLOAD ==================================

download = Download(
    url=URL,
    filename=os.path.basename(URL),
    subfolder='ramp',
    large=False,
    ext=URL.split('.')[-1],
    default_mode='rb',
)

table_names = [
    'analyte',
    #'entity_status_info',
    #'reaction2met',
    #'analytehasontology',
    'metabolite_class',
    #'reaction2protein',
    #'analytehaspathway',
    #'ontology',
    #'reaction_ec_class',
    #'analytesynonym',
    #'pathway',
    #'reaction_protein2met',
    #'catalyzed',
    #'pathway_duplicates',
    'source',
    'chem_props',
    #'pathway_similarity',
    #'version_info',
    #'db_version',
    #'reaction',
]

#opener = download.open()
#path = Path(opener.path.replace('.gz', ''))

#with open(path, 'wb') as f:
#    f.write(opener.result)

#datasets = dict()

#for tbl in table_names:

#    datasets[tbl] = iter_sqlite(opener, table_name=tbl, sqlite_path=path)

def parser(opener, table=None):

    path = Path(opener.path.replace('.gz', ''))

    yield from iter_sqlite(opener, table_name=table, sqlite_path=path)

# =================================== SCHEMA ===================================

MOLECULE_TYPE_TO_ENTITY_TYPE = {
    'gene': EntityTypeCv.GENE,
    'compound': MoleculeSubtypeCv.METABOLITE,
}
SOURCE_TO_ENTITY_TYPE = {
    'hmdb': MoleculeSubtypeCv.METABOLITE,
    'lipidmaps': EntityTypeCv.LIPID,
}
SOURCE_TO_TERM = {
    'hmdb': IdentifierNamespaceCv.HMDB,
    'lipidmaps': IdentifierNamespaceCv.LIPIDMAPS,
    'cas': IdentifierNamespaceCv.CAS,
    'en':  IdentifierNamespaceCv.EC,
    'brenda': IdentifierNamespaceCv.EC,
    'chebi': IdentifierNamespaceCv.CHEBI,
    'chemspider': IdentifierNamespaceCv.CHEMSPIDER,
    'ensembl': IdentifierNamespaceCv.ENSEMBL,
    'entrez': IdentifierNamespaceCv.ENTREZ,
    'gene_symbol': IdentifierNamespaceCv.GENE_NAME_PRIMARY,
}
CLASS_LEVEL_NAME_TO_TERM ={
    'LipidMaps_category': MoleculeAnnotationsCv.LIPID_CATEGORY,
    'LipidMaps_main_class': MoleculeAnnotationsCv.LIPID_MAIN_CLASS,
    'LipidMaps_sub_class': MoleculeAnnotationsCv.LIPID_SUB_CLASS,
    'ClassyFire_class': MoleculeAnnotationsCv.COMPOUND_CLASS,
    'ClassyFire_sub_class': MoleculeAnnotationsCv.COMPOUND_SUBCLASS,
    'ClassyFire_super_class': MoleculeAnnotationsCv.COMPOUND_SUPERCLASS,
}

f = FieldConfig(
    extract={
        'rampID':  r'^(RAMP_[A-Z]+_\d+)$',
        #'hmdbID': r'^(?:hmdb:)?(HMDB\d+)$',
        #'lipidmapsID': r'^(?:LIPIDMAPS:)?(LMSP\d+)$',
        #'sourceID': r'^([a-zA-Z]*):(?:[a-zA-Z]+\d+)$'
    },
    map={
        'type_to_entity': MOLECULE_TYPE_TO_ENTITY_TYPE,
        'source_to_entity': SOURCE_TO_ENTITY_TYPE,
        'source_to_term': SOURCE_TO_TERM,
        'class_level_to_term': CLASS_LEVEL_NAME_TO_TERM,
    },
    transform={
        'lower': lambda x: x.lower(),
        'upper': lambda x: x.upper(),
        'postcolon': lambda x: x.split(':')[-1]
    }
)

# X = Info added in some way
# ? = No clue what is this
# - = Skipped

# analyte
# X  0|rampId|VARCHAR(30)|1||1
# X  1|type|VARCHAR(30)|0||0
# X  2|common_name|TEXT|0||0
#   e.g. RAMP_C_000225217|compound|(10S,14S,16R)-10-F4-NeuroP[13R,17R]

analyte_schema = EntityBuilder(
    entity_type=f('type', map='type_to_entity'),
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.RAMP_ID,
            value=f('rampId', extract='rampID')
        ),
        CV(term=IdentifierNamespaceCv.SYSTEMATIC_NAME, value=f('common_name')),
    ),
)

# metabolite_class
# X  0|ramp_id|VARCHAR(32)|1||2
# X  1|class_source_id|VARCHAR(32)|1||1
# X  2|class_level_name|VARCHAR(128)|1||3
# X  3|class_name|VARCHAR(128)|1||0
# X  4|source|VARCHAR(32)|1||0
#   e.g. RAMP_C_000000001|hmdb:HMDB0000001|ClassyFire_super_class|Organic acids and derivatives|hmdb

metabolite_class_schema = EntityBuilder(
    entity_type=f('source', map='source_to_entity'),
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.RAMP_ID,
            value=f('ramp_id', extract='rampID')
        ),
        CV(
            term=f('source', map='source_to_term', transform='lower'),
            value=f('class_source_id', transform='postcolon'),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=f('class_level_name', map='class_level_to_term'),
            value=f('class_name')
        ),
    ),
)

# source
# X  0|sourceId|VARCHAR(30)|1||0
# X  1|rampId|VARCHAR(30)|0||0
# X  2|IDtype|VARCHAR(30)|0||0
# X  3|geneOrCompound|VARCHAR(30)|0||0
# X  4|commonName|VARCHAR(256)|0||0
# X  5|priorityHMDBStatus|VARCHAR(32)|0||0
# -  6|dataSource|VARCHAR(32)|0||0
# ?  7|pathwayCount|INTEGER|1|'0'|0
#   e.g. LIPIDMAPS:LMFA04010056|RAMP_C_000225217|LIPIDMAPS|compound|(10S,14S,16R)-10-F4-NeuroP[13R,17R]|no_HMDB_status|lipidmaps|0

source_schema = EntityBuilder(
    entity_type=f('geneOrCompound', map='type_to_entity'),
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.RAMP_ID,
            value=f('rampId', extract='rampID')
        ),
        CV(
            term=f('IDtype', map='source_to_term', transform='lower'),
            value=f('sourceId', transform='postcolon'),
        ),
        CV(term=IdentifierNamespaceCv.NAME, value=f('commonName')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=AssayAnnotationsCv.ASSAY_CATEGORY,
            value=f('priorityHMDBStatus')
        ),
    ),
)

# chem_props
# X  0|ramp_id|VARCHAR(30)|1||0
# X  1|chem_data_source|VARCHAR(32)|0||0
# X  2|chem_source_id|VARCHAR(45)|0||0
# X  3|iso_smiles|VARCHAR(256)|0||0
# -  4|inchi_key_prefix|VARCHAR(32)|0||0
# X  5|inchi_key|VARCHAR(32)|0||0
# X  6|inchi|VARCHAR(4096)|0||0
# X  7|mw|FLOAT|0||0
# X  8|monoisotop_mass|FLOAT|0||0
# X  9|common_name|VARCHAR(1024)|0||0
# X  10|mol_formula|VARCHAR(64)|0||0
#   e.g. RAMP_C_000000001|hmdb|hmdb:HMDB0000001|[H]OC(=O)[C@@]([H])(N([H])[H])C([H])([H])C1=C([H])N(C([H])=N1)C([H])([H])[H]|BRMWTNUJHUMWMS|BRMWTNUJHUMWMS-LURJTMIESA-N|InChI=1S/C7H11N3O2/c1-10-3-5(9-4-10)2-6(8)7(11)12/h3-4,6H,2,8H2,1H3,(H,11,12)/t6-/m0/s1|169.1811|169.085126611|1-Methylhistidine|C7H11N3O2

chem_props_schema = EntityBuilder(
    entity_type=EntityTypeCv.SMALL_MOLECULE,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.RAMP_ID,
            value=f('ramp_id', extract='rampID')
        ),
        CV(
            term=f('chem_data_source', map='source_to_term', transform='lower'),
            value=f('chem_source_id', transform='postcolon'),
        ),
        CV(term=IdentifierNamespaceCv.SMILES, value=f('iso_smiles')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('inchi_key')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=f('inchi')),
        CV(
            term=IdentifierNamespaceCv.MOLECULAR_FORMULA,
            value=f('mol_formula')
        ),
        CV(term=IdentifierNamespaceCv.NAME, value=f('common_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('mw')),
        CV(
            term=MoleculeAnnotationsCv.MW_MONOISOTOPIC,
            value=f('monoisotop_mass')
        ),
    )
)

# ================================= RESOURCE ===================================

kwargs = {
    t: Dataset(
        download=download,
        mapper=locals().get('%s_schema' % t),
        raw_parser=partial(parser, table=t)
    )
    for t in table_names
}

resource = Resource(config=config, **kwargs)

# ================================= REFERENCE ==================================

# SQL tables and content:

# analyte
#   0|rampId|VARCHAR(30)|1||1
#   1|type|VARCHAR(30)|0||0
#   2|common_name|TEXT|0||0
#   e.g. RAMP_C_000225217|compound|(10S,14S,16R)-10-F4-NeuroP[13R,17R]

# entity_status_info
#   0|status_category|VARCHAR(64)|1||0
#   1|entity_source_id|VARCHAR(32)|1||0
#   2|entity_source_name|VARCHAR(45)|1||0
#   3|entity_count|INTEGER|1||0
#   e.g. Metabolite-Reaction Associations|rhea|Rhea|325728

# reaction2met
#   0|ramp_rxn_id|VARCHAR(16)|1||0
#   1|rxn_source_id|VARCHAR(16)|1||0
#   2|ramp_cmpd_id|VARCHAR(16)|1||0
#   3|substrate_product|INTEGER|1||0
#   4|met_source_id|VARCHAR(32)|1||0
#   5|met_name|VARCHAR(256)|0||0
#   6|is_cofactor|INTEGER|1|'0'|0
#   e.g. RAMP_R_000000001|rhea:10000|RAMP_C_000000802|0|chebi:15377|H2O|0

# analytehasontology
#   0|rampCompoundId|VARCHAR(30)|0||1
#   1|rampOntologyId|VARCHAR(30)|0||2
#   e.g. RAMP_C_000000001|RAMP_OL_000000001

# metabolite_class
#   0|ramp_id|VARCHAR(32)|1||2
#   1|class_source_id|VARCHAR(32)|1||1
#   2|class_level_name|VARCHAR(128)|1||3
#   3|class_name|VARCHAR(128)|1||0
#   4|source|VARCHAR(32)|1||0
#   e.g. RAMP_C_000000001|hmdb:HMDB0000001|ClassyFire_super_class|Organic acids and derivatives|hmdb

# reaction2protein
#   0|ramp_rxn_id|VARCHAR(16)|1||0
#   1|rxn_source_id|VARCHAR(16)|1||0
#   2|ramp_gene_id|VARCHAR(16)|1||0
#   3|uniprot|VARCHAR(16)|1||0
#   4|protein_name|VARCHAR(16)|1||0
#   5|is_reviewed|INTEGER|1|0|0
#   e.g. RAMP_R_000000007|rhea:10024|RAMP_G_000013905|uniprot:Q9Y5N5|N6AMT1|1

# analytehaspathway
#   0|rampId|VARCHAR(30)|0||0
#   1|pathwayRampId|VARCHAR(30)|0||0
#   2|pathwaySource|VARCHAR(30)|0||0
#   e.g. RAMP_C_000000001|RAMP_P_000000001|hmdb

# ontology
#   0|rampOntologyId|varchar(30)|0||1
#   1|commonName|varchar(64)|0||0
#   2|HMDBOntologyType|varchar(30)|0||0
#   3|metCount|INT|0||0
#   e.g. RAMP_OL_000000001|Saliva|Biofluid and excreta|1216

# reaction_ec_class
#   0|ramp_rxn_id|varchar(16)|1||0
#   1|rxn_source_id|varchar(16)|1||0
#   2|rxn_class_ec|varchar(16)|1||0
#   3|ec_level|INT|1||0
#   4|rxn_class|varchar(256)|1||0
#   5|rxn_class_hierarchy|varchar(512)|1||0
#   e.g. RAMP_R_000000001|rhea:10000|3.-.-.-|1|Hydrolases.|Hydrolases.

# analytesynonym
#   0|Synonym|varchar(500)|0||0
#   1|rampId|varchar(30)|0||0
#   2|geneOrCompound|varchar(30)|0||0
#   3|source|varchar(30)|0||0
#   e.g. 1-Methylhistidine|RAMP_C_000000001|compound|hmdb

# pathway
#   0|pathwayRampId|varchar(30)|0||1
#   1|sourceId|varchar(30)|0||0
#   2|type|varchar(30)|0||0
#   3|pathwayCategory|varchar(30)|0||0
#   4|pathwayName|varchar(250)|0||0
#   e.g. RAMP_P_000000001|SMP0124716|hmdb|smpdb3|1-Methylhistidine Metabolism

# reaction_protein2met
#   0|ramp_rxn_id|VARCHAR(16)|1||0
#   1|rxn_source_id|VARCHAR(16)|1||0
#   2|ramp_gene_id|VARCHAR(16)|1||0
#   3|gene_source_id|VARCHAR(16)|1||0
#   4|substrate_product|INTEGER|1||0
#   5|ramp_cmpd_id|VARCHAR(16)|1||0
#   6|cmpd_source_id|VARCHAR(45)|1||0
#   7|cmpd_name|VARCHAR(256)|0||0
#   8|is_cofactor|INTEGER|1|'0'|0
#   e.g RAMP_R_000000007|rhea:10024|RAMP_G_000013905|uniprot:Q9Y5N5|0|RAMP_C_000261687|rhea-comp:9845|L-lysyl-[histone]|1

# catalyzed
#   0|rampCompoundId|VARCHAR(30)|0||1
#   1|rampGeneId|VARCHAR(30)|0||2
#   2|proteinType|VARCHAR(32)|0||0
#   e.g. RAMP_C_000000001|RAMP_G_000000001|Unknown

# pathway_duplicates
#   0|pathwayRampId1|varchar(30)|1||0
#   1|pathwayRampId2|varchar(30)|1||0
#   e.g. RAMP_P_000050053|RAMP_P_000050055

# source
#   0|sourceId|VARCHAR(30)|1||0
#   1|rampId|VARCHAR(30)|0||0
#   2|IDtype|VARCHAR(30)|0||0
#   3|geneOrCompound|VARCHAR(30)|0||0
#   4|commonName|VARCHAR(256)|0||0
#   5|priorityHMDBStatus|VARCHAR(32)|0||0
#   6|dataSource|VARCHAR(32)|0||0
#   7|pathwayCount|INTEGER|1|'0'|0
#   e.g. LIPIDMAPS:LMFA04010056|RAMP_C_000225217|LIPIDMAPS|compound|(10S,14S,16R)-10-F4-NeuroP[13R,17R]|no_HMDB_status|lipidmaps|0

# chem_props
#   0|ramp_id|VARCHAR(30)|1||0
#   1|chem_data_source|VARCHAR(32)|0||0
#   2|chem_source_id|VARCHAR(45)|0||0
#   3|iso_smiles|VARCHAR(256)|0||0
#   4|inchi_key_prefix|VARCHAR(32)|0||0
#   5|inchi_key|VARCHAR(32)|0||0
#   6|inchi|VARCHAR(4096)|0||0
#   7|mw|FLOAT|0||0
#   8|monoisotop_mass|FLOAT|0||0
#   9|common_name|VARCHAR(1024)|0||0
#   10|mol_formula|VARCHAR(64)|0||0
#   e.g. RAMP_C_000000001|hmdb|hmdb:HMDB0000001|[H]OC(=O)[C@@]([H])(N([H])[H])C([H])([H])C1=C([H])N(C([H])=N1)C([H])([H])[H]|BRMWTNUJHUMWMS|BRMWTNUJHUMWMS-LURJTMIESA-N|InChI=1S/C7H11N3O2/c1-10-3-5(9-4-10)2-6(8)7(11)12/h3-4,6H,2,8H2,1H3,(H,11,12)/t6-/m0/s1|169.1811|169.085126611|1-Methylhistidine|C7H11N3O2

# pathway_similarity
#   0|pathwayRampId|varchar(30)|1||1
#   1|analyte_blob|BLOB|0||0
#   2|metabolite_blob|BLOB|0||0
#   3|gene_blob|BLOB|0||0
#   4|metabolite_count|INTEGER|0|0|0
#   5|gene_count|INTEGER|0|0|0
#   e.g. RAMP_P_000000003|[blob]|[blob]|[blob]|35|33


# version_info
#   0|ramp_db_version|VARCHAR(16)|1||0
#   1|db_mod_date|DATE|1||0
#   2|status|VARCHAR(16)|1||0
#   3|data_source_id|VARCHAR(32)|1||0
#   4|data_source_name|VARCHAR(128)|1||0
#   5|data_source_url|VARCHAR(128)|1||0
#   6|data_source_version|VARCHAR(128)|1||0
#   e.g. 3.0.7|2025-03-06 00:00:39.153317|current|hmdb|HMDB|https://hmdb.ca/|v5.0 (2021-11-17)

# db_version
#   0|ramp_version|VARCHAR(20)|1||0
#   1|load_timestamp|DATETIME|1|'CURRENT_TIMESTAMP'|0
#   2|version_notes|VARCHAR(256)|0||0
#   3|met_intersects_json|VARCHAR(10000)|0||0
#   4|gene_intersects_json|VARCHAR(10000)|0||0
#   5|met_intersects_json_pw_mapped|VARCHAR(10000)|0||0
#   6|gene_intersects_json_pw_mapped|VARCHAR(10000)|0||0
#   7|db_sql_url|VARCHAR(256)|0||0
#   e.g. 3.0.7|2025-03-06 00:00:39.153317|20240822 data update, new metabolite data from refmet, new datasource for pathways from PFOCR, updated MW check for mismerged metabolites, added field for best analyte name|[{"id": "cmpd_src_set_1", "sets": ["WikiPathways"], "size": 655}, {"id": "cmpd_src_set_2", "sets": ["LIPIDMAPS"], "size": 36815}, {"id": "cmpd_src_set_3", "sets": ["Reactome"], "size": 598}, {"id": "cmpd_src_set_4", "sets": ["RefMet"], "size": 176740}, {"id": "cmpd_src_set_5", "sets": ["HMDB"], "size": 184171}, {"id": "cmpd_src_set_6", "sets": ["Rhea"], "size": 8505}, {"id": "cmpd_src_set_7", "sets": ["PFOCR"], "size": 1915}, {"id": "cmpd_src_set_8", "sets": ["KEGG", "HMDB"], "size": 1270}, {"id": "cmpd_src_set_9", "sets": ["WikiPathways", "LIPIDMAPS"], "size": 198}, {"id": "cmpd_src_set_10", "sets": ["WikiPathways", "Reactome"], "size": 54}, {"id": "cmpd_src_set_11", "sets": ["WikiPathways", "RefMet"], "size": 7}, {"id": "cmpd_src_set_12", "sets": ["WikiPathways", "HMDB"], "size": 126}, {"id": "cmpd_src_set_13", "sets": ["WikiPathways", "Rhea"], "size": 74}, {"id": "cmpd_src_set_14", "sets": ["WikiPathways", "KEGG"], "size": 45}, {"id": "cmpd_src_set_15", "sets": ["WikiPathways", "PFOCR"], "size": 152}, {"id": "cmpd_src_set_16", "sets": ["LIPIDMAPS", "Reactome"], "size": 20}, {"id": "cmpd_src_set_17", "sets": ["LIPIDMAPS", "RefMet"], "size": 1625}, {"id": "cmpd_src_set_18", "sets": ["LIPIDMAPS", "HMDB"], "size": 4607}, {"id": "cmpd_src_set_19", "sets": ["LIPIDMAPS", "Rhea"], "size": 479}, {"id": "cmpd_src_set_20", "sets": ["LIPIDMAPS", "PFOCR"], "size": 77}, {"id": "cmpd_src_set_21", "sets": ["Reactome", "RefMet"], "size": 8}, {"id": "cmpd_src_set_22", "sets": ["Reactome", "HMDB"], "size": 77}, {"id": "cmpd_src_set_23", "sets": ["Reactome", "Rhea"], "size": 149}, {"id": "cmpd_src_set_24", "sets": ["Reactome", "PFOCR"], "size": 78}, {"id": "cmpd_src_set_25", "sets": ["RefMet", "HMDB"], "size": 1620}, {"id": "cmpd_src_set_26", "sets": ["RefMet", "Rhea"], "size": 401}, {"id": "cmpd_src_set_27", "sets": ["RefMet", "PFOCR"], "size": 126}, {"id": "cmpd_src_set_28", "sets": ["HMDB", "Rhea"], "size": 175}, {"id": "cmpd_src_set_29", "sets": ["HMDB", "PFOCR"], "size": 302}, {"id": "cmpd_src_set_30", "sets": ["Rhea", "PFOCR"], "size": 317}, {"id": "cmpd_src_set_31", "sets": ["KEGG", "WikiPathways", "HMDB"], "size": 47}, {"id": "cmpd_src_set_32", "sets": ["KEGG", "LIPIDMAPS", "HMDB"], "size": 170}, {"id": "cmpd_src_set_33", "sets": ["KEGG", "Reactome", "HMDB"], "size": 21}, {"id": "cmpd_src_set_34", "sets": ["KEGG", "RefMet", "HMDB"], "size": 884}, {"id": "cmpd_src_set_35", "sets": ["KEGG", "HMDB", "Rhea"], "size": 159}, {"id": "cmpd_src_set_36", "sets": ["KEGG", "HMDB", "PFOCR"], "size": 133}, {"id": "cmpd_src_set_37", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome"], "size": 5}, {"id": "cmpd_src_set_38", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet"], "size": 33}, {"id": "cmpd_src_set_39", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB"], "size": 21}, {"id": "cmpd_src_set_40", "sets": ["WikiPathways", "LIPIDMAPS", "Rhea"], "size": 44}, {"id": "cmpd_src_set_41", "sets": ["WikiPathways", "LIPIDMAPS", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_42", "sets": ["WikiPathways", "Reactome", "HMDB"], "size": 9}, {"id": "cmpd_src_set_43", "sets": ["WikiPathways", "Reactome", "Rhea"], "size": 48}, {"id": "cmpd_src_set_44", "sets": ["WikiPathways", "Reactome", "PFOCR"], "size": 38}, {"id": "cmpd_src_set_45", "sets": ["WikiPathways", "RefMet", "HMDB"], "size": 50}, {"id": "cmpd_src_set_46", "sets": ["WikiPathways", "RefMet", "Rhea"], "size": 2}, {"id": "cmpd_src_set_47", "sets": ["WikiPathways", "RefMet", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_48", "sets": ["WikiPathways", "HMDB", "Rhea"], "size": 5}, {"id": "cmpd_src_set_49", "sets": ["WikiPathways", "HMDB", "PFOCR"], "size": 32}, {"id": "cmpd_src_set_50", "sets": ["WikiPathways", "Rhea", "PFOCR"], "size": 19}, {"id": "cmpd_src_set_51", "sets": ["WikiPathways", "KEGG", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_52", "sets": ["LIPIDMAPS", "Reactome", "RefMet"], "size": 12}, {"id": "cmpd_src_set_53", "sets": ["LIPIDMAPS", "Reactome", "HMDB"], "size": 3}, {"id": "cmpd_src_set_54", "sets": ["LIPIDMAPS", "Reactome", "Rhea"], "size": 19}, {"id": "cmpd_src_set_55", "sets": ["LIPIDMAPS", "Reactome", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_56", "sets": ["LIPIDMAPS", "RefMet", "HMDB"], "size": 1003}, {"id": "cmpd_src_set_57", "sets": ["LIPIDMAPS", "RefMet", "Rhea"], "size": 219}, {"id": "cmpd_src_set_58", "sets": ["LIPIDMAPS", "RefMet", "PFOCR"], "size": 20}, {"id": "cmpd_src_set_59", "sets": ["LIPIDMAPS", "HMDB", "Rhea"], "size": 82}, {"id": "cmpd_src_set_60", "sets": ["LIPIDMAPS", "HMDB", "PFOCR"], "size": 13}, {"id": "cmpd_src_set_61", "sets": ["LIPIDMAPS", "Rhea", "PFOCR"], "size": 55}, {"id": "cmpd_src_set_62", "sets": ["Reactome", "RefMet", "HMDB"], "size": 29}, {"id": "cmpd_src_set_63", "sets": ["Reactome", "RefMet", "Rhea"], "size": 6}, {"id": "cmpd_src_set_64", "sets": ["Reactome", "RefMet", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_65", "sets": ["Reactome", "HMDB", "Rhea"], "size": 4}, {"id": "cmpd_src_set_66", "sets": ["Reactome", "HMDB", "PFOCR"], "size": 27}, {"id": "cmpd_src_set_67", "sets": ["Reactome", "Rhea", "PFOCR"], "size": 24}, {"id": "cmpd_src_set_68", "sets": ["RefMet", "HMDB", "Rhea"], "size": 167}, {"id": "cmpd_src_set_69", "sets": ["RefMet", "HMDB", "PFOCR"], "size": 171}, {"id": "cmpd_src_set_70", "sets": ["RefMet", "Rhea", "PFOCR"], "size": 75}, {"id": "cmpd_src_set_71", "sets": ["HMDB", "Rhea", "PFOCR"], "size": 28}, {"id": "cmpd_src_set_72", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "HMDB"], "size": 8}, {"id": "cmpd_src_set_73", "sets": ["KEGG", "WikiPathways", "Reactome", "HMDB"], "size": 7}, {"id": "cmpd_src_set_74", "sets": ["KEGG", "WikiPathways", "RefMet", "HMDB"], "size": 52}, {"id": "cmpd_src_set_75", "sets": ["KEGG", "WikiPathways", "HMDB", "Rhea"], "size": 11}, {"id": "cmpd_src_set_76", "sets": ["WikiPathways", "HMDB", "KEGG"], "size": 2}, {"id": "cmpd_src_set_77", "sets": ["KEGG", "WikiPathways", "HMDB", "PFOCR"], "size": 31}, {"id": "cmpd_src_set_78", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "HMDB"], "size": 4}, {"id": "cmpd_src_set_79", "sets": ["KEGG", "LIPIDMAPS", "RefMet", "HMDB"], "size": 179}, {"id": "cmpd_src_set_80", "sets": ["KEGG", "LIPIDMAPS", "HMDB", "Rhea"], "size": 68}, {"id": "cmpd_src_set_81", "sets": ["KEGG", "LIPIDMAPS", "HMDB", "PFOCR"], "size": 20}, {"id": "cmpd_src_set_82", "sets": ["KEGG", "Reactome", "RefMet", "HMDB"], "size": 52}, {"id": "cmpd_src_set_83", "sets": ["KEGG", "Reactome", "HMDB", "Rhea"], "size": 13}, {"id": "cmpd_src_set_84", "sets": ["KEGG", "Reactome", "HMDB", "PFOCR"], "size": 16}, {"id": "cmpd_src_set_85", "sets": ["KEGG", "RefMet", "HMDB", "Rhea"], "size": 319}, {"id": "cmpd_src_set_86", "sets": ["KEGG", "RefMet", "HMDB", "PFOCR"], "size": 316}, {"id": "cmpd_src_set_87", "sets": ["KEGG", "HMDB", "Rhea", "PFOCR"], "size": 53}, {"id": "cmpd_src_set_88", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet"], "size": 7}, {"id": "cmpd_src_set_89", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "HMDB"], "size": 3}, {"id": "cmpd_src_set_90", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "Rhea"], "size": 5}, {"id": "cmpd_src_set_91", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB"], "size": 78}, {"id": "cmpd_src_set_92", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "Rhea"], "size": 8}, {"id": "cmpd_src_set_93", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "KEGG"], "size": 1}, {"id": "cmpd_src_set_94", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_95", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB", "Rhea"], "size": 6}, {"id": "cmpd_src_set_96", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_97", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB"], "size": 4}, {"id": "cmpd_src_set_98", "sets": ["WikiPathways", "Reactome", "RefMet", "Rhea"], "size": 2}, {"id": "cmpd_src_set_99", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea"], "size": 6}, {"id": "cmpd_src_set_100", "sets": ["WikiPathways", "Reactome", "HMDB", "PFOCR"], "size": 10}, {"id": "cmpd_src_set_101", "sets": ["WikiPathways", "Reactome", "Rhea", "KEGG"], "size": 1}, {"id": "cmpd_src_set_102", "sets": ["WikiPathways", "Reactome", "Rhea", "PFOCR"], "size": 27}, {"id": "cmpd_src_set_103", "sets": ["WikiPathways", "Reactome", "KEGG", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_104", "sets": ["WikiPathways", "RefMet", "HMDB", "Rhea"], "size": 10}, {"id": "cmpd_src_set_105", "sets": ["WikiPathways", "RefMet", "HMDB", "PFOCR"], "size": 28}, {"id": "cmpd_src_set_106", "sets": ["WikiPathways", "RefMet", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_107", "sets": ["WikiPathways", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_108", "sets": ["WikiPathways", "Rhea", "KEGG", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_109", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "HMDB"], "size": 8}, {"id": "cmpd_src_set_110", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "Rhea"], "size": 10}, {"id": "cmpd_src_set_111", "sets": ["LIPIDMAPS", "Reactome", "HMDB", "Rhea"], "size": 3}, {"id": "cmpd_src_set_112", "sets": ["LIPIDMAPS", "Reactome", "HMDB", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_113", "sets": ["LIPIDMAPS", "Reactome", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_114", "sets": ["LIPIDMAPS", "RefMet", "HMDB", "Rhea"], "size": 294}, {"id": "cmpd_src_set_115", "sets": ["LIPIDMAPS", "RefMet", "HMDB", "PFOCR"], "size": 40}, {"id": "cmpd_src_set_116", "sets": ["LIPIDMAPS", "RefMet", "Rhea", "PFOCR"], "size": 35}, {"id": "cmpd_src_set_117", "sets": ["LIPIDMAPS", "HMDB", "Rhea", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_118", "sets": ["Reactome", "RefMet", "HMDB", "Rhea"], "size": 4}, {"id": "cmpd_src_set_119", "sets": ["Reactome", "RefMet", "HMDB", "PFOCR"], "size": 9}, {"id": "cmpd_src_set_120", "sets": ["Reactome", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_121", "sets": ["RefMet", "HMDB", "Rhea", "PFOCR"], "size": 48}, {"id": "cmpd_src_set_122", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "HMDB"], "size": 10}, {"id": "cmpd_src_set_123", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "RefMet", "HMDB"], "size": 27}, {"id": "cmpd_src_set_124", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "HMDB", "Rhea"], "size": 7}, {"id": "cmpd_src_set_125", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB", "KEGG"], "size": 1}, {"id": "cmpd_src_set_126", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "HMDB", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_127", "sets": ["KEGG", "WikiPathways", "Reactome", "RefMet", "HMDB"], "size": 8}, {"id": "cmpd_src_set_128", "sets": ["KEGG", "WikiPathways", "Reactome", "HMDB", "Rhea"], "size": 16}, {"id": "cmpd_src_set_129", "sets": ["KEGG", "WikiPathways", "Reactome", "HMDB", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_130", "sets": ["KEGG", "WikiPathways", "RefMet", "HMDB", "Rhea"], "size": 47}, {"id": "cmpd_src_set_131", "sets": ["WikiPathways", "RefMet", "HMDB", "KEGG"], "size": 5}, {"id": "cmpd_src_set_132", "sets": ["KEGG", "WikiPathways", "RefMet", "HMDB", "PFOCR"], "size": 61}, {"id": "cmpd_src_set_133", "sets": ["KEGG", "WikiPathways", "HMDB", "Rhea", "PFOCR"], "size": 5}, {"id": "cmpd_src_set_134", "sets": ["WikiPathways", "HMDB", "KEGG", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_135", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "RefMet", "HMDB"], "size": 5}, {"id": "cmpd_src_set_136", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "HMDB", "Rhea"], "size": 3}, {"id": "cmpd_src_set_137", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "HMDB", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_138", "sets": ["KEGG", "LIPIDMAPS", "RefMet", "HMDB", "Rhea"], "size": 135}, {"id": "cmpd_src_set_139", "sets": ["KEGG", "LIPIDMAPS", "RefMet", "HMDB", "PFOCR"], "size": 77}, {"id": "cmpd_src_set_140", "sets": ["KEGG", "LIPIDMAPS", "HMDB", "Rhea", "PFOCR"], "size": 24}, {"id": "cmpd_src_set_141", "sets": ["KEGG", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 36}, {"id": "cmpd_src_set_142", "sets": ["KEGG", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 60}, {"id": "cmpd_src_set_143", "sets": ["KEGG", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 11}, {"id": "cmpd_src_set_144", "sets": ["KEGG", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 253}, {"id": "cmpd_src_set_145", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB"], "size": 12}, {"id": "cmpd_src_set_146", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "Rhea"], "size": 7}, {"id": "cmpd_src_set_147", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_148", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea"], "size": 6}, {"id": "cmpd_src_set_149", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_150", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea"], "size": 60}, {"id": "cmpd_src_set_151", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "PFOCR"], "size": 17}, {"id": "cmpd_src_set_152", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "Rhea", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_153", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_154", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 1}, {"id": "cmpd_src_set_155", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 5}, {"id": "cmpd_src_set_156", "sets": ["WikiPathways", "Reactome", "RefMet", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_157", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_158", "sets": ["WikiPathways", "Reactome", "Rhea", "KEGG", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_159", "sets": ["WikiPathways", "RefMet", "HMDB", "Rhea", "KEGG"], "size": 5}, {"id": "cmpd_src_set_160", "sets": ["WikiPathways", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_161", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 12}, {"id": "cmpd_src_set_162", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_163", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "Rhea", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_164", "sets": ["LIPIDMAPS", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_165", "sets": ["LIPIDMAPS", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 54}, {"id": "cmpd_src_set_166", "sets": ["Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_167", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB"], "size": 14}, {"id": "cmpd_src_set_168", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea"], "size": 7}, {"id": "cmpd_src_set_169", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea"], "size": 41}, {"id": "cmpd_src_set_170", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "KEGG"], "size": 4}, {"id": "cmpd_src_set_171", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "PFOCR"], "size": 35}, {"id": "cmpd_src_set_172", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "HMDB", "Rhea", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_173", "sets": ["KEGG", "WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 40}, {"id": "cmpd_src_set_174", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "KEGG"], "size": 1}, {"id": "cmpd_src_set_175", "sets": ["KEGG", "WikiPathways", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 24}, {"id": "cmpd_src_set_176", "sets": ["KEGG", "WikiPathways", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 28}, {"id": "cmpd_src_set_177", "sets": ["WikiPathways", "Reactome", "HMDB", "KEGG", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_178", "sets": ["KEGG", "WikiPathways", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 76}, {"id": "cmpd_src_set_179", "sets": ["WikiPathways", "RefMet", "HMDB", "KEGG", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_180", "sets": ["WikiPathways", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_181", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 8}, {"id": "cmpd_src_set_182", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_183", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_184", "sets": ["KEGG", "LIPIDMAPS", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 160}, {"id": "cmpd_src_set_185", "sets": ["KEGG", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 52}, {"id": "cmpd_src_set_186", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 10}, {"id": "cmpd_src_set_187", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "KEGG"], "size": 13}, {"id": "cmpd_src_set_188", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_189", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_190", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 13}, {"id": "cmpd_src_set_191", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "KEGG", "PFOCR"], "size": 11}, {"id": "cmpd_src_set_192", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_193", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_194", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 46}, {"id": "cmpd_src_set_195", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 8}, {"id": "cmpd_src_set_196", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea", "KEGG"], "size": 1}, {"id": "cmpd_src_set_197", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_198", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 53}, {"id": "cmpd_src_set_199", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea", "KEGG"], "size": 7}, {"id": "cmpd_src_set_200", "sets": ["KEGG", "WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 167}, {"id": "cmpd_src_set_201", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "KEGG", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_202", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 5}, {"id": "cmpd_src_set_203", "sets": ["WikiPathways", "RefMet", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_204", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 16}, {"id": "cmpd_src_set_205", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_206", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 97}, {"id": "cmpd_src_set_207", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "KEGG"], "size": 4}, {"id": "cmpd_src_set_208", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 123}, {"id": "cmpd_src_set_209", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "KEGG", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_210", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_211", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 16}]|[{"id": "gene_src_set_1", "sets": ["WikiPathways"], "size": 1407}, {"id": "gene_src_set_2", "sets": ["Reactome"], "size": 1414}, {"id": "gene_src_set_3", "sets": ["HMDB"], "size": 1676}, {"id": "gene_src_set_4", "sets": ["Rhea"], "size": 13239}, {"id": "gene_src_set_5", "sets": ["PFOCR"], "size": 4564}, {"id": "gene_src_set_6", "sets": ["WikiPathways", "Reactome"], "size": 223}, {"id": "gene_src_set_7", "sets": ["WikiPathways", "HMDB"], "size": 28}, {"id": "gene_src_set_8", "sets": ["WikiPathways", "Rhea"], "size": 15}, {"id": "gene_src_set_9", "sets": ["WikiPathways", "PFOCR"], "size": 1026}, {"id": "gene_src_set_10", "sets": ["Reactome", "HMDB"], "size": 102}, {"id": "gene_src_set_11", "sets": ["Reactome", "Rhea"], "size": 38}, {"id": "gene_src_set_12", "sets": ["Reactome", "PFOCR"], "size": 2058}, {"id": "gene_src_set_13", "sets": ["HMDB", "Rhea"], "size": 353}, {"id": "gene_src_set_14", "sets": ["WikiPathways", "Reactome", "HMDB"], "size": 39}, {"id": "gene_src_set_15", "sets": ["WikiPathways", "Reactome", "Rhea"], "size": 6}, {"id": "gene_src_set_16", "sets": ["WikiPathways", "Reactome", "PFOCR"], "size": 2789}, {"id": "gene_src_set_17", "sets": ["WikiPathways", "HMDB", "Rhea"], "size": 16}, {"id": "gene_src_set_18", "sets": ["WikiPathways", "HMDB", "PFOCR"], "size": 178}, {"id": "gene_src_set_19", "sets": ["WikiPathways", "Rhea", "PFOCR"], "size": 33}, {"id": "gene_src_set_20", "sets": ["Reactome", "HMDB", "Rhea"], "size": 89}, {"id": "gene_src_set_21", "sets": ["Reactome", "HMDB", "PFOCR"], "size": 499}, {"id": "gene_src_set_22", "sets": ["Reactome", "Rhea", "PFOCR"], "size": 119}, {"id": "gene_src_set_23", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea"], "size": 64}, {"id": "gene_src_set_24", "sets": ["WikiPathways", "Reactome", "HMDB", "PFOCR"], "size": 1364}, {"id": "gene_src_set_25", "sets": ["WikiPathways", "Reactome", "Rhea", "PFOCR"], "size": 164}, {"id": "gene_src_set_26", "sets": ["WikiPathways", "HMDB", "Rhea", "PFOCR"], "size": 153}, {"id": "gene_src_set_27", "sets": ["Reactome", "HMDB", "Rhea", "PFOCR"], "size": 560}, {"id": "gene_src_set_28", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 1989}]|[{"id": "cmpd_src_set_1", "sets": ["WikiPathways"], "size": 655}, {"id": "cmpd_src_set_2", "sets": ["Reactome"], "size": 598}, {"id": "cmpd_src_set_3", "sets": ["HMDB"], "size": 32}, {"id": "cmpd_src_set_4", "sets": ["PFOCR"], "size": 1915}, {"id": "cmpd_src_set_5", "sets": ["KEGG", "HMDB"], "size": 17}, {"id": "cmpd_src_set_6", "sets": ["WikiPathways", "LIPIDMAPS"], "size": 198}, {"id": "cmpd_src_set_7", "sets": ["WikiPathways", "Reactome"], "size": 54}, {"id": "cmpd_src_set_8", "sets": ["WikiPathways", "RefMet"], "size": 7}, {"id": "cmpd_src_set_9", "sets": ["WikiPathways", "HMDB"], "size": 126}, {"id": "cmpd_src_set_10", "sets": ["WikiPathways", "Rhea"], "size": 74}, {"id": "cmpd_src_set_11", "sets": ["WikiPathways", "KEGG"], "size": 45}, {"id": "cmpd_src_set_12", "sets": ["WikiPathways", "PFOCR"], "size": 152}, {"id": "cmpd_src_set_13", "sets": ["LIPIDMAPS", "Reactome"], "size": 20}, {"id": "cmpd_src_set_14", "sets": ["LIPIDMAPS", "HMDB"], "size": 10}, {"id": "cmpd_src_set_15", "sets": ["LIPIDMAPS", "PFOCR"], "size": 77}, {"id": "cmpd_src_set_16", "sets": ["Reactome", "RefMet"], "size": 8}, {"id": "cmpd_src_set_17", "sets": ["Reactome", "HMDB"], "size": 77}, {"id": "cmpd_src_set_18", "sets": ["Reactome", "Rhea"], "size": 149}, {"id": "cmpd_src_set_19", "sets": ["Reactome", "PFOCR"], "size": 78}, {"id": "cmpd_src_set_20", "sets": ["RefMet", "HMDB"], "size": 3}, {"id": "cmpd_src_set_21", "sets": ["RefMet", "PFOCR"], "size": 126}, {"id": "cmpd_src_set_22", "sets": ["HMDB", "PFOCR"], "size": 302}, {"id": "cmpd_src_set_23", "sets": ["Rhea", "PFOCR"], "size": 317}, {"id": "cmpd_src_set_24", "sets": ["KEGG", "WikiPathways", "HMDB"], "size": 47}, {"id": "cmpd_src_set_25", "sets": ["KEGG", "LIPIDMAPS", "HMDB"], "size": 4}, {"id": "cmpd_src_set_26", "sets": ["KEGG", "Reactome", "HMDB"], "size": 21}, {"id": "cmpd_src_set_27", "sets": ["KEGG", "RefMet", "HMDB"], "size": 15}, {"id": "cmpd_src_set_28", "sets": ["KEGG", "HMDB", "Rhea"], "size": 5}, {"id": "cmpd_src_set_29", "sets": ["KEGG", "HMDB", "PFOCR"], "size": 133}, {"id": "cmpd_src_set_30", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome"], "size": 5}, {"id": "cmpd_src_set_31", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet"], "size": 33}, {"id": "cmpd_src_set_32", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB"], "size": 21}, {"id": "cmpd_src_set_33", "sets": ["WikiPathways", "LIPIDMAPS", "Rhea"], "size": 44}, {"id": "cmpd_src_set_34", "sets": ["WikiPathways", "LIPIDMAPS", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_35", "sets": ["WikiPathways", "Reactome", "HMDB"], "size": 9}, {"id": "cmpd_src_set_36", "sets": ["WikiPathways", "Reactome", "Rhea"], "size": 48}, {"id": "cmpd_src_set_37", "sets": ["WikiPathways", "Reactome", "PFOCR"], "size": 38}, {"id": "cmpd_src_set_38", "sets": ["WikiPathways", "RefMet", "HMDB"], "size": 50}, {"id": "cmpd_src_set_39", "sets": ["WikiPathways", "RefMet", "Rhea"], "size": 2}, {"id": "cmpd_src_set_40", "sets": ["WikiPathways", "RefMet", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_41", "sets": ["WikiPathways", "HMDB", "Rhea"], "size": 5}, {"id": "cmpd_src_set_42", "sets": ["WikiPathways", "HMDB", "PFOCR"], "size": 32}, {"id": "cmpd_src_set_43", "sets": ["WikiPathways", "Rhea", "PFOCR"], "size": 19}, {"id": "cmpd_src_set_44", "sets": ["WikiPathways", "KEGG", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_45", "sets": ["LIPIDMAPS", "Reactome", "RefMet"], "size": 12}, {"id": "cmpd_src_set_46", "sets": ["LIPIDMAPS", "Reactome", "HMDB"], "size": 3}, {"id": "cmpd_src_set_47", "sets": ["LIPIDMAPS", "Reactome", "Rhea"], "size": 19}, {"id": "cmpd_src_set_48", "sets": ["LIPIDMAPS", "Reactome", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_49", "sets": ["LIPIDMAPS", "RefMet", "HMDB"], "size": 11}, {"id": "cmpd_src_set_50", "sets": ["LIPIDMAPS", "RefMet", "PFOCR"], "size": 20}, {"id": "cmpd_src_set_51", "sets": ["LIPIDMAPS", "HMDB", "PFOCR"], "size": 13}, {"id": "cmpd_src_set_52", "sets": ["LIPIDMAPS", "Rhea", "PFOCR"], "size": 55}, {"id": "cmpd_src_set_53", "sets": ["Reactome", "RefMet", "HMDB"], "size": 29}, {"id": "cmpd_src_set_54", "sets": ["Reactome", "RefMet", "Rhea"], "size": 6}, {"id": "cmpd_src_set_55", "sets": ["Reactome", "RefMet", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_56", "sets": ["Reactome", "HMDB", "Rhea"], "size": 4}, {"id": "cmpd_src_set_57", "sets": ["Reactome", "HMDB", "PFOCR"], "size": 27}, {"id": "cmpd_src_set_58", "sets": ["Reactome", "Rhea", "PFOCR"], "size": 24}, {"id": "cmpd_src_set_59", "sets": ["RefMet", "HMDB", "Rhea"], "size": 2}, {"id": "cmpd_src_set_60", "sets": ["RefMet", "HMDB", "PFOCR"], "size": 171}, {"id": "cmpd_src_set_61", "sets": ["RefMet", "Rhea", "PFOCR"], "size": 75}, {"id": "cmpd_src_set_62", "sets": ["HMDB", "Rhea", "PFOCR"], "size": 28}, {"id": "cmpd_src_set_63", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "HMDB"], "size": 8}, {"id": "cmpd_src_set_64", "sets": ["KEGG", "WikiPathways", "Reactome", "HMDB"], "size": 7}, {"id": "cmpd_src_set_65", "sets": ["KEGG", "WikiPathways", "RefMet", "HMDB"], "size": 52}, {"id": "cmpd_src_set_66", "sets": ["KEGG", "WikiPathways", "HMDB", "Rhea"], "size": 11}, {"id": "cmpd_src_set_67", "sets": ["WikiPathways", "HMDB", "KEGG"], "size": 2}, {"id": "cmpd_src_set_68", "sets": ["KEGG", "WikiPathways", "HMDB", "PFOCR"], "size": 31}, {"id": "cmpd_src_set_69", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "HMDB"], "size": 4}, {"id": "cmpd_src_set_70", "sets": ["KEGG", "LIPIDMAPS", "RefMet", "HMDB"], "size": 8}, {"id": "cmpd_src_set_71", "sets": ["KEGG", "LIPIDMAPS", "HMDB", "Rhea"], "size": 2}, {"id": "cmpd_src_set_72", "sets": ["KEGG", "LIPIDMAPS", "HMDB", "PFOCR"], "size": 20}, {"id": "cmpd_src_set_73", "sets": ["KEGG", "Reactome", "RefMet", "HMDB"], "size": 52}, {"id": "cmpd_src_set_74", "sets": ["KEGG", "Reactome", "HMDB", "Rhea"], "size": 13}, {"id": "cmpd_src_set_75", "sets": ["KEGG", "Reactome", "HMDB", "PFOCR"], "size": 16}, {"id": "cmpd_src_set_76", "sets": ["KEGG", "RefMet", "HMDB", "Rhea"], "size": 10}, {"id": "cmpd_src_set_77", "sets": ["KEGG", "RefMet", "HMDB", "PFOCR"], "size": 316}, {"id": "cmpd_src_set_78", "sets": ["KEGG", "HMDB", "Rhea", "PFOCR"], "size": 53}, {"id": "cmpd_src_set_79", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet"], "size": 7}, {"id": "cmpd_src_set_80", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "HMDB"], "size": 3}, {"id": "cmpd_src_set_81", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "Rhea"], "size": 5}, {"id": "cmpd_src_set_82", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB"], "size": 78}, {"id": "cmpd_src_set_83", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "Rhea"], "size": 8}, {"id": "cmpd_src_set_84", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "KEGG"], "size": 1}, {"id": "cmpd_src_set_85", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_86", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB", "Rhea"], "size": 6}, {"id": "cmpd_src_set_87", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_88", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB"], "size": 4}, {"id": "cmpd_src_set_89", "sets": ["WikiPathways", "Reactome", "RefMet", "Rhea"], "size": 2}, {"id": "cmpd_src_set_90", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea"], "size": 6}, {"id": "cmpd_src_set_91", "sets": ["WikiPathways", "Reactome", "HMDB", "PFOCR"], "size": 10}, {"id": "cmpd_src_set_92", "sets": ["WikiPathways", "Reactome", "Rhea", "KEGG"], "size": 1}, {"id": "cmpd_src_set_93", "sets": ["WikiPathways", "Reactome", "Rhea", "PFOCR"], "size": 27}, {"id": "cmpd_src_set_94", "sets": ["WikiPathways", "Reactome", "KEGG", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_95", "sets": ["WikiPathways", "RefMet", "HMDB", "Rhea"], "size": 10}, {"id": "cmpd_src_set_96", "sets": ["WikiPathways", "RefMet", "HMDB", "PFOCR"], "size": 28}, {"id": "cmpd_src_set_97", "sets": ["WikiPathways", "RefMet", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_98", "sets": ["WikiPathways", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_99", "sets": ["WikiPathways", "Rhea", "KEGG", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_100", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "HMDB"], "size": 8}, {"id": "cmpd_src_set_101", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "Rhea"], "size": 10}, {"id": "cmpd_src_set_102", "sets": ["LIPIDMAPS", "Reactome", "HMDB", "Rhea"], "size": 3}, {"id": "cmpd_src_set_103", "sets": ["LIPIDMAPS", "Reactome", "HMDB", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_104", "sets": ["LIPIDMAPS", "Reactome", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_105", "sets": ["LIPIDMAPS", "RefMet", "HMDB", "Rhea"], "size": 7}, {"id": "cmpd_src_set_106", "sets": ["LIPIDMAPS", "RefMet", "HMDB", "PFOCR"], "size": 40}, {"id": "cmpd_src_set_107", "sets": ["LIPIDMAPS", "RefMet", "Rhea", "PFOCR"], "size": 35}, {"id": "cmpd_src_set_108", "sets": ["LIPIDMAPS", "HMDB", "Rhea", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_109", "sets": ["Reactome", "RefMet", "HMDB", "Rhea"], "size": 4}, {"id": "cmpd_src_set_110", "sets": ["Reactome", "RefMet", "HMDB", "PFOCR"], "size": 9}, {"id": "cmpd_src_set_111", "sets": ["Reactome", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_112", "sets": ["RefMet", "HMDB", "Rhea", "PFOCR"], "size": 48}, {"id": "cmpd_src_set_113", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "HMDB"], "size": 10}, {"id": "cmpd_src_set_114", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "RefMet", "HMDB"], "size": 27}, {"id": "cmpd_src_set_115", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "HMDB", "Rhea"], "size": 7}, {"id": "cmpd_src_set_116", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB", "KEGG"], "size": 1}, {"id": "cmpd_src_set_117", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "HMDB", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_118", "sets": ["KEGG", "WikiPathways", "Reactome", "RefMet", "HMDB"], "size": 8}, {"id": "cmpd_src_set_119", "sets": ["KEGG", "WikiPathways", "Reactome", "HMDB", "Rhea"], "size": 16}, {"id": "cmpd_src_set_120", "sets": ["KEGG", "WikiPathways", "Reactome", "HMDB", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_121", "sets": ["KEGG", "WikiPathways", "RefMet", "HMDB", "Rhea"], "size": 47}, {"id": "cmpd_src_set_122", "sets": ["WikiPathways", "RefMet", "HMDB", "KEGG"], "size": 5}, {"id": "cmpd_src_set_123", "sets": ["KEGG", "WikiPathways", "RefMet", "HMDB", "PFOCR"], "size": 61}, {"id": "cmpd_src_set_124", "sets": ["KEGG", "WikiPathways", "HMDB", "Rhea", "PFOCR"], "size": 5}, {"id": "cmpd_src_set_125", "sets": ["WikiPathways", "HMDB", "KEGG", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_126", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "RefMet", "HMDB"], "size": 5}, {"id": "cmpd_src_set_127", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "HMDB", "Rhea"], "size": 3}, {"id": "cmpd_src_set_128", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "HMDB", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_129", "sets": ["KEGG", "LIPIDMAPS", "RefMet", "HMDB", "Rhea"], "size": 6}, {"id": "cmpd_src_set_130", "sets": ["KEGG", "LIPIDMAPS", "RefMet", "HMDB", "PFOCR"], "size": 77}, {"id": "cmpd_src_set_131", "sets": ["KEGG", "LIPIDMAPS", "HMDB", "Rhea", "PFOCR"], "size": 24}, {"id": "cmpd_src_set_132", "sets": ["KEGG", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 36}, {"id": "cmpd_src_set_133", "sets": ["KEGG", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 60}, {"id": "cmpd_src_set_134", "sets": ["KEGG", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 11}, {"id": "cmpd_src_set_135", "sets": ["KEGG", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 253}, {"id": "cmpd_src_set_136", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB"], "size": 12}, {"id": "cmpd_src_set_137", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "Rhea"], "size": 7}, {"id": "cmpd_src_set_138", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_139", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea"], "size": 6}, {"id": "cmpd_src_set_140", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_141", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea"], "size": 60}, {"id": "cmpd_src_set_142", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "PFOCR"], "size": 17}, {"id": "cmpd_src_set_143", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "Rhea", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_144", "sets": ["WikiPathways", "LIPIDMAPS", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_145", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 1}, {"id": "cmpd_src_set_146", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 5}, {"id": "cmpd_src_set_147", "sets": ["WikiPathways", "Reactome", "RefMet", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_148", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_149", "sets": ["WikiPathways", "Reactome", "Rhea", "KEGG", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_150", "sets": ["WikiPathways", "RefMet", "HMDB", "Rhea", "KEGG"], "size": 5}, {"id": "cmpd_src_set_151", "sets": ["WikiPathways", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_152", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 12}, {"id": "cmpd_src_set_153", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_154", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "Rhea", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_155", "sets": ["LIPIDMAPS", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_156", "sets": ["LIPIDMAPS", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 54}, {"id": "cmpd_src_set_157", "sets": ["Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_158", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB"], "size": 14}, {"id": "cmpd_src_set_159", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea"], "size": 7}, {"id": "cmpd_src_set_160", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea"], "size": 41}, {"id": "cmpd_src_set_161", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "KEGG"], "size": 4}, {"id": "cmpd_src_set_162", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "PFOCR"], "size": 35}, {"id": "cmpd_src_set_163", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "HMDB", "Rhea", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_164", "sets": ["KEGG", "WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 40}, {"id": "cmpd_src_set_165", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "KEGG"], "size": 1}, {"id": "cmpd_src_set_166", "sets": ["KEGG", "WikiPathways", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 24}, {"id": "cmpd_src_set_167", "sets": ["KEGG", "WikiPathways", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 28}, {"id": "cmpd_src_set_168", "sets": ["WikiPathways", "Reactome", "HMDB", "KEGG", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_169", "sets": ["KEGG", "WikiPathways", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 76}, {"id": "cmpd_src_set_170", "sets": ["WikiPathways", "RefMet", "HMDB", "KEGG", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_171", "sets": ["WikiPathways", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_172", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 8}, {"id": "cmpd_src_set_173", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_174", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_175", "sets": ["KEGG", "LIPIDMAPS", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 160}, {"id": "cmpd_src_set_176", "sets": ["KEGG", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 52}, {"id": "cmpd_src_set_177", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 10}, {"id": "cmpd_src_set_178", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "KEGG"], "size": 13}, {"id": "cmpd_src_set_179", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 3}, {"id": "cmpd_src_set_180", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_181", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 13}, {"id": "cmpd_src_set_182", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "KEGG", "PFOCR"], "size": 11}, {"id": "cmpd_src_set_183", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_184", "sets": ["LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_185", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea"], "size": 46}, {"id": "cmpd_src_set_186", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "PFOCR"], "size": 8}, {"id": "cmpd_src_set_187", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea", "KEGG"], "size": 1}, {"id": "cmpd_src_set_188", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 2}, {"id": "cmpd_src_set_189", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 53}, {"id": "cmpd_src_set_190", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea", "KEGG"], "size": 7}, {"id": "cmpd_src_set_191", "sets": ["KEGG", "WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 167}, {"id": "cmpd_src_set_192", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "KEGG", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_193", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 5}, {"id": "cmpd_src_set_194", "sets": ["WikiPathways", "RefMet", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 6}, {"id": "cmpd_src_set_195", "sets": ["KEGG", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 16}, {"id": "cmpd_src_set_196", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 7}, {"id": "cmpd_src_set_197", "sets": ["WikiPathways", "Reactome", "RefMet", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 97}, {"id": "cmpd_src_set_198", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "KEGG"], "size": 4}, {"id": "cmpd_src_set_199", "sets": ["KEGG", "WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "PFOCR"], "size": 123}, {"id": "cmpd_src_set_200", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "KEGG", "PFOCR"], "size": 4}, {"id": "cmpd_src_set_201", "sets": ["WikiPathways", "LIPIDMAPS", "RefMet", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 1}, {"id": "cmpd_src_set_202", "sets": ["WikiPathways", "LIPIDMAPS", "Reactome", "RefMet", "HMDB", "Rhea", "KEGG", "PFOCR"], "size": 16}]|[{"id": "gene_src_set_1", "sets": ["WikiPathways"], "size": 1407}, {"id": "gene_src_set_2", "sets": ["Reactome"], "size": 1414}, {"id": "gene_src_set_3", "sets": ["HMDB"], "size": 530}, {"id": "gene_src_set_4", "sets": ["PFOCR"], "size": 4564}, {"id": "gene_src_set_5", "sets": ["WikiPathways", "Reactome"], "size": 223}, {"id": "gene_src_set_6", "sets": ["WikiPathways", "HMDB"], "size": 28}, {"id": "gene_src_set_7", "sets": ["WikiPathways", "Rhea"], "size": 15}, {"id": "gene_src_set_8", "sets": ["WikiPathways", "PFOCR"], "size": 1026}, {"id": "gene_src_set_9", "sets": ["Reactome", "HMDB"], "size": 102}, {"id": "gene_src_set_10", "sets": ["Reactome", "Rhea"], "size": 38}, {"id": "gene_src_set_11", "sets": ["Reactome", "PFOCR"], "size": 2058}, {"id": "gene_src_set_12", "sets": ["HMDB", "Rhea"], "size": 37}, {"id": "gene_src_set_13", "sets": ["WikiPathways", "Reactome", "HMDB"], "size": 39}, {"id": "gene_src_set_14", "sets": ["WikiPathways", "Reactome", "Rhea"], "size": 6}, {"id": "gene_src_set_15", "sets": ["WikiPathways", "Reactome", "PFOCR"], "size": 2789}, {"id": "gene_src_set_16", "sets": ["WikiPathways", "HMDB", "Rhea"], "size": 16}, {"id": "gene_src_set_17", "sets": ["WikiPathways", "HMDB", "PFOCR"], "size": 178}, {"id": "gene_src_set_18", "sets": ["WikiPathways", "Rhea", "PFOCR"], "size": 33}, {"id": "gene_src_set_19", "sets": ["Reactome", "HMDB", "Rhea"], "size": 89}, {"id": "gene_src_set_20", "sets": ["Reactome", "HMDB", "PFOCR"], "size": 499}, {"id": "gene_src_set_21", "sets": ["Reactome", "Rhea", "PFOCR"], "size": 119}, {"id": "gene_src_set_22", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea"], "size": 64}, {"id": "gene_src_set_23", "sets": ["WikiPathways", "Reactome", "HMDB", "PFOCR"], "size": 1364}, {"id": "gene_src_set_24", "sets": ["WikiPathways", "Reactome", "Rhea", "PFOCR"], "size": 164}, {"id": "gene_src_set_25", "sets": ["WikiPathways", "HMDB", "Rhea", "PFOCR"], "size": 153}, {"id": "gene_src_set_26", "sets": ["Reactome", "HMDB", "Rhea", "PFOCR"], "size": 560}, {"id": "gene_src_set_27", "sets": ["WikiPathways", "Reactome", "HMDB", "Rhea", "PFOCR"], "size": 1989}]|

# reaction
#   0|ramp_rxn_id|VARCHAR(16)|1||1
#   1|rxn_source_id|VARCHAR(16)|1||0
#   2|status|INTEGER|1||0
#   3|is_transport|INTEGER|1||0
#   4|direction|VARCHAR(8)|1||0
#   5|label|VARCHAR(256)|1||0
#   6|equation|VARCHAR(256)|1||0
#   7|html_equation|VARCHAR(256)|1||0
#   8|ec_num|VARCHAR(256)|0||0
#   9|has_human_prot|INTEGER|1||0
#   10|only_human_mets|INTEGER|1||0
#   e.g. RAMP_R_000000001|rhea:10000|1|0|UN|H2O + pentanamide = NH4(+) + pentanoate|H2O + pentanamide = NH4(+) + pentanoate|H<small><sub>2</sub></small>O + pentanamide = NH<small><sub>4</sub></small><small><sup>+</sup></small> + pentanoate|3.5.1.50|0|0
