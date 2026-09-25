"""
Parse MetaChatDB and emit Entity records.

This module converts MetaChatDB entries into Entity records using the
declarative schema pattern.
"""

from __future__ import annotations
import os

from pypath.inputs_v2.base import iter_tsv
from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
)
from pypath.internals.tabular_builder import (
    MembershipBuilder,
    Member,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    MoleculeAnnotationsCv,
    BiologicalRoleCv,
    MoleculeSubtypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    ProteinFunctionalClassCv,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

BASE_URL = (
    'https://raw.githubusercontent.com/SonghaoLuo/MetaChat/refs/heads/main/'
    'metachat/preprocessing/_data/MetaChatDB/MetaChatDB_%s.tsv'
)

TAXON_ID = {
    'human': '9606',
    'mouse': '10090',
}

files = ['human', 'mouse']

config = ResourceConfig(
    id=ResourceCv.METACHATDB,
    name='MetaChatDB',
    url='https://github.com/SonghaoLuo/MetaChat',
    license=LicenseCV.MIT,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='',
    primary_category='interactions',
    description=(
        'MetaChatDB is a literature-supported database for metabolite-sensor '
        'interactions for both human and mouse.'
    ),
)

download = {
    k: Download(
        url=BASE_URL % k,
        filename=os.path.basename(BASE_URL % k),
        subfolder='MetaChatDB',
        large=True,
        ext='.tsv',
        default_mode='r',
    )
    for k in files
}

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={
        'wp_id': r'^(WP\d+):',
        'notna': r'\b(.*)\b(?<!NA)',
        'pmid': r'(\d+)',
    },
    map={
        'type_to_entity': {
            'nuclear receptor': ProteinFunctionalClassCv.NUCLEAR_HORMONE_RECEPTOR,
            'receptor': EntityTypeCv.RECEPTOR,
            'transporter': ProteinFunctionalClassCv.TRANSPORTER,
        },
        'class_to_entity': {
            'Benzenoids': MoleculeSubtypeCv.SYNTHETIC_ORGANIC,
            'Homogeneous metal compounds': MoleculeSubtypeCv.INORGANIC,
            'Homogeneous non-metal compounds': MoleculeSubtypeCv.INORGANIC,
            'Lipids and lipid-like molecules': MoleculeSubtypeCv.LIPID,
            'Nucleosides, nucleotides, and analogues': MoleculeSubtypeCv.NUCLEIC_ACID,
            'Organic acids and derivatives': MoleculeSubtypeCv.METABOLITE,
            'Organic nitrogen compounds': MoleculeSubtypeCv.METABOLITE,
            'Organic oxygen compounds': MoleculeSubtypeCv.METABOLITE,
            'Organoheterocyclic compounds': MoleculeSubtypeCv.METABOLITE,
        }
    },
)

def schema(key):


    return EntityBuilder(
        entity_type=EntityTypeCv.INTERACTION,
        annotations=AnnotationsBuilder(
            CV(
                term=IdentifierNamespaceCv.PUBMED,
                value=f('Evidences', delimiter='; ', extract='pmid')
            ),
            CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=TAXON_ID[key])
        ),
        membership=MembershipBuilder(
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.NAME, value=f('Sensor.Name')),
                        CV(
                            term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                            value=f('Sensor.Gene')
                        ),
                        CV()
                    )
                ),
                annotations=AnnotationsBuilder(
                    CV(term=f('Sensor.Type', delimiter=',', map='type_to_entity')),
                    CV(
                        term=IdentifierNamespaceCv.WIKIPATHWAYS,
                        value=f('Sensor.Pathway', delimiter='; ', extract='wp_id')
                    ),
                ),
            ),
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.SMALL_MOLECULE,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.HMDB, value=f('HMDB.ID')),
                        CV(
                            term=IdentifierNamespaceCv.NAME,
                            value=f('Metabolite.Name')
                        ),
                    )
                ),
                annotations=AnnotationsBuilder(
                    CV(term=f('Metabolite.Class', map='class_to_entity')),
                    CV(
                        term=BiologicalRoleCv.PATHWAY,
                        value=f('Metabolite.Pathway', delimiter='; ', extract='notna')
                    ),
                    CV(
                        term=MoleculeAnnotationsCv.BIOSPECIMEN_LOCATION,
                        value=f('Long.Range.Channel', delimiter=',')
                    )
                ),
            ),
        ),
    )

# ================================= RESOURCE ===================================

resource = Resource(
    config=config,
    **{
        k: Dataset(
            download=download[k],
            mapper=schema(k),
            raw_parser=iter_tsv
        )
        for k in files
    }
)

# ================================= REFERENCE ==================================
# X             X                                                                                                                           X               X                                                                                                           X           X                               NO                      X                                   X                                                                       X                                                           X
# Sensor.Gene	Sensor.Name	                                                                                                                Sensor.Type	    Sensor.Pathway	                                                                                            HMDB.ID	    Metabolite.Name	                Metabolite.transporter	Metabolite.Class	                Metabolite.Pathway	                                                    Long.Range.Channel	                                        Evidences	                    Sources
# ABCA1	        Phospholipid-transporting ATPase ABCA1 (EC 7.6.2.1) (ATP-binding cassette sub-family A member 1) (ATP-binding casset...	    transporter	    WP2011: SREBF and miR33 in cholesterol and lipid homeostasis; WP299: Nuclear receptors in lipid metab...	HMDB0006247	25-Hydroxycholesterol		                            Lipids and lipid-like molecules	    Primary bile acid biosynthesis	                                        Blood,Urine	                                                16611739	                    Mebocost; MRCLinkdb
# ABCA1	        Phospholipid-transporting ATPase ABCA1 (EC 7.6.2.1) (ATP-binding cassette sub-family A member 1) (ATP-binding casset...	    transporter	    WP2011: SREBF and miR33 in cholesterol and lipid homeostasis; WP299: Nuclear receptors in lipid metab...	HMDB0000067	Cholesterol		                                        Lipids and lipid-like molecules	    Primary bile acid biosynthesis; Steroid Biosynthesis; Steroidogenesis	Bile,Blood,Cerebrospinal Fluid (CSF),Feces,Saliva,Urine	    39978466; 36521550; 37181814	MetalinkDB
# ABCA13	    ATP-binding cassette sub-family A member 13 (EC 7.6.2.-)	                                                                transporter		                                                                                                            HMDB0000067	Cholesterol		                                        Lipids and lipid-like molecules	    Primary bile acid biosynthesis; Steroid Biosynthesis; Steroidogenesis	Bile,Blood,Cerebrospinal Fluid (CSF),Feces,Saliva,Urine	    33478937	                    MetaChatDB
# ABCB4	        Phosphatidylcholine translocator ABCB4 (EC 7.6.2.1) (ATP-binding cassette sub-family B member 4) (Multidrug resistan...	    transporter	    WP2879: Farnesoid X receptor pathway; WP2882: Nuclear receptors meta pathway; WP299: Nuclear receptor...	HMDB0008138	PC(18:2(9Z,12Z)/18:2(9Z,12Z))		                    Lipids and lipid-like molecules	    NA	                                                                    Blood,Feces,Saliva,Urine	                                32917728	                    MetaChatDB
# ABCB6	        ATP-binding cassette sub-family B member 6 (ABC-type heme transporter ABCB6) (EC 7.6.2.5) (Mitochondrial ABC transpo...     transporter		                                                                                                            HMDB0001261	Coproporphyrinogen III		                            Organoheterocyclic compounds	    Porphyrin Metabolism	                                                Blood	                                                    23792964	                    MetaChatDB
