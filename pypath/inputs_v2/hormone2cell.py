"""
Parse Hormone2Cell data and emit Entity records.

This module converts hormone-cell interaction information from Hormone2Cell into
Entity records using the declarative schema pattern.
"""

from __future__ import annotations

import pandas as pd

from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers import brenda as _parsers
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    OntologyCv,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

URL = 'https://hormonecellatlas.org.uk/accounts/downloads/file/table-s2b/'

config = ResourceConfig(
    id=ResourceCv.HORMONE2CELL,
    name='Hormone2Cell',
    url='https://hormonecellatlas.org.uk/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='42207862',
    primary_category='intercell',
    description=(
        'The Hormone Cell Atlas uses single cell transcriptomics to map '
        'hormone production and action across 47 disease-free human tissues at '
        'cellular resolution.'
    ),
)


download = Download(
    url=URL,
    filename='hormone2cell.xlsx',
    subfolder='hormone2cell',
    large=True,
    ext='xlsx',
    default_mode='rb',
)

def parser(opener):

    df = pd.read_excel(opener.path, skiprows=3)

    yield from df.to_dict(orient='records')

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={},
    map={},
    transform={},
)

schema = EntityBuilder(
    #entity_type=EntityTypeCv.
)

# ================================= RESOURCE ===================================

#resource = Resource(
#    config=config,
#    data=Dataset(
#        download=download,
#        mapper=schema,
#        raw_parser=parser,
#    ),
#)

# ================================= REFERENCE ==================================
# S1B
# Tissue	        System	        Classical endocrine cell type	Neuroendocrine cell type	Lineage	    Tissue level broad-grained cell type (celltype_level1)	    Tissue level fine-grained cell type (celltype_level2)	    Tissue level Cell Type Broad 1	    Tissue level Cell Type Broad 2	    Tissue level Cell Type Broad 3	    Celltype_count
# Fallopian Tube	Reproductive	not	                            not	                        epithelial  epithelial_cell_nonciliated                             	epithelial_cell_nonciliated_1	                            epithelial	                        epithelial_cell_nonciliated	        epithelial_cell_nonciliated	        7477
# Fallopian Tube	Reproductive	not	                            not	                        immune      t_cell                              	                    t_cell_cd4	                                                T_cell	                            T_cell	                            CD4_T_cell	                        5782
# Fallopian Tube	Reproductive	not	                            not	                        mesenchymal smooth_muscle_cell                              	        smooth_muscle_cell	                                        smooth_muscle_cell              	smooth_muscle_cell	                smooth_muscle_cell	                6438
# Fallopian Tube	Reproductive	not	                            not	                        epithelial  epithelial_cell_ciliated                                	epithelial_cell_ciliated	                                epithelial	                        epithelial_cell_ciliated	        epithelial_cell_ciliated	        4574
# Fallopian Tube	Reproductive	not	                            not	                        mesenchymal myofibroblast                               	            myofibroblast_3	                                            myofibroblast	                    myofibroblast	                    myofibroblast	                    4863
# Fallopian Tube	Reproductive	not	                            not	                        epithelial  epithelial_cell_nonciliated                             	epithelial_cell_nonciliated_2	                            epithelial	                        epithelial_cell_nonciliated     	epithelial_cell_nonciliated	        7786

# S2C
# gene	    category	ID	        status_in_h2c
# POMC	    hormone	    hormone_1	h2c_included
# INHBA 	hormone	    hormone_2	h2c_included
# INHBB	    hormone	    hormone_3	h2c_included
# INHBE 	hormone	    hormone_4	h2c_included
# C1QTNF12	hormone	    hormone_5	h2c_included

# S2D
# hormone_short	hpc_include1	hpc_include2	hpc_include3	hpc_include4	hpc_include5	hpc_include6	hpc_include7	hpc_include8	hpc_exclude1	hpc_exclude2	receptor_known	hormone_type_broad	hormone_type_fine
# adm2	        ADM2										                                                                                                                    yes	            peptide_protein 	peptide
# adrenaline	TH	            DBH	            PNMT	        SLC18A1							                                                                                yes	            amine_derived   	amine_derived
# adrenaline	TH	            DBH	            PNMT	        SLC18A2							                                                                                yes	            amine_derived   	amine_derived
# adropin	    ENHO										                                                                                                                    yes	            peptide_protein 	peptide
# agrp	        AGRP										                                                                                                                    yes	            peptide_protein 	peptide
# aldosterone	STAR	        CYP11A1	        FDX1	        FDXR	        HSD3B1	        CYP21A2	        CYP11B2		                    CYP17A1		                    yes	            steroid 	        steroid
# aldosterone	STAR	        CYP11A1	        FDX1	        FDXR	        HSD3B2	        CYP21A2	        CYP11B2		                    CYP17A1		                    yes	            steroid 	        steroid
# amh         	AMH									                                        	                                                                                yes	            peptide_protein 	protein_multimer
# amylin_iapp	IAPP	        PCSK1	        PCSK2	        SCG5							                                                                                yes	            peptide_protein 	peptide
# androgens	    STAR	        CYP11A1	        FDX1	        FDXR	        CYP17A1	        CYB5A			                                CYP19A1	        CYP21A2	        yes	            steroid 	        steroid

# S2E
# hormone_short	    receptorgene1	    receptorgene2	    gene1_name	                gene2_name	                receptor_type_broad	    receptor_status
# acth      	    MC2R		                            melanocortin 2 receptor		                            gpcr                	established
# activin_a     	ACVR2A	            ACVR1B	            Activin A receptor_type2A	Activin A receptor_type1B	enzyme_linked	        established
# activin_a     	ACVR2B	            ACVR1B	            Activin A receptor_type2B	Activin A receptor_type1B	enzyme_linked	        established
# activin_ab        ACVR2A	            ACVR1B	            Activin A receptor_type2A	Activin A receptor_type1B	enzyme_linked	        established
# activin_ab        ACVR2A	            ACVR1C	            Activin A receptor_type2A	Activin A receptor_type1C	enzyme_linked	        established
# activin_ab        ACVR2B	            ACVR1B	            Activin A receptor_type2B	Activin A receptor_type1B	enzyme_linked	        established
# activin_ab        ACVR2B	            ACVR1C	            Activin A receptor_type2B	Activin A receptor_type1C	enzyme_linked	        established
# activin_b     	ACVR2A	            ACVR1B	            Activin A receptor_type2A	Activin A receptor_type1B	enzyme_linked	        established
# activin_b     	ACVR2A	            ACVR1C	            Activin A receptor_type2A	Activin A receptor_type1C	enzyme_linked	        established
# activin_b     	ACVR2B	            ACVR1B	            Activin A receptor_type2B	Activin A receptor_type1B	enzyme_linked	        established
# activin_b     	ACVR2B	            ACVR1C	            Activin A receptor_type2B	Activin A receptor_type1C	enzyme_linked	        established
# activin_e     	ACVR2A	            ACVR1C	            Activin A receptor_type2A	Activin A receptor_type1C	enzyme_linked	        established
# activin_e     	ACVR2B	            ACVR1C	            Activin A receptor_type2B	Activin A receptor_type1C	enzyme_linked	        established

# S2F
# gene_id	gene_name	                                                uniprot_id	h2c_present	type	            analysis_cat	reaction	                                        RHEA	    other genes	    other_name
# CYP11A1	Cholesterol side-chain cleavage enzyme, mitochondrial	    P05108	    1	        enzyme_synthesis	synthesis	    Pregnenolone_from_cholesterol	                    Rhea:35739	FDX1, FDXR	    P450scc
# CYP11B1	Cytochrome P450 11B1, mitochondrial	                        P15538	    1	        enzyme_synthesis	synthesis	    Corticosterone_from_11deoxycorticosterone	        Rhea:46104
# CYP11B1	Cytochrome P450 11B1, mitochondrial                 	    P15538	    1	        enzyme_synthesis	synthesis	    Cortisol_from_11deoxycortisol	                    Rhea:46100
# CYP11B2	Cytochrome P450 11B2, mitochondrial (aldosterone synthase)	P19099	    1	        enzyme_synthesis	synthesis	    Aldosterone_from_18-hydroxycorticosterone	        Rhea:50792
# CYP11B2	Cytochrome P450 11B2, mitochondrial (aldosterone synthase)	P19099	    1	        enzyme_synthesis	synthesis	    Corticosterone_from_11deoxycorticosterone	        Rhea:46104
# CYP11B2	Cytochrome P450 11B2, mitochondrial (aldosterone synthase)	P19099	    1	        enzyme_synthesis	synthesis	    Cortisol_from_11deoxycortisol	                    Rhea:46100
# CYP11B2	Cytochrome P450 11B2, mitochondrial (aldosterone synthase)	P19099	    1	        enzyme_synthesis	synthesis	    18-hydroxycorticosterone_from_hydroxycorticosterone	Rhea:11872
# CYP17A1	Steroid 17-alpha-hydroxylase/17,20 lyase	                enzyme	    1	        enzyme_synthesis	synthesis	    17alpha-hydroxyprogesterone_from_progesterone	    Rhea:50236
# CYP17A1	Steroid 17-alpha-hydroxylase/17,20 lyase	                P05093	    1	        enzyme_synthesis	synthesis	    Androstenedione_from_17alpha-hydroxyprogesterone	Rhea:14753	CYB5A
# CYP17A1	Steroid 17-alpha-hydroxylase/17,20 lyase	                enzyme	    1	        enzyme_synthesis	synthesis	    17alpha-hydroxyprogesterone_from_progesterone	    Rhea:46308


# interactions + cell type