"""
Parse Hormone2Cell data and emit Entity records.

This module converts hormone-cell interaction information from Hormone2Cell into
Entity records using the declarative schema pattern.
"""

from __future__ import annotations
import re
from functools import partial

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
    MembershipBuilder,
    MembersFromList,
    Member,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    InterCellAnnotations,
    MoleculeAnnotationsCv,
    AssayAnnotationsCv,
    BiologicalRoleCv,
    MoleculeSubtypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    ParticipantMetadataCv,
    LicenseCV,
    ProteinFunctionalClassCv,
    OntologyCv,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

URL = 'https://rescued.omnipathdb.org/hormone2cell_s1_to_s6.xlsx'

GROUP_OR_PROHORMONE = re.compile(r'(?:.*_all$)|(?:^pro_)')
# AGG_STR joins all unique values by comma and removes starting/trailing commas
AGG_STR = lambda x: ','.join({str(i) for i in x if i and pd.notna(i)}).strip(',')
TABLES_HORMONES = ['S2B', 'S2D', 'S2E', 'S3A', 'S3B', 'S3C']
TABLES_RECEPTORS = ['S6A', 'S6B', 'S6C']
PARAMS = {
    'S2B': {
        'skiprows': 3,
        'filters': {
            'hormone_ID_unique': ['group_name', 'prohormone'],
            'include': 'NA'
        },
        'pk': 'hormone_short',
        'usecols': [
            'hormone_figures',
            'hormone_display',
            'hormone_other',
            'ref_hormone',
            'ref_receptor'
        ],
    },

# NOTE: Not sure how to parse Table S2D further...
#       The values in the columns named r"hpc_include\d" have different meaning
#       depending on the type of hormone.
#       See map['hormone_type_to_entity'] above for reference:
#       - amine_derived: Enzyme(s) synthesizing the amine-derived hormone
#       - peptide: Enzyme(s) + protein whose cleavage generates the peptides
#       - prostaglandin: Enzyme(s) synthesizing the prostaglandin hormone
#       - protein_monomer: The protein gene
#       - protein_multimer: The genes for the protein monomers of the complex
#       - steroid: Enzyme(s) synthesizing the steroid hormone

    'S2D': {
        'skiprows': 2,
        'filters': {
                'hormone_short': GROUP_OR_PROHORMONE,
                'receptor_known': 'no'
        },
        'pk': 'hormone_short',
        'usecols': ['hormone_type_fine'],
    },
    'S2E': { # Do not aggregate this table contains individual interactions
        'skiprows': 3,
        'filters': {
            'receptor_status': [
                'low_confidence',
                'other_target',
                'receptor_unknown',
                'prohormone'
            ],
        },
    },
    'S3A': {
        'skiprows': 2,
        'filters': {'Hormone_short': GROUP_OR_PROHORMONE},
        'pk': 'Hormone_short',
        'usecols': [
            'Tissue',
            'Tissue level fine-grained cell type (celltype_level2)'
        ],
    },
    'S3B': {
        'skiprows': 3,
        'filters': {'Hormone_short': GROUP_OR_PROHORMONE},
        'pk': 'Hormone_short',
        'usecols': [
            'Tissue',
            'Tissue level broad-grained cell type (celltype_level1)'
        ],
    },
    'S3C': { # Similar to A/B but specific for brain
        'skiprows': 2,
        'filters': {'Hormone_short': GROUP_OR_PROHORMONE},
        'pk': 'Hormone_short',
        'usecols': ['Tissue'],
    },
    'S6A': {
        'skiprows': 2,
        'filters': {'Hormone_short': GROUP_OR_PROHORMONE},
        'pk': 'Hormone_short',
        'usecols': [
            'Tissue',
            'Tissue level fine-grained cell type (celltype_level2)'
        ],
    },
    'S6B': {
        'skiprows': 3,
        'filters': {'Hormone_short': GROUP_OR_PROHORMONE},
        'pk': 'Hormone_short',
        'usecols': [
            'Tissue',
            'Tissue level broad-grained cell type (celltype_level1)'
        ],
    },
    'S6C': { # Similar to A/B but specific for brain
        'skiprows': 3,
        'filters': {'Hormone_short': GROUP_OR_PROHORMONE},
        'pk': 'Hormone_short',
        'usecols': ['Tissue'],
    },
}

config = ResourceConfig(
    id=ResourceCv.HORMONE2CELL,
    name='Hormone2Cell',
    url='https://hormonecellatlas.org.uk/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.STATIC,
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
    filename='hormone2cell_s1_to_s6.xlsx',
    subfolder='hormone2cell',
    large=True,
    ext='.xlsx',
    default_mode='r',
)

def _process_table(
        opener,
        key='',
        skiprows=0,
        filters={},
        pk='hormone_short',
        usecols=[]
    ):
    '''
    Processes a given table and returns an aggregated DataFrame merging rows
    with the same given `pk` (primary key) and values of the other given columns
    merged as defined in `AGG_STR` function. If no `usecols` are provided,
    returns the DataFrame without aggregation (still applies filters if any).

    Applies lowercase to all column names a the end (eases process downstream).

    * Arguments:
        - *opener* [Opener]: The file opener instance
        - *key* [str]: Number of the spreadsheet to process, i.e. "Table [key]".
        - *skiprows* [int]: Number of rows to skip at the start of the
          spreadsheet (the table caption/empty rows)
        - *filters* [dict]: Key value pairs indicating the column (key) and
          values to filter **out** (value) from the table. Values can be passed
          as string (single value to filter out), list of strings or regex
          pattern (for multiple values)
        - *pk* [str]: Primary key to merge rows based on. This is, the column
          name defining the keys with same value that will be merged into a
          single row.
        - *usecols* [list]: List of columns to merge their values based on their
          `pk` using the AGG_STR function.
    '''

    df = pd.read_excel(
        opener.path,
        sheet_name=f'Table {key}',
        skiprows=skiprows
    )

    # Skipping rows according to filter(s) if any
    if filters and isinstance(filters, dict):

        for k, vals in filters.items():

            # Process regex-based filter
            if isinstance(vals, re.Pattern):

                vals = [v for v in df[k].values if vals.match(v)]

            vals = [vals] if not isinstance(vals, list) else vals

            for v in vals:

                df = df.loc[df[k] != v]

    if usecols:

        df = df.groupby(pk).aggregate({k: AGG_STR for k in usecols})
        df.reset_index(inplace=True)

    df.columns = [c.lower() for c in df.columns]

    return df


def _merge_cols(df, newcol, cols):

    df[newcol] = df[cols].agg(AGG_STR, axis=1)
    df.drop(cols, axis=1, inplace=True)

    return df


def _merge_tables(tables):

    result = tables[0]

    for t in tables[1:]:

        commoncols = set(result.columns).intersection(t.columns)
        commoncols.remove('hormone_short')

        result = result.merge(t, how='outer', on='hormone_short')

        # Merging common columns if any
        for c in commoncols:

            result = _merge_cols(result, c, [f'{c}_x', f'{c}_y'])

    return result


def parser(opener):

    tables = {k: _process_table(opener, k, **v) for k, v in PARAMS.items()}
    hormones = _merge_tables([tables[i] for i in TABLES_HORMONES])
    hormones = _merge_cols(
        hormones,
        'cell_types',
        [
            'tissue level fine-grained cell type (celltype_level2)',
            'tissue level broad-grained cell type (celltype_level1)'
        ]
    )
    receptors = _merge_tables([tables[i] for i in TABLES_RECEPTORS])
    receptors = _merge_cols(
        receptors,
        'cell_types',
        [
            'tissue level fine-grained cell type (celltype_level2)',
            'tissue level broad-grained cell type (celltype_level1)'
        ]
    )

    result = hormones.merge(
        receptors,
        how='outer',
        on='hormone_short',
        suffixes=('', '_receptor')
    )

    yield from result.to_dict(orient='records')


# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={
        'pmid': r'/(\d+)/?$'
    },
    map={
        'hormone_type_to_molecule': {
            'amine_derived': EntityTypeCv.SMALL_MOLECULE,
            'peptide': MoleculeSubtypeCv.PEPTIDE,
            'prostaglandin': MoleculeSubtypeCv.LIPID,
            'protein_monomer': EntityTypeCv.PROTEIN,
            'protein_multimer': EntityTypeCv.COMPLEX,
            'steroid': EntityTypeCv.SMALL_MOLECULE,
        },
        'receptor_type': {
            'enzyme_linked': ProteinFunctionalClassCv.CATALYTIC_RECEPTOR,
            'gpcr': ProteinFunctionalClassCv.GPCR,
            'nhr': ProteinFunctionalClassCv.NUCLEAR_HORMONE_RECEPTOR,
            'other': ProteinFunctionalClassCv.OTHER_PROTEIN,
            'other ': ProteinFunctionalClassCv.OTHER_PROTEIN,
        },
    },
    transform={},
)

# NOTE: One of the major challenges for integrating this resource is mapping the
#       hormones to existing entities. For instance, protein-based ones are
#       fine since can map to UniProt but other molecules are not
#       cross-referenced to other resources identifiers or molecular formula,
#       structure, etc... We just have their name which won't be 100% mappable
#       Any other ideas apart from doing it manually and hard-coded, feel free
#       to suggest or implement

schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.HORMONE,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.ABBREVIATED_NAME,
                        value=f('hormone_short')
                    ),
                    CV(
                        term=IdentifierNamespaceCv.NAME,
                        value=f('hormone_display')
                    ),
                    CV(
                        term=IdentifierNamespaceCv.ABBREVIATED_NAME,
                        value=f('hormone_figures')
                    ),
                    CV(
                        term=IdentifierNamespaceCv.SYNONYM,
                        value=f('hormone_other', delimiter=', ')
                    ),
                )
            ),
            annotations=AnnotationsBuilder(
                CV(term=InterCellAnnotations.LIGAND),
                CV(
                    term=ParticipantMetadataCv.PARTICIPANT_XREF,
                    value=f('ref_hormone', delimiter='\n', extract='pmid')
                ),
                CV(term=f('hormone_type_fine', map='hormone_type_to_molecule')),
                CV(term=AssayAnnotationsCv.TISSUE, value=f('tissue')),
                CV(term=AssayAnnotationsCv.CELL_TYPE, value=f('cell_type')),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.RECEPTOR,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                        value=f('receptorgene1')
                    ),
                    CV(
                        term=IdentifierNamespaceCv.NAME,
                        value=f('gene1_name')
                    )
                )
            ),
            annotations=AnnotationsBuilder(
                CV(term=EntityTypeCv.PROTEIN),
                CV(term=InterCellAnnotations.RECEPTOR),
                CV(term=f('receptor_type_broad', map='receptor_type')),
                CV(
                    term=InteractionMetadataCv.PARTICIPANT_XREF,
                    value=f('ref_receptor', delimiter='\n', extract='pmid')
                ),
                CV(term=AssayAnnotationsCv.TISSUE, value=f('tissue_receptor')),
                CV(
                    term=AssayAnnotationsCv.CELL_TYPE,
                    value=f('cell_type_receptor')
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.RECEPTOR,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                        value=f('receptorgene2')
                    ),
                    CV(
                        term=IdentifierNamespaceCv.NAME,
                        value=f('gene2_name')
                    )
                )
            ),
            annotations=AnnotationsBuilder(
                CV(term=EntityTypeCv.PROTEIN),
                CV(term=InterCellAnnotations.RECEPTOR),
                CV(term=f('receptor_type_broad', map='receptor_type')),
                CV(
                    term=InteractionMetadataCv.PARTICIPANT_XREF,
                    value=f('ref_receptor', delimiter='\n', extract='pmid')
                ),
                CV(term=AssayAnnotationsCv.TISSUE, value=f('tissue_receptor')),
                CV(
                    term=AssayAnnotationsCv.CELL_TYPE,
                    value=f('cell_type_receptor')
                ),
            ),
        ),
    )
)

# ================================= RESOURCE ===================================

resource = Resource(
    config=config,
    hormone2cell=Dataset(
        download=download,
        mapper=schema,
        raw_parser=parser
    )
)

# ================================= REFERENCE ==================================
# S2B
#               X                   Y
# hormone_ID	hormone_ID_unique 	hormone_short	hormone_figures	hormone_display	    hormone_other	            Tier	include	exclude	classical	regulatory	gut_pancreatic hormones	neuropeptides	prostaglandins	reprodutive	cardiovascular	other	ref_hormone	                                ref_receptor
# entry_001	    hormone_001	        acth	        ACTH	        ACTH	            Adrenocorticotropic hormone	1	    1	    NA	    1	        1	        NA	                    1	            NA	            NA	        NA	            NA	    https://pubmed.ncbi.nlm.nih.gov/30156493/	https://pubmed.ncbi.nlm.nih.gov/26793988/ https://pubmed.ncbi.nlm.nih.gov/28220105/
# entry_002	    hormone_002	        activin_a	    Activin A	    Activin A	 	                                2	    1	    NA	    NA	        NA	        NA	                    NA	            NA	            1	        NA	            NA	    https://pubmed.ncbi.nlm.nih.gov/27328872/   https://pubmed.ncbi.nlm.nih.gov/3012369/ https://pubmed.ncbi.nlm.nih.gov/30540228/	https://pubmed.ncbi.nlm.nih.gov/1646080/ https://pubmed.ncbi.nlm.nih.gov/15196700/ https://pubmed.ncbi.nlm.nih.gov/19273500/ https://pubmed.ncbi.nlm.nih.gov/26884470/
# entry_003	    hormone_003	        activin_ab	    Activin AB	    Activin AB	 	                                2	    1	    NA	    NA	        NA	        NA	                    NA	            NA	            NA	        NA	            1	    https://pubmed.ncbi.nlm.nih.gov/3086749/    https://pubmed.ncbi.nlm.nih.gov/27328872/	https://pubmed.ncbi.nlm.nih.gov/15196700/ https://pubmed.ncbi.nlm.nih.gov/26884470/
# entry_004	    hormone_004	        activin_b	    Activin B	    Activin B	 	                                2	    1	    NA	    NA	        NA	        NA	                    NA	            NA	            NA	        NA	            1	    https://pubmed.ncbi.nlm.nih.gov/1644823/    https://pubmed.ncbi.nlm.nih.gov/27328872/	https://pubmed.ncbi.nlm.nih.gov/15196700/ https://pubmed.ncbi.nlm.nih.gov/26884470/
# entry_005	    hormone_005	        activin_e	    Activin E	    Activin E	 	                                3	    1	    NA	    NA	        NA	        NA	                    NA	            NA	            NA	        NA	            NA	    https://pubmed.ncbi.nlm.nih.gov/8941337/    https://pubmed.ncbi.nlm.nih.gov/38533769/	https://pubmed.ncbi.nlm.nih.gov/38533769/

# S2D
# Y
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
# Y
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

# S3A
# Y
# Hormone_short	Tier	Tissue	                        Tissue level fine-grained cell type (celltype_level2)	Source Assay	Strength	        Cell count from scRNA-seq	Cell count from snRNA-seq
# cart	        N	    Stomach	                        enteroendocrine_cell	                                cell_only	    1.6101311	        162	                        0
# cart	        N	    Pituitary	                    gonadotroph	                                            cell_only	    0.636176012788019	1055	                    936
# cart	        N	    Cerebral Cortex	                neuron_excitatory_L2negative3_CUX2_[...]	            cell_only	    0.213208141785439	14662	                    146995
# glp2	        2	    Pancreas	                    gamma_cell	                                            cell_only	    5.23447505055456	118	                        797
# glp2	        2	    Large Intestine	                enteroendocrine_cell	                                cell_only	    4.133497	        396	                        0

# S3B
# Y
# Hormone_short	Tier	Tissue	    Tissue level broad-grained cell type (celltype_level1)	Strength_max	    Celltype_Expressing_Count	Celltype_Count
# activin_a	    2	    Adipose	    dendritic_cell	                                        0.287404055017379	1	                        3
# activin_a	    2	    Adipose	    macrophage	                                            0.416829472890448	1	                        4
# activin_a	    2	    Adipose	    mast_cell	                                            0.558139440494919	1	                        1
# activin_a	    2	    Adipose	    monocyte	                                            0.900293639976023	2	                        3
# activin_a	    2	    Adipose	    vascular_smooth_muscle_cell	                            0.265987057185408	1	                        1
# activin_b	    2	    Adipose	    adipocyte	                                            0.126652316578524	3	                        8
# activin_b	    2	    Adipose	    endothelial_cell	                                    1.06644160266406	2	                        2

# S3C
# Y
# Hormone_short	Tier	Tissue	        Anatomical-region–resolved tissue level fine-grained cell type (celltype_level2)	Strength	        Source Assay
# somatostatin	2	    CerebralCortex	neuron_inhibitory_SST_NPY_____FRO	                                                1.24817403520624	cell_only
# somatostatin	2	    CerebralCortex	neuron_inhibitory_SST_ADAMTSL1_____PCx	                                            1.11211415255699	cell_only
# somatostatin	2	    CerebralCortex	neuron_inhibitory_SST_HTR2C_____PCx	                                                1.06945688030573	cell_only
# somatostatin	2	    CerebralCortex	neuron_inhibitory_SST_HS3ST5_____PCx	                                            1.0405526595268	    cell_only
# somatostatin	2	    CerebralCortex	neuron_inhibitory_SST_ADAMTSL1_____TCx	                                            1.02391851305284	cell_only
