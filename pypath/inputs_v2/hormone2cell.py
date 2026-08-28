"""
Parse Hormone2Cell data and emit Entity records.

This module converts hormone-cell interaction information from Hormone2Cell into
Entity records using the declarative schema pattern.
"""

from __future__ import annotations
import re

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

SPLITTABLE = re.compile(r'^Table (S\d\D)(?:-(.*):(.*$))?')
GROUP_OR_PROHORMONE = re.compile(r'(?:.*_all$)|(?:^pro_)')

# NOTE: Some tables need to be processed separately in the parser as they have
#       different schemas depending in the type of entity. Therefore we define
#       the format [tablename]-[columnname]:[value]
#       See parser() for more info on how these parameters below are used
sheet_skiprows_filter = {
    'Table S2B': {
        'skiprows': 3,
        'filter': {'hormone_ID_unique': 'group_name', 'include': 'NA'}
    },
    'Table S2D': {
        'skiprows': 2,
        'filter': {
                'hormone_short': GROUP_OR_PROHORMONE,
                'receptor_known': 'no'
        }
    },
    'Table S2E': {
        'skiprows': 3,
        'filter': {
            'receptor_status': [
                'low_confidence',
                'other_target',
                'receptor_unknown',
                'prohormone'
            ],
        }
    },
    'Table S3A': {
        'skiprows': 2,
        'filter': {'Hormone_short': GROUP_OR_PROHORMONE}
    },
    'Table S3B': {'skiprows': 3, 'filter': {}},
    'Table S3C': {'skiprows': 2, 'filter': {}},
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

def parser(opener, key, skiprows=0, filters={}, merge={}, sep=','):
    '''
    Parses a given table and returns dict of entries for the schema to process

    * Arguments:
        - *opener* [Opener]: The file opener instance
        - *key* [str]: Name of the spreadsheet to process. May include a column
          and value to subset the table like "[tablename]-[columnname]:[value]"
          This will generate a subtable where only the rows whose values are
          [value] in the given [columnname] are returned
        - *skiprows* [int]: Number of rows to skip at the start of the
          spreadsheet (the table caption/empty rows)
        - *filters* [dict]: Key value pairs indicating the column (key) and
          values to filter **out** (value) from the table. Values can be passed
          as string (single value to filter out), list of strings or regex
          pattern (for multiple values)
        - *merge* [dict]: Key value pairs indicating the name of the resulting
          merged column (key) and list of columns to merge (value). The
          resulting merged values are the result of joining them by `sep`
        - *sep* [str]: Separator to join the merged columns, defaults to ","
    '''

    table, column, value = SPLITTABLE.findall(key)[0]

    df = pd.read_excel(
        opener.path,
        sheet_name=f'Table {table}',
        skiprows=skiprows
    )

    # Subsetting into subtable if any
    if column and value:

        df.query(f'{column} == "{value}"', inplace=True)

    # Skipping rows according to filter(s) if any
    if filters and isinstance(filters, dict):

        for k, vals in filters.items():

            # Process regex-based filter
            if isinstance(vals, re.Pattern):

                vals = [v for v in df[k].values if vals.match(v)]

            vals = [vals] if not isinstance(vals, list) else vals

            for v in vals:

                df.query(f'{k} != "{v}"', inplace=True)

    # Merging columns if any
    if merge and isinstance(merge, dict):

        for newcol, cols in merge.items():

            df[newcol] = df[cols].agg(sep.join, axis=1)
            df.drop(cols, axis=1, inplace=True)

    yield from df.to_dict(orient='records')

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={
        'pmid': r'/(\d+)/?$'
    },
    map={
        'hormone_type_to_entity': {
            'amine_derived': EntityTypeCv.REACTION,
            'peptide': EntityTypeCv.REACTION,
            'prostaglandin': EntityTypeCv.REACTION,
            'protein_monomer': EntityTypeCv.PROTEIN,
            'protein_multimer': EntityTypeCv.COMPLEX,
            'steroid': EntityTypeCv.REACTION,
        },
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
schema_S2B = EntityBuilder(
    entity_type=EntityTypeCv.HORMONE,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.H2C_ID, value=f('hormone_short')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('hormone_display')),
        CV(
            term=IdentifierNamespaceCv.ABBREVIATED_NAME,
            value=f('hormone_figures')
        ),
        CV(
            term=IdentifierNamespaceCv.SYNONYM,
            value=f('hormone_other', delimiter=', ')
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=ParticipantMetadataCv.PARTICIPANT_XREF,
            value=f('ref_hormone', delimiter='\n', extract='pmid')
        ),
        CV(
            term=InteractionMetadataCv.INTERACTION_XREF,
            value=f('ref_receptor', delimiter='\n', extract='pmid')
        )
    ),
)

# NOTE: Not sure how to parse this one further...
#       The values in the columns named r"hpc_include\d" have different meaning
#       depending on the type of hormone.
#       See map['hormone_type_to_entity'] above for reference:
#       - amine_derived: Enzyme(s) synthesizing the amine-derived hormone
#       - peptide: Enzyme(s) + protein whose cleavage generates the peptides
#       - prostaglandin: Enzyme(s) synthesizing the prostaglandin hormone
#       - protein_monomer: The protein gene
#       - protein_multimer: The genes for the protein monomers of the complex
#       - steroid: Enzyme(s) synthesizing the steroid hormone
schema_S2D = EntityBuilder(
    entity_type=f('hormone_type_fine', map='hormone_type_to_molecule'),
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.H2C_ID, value=f('hormone_short')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=EntityTypeCv.HORMONE)
    ),
)

schema_S2E = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.HORMONE,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.H2C_ID,
                        value=f('hormone_short')
                    )
                )
            ),
            annotations=AnnotationsBuilder(
                CV(term=InterCellAnnotations.LIGAND),
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
                CV(term=InterCellAnnotations.RECEPTOR),
                CV(term=f('receptor_type_broad', map='receptor_type'))
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
                CV(term=InterCellAnnotations.RECEPTOR),
                CV(term=f('receptor_type_broad', map='receptor_type'))
            ),
        ),
    )
)

schema_S3A = EntityBuilder(
    entity_type=
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
