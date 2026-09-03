"""
Parse PMI-DB data and emit Entity records.

This module converts protein-metabolite interaction information from PMI-DB into
Entity records using the declarative schema pattern.
"""

from __future__ import annotations
from functools import partial

from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
)
from pypath.inputs_v2.parsers.base import iter_tsv, _first_handle
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    AssociationBuilder,
    AssociationsBuilder,
    CV,
    Member,
    MembershipBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    OntologyCv,
    AssayAnnotationsCv,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

BASE_URL = 'http://easybioai.com/PMIDB/static/%s.txt'
files_header = {
    'interaction': None, # Already present
    'protein_infor': ['UniProt', 'gene', 'GO', 'link'],
    'metabolites_infor': ['kegg', 'hmdb', 'name', 'link', 'formula'],
}

config = ResourceConfig(
    id=ResourceCv.PMIDB,
    name='PMI-DB',
    url='http://easybioai.com/PMIDB/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='33554247',
    primary_category='interactions',
    annotation_ontologies=(OntologyCv.GENE_ONTOLOGY,),
    description=(
        'PMIDB provides both interactions and non-interactions between protein '
        'and metabolite, which not only reduces the experimental cost for '
        'biological experimenters, but also facilitates the construction of '
        'more accurate algorithms for researchers using machine learning.'
    ),
)

download = {
    f: Download(
        url=BASE_URL % f,
        filename=f'{f}.txt',
        subfolder='PMIDB',
        large=True,
        ext='.txt',
        default_mode='r',
    )
    for f in files_header.keys()
}

def parser(opener, header = None, sep = '\t', **_kwargs):
    '''
    Parses plain text file (TSV by default) to return an iterable of records.
    Accepts custom header for tables without one.
    '''

    if header is not None and iter(header):

        entries = [line.strip().split(sep) for line in _first_handle(opener)]

        if x := len(entries[0]) - len(header):

            raise KeyError(
                'Length of header is %s than the number of entries'
                % ('smaller' if x > 0 else 'larger')
            )

        yield from [
            {k: v for k, v in zip(header, line)}
            for line in entries
        ]

    else: # Interactions file includes header but needs filtering (see below)

        result = iter_tsv(opener)
        # NOTE: Filtering out FALSE-labeled interactions
        yield from [r for r in result if r['interaction'] == 'TRUE']

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={
        'notna': r'\b([a-z0-9]+)\b(?<!NA)',
        'pmid': r'^PMID:(\d+)$',
        'go': r'(GO:\d+)',
        'hmdb': r'(HMDB\d+)'
    },
    map={
        'organism': {
            'human': '9606',
            'mouse': '10090',
            'yeast': '4932', # Presumably
            'Escherichia coli': '562',
        }
    },
    transform={},
)

schema_interaction = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    annotations=AnnotationsBuilder(
        CV(
            term=AssayAnnotationsCv.ASSAY_CATEGORY,
            value=f('probe', extract='notna')
        ),
        CV(
            term=IdentifierNamespaceCv.PUBMED,
            value=f('literature', extract='pmid')
        ),
        CV(term=AssayAnnotationsCv.DESCRIPTION, value=f('Technology')),
        CV(
            term=IdentifierNamespaceCv.NCBI_TAX_ID,
            value=f('Sample source', map='organism')
        ),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.UNIPROT,
                        value=f('Uniprot_KB_id')
                    ),
                ),
            )
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.KEGG_COMPOUND,
                        value=f('KEGG')
                    ),
                ),
            )
        ),
    ),
)

schema_protein_infor = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('UniProt')),
        CV(
            term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
            value=f('gene', delimiter=',')
        ),

    ),
    associations=AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=EntityTypeCv.CV_TERM,
            object_identifier_type=IdentifierNamespaceCv.CV_TERM_ACCESSION,
            object_identifier=f('GO', delimiter=','),
        ),
    ),
)

schema_metabolites_infor = EntityBuilder(
    entity_type=EntityTypeCv.SMALL_MOLECULE,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.KEGG, value=f('kegg')),
        CV(
            term=IdentifierNamespaceCv.HMDB,
            value=f('hmdb', extract='hmdb')
        ),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.SMILES, value=f('formula'))
    ),
)

# ================================= RESOURCE ===================================

resource = Resource(
    config=config,
    **{
        k: Dataset(
            download=download[k],
            mapper=locals().get(f'schema_{k}'),
            raw_parser=partial(parser, header=v),
        )
        for k, v in files_header.items()
    }
)

# ================================= REFERENCE ==================================

# interaction.txt
# KEGG  	Uniprot_KB_id	interaction	probe	mean	literature	    Technology	Sample source	    MS Quality control method	Candidate selection cutoff
# C00002	P00350	        TRUE	    NA	    NA	    PMID:29307493	LiP-SMap	Escherichia coli	DDA and DIA	                FC>2,Q<0.01
# XXX: We discarded those with interaction = FALSE, also not using columns: mean, MS Quality control method and Candidate selection cutoff

# metabolites_infor.txt
# kegg      hmdb        name                        link                                            formula
# C00002	HMDB00538	Adenosine 5'-triphosphate	https://www.genome.jp/dbget-bin/www_bget?C00002	NC1=NC=NC2=C1N=CN2[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O
# XXX: File has no header, added custom one, column link is skipped

# protein_infor.txt
# UniProt   gene    GO  link
# P00350	gnd	    NA	https://www.uniprot.org/uniprot/P00350
# XXX: File has no header, added custom one, column link is skipped
