import collections

__all__ = [
    'MrclinksdbRaw',
    'MrclinksdbInteraction',
    'MrclinksdbMetaboliteCell',
]


MrclinksdbRaw = collections.namedtuple(
    'MrclinksdbRaw',
    [
        'mrid',
        'hmdb_id',
        'metabolite_name',
        'pubchem_cid_sid',
        'molecular_formula',
        'kingdom',
        'super_class_',
        'class_',
        'canonical_smiles',
        'receptor_gene_id',
        'receptor_uniprot_id',
        'receptor_symbol',
        'protein_name',
        'pmid',
        'other_db_',
    ]
)

MrclinksdbInteraction = collections.namedtuple(
    'MrclinksdbInteraction',
    [
        'mrid',
        'hmdb',
        'name',
        'pubchem',
        'pubchem_sid',
        'formula',
        'compound_kingdom',
        'compound_superclass',
        'compound_class',
        'smiles',
        'receptor_entrez',
        'receptor_uniprot',
        'receptor_genesymbol',
        'pmids',
        'resource',
        'receptor_location',
    ],
)

MrclinksdbMetaboliteCell = collections.namedtuple(
    'MrclinksdbMetaboliteCell',
    [
        'hmdb',
        'metabolite',
        'interaction',
        'cell_type',
        'experimental_subject',
        'disease',
        'effect',
        'pmids',
    ],
)
