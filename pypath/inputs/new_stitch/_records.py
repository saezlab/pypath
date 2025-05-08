import collections

StitchActions = collections.namedtuple(
    'StitchActions',
    [
        'chemical_id',
        'chemical_acting',
        'is_stereospecific',
        'protein_id',
        'protein_acting',
        'mode',
        'action',
        'score',
        'ncbi_taxa',
    ],
)

ParsedIds = collections.namedtuple(
    'ParsedIds',
    [
        'chemical_id',
        'stereospecific',
        'protein_id',
        'ncbi_taxa',
        'protein_acting',
        'chemical_acting',
    ],
)

StitchLinks = collections.namedtuple(
    'StitchLinks',
    [
        'chemical_id',
        'protein_id',
        'experimental',
        'prediction',
        'database',
        'textmining',
        'combined_score',
    ],
)

StitchInteractions = collections.namedtuple(
    'StitchInteractions',
    [
        'chemical_id',
        'chemical_acting',
        'is_stereospecific',
        'protein_id',
        'protein_acting',
        'mode',
        'action',
        'ncbi_taxa',
        'experimental',
        'prediction',
        'database',
        'textmining',
        'combined_score',
    ],
)