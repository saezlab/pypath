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