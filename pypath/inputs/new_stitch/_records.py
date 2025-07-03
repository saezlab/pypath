import collections

StitchAction = collections.namedtuple(
    'StitchAction',
    [
        'source',
        'target',
        'directed',
        'mode',
        'activation',
        'inhibition',
        'score',
    ],
)

Entity = collections.namedtuple(
    'Entity',
    [
        'id',
        'type',
        'stereospecific',
        'ncbi_tax_id',
    ],
)
Entity.__new__.__defaults__ = (None, None)


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
        'ncbi_tax_id',
        'stereospecific',
    ],
)

StitchInteractions = collections.namedtuple(
    'StitchInteractions',
    [
        'source',
        'target',
        'directed',
        'mode',
        'activation',
        'inhibition',
        'experimental',
        'prediction',
        'database',
        'textmining',
        'combined_score',
        'final_score',
    ],
)
