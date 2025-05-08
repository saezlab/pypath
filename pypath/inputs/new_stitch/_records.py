import collections

StitchAction = collections.namedtuple(
    'StitchAction',
    [
        'source',
        'target',
        'directed',
        'mode',
        'activation',
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
        'final_score',
    ],
)
