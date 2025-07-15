import collections

all = [
    "BindingdbInteraction",
    "BindingdbLigand",
    "BindingdbTarget",
    "AllostericRegulation",
    "ReactionConstant",
]

BindingdbInteraction = collections.namedtuple(
    "BindingdbInteraction",
    [
        "ligand",
        "target",
    ],
)

BindingdbLigand = collections.namedtuple(
    "BindingdbLigand",
    [
        "name",
        "smiles",
        "inchi",
        "inchi_key",
        "pubchem",
    ]
)

BindingdbTarget = collections.namedtuple(
    "BindingdbTarget",
    [
        "name",
        "organism",
        "ncbi_tax_id",
        "uniprot",
        "regions_mutations",
    ]
)

AllostericRegulation = collections.namedtuple(
    'AllostericRegulation',
    [
        'action',
        'compound',
        'organism',
        'protein',
        'id_type',
        'wrong_ec',
        'pubmeds',
        'reaction_constants',
    ]
)

ReactionConstant = collections.namedtuple(
    'ReactionConstant',
    [
        'type',
        'value',
        'conditions',
        'pubmeds',
    ]
)