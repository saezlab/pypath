import collections

from . import _raw

BindingdbInteraction = colletions.namedtuple(
    "BindingdbInteraction",
    [
        "ligand",
        "target",
        "action",
        "action_type",
        "is_stimulation",
        "is_inhibition",
        "endogenous",
        "affinity_high",
        "affinity_low",
        "affinity_median",
        "affinity_units",
        "primary_target",
        "pubmed",
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

