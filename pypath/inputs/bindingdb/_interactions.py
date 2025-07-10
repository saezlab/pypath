import collections
from collections.abc import Generator

from . import _raw

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
        "uniprot",
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


def interactions(dataset: str = 'All',
                 max_lines: int | None = None) -> Generator[BindingdbInteraction]:
    
    uniprot_mapping = _raw.mapping()

    for record in _raw.table(dataset = dataset, max_lines = max_lines):



        yield BindingdbInteraction(
            ligand = BindingdbLigand(
                name = record['BindingDB Ligand Name'],
                smiles = record['Ligand SMILES'],
                inchi = record['Ligand InChI'],
                inchi_key = record['Ligand InChI Key'],
                pubchem = record['PubChem CID'],
            ),
            target = BindingdbTarget(
                name = record['Target Name'],
                organism = record['Target Source Organism According to Curator or DataSource'],
                uniprot = uniprot_mapping.get(record['Target Name']), #TODO : create regex split to access all uniprots
            ),
        )
