"""
Parse NeuronChat data and emit Entity records.

This module converts NeuronChat interaction data into Entity records using the
declarative schema pattern.
"""

from __future__ import annotations
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    MacromolecularComplex,
    ChemicalEntity,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
    Relation,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.neuronchat import iter_neuronchat

config = ResourceConfig(
    id=ResourceCv.NEURONCHAT,
    name='NeuronChat',
    url='https://github.com/Wei-BioMath/NeuronChat',
    license=LicenseCV.GPL_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='36854676',
    primary_category='interactions',
    description='NeuronChat is a manually curated resource of neural-specific intercellular molecular interactions, designed for inferring neuron-neuron communication from single-cell and spatial transcriptomics data.',
)
SMALL_MOLECULE_LIGANDS = {
    '5HT': 'Serotonin',
    'Ach': 'Acetylcholine',
    'CO': 'CO',
    'DA': 'Dopamine',
    'Epi': 'Epinephrine',
    'GABA': 'GABA',
    'Glu': 'Glutamate',
    'Gly': 'Glycine',
    'NE': 'Noradrenaline',
    'NO': 'NO',
}


def _clean(value: object) -> str:
    return str(value or '').strip()


def _list(value: object) -> list[str]:
    if isinstance(value, list):
        return [_clean(item) for item in value if _clean(item)]
    text = _clean(value)
    return [text] if text else []


def _interaction_ligand_token(row: dict[str, object]) -> str:
    return _clean(row.get('interaction_name')).split('_', 1)[0]


def _identifier(type_: object, value: object) -> Identifier | None:
    value = _clean(value)
    return Identifier(type=type_, value=value) if value else None


def _annotation(term: object, value: object = None) -> Annotation | None:
    if value is None:
        return Annotation(term=term)
    value = _clean(value)
    return Annotation(term=term, value=value) if value else None


def _identifiers(*items: Identifier | None) -> list[Identifier]:
    out: list[Identifier] = []
    seen: set[tuple[object, str]] = set()
    for item in items:
        if item is None:
            continue
        key = (item.type, item.value)
        if key in seen:
            continue
        out.append(item)
        seen.add(key)
    return out


def _annotations(*items: Annotation | None) -> list[Annotation] | None:
    out: list[Annotation] = []
    seen: set[tuple[object, object, object]] = set()
    for item in items:
        if item is None:
            continue
        key = (format_term(item.term), format_term(item.value), item.units)
        if key in seen:
            continue
        out.append(item)
        seen.add(key)
    return out or None


def _protein(
    name: object, taxon_id: str, *annotations: Annotation | None
) -> Entity:
    return Entity(
        type=Protein,
        identifiers=_identifiers(_identifier(Namespace.GENESYMBOL, name)),
        annotations=_annotations(
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
            *annotations,
        ),
    )


def _source_entity(row: dict[str, object], taxon_id: str) -> Entity:
    ligand_token = _interaction_ligand_token(row)
    if ligand_token in SMALL_MOLECULE_LIGANDS:
        name = SMALL_MOLECULE_LIGANDS[ligand_token]
        return Entity(
            type=ChemicalEntity,
            identifiers=[
                Identifier(type=Namespace.NAME, value=name),
                Identifier(type=Namespace.SYNONYM, value=ligand_token),
            ],
            annotations=[],
        )
    return _protein(ligand_token, taxon_id)


def _target_entity(row: dict[str, object], taxon_id: str) -> Entity | None:
    subunits = _list(row.get('receptor_subunit'))
    if not subunits:
        return None
    if len(subunits) == 1:
        return _protein(subunits[0], taxon_id)
    name = '+'.join(subunits)
    return Entity(
        type=MacromolecularComplex,
        identifiers=[Identifier(type=Namespace.NAME, value=name)],
        annotations=[
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}')
        ],
        membership=[
            Membership(member=_protein(subunit, taxon_id))
            for subunit in subunits
        ],
    )


def interactions_schema(taxon_id: str):
    def mapper(row: dict[str, object]) -> Relation | None:
        target = _target_entity(row, taxon_id)
        if target is None or not _interaction_ligand_token(row):
            return None
        return Relation(
            subject=_source_entity(row, taxon_id),
            predicate=slots.interacts_with,
            object=target,
            annotations=[*[]],
            identifiers=[
                Identifier(
                    type=Namespace.NAME,
                    value=_clean(row.get('interaction_name')),
                )
            ],
        )

    return mapper


resource = Resource(
    config,
    human_interactions=Dataset(
        download=Download(
            url='https://github.com/Wei-BioMath/NeuronChat/raw/main/data/interactionDB_human.rda',
            filename='neuronchat_interactions_human.rda',
            subfolder='neuronchat',
            default_mode='rb',
            ext='rda',
            encoding=None,
        ),
        mapper=interactions_schema('9606'),
        raw_parser=iter_neuronchat,
    ),
    mouse_interactions=Dataset(
        download=Download(
            url='https://github.com/Wei-BioMath/NeuronChat/raw/main/data/interactionDB_mouse.rda',
            filename='neuronchat_interactions_mouse.rda',
            subfolder='neuronchat',
            default_mode='rb',
            ext='rda',
            encoding=None,
        ),
        mapper=interactions_schema('10090'),
        raw_parser=iter_neuronchat,
    ),
)
