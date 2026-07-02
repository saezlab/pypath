"""
Parse NeuronChat data and emit Entity records.

This module converts NeuronChat interaction data into Entity records using the
declarative schema pattern.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    MoleculeSubtypeCv,
    UpdateCategoryCV,
    ResourceCv,
    InteractionMetadataCv,
    ParticipantMetadataCv,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.neuronchat import iter_neuronchat

config = ResourceConfig(
    id=ResourceCv.NEURONCHAT,
    name='NeuronChat',
    url='https://github.com/Wei-BioMath/NeuronChat',
    license=LicenseCV.GPL_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='36854676', # Zhao et al. 2023 DOI: 10.1038/s41467-023-36800-w
    primary_category='interactions',
    description=(
        'NeuronChat is a manually curated resource of neural-specific '
        'intercellular molecular interactions, designed for inferring '
        'neuron-neuron communication from single-cell and spatial '
        'transcriptomics data.'
    ),
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


def _ligand_name(row: dict[str, object]) -> str:
    token = _interaction_ligand_token(row)
    return SMALL_MOLECULE_LIGANDS.get(token, token)


def _source_label(row: dict[str, object]) -> str:
    return f'{_clean(row.get("interaction_name"))}_source'


def _target_label(row: dict[str, object]) -> str:
    return f'{_clean(row.get("interaction_name"))}_target'


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
        key = (item.term, item.value, item.units)
        if key in seen:
            continue
        out.append(item)
        seen.add(key)
    return out or None


def _protein(name: object, taxon_id: str, *annotations: Annotation | None) -> Entity:
    return Entity(
        type=EntityTypeCv.PROTEIN,
        identifiers=_identifiers(
            _identifier(IdentifierNamespaceCv.GENE_NAME_PRIMARY, name),
        ),
        annotations=_annotations(
            Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
            *annotations,
        ),
    )


def _source_entity(row: dict[str, object], taxon_id: str) -> Entity:
    ligand_token = _interaction_ligand_token(row)
    if ligand_token in SMALL_MOLECULE_LIGANDS:
        return Entity(
            type=EntityTypeCv.CHEMICAL,
            identifiers=_identifiers(
                _identifier(IdentifierNamespaceCv.NAME, _ligand_name(row)),
                _identifier(IdentifierNamespaceCv.ABBREVIATED_NAME, ligand_token),
            ),
            annotations=_annotations(
                Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
                Annotation(
                    term=MoleculeAnnotationsCv.MOLECULE_SUBTYPE,
                    value=MoleculeSubtypeCv.METABOLITE,
                ),
                _annotation(MoleculeAnnotationsCv.SOURCE_STATUS, _source_label(row)),
            ),
        )

    return _protein(
        ligand_token,
        taxon_id,
        _annotation(MoleculeAnnotationsCv.SOURCE_STATUS, _source_label(row)),
    )


def _target_entity(row: dict[str, object], taxon_id: str) -> Entity:
    subunits = _list(row.get('receptor_subunit'))
    source_status = _annotation(MoleculeAnnotationsCv.SOURCE_STATUS, _target_label(row))
    if len(subunits) == 1:
        return _protein(subunits[0], taxon_id, source_status)

    return Entity(
        type=EntityTypeCv.COMPLEX,
        identifiers=_identifiers(
            _identifier(IdentifierNamespaceCv.NAME, _target_label(row)),
        ),
        annotations=_annotations(source_status),
        membership=[
            Membership(member=_protein(subunit, taxon_id))
            for subunit in subunits
        ] or None,
    )


def interactions_schema(taxon_id: str):
    def mapper(row: dict[str, object]) -> Entity:
        return Entity(
            type=EntityTypeCv.INTERACTION,
            identifiers=_identifiers(
                _identifier(IdentifierNamespaceCv.NAME, row.get('interaction_name')),
            ),
            annotations=_annotations(
                _annotation(InteractionMetadataCv.LIGAND_TYPE, row.get('ligand_type')),
                _annotation(InteractionMetadataCv.INTERACTION_TYPE, row.get('interaction_type')),
            ),
            membership=[
                Membership(
                    member=_source_entity(row, taxon_id),
                    annotations=[Annotation(term=ParticipantMetadataCv.SOURCE)],
                ),
                Membership(
                    member=_target_entity(row, taxon_id),
                    annotations=[Annotation(term=ParticipantMetadataCv.TARGET)],
                ),
            ],
        )

    return mapper


resource = Resource(
    config,
    human_interactions=Dataset(
        download=Download(
            url=f'https://github.com/Wei-BioMath/NeuronChat/raw/main/data/interactionDB_human.rda',
            filename=f'neuronchat_interactions_human.rda',
            subfolder='neuronchat',
            default_mode='rb',
            ext='rda',
            encoding=None, # avoid encoding in binary mode
        ),
        mapper=interactions_schema('9606'),
        raw_parser=iter_neuronchat,
    ),
    mouse_interactions=Dataset(
        download=Download(
            url=f'https://github.com/Wei-BioMath/NeuronChat/raw/main/data/interactionDB_mouse.rda',
            filename=f'neuronchat_interactions_mouse.rda',
            subfolder='neuronchat',
            default_mode='rb',
            ext='rda',
            encoding=None, # avoid encoding in binary mode
        ),
        mapper=interactions_schema('10090'),
        raw_parser=iter_neuronchat,
    ),
)
