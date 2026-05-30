"""Parse DrugCentral drug-target interaction data and emit Entity records."""

from __future__ import annotations

import csv
from functools import cache

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.share.downloads import download_and_open
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)

DRUGCENTRAL_STRUCTURES_URL = (
    'https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/'
    'structures.smiles.tsv'
)

config = ResourceConfig(
    id=ResourceCv.DRUGCENTRAL,
    name='DrugCentral',
    url='https://drugcentral.org/',
    license=LicenseCV.CC_BY_SA_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33151290',
    primary_category='interactions',
    description=(
        'DrugCentral is an online drug information resource with active '
        'ingredients, chemical entities, pharmacologic action, indications, '
        'mechanism of action and drug-target interaction data.'
    ),
)


def _clean(value: object) -> str:
    return str(value or '').strip().strip('"')


def _values(value: object, delimiter: str = '|') -> list[str]:
    return [
        item.strip()
        for item in _clean(value).split(delimiter)
        if item.strip()
    ]


def _annotation(
    term: object,
    value: object = None,
    units: object = None,
) -> Annotation | None:
    value = _clean(value)
    units = _clean(units)
    if not value:
        return None
    return Annotation(term=term, value=value, units=units or None)


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


def _identifiers(*items: Identifier | None) -> list[Identifier]:
    out: list[Identifier] = []
    seen: set[tuple[object, str]] = set()
    for item in items:
        if item is None or not _clean(item.value):
            continue
        key = (item.type, _clean(item.value))
        if key in seen:
            continue
        out.append(item)
        seen.add(key)
    return out


def _identifier(term: object, value: object) -> Identifier | None:
    value = _clean(value)
    return Identifier(type=term, value=value) if value else None


@cache
def _structure_inchikeys() -> dict[str, str]:
    opener = download_and_open(
        url=DRUGCENTRAL_STRUCTURES_URL,
        filename='structures.smiles.tsv',
        subfolder='drugcentral',
        large=True,
        default_mode='r',
    )
    try:
        return {
            _clean(row.get('ID')): _clean(row.get('InChIKey'))
            for row in csv.DictReader(opener.result, delimiter='\t')
            if _clean(row.get('ID')) and _clean(row.get('InChIKey'))
        }
    finally:
        opener.close()


def _interaction_id(row: dict[str, object]) -> str:
    fields = (
        row.get('STRUCT_ID'),
        row.get('ACCESSION'),
        row.get('ACT_TYPE'),
        row.get('ACTION_TYPE'),
        row.get('ACT_SOURCE'),
    )
    return ':'.join(_clean(value) or '-' for value in fields)


def _taxon_id(row: dict[str, object]) -> str:
    return {
        'Homo sapiens': '9606',
        'Mus musculus': '10090',
        'Rattus norvegicus': '10116',
    }.get(_clean(row.get('ORGANISM')), _clean(row.get('ORGANISM')))


def _organism_annotation(row: dict[str, object]) -> Annotation | None:
    taxon_id = _taxon_id(row)
    if taxon_id.isdigit():
        return Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id)
    return _annotation(IdentifierNamespaceCv.SCIENTIFIC_NAME, taxon_id)


def _small_molecule(row: dict[str, object]) -> Entity:
    struct_id = _clean(row.get('STRUCT_ID'))
    return Entity(
        type=EntityTypeCv.CHEMICAL,
        identifiers=_identifiers(
            _identifier(IdentifierNamespaceCv.DRUGCENTRAL, struct_id),
            _identifier(
                IdentifierNamespaceCv.STANDARD_INCHI_KEY,
                _structure_inchikeys().get(struct_id),
            ),
            _identifier(IdentifierNamespaceCv.NAME, row.get('DRUG_NAME')),
        ),
    )


def _protein(
    row: dict[str, object],
    *,
    accession: object,
    gene: object = None,
    swissprot: object = None,
    name: object = None,
) -> Entity:
    return Entity(
        type=EntityTypeCv.PROTEIN,
        identifiers=_identifiers(
            _identifier(IdentifierNamespaceCv.UNIPROT, accession),
            _identifier(IdentifierNamespaceCv.GENE_NAME_PRIMARY, gene),
            _identifier(IdentifierNamespaceCv.UNIPROT_ENTRY_NAME, swissprot),
            _identifier(IdentifierNamespaceCv.NAME, name),
        ),
        annotations=_annotations(
            _organism_annotation(row),
            _annotation(InteractionMetadataCv.INTERACTION_ANNOTATION, row.get('TARGET_CLASS')),
            _annotation(InteractionMetadataCv.INTERACTION_ANNOTATION, row.get('TDL')),
        ),
    )


def _aligned(values: list[str], index: int) -> str:
    return values[index] if index < len(values) else ''


def _target_group_type(row: dict[str, object]) -> EntityTypeCv:
    name = _clean(row.get('TARGET_NAME')).lower()
    if len(_values(row.get('ACCESSION'))) < 2:
        return EntityTypeCv.PROTEIN
    if (
        'complex' in name
        or '/' in name
        or ('receptor' in name and ',' in name)
        or 'urease' in name
        or 'transporting atpase' in name
    ):
        return EntityTypeCv.COMPLEX
    return EntityTypeCv.PROTEIN_FAMILY


def _target_entity(row: dict[str, object]) -> Entity:
    accessions = _values(row.get('ACCESSION'))
    genes = _values(row.get('GENE'))
    swissprots = _values(row.get('SWISSPROT'))
    target_name = _clean(row.get('TARGET_NAME'))

    if len(accessions) < 2:
        return _protein(
            row,
            accession=accessions[0] if accessions else '',
            gene=genes[0] if genes else '',
            swissprot=swissprots[0] if swissprots else '',
            name=target_name,
        )

    return Entity(
        type=_target_group_type(row),
        identifiers=_identifiers(
            _identifier(IdentifierNamespaceCv.NAME, target_name),
        ),
        annotations=_annotations(
            _organism_annotation(row),
            _annotation(InteractionMetadataCv.INTERACTION_ANNOTATION, row.get('TARGET_CLASS')),
            _annotation(InteractionMetadataCv.INTERACTION_ANNOTATION, row.get('TDL')),
        ),
        membership=[
            Membership(
                member=_protein(
                    row,
                    accession=accession,
                    gene=_aligned(genes, index),
                    swissprot=_aligned(swissprots, index),
                ),
            )
            for index, accession in enumerate(accessions)
        ],
    )


def map_drugcentral_interaction(row: dict[str, object]) -> Entity:
    """Map a DrugCentral raw row to a drug-target interaction entity."""
    act_type = _clean(row.get('ACT_TYPE'))
    act_value = _clean(row.get('ACT_VALUE'))
    act_unit = _clean(row.get('ACT_UNIT'))
    activity = '='.join(item for item in (act_type, act_value) if item)

    return Entity(
        type=EntityTypeCv.INTERACTION,
        identifiers=_identifiers(
            _identifier(IdentifierNamespaceCv.DRUGCENTRAL, _interaction_id(row)),
            _identifier(
                IdentifierNamespaceCv.NAME,
                f'{_clean(row.get("DRUG_NAME"))} - {_clean(row.get("TARGET_NAME"))}',
            ),
        ),
        annotations=_annotations(
            _organism_annotation(row),
            _annotation(InteractionMetadataCv.INTERACTION_PARAMETER, activity, act_unit),
            _annotation(InteractionMetadataCv.INTERACTION_ANNOTATION, row.get('ACT_COMMENT')),
            _annotation(InteractionMetadataCv.INTERACTION_ANNOTATION, row.get('ACT_SOURCE')),
            _annotation(InteractionMetadataCv.INTERACTION_ANNOTATION, row.get('RELATION')),
            _annotation(InteractionMetadataCv.INTERACTION_ANNOTATION, row.get('MOA_SOURCE')),
            _annotation(InteractionMetadataCv.INTERACTION_XREF, row.get('ACT_SOURCE_URL')),
            _annotation(InteractionMetadataCv.INTERACTION_XREF, row.get('MOA_SOURCE_URL')),
            _annotation(InteractionMetadataCv.CONTROL_TYPE, row.get('ACTION_TYPE')),
            _annotation(
                InteractionMetadataCv.INTERACTION_ANNOTATION,
                'mechanism_of_action' if _clean(row.get('MOA')) == '1' else None,
            ),
        ),
        membership=[
            Membership(member=_small_molecule(row)),
            Membership(member=_target_entity(row)),
        ],
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz',
            filename='drug.target.interaction.tsv.gz',
            subfolder='drugcentral',
            ext='gz',
        ),
        mapper=map_drugcentral_interaction,
        raw_parser=iter_tsv,
    ),
)
