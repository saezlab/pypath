"""Parse Cellinker ligand-receptor and metabolite-protein data."""

from __future__ import annotations

from collections.abc import Callable
from urllib.parse import quote

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    InterCellAnnotations,
    LicenseCV,
    MoleculeAnnotationsCv,
    ProteinFunctionalClassCv,
    ReactionAnnotationsCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)

BASE_URL = 'https://www.cellknowledge.com.cn/cellinker/download/'

SPECIES = {
    'human': {
        'label': 'Homo sapiens',
        'taxon_id': '9606',
        'protein_prefix': 'LRH',
        'metabolite_prefix': 'MRH',
        'gene_id': 'Human_geneID',
        'gene_name': 'GENE_NAME',
        'uniprot': 'UNIPROT_ID',
    },
    'mouse': {
        'label': 'Mus musculus',
        'taxon_id': '10090',
        'protein_prefix': 'LRM',
        'metabolite_prefix': 'MRM',
        'gene_id': 'Mouse_geneID',
        'gene_name': 'Mouse_gene symbol',
        'uniprot': 'Mouse_uniprot',
    },
}


config = ResourceConfig(
    id=ResourceCv.CELLINKER,
    name='Cellinker',
    url='https://www.cellknowledge.com.cn/cellinker/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33471060',
    primary_category='interactions',
    description=(
        'Cellinker is a literature-curated repository of ligand-receptor '
        'interactions involved in cell-cell communication, including protein '
        'ligands, metabolite ligands, receptors, enzymes and transporters.'
    ),
)


def _download(species: str, dataset: str) -> Download:
    label = SPECIES[species]['label']
    title = {
        'protein_interactions': f'{label} protein L-R interaction',
        'metabolite_interactions': f'{label} metabolite L-R interaction',
        'enzymes': f'{label} enzyme',
        'transporters': f'{label} transporter protein',
    }[dataset]
    return Download(
        url=f'{BASE_URL}{quote(f"{title}.txt")}',
        filename=f'{title}.txt',
        subfolder='cellinker',
        encoding='latin-1',
    )


def _clean(value: object) -> str:
    return str(value or '').strip().strip('"')


def _row_value(row: dict[str, object], *keys: str) -> str:
    for key in keys:
        value = _clean(row.get(key))
        if value:
            return value
    return ''


def _split(value: object, delimiter: str = ';') -> list[str]:
    return [
        item.strip()
        for item in _clean(value).split(delimiter)
        if item.strip()
    ]


def _annotation(term: object, value: object = None) -> Annotation | None:
    value = _clean(value)
    return Annotation(term=term, value=value) if value else None


def _flag(term: object, enabled: bool) -> Annotation | None:
    return Annotation(term=term) if enabled else None


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


def _pubmed_annotations(value: object) -> list[Annotation]:
    return [
        Annotation(term=IdentifierNamespaceCv.PUBMED, value=pmid)
        for pmid in _split(value)
    ]


def _database_annotations(value: object) -> list[Annotation]:
    return [
        Annotation(term=InteractionMetadataCv.INTERACTION_XREF, value=database)
        for database in _split(value)
    ]


def _pubchem_identifiers(value: object) -> list[Identifier]:
    identifiers: list[Identifier] = []
    for token in _split(value):
        if token.startswith('CID:'):
            identifiers.append(
                Identifier(
                    type=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                    value=token.removeprefix('CID:'),
                )
            )
        elif token.startswith('SID:'):
            identifiers.append(
                Identifier(
                    type=IdentifierNamespaceCv.PUBCHEM_SUBSTANCE,
                    value=token.removeprefix('SID:'),
                )
            )
    return identifiers


def _protein(
    *,
    uniprot: object,
    gene_id: object,
    gene_name: object,
    taxon_id: str,
    name: object = None,
    annotations: list[Annotation] | None = None,
) -> Entity:
    return Entity(
        type=EntityTypeCv.PROTEIN,
        identifiers=_identifiers(
            _identifier(IdentifierNamespaceCv.UNIPROT, uniprot),
            _identifier(IdentifierNamespaceCv.ENTREZ, gene_id),
            _identifier(IdentifierNamespaceCv.GENE_NAME_PRIMARY, gene_name),
            _identifier(IdentifierNamespaceCv.NAME, name),
        ),
        annotations=_annotations(
            Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
            *(annotations or []),
        ),
    )


def _metabolite(row: dict[str, object]) -> Entity:
    return Entity(
        type=EntityTypeCv.CHEMICAL,
        identifiers=_identifiers(
            _identifier(IdentifierNamespaceCv.HMDB, _row_value(row, 'HMDB ID', 'HMDB_ID')),
            _identifier(IdentifierNamespaceCv.NAME, _row_value(row, 'Metabolite name', 'METABOLITE_NAME', 'mETABOLITE_NAME')),
            _identifier(IdentifierNamespaceCv.SMILES, row.get('Canonical SMILES')),
            *_pubchem_identifiers(row.get('PubChem CID/SID')),
            _identifier(IdentifierNamespaceCv.MOLECULAR_FORMULA, row.get('Molecular Formula')),
        ),
        annotations=_annotations(
            _annotation(MoleculeAnnotationsCv.COMPOUND_KINGDOM, row.get('Kingdom')),
            _annotation(MoleculeAnnotationsCv.COMPOUND_SUPERCLASS, row.get('Super Class')),
            _annotation(MoleculeAnnotationsCv.COMPOUND_CLASS, row.get('Class')),
        ),
    )


def _cellinker_id(row: dict[str, object], fallback_parts: tuple[object, ...]) -> str:
    direct = _row_value(row, 'LRID', 'MRID')
    if direct:
        return direct
    return ':'.join(_clean(part) or '-' for part in fallback_parts)


def _make_protein_interaction_mapper(
    species: str,
) -> Callable[[dict[str, object]], Entity]:
    taxon_id = SPECIES[species]['taxon_id']

    def mapper(row: dict[str, object]) -> Entity:
        return Entity(
            type=EntityTypeCv.INTERACTION,
            identifiers=_identifiers(
                _identifier(IdentifierNamespaceCv.CELLINKER, row.get('LRID')),
                _identifier(
                    IdentifierNamespaceCv.NAME,
                    f'{_clean(row.get("ligand_symbol"))} - {_clean(row.get("Receptor_symbol"))}',
                ),
            ),
            annotations=_annotations(
                Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
                *_pubmed_annotations(row.get('PMID')),
                *_database_annotations(row.get('Database')),
            ),
            membership=[
                Membership(
                    member=_protein(
                        uniprot=row.get('Ligand_Uniprot'),
                        gene_id=row.get('Ligand_geneid'),
                        gene_name=row.get('ligand_symbol'),
                        taxon_id=taxon_id,
                    ),
                    annotations=[Annotation(term=InterCellAnnotations.LIGAND)],
                ),
                Membership(
                    member=_protein(
                        uniprot=row.get('Receptor_Uniprot'),
                        gene_id=row.get('Receptor_geneid'),
                        gene_name=row.get('Receptor_symbol'),
                        taxon_id=taxon_id,
                    ),
                    annotations=[Annotation(term=InterCellAnnotations.RECEPTOR)],
                ),
            ],
        )

    return mapper


def _make_metabolite_interaction_mapper(
    species: str,
) -> Callable[[dict[str, object]], Entity]:
    taxon_id = SPECIES[species]['taxon_id']

    def mapper(row: dict[str, object]) -> Entity:
        return Entity(
            type=EntityTypeCv.INTERACTION,
            identifiers=_identifiers(
                _identifier(IdentifierNamespaceCv.CELLINKER, row.get('MRID')),
                _identifier(
                    IdentifierNamespaceCv.NAME,
                    f'{_row_value(row, "Metabolite name")} - {_row_value(row, "Receptor_symbol")}',
                ),
            ),
            annotations=_annotations(
                Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
                *_pubmed_annotations(row.get('PMID')),
                *_database_annotations(_row_value(row, 'Other.DB', 'Database')),
            ),
            membership=[
                Membership(
                    member=_metabolite(row),
                    annotations=[Annotation(term=InterCellAnnotations.LIGAND)],
                ),
                Membership(
                    member=_protein(
                        uniprot=_row_value(row, 'Receptor uniprot_ id', 'Receptor_uniprot_ id'),
                        gene_id=_row_value(row, 'Receptor_gene ID', 'Receptor_geneID'),
                        gene_name=row.get('Receptor_symbol'),
                        name=row.get('protein name'),
                        taxon_id=taxon_id,
                    ),
                    annotations=[Annotation(term=InterCellAnnotations.RECEPTOR)],
                ),
            ],
        )

    return mapper


def _make_metabolite_protein_mapper(
    species: str,
    role: str,
) -> Callable[[dict[str, object]], Entity]:
    taxon_id = SPECIES[species]['taxon_id']
    protein_keys = SPECIES[species]

    def protein_role() -> BiologicalRoleCv:
        return (
            BiologicalRoleCv.CONTROLLER
            if role == 'transporter'
            else BiologicalRoleCv.ENZYME
        )

    def metabolite_role(row: dict[str, object]) -> BiologicalRoleCv | None:
        value = _clean(row.get('enzyme product/substrate')).lower()
        if value == 'product':
            return BiologicalRoleCv.PRODUCT
        if value == 'substrate':
            return BiologicalRoleCv.SUBSTRATE
        return None

    def mapper(row: dict[str, object]) -> Entity:
        protein_uniprot = _row_value(row, protein_keys['uniprot'])
        metabolite_id = _row_value(row, 'HMDB_ID')
        cellinker_id = _cellinker_id(row, (role, metabolite_id, protein_uniprot))
        return Entity(
            type=EntityTypeCv.INTERACTION,
            identifiers=_identifiers(
                _identifier(IdentifierNamespaceCv.CELLINKER, cellinker_id),
                _identifier(
                    IdentifierNamespaceCv.NAME,
                    f'{_row_value(row, "METABOLITE_NAME", "mETABOLITE_NAME")} - '
                    f'{_row_value(row, "ENZYME_NAME")}',
                ),
            ),
            annotations=_annotations(
                Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=taxon_id),
                _annotation(ReactionAnnotationsCv.XREF, row.get('REACTIONS')),
            ),
            membership=[
                Membership(
                    member=_metabolite(row),
                    annotations=_annotations(
                        Annotation(term=InterCellAnnotations.LIGAND),
                        _flag(metabolite_role(row), metabolite_role(row) is not None),
                    ),
                ),
                Membership(
                    member=_protein(
                        uniprot=protein_uniprot,
                        gene_id=_row_value(row, protein_keys['gene_id']),
                        gene_name=_row_value(row, protein_keys['gene_name']),
                        name=row.get('ENZYME_NAME'),
                        taxon_id=taxon_id,
                        annotations=[
                            item
                            for item in (
                                _annotation(
                                    MoleculeAnnotationsCv.PROTEIN_FUNCTIONAL_CLASS,
                                    ProteinFunctionalClassCv.TRANSPORTER
                                    if role == 'transporter'
                                    else row.get('type'),
                                ),
                            )
                            if item is not None
                        ],
                    ),
                    annotations=_annotations(
                        Annotation(term=protein_role()),
                        _flag(InterCellAnnotations.RECEPTOR, role == 'transporter'),
                    ),
                ),
            ],
        )

    return mapper


resource = Resource(
    config,
    human_protein_interactions=Dataset(
        download=_download('human', 'protein_interactions'),
        mapper=_make_protein_interaction_mapper('human'),
        raw_parser=iter_tsv,
    ),
    mouse_protein_interactions=Dataset(
        download=_download('mouse', 'protein_interactions'),
        mapper=_make_protein_interaction_mapper('mouse'),
        raw_parser=iter_tsv,
    ),
    human_metabolite_interactions=Dataset(
        download=_download('human', 'metabolite_interactions'),
        mapper=_make_metabolite_interaction_mapper('human'),
        raw_parser=iter_tsv,
    ),
    mouse_metabolite_interactions=Dataset(
        download=_download('mouse', 'metabolite_interactions'),
        mapper=_make_metabolite_interaction_mapper('mouse'),
        raw_parser=iter_tsv,
    ),
    human_enzymes=Dataset(
        download=_download('human', 'enzymes'),
        mapper=_make_metabolite_protein_mapper('human', 'enzyme'),
        raw_parser=iter_tsv,
    ),
    mouse_enzymes=Dataset(
        download=_download('mouse', 'enzymes'),
        mapper=_make_metabolite_protein_mapper('mouse', 'enzyme'),
        raw_parser=iter_tsv,
    ),
    human_transporters=Dataset(
        download=_download('human', 'transporters'),
        mapper=_make_metabolite_protein_mapper('human', 'transporter'),
        raw_parser=iter_tsv,
    ),
    mouse_transporters=Dataset(
        download=_download('mouse', 'transporters'),
        mapper=_make_metabolite_protein_mapper('mouse', 'transporter'),
        raw_parser=iter_tsv,
    ),
)
