"""Parse Cellinker ligand-receptor and metabolite-protein data."""

from __future__ import annotations
from collections.abc import Callable
from urllib.parse import quote
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import Protein, ChemicalEntity, slots
from omnipath_core.naming import Namespace
from omnipath_core.silver_schema import format_term
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
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
    description='Cellinker is a literature-curated repository of ligand-receptor interactions involved in cell-cell communication, including protein ligands, metabolite ligands, receptors, enzymes and transporters.',
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
        item.strip() for item in _clean(value).split(delimiter) if item.strip()
    ]


def _annotation(term: object, value: object = None) -> Annotation | None:
    value = _clean(value)
    return Annotation(term=term, value=value) if value else None


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
        Annotation(term=slots.publications, value=f'PMID:{pmid}')
        for pmid in _split(value)
    ]


def _pubchem_identifiers(value: object) -> list[Identifier]:
    identifiers: list[Identifier] = []
    for token in _split(value):
        if token.startswith('CID:'):
            identifiers.append(
                Identifier(
                    type=Namespace.PUBCHEM, value=token.removeprefix('CID:')
                )
            )
        elif token.startswith('SID:'):
            identifiers.append(
                Identifier(
                    type=Namespace.PUBCHEM_SUBSTANCE,
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
) -> Entity | None:
    identifiers = _identifiers(
        _identifier(Namespace.UNIPROT, uniprot),
        _identifier(Namespace.ENTREZ, gene_id),
        _identifier(Namespace.GENESYMBOL, gene_name),
    )
    if not identifiers:
        identifiers = _identifiers(_identifier(Namespace.NAME, name))
    if not identifiers:
        return None
    return Entity(
        type=Protein,
        identifiers=_identifiers(
            *identifiers, _identifier(Namespace.NAME, name)
        ),
        annotations=_annotations(
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
            *(annotations or []),
        ),
    )


def _metabolite(row: dict[str, object]) -> Entity:
    return Entity(
        type=ChemicalEntity,
        identifiers=_identifiers(
            _identifier(Namespace.HMDB, _row_value(row, 'HMDB ID', 'HMDB_ID')),
            _identifier(Namespace.SMILES, row.get('Canonical SMILES')),
            *_pubchem_identifiers(row.get('PubChem CID/SID')),
            _identifier(
                Namespace.NAME,
                _row_value(
                    row, 'Metabolite name', 'METABOLITE_NAME', 'mETABOLITE_NAME'
                ),
            ),
        ),
        annotations=_annotations(
            _annotation(
                slots.has_chemical_formula, row.get('Molecular Formula')
            )
        ),
    )


def _make_protein_interaction_mapper(
    species: str,
) -> Callable[[dict[str, object]], Relation | None]:
    taxon_id = SPECIES[species]['taxon_id']

    def mapper(row: dict[str, object]) -> Relation | None:
        ligand = _protein(
            uniprot=row.get('Ligand_Uniprot'),
            gene_id=row.get('Ligand_geneid'),
            gene_name=row.get('ligand_symbol'),
            taxon_id=taxon_id,
        )
        receptor = _protein(
            uniprot=row.get('Receptor_Uniprot'),
            gene_id=row.get('Receptor_geneid'),
            gene_name=row.get('Receptor_symbol'),
            taxon_id=taxon_id,
        )
        if ligand is None or receptor is None:
            return None
        return Relation(
            subject=ligand,
            predicate=slots.interacts_with,
            object=receptor,
            identifiers=_identifiers(
                _identifier(Namespace.CELLINKER, row.get('LRID')),
                _identifier(
                    Namespace.NAME,
                    f'{_clean(row.get("ligand_symbol"))} - {_clean(row.get("Receptor_symbol"))}',
                ),
            ),
            annotations=_annotations(
                Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
                *_pubmed_annotations(row.get('PMID')),
                *[
                    Annotation(term=slots.supporting_data_source, value=v)
                    for v in _split(_row_value(row, 'Other.DB', 'Database'))
                ],
            ),
        )

    return mapper


def _make_metabolite_interaction_mapper(
    species: str,
) -> Callable[[dict[str, object]], Relation | None]:
    taxon_id = SPECIES[species]['taxon_id']

    def mapper(row: dict[str, object]) -> Relation | None:
        metabolite = _metabolite(row)
        receptor = _protein(
            uniprot=_row_value(
                row, 'Receptor uniprot_ id', 'Receptor_uniprot_ id'
            ),
            gene_id=_row_value(row, 'Receptor_gene ID', 'Receptor_geneID'),
            gene_name=row.get('Receptor_symbol'),
            name=row.get('protein name'),
            taxon_id=taxon_id,
        )
        if receptor is None or not metabolite.identifiers:
            return None
        return Relation(
            subject=metabolite,
            predicate=slots.interacts_with,
            object=receptor,
            identifiers=_identifiers(
                _identifier(Namespace.CELLINKER, row.get('MRID')),
                _identifier(
                    Namespace.NAME,
                    f'{_row_value(row, "Metabolite name")} - {_row_value(row, "Receptor_symbol")}',
                ),
            ),
            annotations=_annotations(
                Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
                *_pubmed_annotations(row.get('PMID')),
                *[
                    Annotation(term=slots.supporting_data_source, value=v)
                    for v in _split(_row_value(row, 'Other.DB', 'Database'))
                ],
            ),
        )

    return mapper


def _make_metabolite_protein_mapper(
    species: str, role: str
) -> Callable[[dict[str, object]], Relation | None]:
    taxon_id = SPECIES[species]['taxon_id']
    protein_keys = SPECIES[species]

    def mapper(row: dict[str, object]) -> Relation | None:
        protein_uniprot = _row_value(row, protein_keys['uniprot'])
        metabolite = _metabolite(row)
        protein = _protein(
            uniprot=protein_uniprot,
            gene_id=_row_value(row, protein_keys['gene_id']),
            gene_name=_row_value(row, protein_keys['gene_name']),
            name=row.get('ENZYME_NAME'),
            taxon_id=taxon_id,
        )
        if protein is None or not metabolite.identifiers:
            return None
        predicate = slots.associated_with
        return Relation(
            subject=metabolite,
            predicate=predicate,
            object=protein,
            identifiers=_identifiers(
                _identifier(
                    Namespace.CELLINKER, _row_value(row, 'LRID', 'MRID')
                ),
                _identifier(
                    Namespace.NAME,
                    f'{_row_value(row, "METABOLITE_NAME", "mETABOLITE_NAME")} - {_row_value(row, "ENZYME_NAME")}',
                ),
            ),
            annotations=_annotations(
                Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}')
            ),
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
