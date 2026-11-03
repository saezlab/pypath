"""
Parse CellChatDB data and emit Entity records.

This module converts CellChatDB ligand-receptor interactions, complexes and
cofactor groups into Entity records.
"""

from __future__ import annotations
from collections.abc import Generator
import functools
import re
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    ChemicalEntity,
    Protein,
    NamedThing,
    MacromolecularComplex,
    slots,
    DirectionQualifierEnum,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from omnipath_core.source_attributes import PARTICIPANT_ROLE
from omnipath_core.silver_schema import format_term
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
    Relation,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cellchat import (
    NON_PROTEIN,
    iter_cellchat_cofactors,
    iter_cellchat_cofactor_interactions,
    iter_cellchat_complexes,
    iter_cellchat_interactions,
)

config = ResourceConfig(
    id=ResourceCv.CELLCHAT,
    name='CellChat',
    full='CellChatDB',
    url='https://github.com/jinworks/CellChat',
    license=LicenseCV.GPL_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33597522',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS, Content.PROTEIN_COMPLEXES),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Ligand–receptor interactions, protein complexes and signaling cofactors from CellChatDB.',
    description='CellChatDB is a manually curated database of literature-supported cell-cell communication interactions for CellChat, including ligand-receptor pairs, multi-subunit complexes and signaling cofactors across multiple species.',
)
BASE_URL = 'https://raw.githubusercontent.com/jinworks/CellChat/main/data/'
SPECIES = {'human': '9606', 'mouse': '10090', 'zebrafish': '7955'}


def _download(species: str) -> Download:
    return Download(
        url=f'{BASE_URL}CellChatDB.{species}.rda',
        filename=f'cellchatdb_{species}.rda',
        subfolder='cellchat',
        default_mode='rb',
        encoding=None,
    )


def _parser(parser: Any, species: str) -> Any:
    return functools.partial(parser, taxon_id=SPECIES[species])


# CellChat names a metabolite ligand by its abbreviation followed by the genes it
# uses as a proxy for the metabolite: the enzymes that make it and the
# transporters that release it ("5-HT-TPH1_SLC18A1", "PGE2-PTGES"), and stores
# those genes as the ligand's subunits. The ligand is the metabolite. Its ChEBI
# term is the one whose name or synonym is the metabolite's name; where the
# abbreviation names stereoisomers, the one the enzyme makes. None: the
# abbreviation does not tell which metabolite, and the ligand keeps it as its
# name only.
METABOLITES: dict[str, tuple[str, str] | None] = {
    '12oxoLTB4': ('CHEBI:27814', '12-dehydro-leukotriene B4'),
    '2AG': ('CHEBI:52392', '2-arachidonoylglycerol'),
    '5-HT': ('CHEBI:28790', 'serotonin'),
    # "5-HT-DOPA-DDC": pathway "SerotoninDopamin", serotonin and dopamine
    # receptors; DDC makes both.
    '5-HT-DOPA': None,
    # A 5alpha-reduced progesterone made by DHRS9, which reduces
    # 5alpha-pregnane-3,20-dione further: which one is not certain.
    '5alphaP': None,
    # Pathway "27HC" and NR1H4, but made by CYP11A1 and CYP3A4.
    '22HC': None,
    'ADO': ('CHEBI:16335', 'adenosine'),
    'ADR': ('CHEBI:28918', '(R)-adrenaline'),
    'Ach': ('CHEBI:15355', 'acetylcholine'),
    'Aldo': ('CHEBI:27584', 'aldosterone'),
    'Androstenedione': ('CHEBI:16422', 'androst-4-ene-3,17-dione'),
    'Androsterone': ('CHEBI:16032', 'androsterone'),
    'Calcitriol': ('CHEBI:17823', 'calcitriol'),
    'Cholesterol': ('CHEBI:16113', 'cholesterol'),
    # CYP11B1 makes cortisol in human and corticosterone in rodents.
    'Cort': None,
    'DHEA': ('CHEBI:28689', 'dehydroepiandrosterone'),
    'DHEAS': ('CHEBI:16814', 'dehydroepiandrosterone sulfate'),
    'DHT': ('CHEBI:16330', '17beta-hydroxy-5alpha-androstan-3-one'),
    'DOC': ('CHEBI:16973', '11-deoxycorticosterone'),
    # Pathway "Dopamine", dopamine receptors, made by DDC.
    'DOPA': ('CHEBI:18243', 'dopamine'),
    'Desmosterol': ('CHEBI:17737', 'desmosterol'),
    'E2': ('CHEBI:16469', '17beta-estradiol'),
    'GABA': ('CHEBI:16865', 'gamma-aminobutyric acid'),
    'Glu': ('CHEBI:16015', 'L-glutamic acid'),
    'Gly': ('CHEBI:15428', 'glycine'),
    'Hist': ('CHEBI:18295', 'histamine'),
    'LTB4': ('CHEBI:15647', 'leukotriene B4'),
    'LTC4': ('CHEBI:16978', 'leukotriene C4'),
    'LTD4': ('CHEBI:28666', 'leukotriene D4'),
    'LTE4': ('CHEBI:15650', 'leukotriene E4'),
    'LXA4': ('CHEBI:6498', 'lipoxin A4'),
    'Melatonin': ('CHEBI:16796', 'melatonin'),
    'NE': ('CHEBI:18357', '(R)-noradrenaline'),
    'PGD2': ('CHEBI:15555', 'prostaglandin D2'),
    'PGE2': ('CHEBI:15551', 'prostaglandin E2'),
    'PGF2a': ('CHEBI:15553', 'prostaglandin F2alpha'),
    'PGI2': ('CHEBI:15552', 'prostaglandin I2'),
    'Pregnenolone': ('CHEBI:16581', 'pregnenolone'),
    'Progesterone': ('CHEBI:17026', 'progesterone'),
    'RA': ('CHEBI:15367', 'all-trans-retinoic acid'),
    'T3': ('CHEBI:18258', "3,3',5-triiodo-L-thyronine"),
    'TXA2': ('CHEBI:15627', 'thromboxane A2'),
    'Testosterone': ('CHEBI:17347', 'testosterone'),
}


def metabolite(label: str) -> str | None:
    """The metabolite abbreviation a CellChat ligand label starts with.

    Peptide ligands of the same kind ("b-Endorphin-POMC") and protein hormones
    ("CGA_TSHB") are not in the table and stay proteins.
    """
    matches = [a for a in METABOLITES if label == a or label.startswith(f'{a}-')]
    return max(matches, key=len, default=None)


def _metabolite_entity(abbreviation: str) -> Entity:
    term = METABOLITES[abbreviation]
    if term is None:
        identifiers = [Identifier(type=Namespace.NAME, value=abbreviation)]
    else:
        identifiers = [
            Identifier(type=Namespace.CHEBI, value=term[0]),
            Identifier(type=Namespace.NAME, value=term[1]),
            Identifier(type=Namespace.SYNONYM, value=abbreviation),
        ]
    return Entity(type=ChemicalEntity, identifiers=identifiers)


def _row_metabolite(row: dict[str, Any], label_key: str) -> str | None:
    if _first(row.get('annotation')) != NON_PROTEIN:
        return None
    return metabolite(_first(row.get(label_key)))


_pmid_re = re.compile('PMID:?\\s*(\\d+)', re.IGNORECASE)
_pmc_re = re.compile('\\bPMC(?:ID)?\\s*:?\\s*(\\d+)', re.IGNORECASE)


def _values(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value).strip()
    return [text] if text else []


def _first(value: Any) -> str:
    vals = _values(value)
    return vals[0] if vals else ''


def _annotations(*items: Annotation | None) -> list[Annotation] | None:
    out: list[Annotation] = []
    seen: set[tuple[Any, Any, Any]] = set()
    for item in items:
        if item is None:
            continue
        key = (format_term(item.term), format_term(item.value), item.units)
        if key not in seen:
            out.append(item)
            seen.add(key)
    return out or None


def _annotation(term: Any, value: Any = None) -> Annotation | None:
    if value is None:
        return None
    value = str(value).strip()
    return Annotation(term=term, value=value) if value else None


def _pubmed_annotations(evidence: str) -> list[Annotation]:
    return [
        Annotation(term=slots.publications, value=f'PMID:{pmid}')
        for pmid in sorted(set(_pmid_re.findall(evidence or '')))
    ]


def _pmc_annotations(evidence: str) -> list[Annotation]:
    return [
        Annotation(term=slots.publications, value=f'PMC:PMC{pmc}')
        for pmc in sorted(set(_pmc_re.findall(evidence or '')))
    ]


def _protein_entity(gene: str, taxon_id: str) -> Entity:
    return Entity(
        type=Protein,
        identifiers=[Identifier(type=Namespace.GENESYMBOL, value=gene)],
        annotations=[
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}')
        ],
    )


def _protein_group_entity(
    *,
    name: str,
    genes: list[str],
    taxon_id: str,
    annotations: list[Annotation] | None = None,
    group_type: type = MacromolecularComplex,
) -> Entity:
    genes = [gene for gene in genes if gene]
    annotations = _annotations(
        Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
        *(annotations or []),
    )
    if len(genes) == 1:
        return Entity(
            type=Protein,
            identifiers=[
                Identifier(type=Namespace.GENESYMBOL, value=genes[0]),
                Identifier(type=Namespace.NAME, value=genes[0]),
            ],
            annotations=annotations,
        )
    return Entity(
        type=group_type if genes else NamedThing,
        identifiers=[Identifier(type=Namespace.NAME, value=name)],
        annotations=annotations,
        membership=[
            Membership(member=_protein_entity(gene, taxon_id)) for gene in genes
        ]
        or None,
    )


# Attributes CellChatDB gives an interaction ("Secreted Signaling", "Cell-Cell
# Contact", "ECM-Receptor", "Non-protein Signaling") and each of its partners
# (subcellular location, whether it is transmembrane, its secreted type), as it
# writes them. The partner attributes are participant-scoped relation annotations,
# with the partner's role: the predicate is symmetric, so the stored subject and
# object need not be the ligand and the receptor.
INTERACTION_CLASS = 'cellchat:annotation'
LOCATION = 'cellchat:location'
TRANSMEMBRANE = 'cellchat:transmembrane'
SECRETED_TYPE = 'cellchat:secreted_type'


def _participant_annotations(row: dict[str, Any], role: str) -> list[Annotation]:
    """The role of a partner (``ligand`` or ``receptor``) and its attributes."""
    return [
        annotation
        for annotation in (
            Annotation(term=PARTICIPANT_ROLE, value=role),
            _annotation(LOCATION, row.get(f'{role}_location')),
            _annotation(TRANSMEMBRANE, row.get(f'{role}_transmembrane')),
            _annotation(SECRETED_TYPE, row.get(f'{role}_secreted_type')),
        )
        if annotation is not None
    ]


def map_cellchat_interaction(row: dict[str, Any]) -> Relation:
    taxon_id = _first(row.get('taxon_id'))
    evidence = _first(row.get('evidence'))
    annotations = _annotations(
        Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
        _annotation(INTERACTION_CLASS, row.get('annotation')),
        *_pubmed_annotations(evidence),
        *_pmc_annotations(evidence),
    )
    abbreviation = _row_metabolite(row, 'ligand_name')
    ligand = (
        _metabolite_entity(abbreviation)._replace(
            annotations=_participant_annotations(row, 'ligand')
        )
        if abbreviation
        else _protein_group_entity(
            name=_first(row.get('ligand_name')),
            genes=_values(row.get('ligand_genes')),
            taxon_id=taxon_id,
            annotations=_participant_annotations(row, 'ligand'),
        )
    )
    receptor = _protein_group_entity(
        name=_first(row.get('receptor_name')),
        genes=_values(row.get('receptor_genes')),
        taxon_id=taxon_id,
        annotations=_participant_annotations(row, 'receptor'),
    )
    return Relation(
        subject=ligand,
        predicate=slots.interacts_with,
        object=receptor,
        identifiers=[
            Identifier(
                type=Namespace.CELLCHAT,
                value=_first(row.get('interaction_name')),
            ),
            Identifier(
                type=Namespace.NAME,
                value=_first(row.get('interaction_name_2'))
                or _first(row.get('interaction_name')),
            ),
        ],
        annotations=[*(annotations or [])],
    )


def map_cellchat_cofactor_interaction(row: dict[str, Any]) -> Relation:
    taxon_id = _first(row.get('taxon_id'))
    evidence = _first(row.get('evidence'))
    effect = _first(row.get('effect'))
    direction = {
        'stimulation': DirectionQualifierEnum.increased,
        'inhibition': DirectionQualifierEnum.decreased,
    }.get(effect)
    if direction is None:
        raise ValueError(f'Unknown CellChat cofactor effect: {effect!r}')
    cofactor_gene = _first(row.get('cofactor_gene'))
    annotations = _annotations(
        Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
        *_pubmed_annotations(evidence),
        *_pmc_annotations(evidence),
    )
    cofactor = _protein_group_entity(
        name=cofactor_gene, genes=[cofactor_gene], taxon_id=taxon_id
    )
    abbreviation = (
        _row_metabolite(row, 'target_name')
        if _first(row.get('target_role')) == 'ligand'
        else None
    )
    target = (
        _metabolite_entity(abbreviation)
        if abbreviation
        else _protein_group_entity(
            name=_first(row.get('target_name')),
            genes=_values(row.get('target_genes')),
            taxon_id=taxon_id,
        )
    )
    return Relation(
        subject=cofactor,
        predicate=slots.affects,
        object=target,
        identifiers=[
            Identifier(
                type=Namespace.CELLCHAT,
                value=_first(row.get('cofactor_interaction_name')),
            ),
            Identifier(
                type=Namespace.NAME,
                value=f'{cofactor_gene} {effect} {_first(row.get("target_name"))}',
            ),
        ],
        annotations=[
            *(annotations or []),
            Annotation(term=slots.object_direction_qualifier, value=direction),
            Annotation(term=slots.original_subject, value=_first(row.get("cofactor_group")) or cofactor_gene),
            Annotation(term=slots.original_predicate, value=_first(row.get("cofactor_role")) or effect),
        ],
    )


def map_cellchat_group(row: dict[str, Any]) -> Entity | None:
    # The genes of a metabolite ligand are no complex: they are in
    # the metabolite genes datasets.
    if row.get('non_protein_ligand') and metabolite(_first(row.get('name'))):
        return None
    return _protein_group_entity(
        name=_first(row.get('name')),
        genes=_values(row.get('genes')),
        taxon_id=_first(row.get('taxon_id')),
    )


def iter_metabolite_genes(
    opener, *, taxon_id: str, **kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Each metabolite ligand with each gene CellChat measures it by, once."""
    seen = set()
    for row in iter_cellchat_interactions(opener, taxon_id=taxon_id, **kwargs):
        if not (abbreviation := _row_metabolite(row, 'ligand_name')):
            continue
        for gene in _values(row.get('ligand_genes')):
            if (abbreviation, gene) not in seen:
                seen.add((abbreviation, gene))
                yield {'metabolite': abbreviation, 'gene': gene, 'taxon_id': taxon_id}


def map_metabolite_gene(row: dict[str, Any]) -> Relation:
    """A gene whose expression CellChat takes for the amount of a metabolite
    ligand: an enzyme that makes it or a transporter that releases it. CellChat
    does not say which, so the relation is an association."""
    taxon_id = _first(row.get('taxon_id'))
    return Relation(
        subject=_protein_entity(_first(row.get('gene')), taxon_id),
        predicate=slots.associated_with,
        object=_metabolite_entity(_first(row.get('metabolite'))),
        annotations=[
            Annotation(term=slots.in_taxon, value=f'NCBITaxon:{taxon_id}'),
        ],
    )


def map_cellchat_cofactor_group(row: dict[str, Any]) -> Entity:
    return _protein_group_entity(
        name=_first(row.get('name')),
        genes=_values(row.get('genes')),
        taxon_id=_first(row.get('taxon_id')),
        group_type=NamedThing,
    )


resource = Resource(
    config,
    human_interactions=Dataset(
        download=_download('human'),
        mapper=map_cellchat_interaction,
        raw_parser=_parser(iter_cellchat_interactions, 'human'),
    ),
    mouse_interactions=Dataset(
        download=_download('mouse'),
        mapper=map_cellchat_interaction,
        raw_parser=_parser(iter_cellchat_interactions, 'mouse'),
    ),
    zebrafish_interactions=Dataset(
        download=_download('zebrafish'),
        mapper=map_cellchat_interaction,
        raw_parser=_parser(iter_cellchat_interactions, 'zebrafish'),
    ),
    human_cofactor_interactions=Dataset(
        download=_download('human'),
        mapper=map_cellchat_cofactor_interaction,
        raw_parser=_parser(iter_cellchat_cofactor_interactions, 'human'),
    ),
    mouse_cofactor_interactions=Dataset(
        download=_download('mouse'),
        mapper=map_cellchat_cofactor_interaction,
        raw_parser=_parser(iter_cellchat_cofactor_interactions, 'mouse'),
    ),
    zebrafish_cofactor_interactions=Dataset(
        download=_download('zebrafish'),
        mapper=map_cellchat_cofactor_interaction,
        raw_parser=_parser(iter_cellchat_cofactor_interactions, 'zebrafish'),
    ),
    human_metabolite_genes=Dataset(
        download=_download('human'),
        mapper=map_metabolite_gene,
        raw_parser=_parser(iter_metabolite_genes, 'human'),
    ),
    mouse_metabolite_genes=Dataset(
        download=_download('mouse'),
        mapper=map_metabolite_gene,
        raw_parser=_parser(iter_metabolite_genes, 'mouse'),
    ),
    zebrafish_metabolite_genes=Dataset(
        download=_download('zebrafish'),
        mapper=map_metabolite_gene,
        raw_parser=_parser(iter_metabolite_genes, 'zebrafish'),
        # CellChatDB.zebrafish has no metabolite ligands.
        allow_empty=True,
    ),
    human_complexes=Dataset(
        download=_download('human'),
        mapper=map_cellchat_group,
        raw_parser=_parser(iter_cellchat_complexes, 'human'),
    ),
    mouse_complexes=Dataset(
        download=_download('mouse'),
        mapper=map_cellchat_group,
        raw_parser=_parser(iter_cellchat_complexes, 'mouse'),
    ),
    zebrafish_complexes=Dataset(
        download=_download('zebrafish'),
        mapper=map_cellchat_group,
        raw_parser=_parser(iter_cellchat_complexes, 'zebrafish'),
    ),
    human_cofactors=Dataset(
        download=_download('human'),
        mapper=map_cellchat_cofactor_group,
        raw_parser=_parser(iter_cellchat_cofactors, 'human'),
    ),
    mouse_cofactors=Dataset(
        download=_download('mouse'),
        mapper=map_cellchat_cofactor_group,
        raw_parser=_parser(iter_cellchat_cofactors, 'mouse'),
    ),
    zebrafish_cofactors=Dataset(
        download=_download('zebrafish'),
        mapper=map_cellchat_cofactor_group,
        raw_parser=_parser(iter_cellchat_cofactors, 'zebrafish'),
    ),
)
