"""
Parse the ligand-receptor pairs of Baccin et al. 2019 and emit Relation records.

Supplementary Table 3 of the paper is the ligand-receptor database of the
RNA-Magnet method (version 1.0.0 of the RNAMagnet package). Each row is a pair
of mouse gene symbols with its origin, the authors' verdict after checking the
literature, and the references. A partner is one gene, several genes joined by
`&` that form a complex, or several genes joined by `|` of which any one is
the partner. Only pairs the authors consider correct are emitted.
"""

from __future__ import annotations
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    MacromolecularComplex,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from omnipath_core.source_attributes import PARTICIPANT_ROLE
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
    Relation,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.baccin2019 import (
    REFERENCE,
    iter_pairs,
    pubmed_ids,
    reference_tokens,
)

config = ResourceConfig(
    id=ResourceCv.BACCIN2019,
    name='Baccin2019',
    full='Baccin 2019 ligand-receptor pairs',
    url='https://doi.org/10.1038/s41556-019-0439-6',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='31871321',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Manually checked mouse ligand–receptor pairs, with heteromeric receptors and their publications.',
    description='Baccin et al. 2019 mapped the cellular and spatial organization of the mouse bone marrow niche by single-cell and spatial transcriptomics, and predicted cell-cell interactions with the RNA-Magnet method. Its ligand-receptor database starts from the pairs of Ramilowski et al. 2015 translated to mouse genes; the authors checked the literature for every pair expressed in bone marrow, marked each as correct or incorrect, and added pairs of their own, including heterodimeric receptors such as the integrins.',
)
download = Download(
    url='https://static-content.springer.com/esm/art%3A10.1038%2Fs41556-019-0439-6/MediaObjects/41556_2019_439_MOESM3_ESM.xlsx',
    filename='baccin2019_supplementary_tables.xlsx',
    subfolder='baccin2019',
    default_mode='rb',
    encoding=None,
)
# The table has no species column; its partners are MGI symbols.
TAXON = 'NCBITaxon:10090'
# `ManualAnnotation`, the authors' verdict on a pair: is it emitted?
# RNAMagnet documents `Correct`, `Incorrect`, `Irrelevant` (correct but not
# relevant at homeostasis) and `Not Expressed` (not checked, because a partner
# is not expressed in bone marrow). `Scavenger` and `Likely` are not
# documented. Unchecked pairs are Ramilowski pairs the authors translated to
# mouse; RNAMagnet 2.2.0, which checked them, marks 287 of the 1,072 incorrect.
ANNOTATION = {
    'Correct': True,
    'Irrelevant': True,
    'Incorrect': False,
    'Not Expressed': False,
    'Scavenger': False,
    'Likely': False,
}
# `Source`: the database a pair was taken from. `Baccin` marks a pair the
# authors added themselves; its evidence is the row's references.
SOURCES = {
    'Ramilowski': 'Ramilowski2015',
    'Baccin': None,
}
# Reference cells that name a KEGG pathway map instead of a publication.
REFERENCE_DATABASES = {
    'mmu04514': 'KEGG',
    'mmu:04650': 'KEGG',
}
# Four partners are the stem of a gene family's symbols (laminins,
# thrombospondins, ADAMs and fibrillins), alone or among alternatives. A stem
# is not a gene and identifies no protein, so its pairs are not emitted.
FAMILY_STEMS = frozenset({'Lam', 'Thbs', 'Adam', 'Fbn'})
# Two ligands join the amyloid precursor with the proteases that cleave it.
# These are not subunits of a complex, and the table marks the pairs of App
# alone with the same receptors as incorrect, so the pairs are not emitted.
PROCESSING_GROUPS = frozenset({'App&Bace1', 'App&Bace1&Psen'})


def _pubmeds(row: dict[str, Any]) -> list[str]:
    """PubMed IDs of a pair in source order, each once."""
    return pubmed_ids(row.get(REFERENCE))


def _data_sources(row: dict[str, Any]) -> list[str]:
    """Databases supporting a pair: its origin and those its references name."""
    source = row.get('Source')
    if source not in SOURCES:
        raise ValueError(f'Unknown Baccin 2019 source {source!r}')
    out = [SOURCES[source]] if SOURCES[source] else []
    for token in reference_tokens(row.get(REFERENCE)):
        database = REFERENCE_DATABASES.get(token)
        if database and database not in out:
            out.append(database)
    return out


def _protein(symbol: str) -> Entity:
    return Entity(
        type=Protein,
        identifiers=[Identifier(type=Namespace.GENESYMBOL, value=symbol)],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def _participant(value: str, role: str) -> Entity | None:
    """A protein, or a complex when the source joins its subunits by `&`.

    The table has a ligand column and a receptor column; the role is on the
    protein or on the complex as a whole, not on its subunits.
    """
    symbols = [token for part in value.split('&') if (token := part.strip())]
    if not symbols:
        return None
    role_annotation = Annotation(term=PARTICIPANT_ROLE, value=role)
    if len(symbols) == 1:
        protein = _protein(symbols[0])
        protein.annotations.append(role_annotation)
        return protein
    return Entity(
        type=MacromolecularComplex,
        # The complex has no accession; the cell is the source's name for it.
        identifiers=[Identifier(type=Namespace.NAME, value=value)],
        annotations=[
            Annotation(term=slots.in_taxon, value=TAXON),
            role_annotation,
        ],
        membership=[Membership(member=_protein(symbol)) for symbol in symbols],
    )


def map_baccin_interaction(row: dict[str, Any]) -> Relation | None:
    verdict = row.get('ManualAnnotation')
    if verdict not in ANNOTATION:
        raise ValueError(f'Unknown Baccin 2019 manual annotation {verdict!r}')
    ligand = str(row.get('ligand') or '').strip()
    receptor = str(row.get('receptor') or '').strip()
    if not ANNOTATION[verdict] or not ligand or not receptor:
        return None
    if ligand in FAMILY_STEMS or receptor in FAMILY_STEMS:
        return None
    if ligand in PROCESSING_GROUPS or receptor in PROCESSING_GROUPS:
        return None
    # The ligand is the subject and the receptor the object.
    return Relation(
        subject=_participant(ligand, 'ligand'),
        predicate=slots.interacts_with,
        object=_participant(receptor, 'receptor'),
        identifiers=[
            Identifier(type=Namespace.NAME, value=name)
            for name in [row.get('Pair.Name')]
            if name
        ],
        annotations=[
            *[
                Annotation(term=slots.supporting_data_source, value=database)
                for database in _data_sources(row)
            ],
            *[
                Annotation(term=slots.publications, value=f'PMID:{pmid}')
                for pmid in _pubmeds(row)
            ],
        ],
    )


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=map_baccin_interaction,
        raw_parser=iter_pairs,
    ),
)
