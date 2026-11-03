"""
Parse the Human Plasma Membrane Receptome (HPMR) and emit Relation records.

HPMR lists, on the page of each human plasma membrane receptor, the ligands
and the other receptors it interacts with, each entry with a comment and the
PubMed records that support it. The website is gone; the pages are read from
the copy OmniPath saved in 2021.
"""

from __future__ import annotations
from collections.abc import Callable
import functools
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    ChemicalEntity,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from omnipath_core.source_attributes import PARTICIPANT_ROLE
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Relation,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.hpmr import NON_GENE, iter_interactions, pubmed_ids

config = ResourceConfig(
    id=ResourceCv.HPMR,
    name='HPMR',
    full='Human Plasma Membrane Receptome',
    # receptome.org now belongs to an unrelated site; this is the publication.
    url='https://doi.org/10.1126/stke.2003.187.re9',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='12815191',
    topics=(Topic.CELL_CELL_COMMUNICATION, Topic.SIGNALING),
    content=(
        Content.LIGAND_RECEPTOR_INTERACTIONS,
        Content.PHYSICAL_INTERACTIONS,
    ),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Literature-curated ligands and interacting receptors of human plasma membrane receptors.',
    description='The Human Plasma Membrane Receptome is a database of the human plasma membrane receptors involved in signal transduction, grouped into families by sequence. For each receptor its curators collected from the literature the ligands it binds and the other receptors it interacts with, with a comment and the supporting publications. The website has been offline for years; OmniPath keeps a copy of its receptor pages.',
)
# receptome.org no longer serves the database. This is the copy of the
# receptor pages that OmniPath saved from the site in November 2021.
download = Download(
    url='https://rescued.omnipathdb.org/hpmr/genes.tar.gz',
    filename='hpmr_genes.tar.gz',
    subfolder='hpmr',
    default_mode='rb',
    encoding=None,
    ext='tar.gz',
)
# The pages give no species; every receptor and ligand gene is human.
TAXON = 'NCBITaxon:9606'
# Ligands that are not gene products, as the site spells them: is the ligand
# a chemical? Light and taste are stimuli; the porin and flagellin are
# bacterial proteins and the RNAs are nucleic acids, none identified by an
# accession.
NON_GENE_LIGANDS = {
    'ADP': True,
    'ATP': True,
    'Calcium': True,
    'Fatty_acids': True,
    'Hemozoin': True,
    'Hydroxyeicosanoids': True,
    'L-2-amino-4-phosphonobutyrate': True,
    'L-glutamate': True,
    'L-serine-O-phosphate': True,
    'Medium_chain_fatty_acids': True,
    'N-acetylaspartylglutamate': True,
    'Oleic_acid': True,
    'PGE2': True,
    'Proton': True,
    'Short-chain_fatty_acids': True,
    'Sphingosine_1-phosphate': True,
    'UDP': True,
    'acetate': True,
    'adenosine': True,
    'alpha-ketoglutarate': True,
    'anandamide': True,
    'beta-alanine': True,
    'bile_acids': True,
    'cation': True,
    'dioleoylphosphatidic_acid': True,
    'estrogen': True,
    'folic_acid': True,
    'glutamate': True,
    'hexarelin_a_synthetic_peptide': True,
    'kynurenic_acid': True,
    'leukotriene_B4': True,
    'leukotriene_C4': True,
    'lipopolysaccharide': True,
    'lysophosphatidic_acid': True,
    'niacin_(nicotinic_acid)': True,
    'oxoeicosanoid': True,
    'propionate': True,
    'prostacyclin': True,
    'sphingosine_1_phosphate': True,
    'succinate': True,
    'thrombosin_A2': True,
    'Bacterial_Flagellin': False,
    'dsRNA': False,
    'light': False,
    'meningococcal_porin_PorB': False,
    'ssRNA': False,
    'taste': False,
}
# For some receptors of a small molecule HPMR names as the ligand an enzyme
# that makes or breaks down that molecule (tryptophan hydroxylase for the
# serotonin receptors, hyaluronan synthase for CD44). The enzyme does not bind
# the receptor and the molecule is not named, so these pairs are not emitted.
ENZYMES_FOR_LIGANDS = frozenset(
    {
        'AANAT',
        'ACP6',
        'ASMT',
        'CHAT',
        'CYP11A1',
        'DBH',
        'GAD1',
        'HAS2',
        'HDC',
        'LIPH',
        'PGDS',
        'PTDSS1',
        'PTGS2',
        'TH',
        'TPH1',
    }
)
# The ligand the site titles `LAMININ, GAMMA-1` carries the symbol and gene
# of DUSP18, a phosphatase. Its pairs are not emitted.
MISLABELLED_LIGANDS = frozenset({'DUSP18'})
def _pubmeds(row: dict[str, Any]) -> list[str]:
    """PubMed IDs of an entry in page order, each once."""
    return pubmed_ids(row.get('references'))


def _roles(role: str | None) -> list[Annotation]:
    """The participant role, where the page states one of the two roles."""
    return [Annotation(term=PARTICIPANT_ROLE, value=role)] if role else []


def _protein(
    symbol: str, entrez: Any = None, role: str | None = None
) -> Entity:
    identifiers = [Identifier(type=Namespace.GENESYMBOL, value=symbol)]
    if entrez:
        identifiers.append(Identifier(type=Namespace.ENTREZ, value=str(entrez)))
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=[
            Annotation(term=slots.in_taxon, value=TAXON),
            *_roles(role),
        ],
    )


def _chemical(name: str, role: str | None = None) -> Entity | None:
    if name not in NON_GENE_LIGANDS:
        raise ValueError(f'Unknown HPMR non-gene ligand {name!r}')
    if not NON_GENE_LIGANDS[name]:
        return None
    return Entity(
        type=ChemicalEntity,
        # The site writes the spaces of a name as underscores.
        identifiers=[
            Identifier(type=Namespace.NAME, value=name.replace('_', ' '))
        ],
        annotations=_roles(role),
    )


def _partner(row: dict[str, Any]) -> Entity | None:
    name = str(row.get('partner') or '').strip()
    role = row.get('partner_role')
    if not name:
        return None
    # The page marks a partner as a `Ligand` or a `non-gene ligand`; a partner
    # marked `Receptor` is another receptor, so it has neither role.
    if role == NON_GENE:
        return _chemical(name, 'ligand')
    if role == 'Ligand' and (
        name in ENZYMES_FOR_LIGANDS or name in MISLABELLED_LIGANDS
    ):
        return None
    return _protein(name, role='ligand' if role == 'Ligand' else None)


def _annotations(row: dict[str, Any]) -> list[Annotation]:
    out = [
        Annotation(term=slots.publications, value=f'PMID:{pmid}')
        for pmid in _pubmeds(row)
    ]
    if comment := str(row.get('comment') or '').strip():
        out.append(Annotation(term=slots.description, value=comment))
    return out


def _mapper(roles: tuple[str, ...]) -> Callable[[dict[str, Any]], Relation | None]:
    def mapper(row: dict[str, Any]) -> Relation | None:
        role = row.get('partner_role')
        receptor = str(row.get('receptor') or '').strip()
        if role not in roles or not receptor:
            return None
        partner = _partner(row)
        if partner is None:
            return None
        # The partner is the subject and the page's receptor the object: a
        # ligand and its receptor, or two receptors that interact. The page's
        # receptor is a receptor only against a ligand; between two receptors
        # the source states no ligand-receptor roles.
        return Relation(
            subject=partner,
            predicate=slots.interacts_with,
            object=_protein(
                receptor,
                row.get('receptor_entrez'),
                'receptor' if role in LIGAND_ROLES else None,
            ),
            identifiers=[],
            annotations=_annotations(row),
        )

    return mapper


LIGAND_ROLES = ('Ligand', NON_GENE)
RECEPTOR_ROLES = ('Receptor',)
map_hpmr_interaction = _mapper(LIGAND_ROLES)
map_hpmr_receptor_interaction = _mapper(RECEPTOR_ROLES)
resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=map_hpmr_interaction,
        raw_parser=functools.partial(iter_interactions, roles=LIGAND_ROLES),
    ),
    receptor_interactions=Dataset(
        download=download,
        mapper=map_hpmr_receptor_interaction,
        raw_parser=functools.partial(iter_interactions, roles=RECEPTOR_ROLES),
    ),
)
