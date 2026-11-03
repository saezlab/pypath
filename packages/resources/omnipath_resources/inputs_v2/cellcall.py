"""
Parse CellCall ligand-receptor pairs and emit Relation records.

CellCall distributes four tables: a core and an extended set, each for human
and, inferred by orthology, for mouse. A row is a ligand-receptor-transcription
factor axis: a ligand-receptor pair and one transcription factor downstream of
it in a KEGG pathway. Ligands, receptors and transcription factors are named
by NCBI Gene IDs and gene symbols. This module emits the ligand-receptor
pairs; the transcription factors and pathways stay in the parsed row.
"""

from __future__ import annotations
from collections.abc import Callable
import functools
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    MacromolecularComplex,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from omnipath_core.source_attributes import PARTICIPANT_ROLE
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cellcall import SUBUNIT_DELIMITER, iter_pairs

config = ResourceConfig(
    id=ResourceCv.CELLCALL,
    name='CellCall',
    url='https://github.com/ShellyCoder/cellcall',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.ARTISTIC_2_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='34331449',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary='Ligand–receptor pairs of the core and extended CellCall sets for human and mouse.',
    description='CellCall is a toolkit that infers cell-cell communication from single-cell transcriptomes by combining ligand and receptor expression with the activity of downstream transcription factors. Its reference data are ligand-receptor-transcription factor axes: ligand-receptor interactions with literature or experimental support collected from NATMI, Cellinker, CellTalkDB, CellChat and STRING, including multi-subunit receptors, each linked to the transcription factors downstream of it in the same branch of a KEGG pathway. A core set and an extended set are provided for human, and the mouse sets were derived from them by orthology.',
)
BASE_URL = (
    'https://raw.githubusercontent.com/ShellyCoder/cellcall/master/inst/extdata/'
)
# The tables have no species column. The file name tells the species: the
# ``homology`` files hold the mouse orthologs of the human sets.
FILES = {
    'human_core': ('new_ligand_receptor_TFs.txt', '9606'),
    'human_extended': ('new_ligand_receptor_TFs_extended.txt', '9606'),
    'mouse_core': ('new_ligand_receptor_TFs_homology.txt', '10090'),
    'mouse_extended': ('new_ligand_receptor_TFs_homology_extended.txt', '10090'),
}


def _download(dataset: str) -> Download:
    filename, _ = FILES[dataset]
    return Download(
        url=f'{BASE_URL}{filename}',
        filename=f'cellcall_{filename}',
        subfolder='cellcall',
    )


def _parser(dataset: str) -> Any:
    return functools.partial(iter_pairs, taxon_id=FILES[dataset][1])


f = FieldConfig(
    transform={
        'taxon_curie': lambda v: f'NCBITaxon:{v}',
    },
)
taxon = CV(term=slots.in_taxon, value=f('taxon_id', transform='taxon_curie'))


def _participant(role: str) -> Callable[[dict[str, Any]], Any]:
    symbol, gene_id = f'{role}_Symbol', f'{role}_ID'
    # The columns are named for the role (`Ligand_*`, `Receptor_*`). A
    # multi-subunit receptor is the receptor as a whole, so the role is on the
    # complex and not on its subunits.
    role_cv = CV(term=PARTICIPANT_ROLE, value=role.lower())
    protein = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(symbol)),
            CV(term=Namespace.ENTREZ, value=f(gene_id)),
        ),
        annotations=AnnotationsBuilder(taxon, role_cv),
    )
    # The source gives a multi-subunit receptor no accession or name of its
    # own, only its subunits; the symbol cell as written names the complex.
    complex_ = EntityBuilder(
        entity_type=MacromolecularComplex,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.NAME, value=f(symbol)),
        ),
        annotations=AnnotationsBuilder(taxon, role_cv),
        membership=MembershipBuilder(
            MembersFromList(
                entity_type=Protein,
                identifiers=IdentifiersBuilder(
                    # Symbol and gene ID of a subunit share a position.
                    CV(
                        term=Namespace.GENESYMBOL,
                        value=f(
                            symbol,
                            delimiter=SUBUNIT_DELIMITER,
                            preserve_indices=True,
                        ),
                    ),
                    CV(
                        term=Namespace.ENTREZ,
                        value=f(
                            gene_id,
                            delimiter=SUBUNIT_DELIMITER,
                            preserve_indices=True,
                        ),
                    ),
                ),
                entity_annotations=AnnotationsBuilder(taxon),
            )
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        composite = SUBUNIT_DELIMITER in (row.get(gene_id) or '')
        return (complex_ if composite else protein).build(row)

    return build


# The ligand is the subject and the receptor the object. A pair's downstream
# transcription factors and KEGG pathways are not part of the relation.
interactions_schema = RelationBuilder(
    subject=_participant('Ligand'),
    predicate=slots.interacts_with,
    object=_participant('Receptor'),
)
resource = Resource(
    config,
    human_core_interactions=Dataset(
        download=_download('human_core'),
        mapper=interactions_schema,
        raw_parser=_parser('human_core'),
    ),
    human_extended_interactions=Dataset(
        download=_download('human_extended'),
        mapper=interactions_schema,
        raw_parser=_parser('human_extended'),
    ),
    mouse_core_interactions=Dataset(
        download=_download('mouse_core'),
        mapper=interactions_schema,
        raw_parser=_parser('mouse_core'),
    ),
    mouse_extended_interactions=Dataset(
        download=_download('mouse_extended'),
        mapper=interactions_schema,
        raw_parser=_parser('mouse_extended'),
    ),
)
