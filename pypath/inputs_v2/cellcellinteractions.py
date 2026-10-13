"""
Parse the Cell-Cell Interaction Database of the Bader lab and emit Entity records.

The database's molecule definition files list the genes that the lab classified as
ligands, receptors and extracellular matrix (ECM) proteins, first from Gene Ontology
and UniProt annotations and then by manual curation. This module reads the table
that holds all three lists in one (``protein_types.txt``):

    hgnc_symbol   curated_protein_type
    SBSPON        ECM/Receptor/Ligand
    A2M           Ligand

A gene is in one, two or all three of the lists; a cell with several roles joins
them with ``/``. Each role is a term of the database, and each gene is associated
with the terms of its roles. The ligand-receptor interaction file of the same page
is not read.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cellcellinteractions import iter_genes, iter_roles
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Association,
    Entity,
    EntityRef,
    Identifier,
)

config = ResourceConfig(
    id=ResourceCv.CELLCELLINTERACTIONS,
    name='CellCellInteractions',
    full='Cell-Cell Interaction Database',
    url='https://baderlab.org/CellCellInteractions',
    # The Zenodo record of the database (doi:10.5281/zenodo.7589953), which holds
    # the files of the Bader lab page, gives the license as CC-BY-4.0. The page
    # states none.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    # The page of the database (baderlab.org/CellCellInteractions): "the automated
    # construction of a cell-cell interaction database"; the receptor, ligand and
    # ECM genes "were defined based on the union of the annotations" of Gene
    # Ontology terms, UniProt keywords and the Human Protein Atlas secretome, and
    # then "manually curated" (genes removed or moved to the right list). The page
    # names no publications as the basis of the curation.
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary=(
        'Genes classified by the Bader lab as ligands, receptors and '
        'extracellular matrix proteins.'
    ),
    description=(
        'The Cell-Cell Interaction Database of the Bader lab at the University of '
        'Toronto is an automatically built set of protein-protein interactions '
        'that can mediate cell-cell communication. It rests on lists of ligand, '
        'receptor and extracellular matrix genes, defined from Gene Ontology and '
        'UniProt annotations and the human secretome, and then curated by hand: '
        'genes that were neither receptors nor ligands were removed and '
        'misclassified genes moved to the right list. A gene can be in several '
        'lists.'
    ),
)

# The roles of the table and how the page of the database defines them: the GO
# terms and annotations that each list was built from, before manual curation.
ROLE_NAMESPACE = 'cellcellinteractions_role'
ROLES = {
    'ECM': (
        'Extracellular matrix',
        'Genes annotated with the Gene Ontology terms extracellular matrix '
        '(GO:0031012), proteinaceous extracellular matrix (GO:0005578), '
        'extracellular matrix structural constituent (GO:0005201), extracellular '
        'matrix protein binding (GO:1990430) or extracellular matrix cell '
        'signalling (GO:0035426), after manual curation.',
    ),
    'Ligand': (
        'Ligand',
        'Genes annotated with the Gene Ontology term receptor binding '
        '(GO:0005102) or labelled as secreted in the Human Protein Atlas '
        'secretome, after manual curation.',
    ),
    'Receptor': (
        'Receptor',
        'Genes annotated with the Gene Ontology terms receptor complex '
        '(GO:0043235), integrin complex (GO:0008305), protein localized to '
        'membrane (GO:0072657), receptor clustering (GO:0043113), receptor '
        'activity (GO:0004872) or external side of plasma membrane (GO:0009897), '
        'or with the UniProt keyword Receptor (KW-0675), after manual curation.',
    ),
}

# The symbols are HGNC symbols, so the genes are human.
TAXON = 'NCBITaxon:9606'


def _role(role: str) -> str:
    if role not in ROLES:
        raise ValueError(f'Unknown Cell-Cell Interaction Database role: {role!r}')
    return role


def map_role(row: dict[str, Any]) -> Entity:
    """A role of the table as a term, with the definition of the page."""
    name, definition = ROLES[_role(row['role'])]
    return Entity(
        type=model.OntologyClass,
        identifiers=[
            Identifier(type=ROLE_NAMESPACE, value=row['role']),
            Identifier(type=Namespace.NAME, value=name),
        ],
        annotations=[Annotation(term=slots.description, value=definition)],
    )


def map_gene(row: dict[str, Any]) -> Entity:
    """A gene associated with the term of each of its roles."""
    return Entity(
        type=model.Gene,
        identifiers=[Identifier(type=Namespace.GENESYMBOL, value=row['symbol'])],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        associations=[
            Association(
                object=EntityRef(model.OntologyClass, ROLE_NAMESPACE, _role(role)),
                predicate=slots.associated_with,
            )
            for role in row['roles']
        ],
    )


# The file the Bader lab page links. The Zenodo record of the database has the
# same file (checked byte for byte).
download = Download(
    url=(
        'https://baderlab.org/CellCellInteractions?'
        'action=AttachFile&do=get&target=protein_types.txt'
    ),
    filename='cellcellinteractions_protein_types.txt',
    subfolder='cellcellinteractions',
    large=True,
)

resource = Resource(
    config,
    terms=Dataset(download=download, mapper=map_role, raw_parser=iter_roles),
    genes=Dataset(download=download, mapper=map_gene, raw_parser=iter_genes),
)
