"""
Parse the Cell Surface Protein Atlas (CSPA) and emit protein-to-list and protein-to-cell-type records.

Bausch-Fluck et al. 2015 applied Cell Surface Capture, which labels and
isolates the N-glycoproteins on the outside of living cells, to 41 human and 31
mouse cell types and identified the proteins by mass spectrometry. The atlas has
1,492 human and 1,296 mouse proteins, each with a confidence category.

Two kinds of statement are imported. A protein is related to the term of the
atlas (a plain list: every protein of the lists is one of its members), with the
category and the annotations of the list on the statement. A protein is related
to each cell type it was identified on (the source's own labels, which carry no
ontology accession), with the identification's probability and number of unique
peptides.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.cspa import (
    CELL_TYPE,
    LIST,
    iter_detections,
    iter_human_proteins,
    iter_mouse_proteins,
    iter_terms,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Relation

config = ResourceConfig(
    id=ResourceCv.CSPA,
    name='CSPA',
    full='Cell Surface Protein Atlas',
    url='https://wollscheidlab.org/CSPA/',
    # The paper (PLoS ONE 2015): "This is an open-access article distributed under the
    # terms of the Creative Commons Public Domain declaration ... " (CC0). The website
    # states no license of its own. The legacy registry says CC BY 4.0.
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='25894527',
    topics=(Topic.SUBCELLULAR_LOCALIZATION, Topic.CELL_TYPES),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL,),
    summary=(
        'Cell surface glycoproteins of human and mouse found by Cell Surface '
        'Capture, with a confidence category and the cell types they were '
        'identified on.'
    ),
    description=(
        'The Cell Surface Protein Atlas is a mass spectrometry derived atlas of '
        'cell surface glycoproteins. Cell Surface Capture was applied to 41 '
        'human and 31 mouse cell types, giving 1,492 human and 1,296 mouse '
        'proteins with experimental evidence for their expression at the cell '
        'surface, 136 G-protein coupled receptors and 75 membrane receptor '
        'tyrosine kinases among them. Each protein has a category: high '
        'confidence, putative or unspecific.'
    ),
)

# The web site (the old address, wlab.ethz.ch/cspa, no longer serves the files).
S1_URL = 'https://wollscheidlab.org/CSPA/data/S1_File.xlsx'
S2_URL = 'https://wollscheidlab.org/CSPA/data/S2_File.xlsx'
# Proteins of the atlas with their category and annotations (S2).
proteins_download = Download(
    url=S2_URL,
    filename='cspa_S2_File.xlsx',
    subfolder='cspa',
    large=True,
    default_mode='rb',
    encoding=None,
)
# Identifications of proteins in cell types (S1).
detections_download = Download(
    url=S1_URL,
    filename='cspa_S1_File.xlsx',
    subfolder='cspa',
    large=True,
    default_mode='rb',
    encoding=None,
)

# The atlas as a plain list (one term), and the source's own cell types.
LIST_NAMESPACE = 'cspa_list'
LIST_ID = 'cell_surface_protein_atlas'
CELL_TYPE_NAMESPACE = 'cspa_cell_type'
TAXA = {'Human': 'NCBITaxon:9606', 'Mouse': 'NCBITaxon:10090'}

# Source attributes of the statement that puts a protein in the list. The
# probability, the count of cell types and the number of predicted
# transmembrane domains are those of the table, Measurements without a unit; the
# flags are the table's columns (``yes``/``no``, ``1``/``0``), written
# ``'true'``/``'false'``.
_N_CELL_TYPES = 'cspa:n_cell_types'
_PHOBIUS_TM = 'cspa:phobius_transmembrane_predicted'
_TM_DOMAINS = 'cspa:predicted_transmembrane_domains'
_GPI = 'cspa:gpi'
_UNIPROT_CELL_SURFACE = 'cspa:uniprot_cell_surface'
# Source attribute of the identification of a protein in a cell type.
_UNIQUE_PEPTIDES = 'cspa:unique_peptides'


def _flag(value: bool) -> str:
    """A flag of the table, read by the parser into a boolean, as 'true' or 'false'."""
    if value is True:
        return 'true'
    if value is False:
        return 'false'
    raise ValueError(f'Not a flag: {value!r}')


def _count(row: dict[str, Any], field: str, column: str):
    """A whole number of the table as a Measurement; the unit is the thing counted."""
    value = measurement(row[field], source_field=column)
    if value is None:
        raise ValueError(f'CSPA {row["uniprot"]}: {column} is {row[field]!r}')
    return value


def _protein(row: dict[str, Any]) -> Entity:
    identifiers = [
        Identifier(type=Namespace.UNIPROT, value=row['uniprot']),
        Identifier(type=Namespace.UNIPROT_ENTRY, value=row['entry_name']),
    ]
    if row.get('gene'):
        identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=row['gene']))
    if row.get('entrez'):
        identifiers.append(Identifier(type=Namespace.ENTREZ, value=row['entrez']))
    if row.get('cd'):
        identifiers.append(Identifier(type=Namespace.SYNONYM, value=row['cd']))
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXA[row['organism']])],
    )


def _probability(row: dict[str, Any]) -> Annotation:
    # The probability of the protein identification, between 0 and 1.
    return Annotation(
        term=slots.has_confidence_score,
        value=measurement(row['probability'], source_field='protein probability'),
    )


def map_surface_protein(row: dict[str, Any]) -> Relation:
    """A protein of the atlas, in the list with its category."""
    annotations = [
        # "1 - high confidence", "2 - putative" or "3 - unspecific", as written.
        Annotation(term=slots.has_confidence_level, value=row['category']),
        _probability(row),
        Annotation(
            term=_N_CELL_TYPES,
            value=_count(
                row, 'n_cell_types', 'count detection in different cell types'
            ),
        ),
        Annotation(term=_PHOBIUS_TM, value=_flag(row['phobius_tm'])),
        Annotation(
            term=_TM_DOMAINS,
            value=_count(row, 'tm_domains', '# predicted TM domains'),
        ),
        Annotation(term=_GPI, value=_flag(row['gpi'])),
    ]
    # The table writes "#N/A" where UniProt does not annotate the protein.
    if row['uniprot_cell_surface']:
        annotations.append(
            Annotation(term=_UNIPROT_CELL_SURFACE, value=_flag(True))
        )
    return Relation(
        subject=_protein(row),
        # The third category lists proteins the atlas calls unspecific, so the
        # list does not state that its members are at the cell surface, and the
        # statement is an association, not `located_in`.
        predicate=slots.associated_with,
        object=Entity(
            type=OntologyClass,
            identifiers=[Identifier(type=LIST_NAMESPACE, value=LIST_ID)],
        ),
        annotations=annotations,
    )


def map_detection(row: dict[str, Any]) -> Relation:
    """A protein identified on the surface of a cell type."""
    return Relation(
        subject=_protein(row),
        # The protein is found on the surface of the cell type's cells. The
        # cell types are the source's own labels, the cells of a line or a
        # primary culture, an anatomical entity as `expressed_in` requires.
        predicate=slots.expressed_in,
        object=Entity(
            type=OntologyClass,
            identifiers=[
                Identifier(type=CELL_TYPE_NAMESPACE, value=row['cell_type'])
            ],
        ),
        annotations=[
            _probability(row),
            Annotation(
                term=_UNIQUE_PEPTIDES,
                value=_count(row, 'peptides', 'num unique peps'),
            ),
        ],
    )


def map_term(row: dict[str, Any]) -> Entity:
    """The list of the atlas, or a cell type."""
    if row['kind'] == LIST:
        return Entity(
            type=OntologyClass,
            identifiers=[
                Identifier(type=LIST_NAMESPACE, value=LIST_ID),
                Identifier(type=Namespace.NAME, value='Cell Surface Protein Atlas'),
            ],
            annotations=[
                Annotation(
                    term=slots.description,
                    value=(
                        'Proteins identified as cell surface glycoproteins by '
                        'Cell Surface Capture and mass spectrometry on human '
                        'and mouse cell types.'
                    ),
                )
            ],
        )
    if row['kind'] == CELL_TYPE:
        return Entity(
            type=OntologyClass,
            identifiers=[
                Identifier(type=CELL_TYPE_NAMESPACE, value=row['cell_type']),
                Identifier(type=Namespace.NAME, value=row['cell_type']),
            ],
        )
    raise ValueError(f'Unknown kind of CSPA term: {row["kind"]!r}')


resource = Resource(
    config,
    terms=Dataset(
        download=detections_download,
        mapper=map_term,
        raw_parser=iter_terms,
    ),
    human_proteins=Dataset(
        download=proteins_download,
        mapper=map_surface_protein,
        raw_parser=iter_human_proteins,
    ),
    mouse_proteins=Dataset(
        download=proteins_download,
        mapper=map_surface_protein,
        raw_parser=iter_mouse_proteins,
    ),
    detections=Dataset(
        download=detections_download,
        mapper=map_detection,
        raw_parser=iter_detections,
    ),
)
