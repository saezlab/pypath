"""
Cell type marker genes of PanglaoDB.

PanglaoDB is a database of single cell RNA sequencing experiments with a
compendium of cell type markers. Each row of its marker table says that a gene
is a marker of a cell type, for human, mouse or both. The module imports:

* ``markers``: a statement that a gene is expressed in a cell type, one for each
  species the row states, with the sensitivity and specificity of that species;
* ``terms``: the cell types, organs and germ layers of the table as entities
  with names, each cell type ``part_of`` its organ.

The cell types, organs and germ layers are PanglaoDB's own labels (the table has
no ontology identifiers), so each is a source term in a namespace of its own.
"""

from __future__ import annotations

from biolink_model.datamodel.model import Gene, OntologyClass, slots
from omnipath_core.measurements import Measurement
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.panglaodb import iter_markers, iter_terms
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import EntityRef, OntologyRelation, Relation
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.PANGLAODB,
    name='PanglaoDB',
    url='https://panglaodb.se/',
    # The site states no terms of use for the marker table; it asks to cite the
    # paper ("How to cite", checked 2026-10-12). The paper is open access (CC BY
    # 4.0), which is the paper's license, not the data's. The legacy registry
    # (resources.json) has the placeholder "Informal".
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.STATIC,  # the table is dated 27 March 2020
    pubmed='30951143',
    topics=(Topic.CELL_TYPES,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED, KnowledgeOrigin.EXPERIMENTAL),
    summary='Marker genes of human and mouse cell types, with their organs and germ layers.',
    description=(
        'PanglaoDB is a web server and database of single cell RNA sequencing '
        'datasets of mouse and human. Its marker table lists, for each cell '
        'type, genes that mark it: a mix of canonical markers from the '
        'literature and genes found expressed in the cell clusters of the '
        'datasets, with the sensitivity and specificity of each marker in human '
        'and mouse and the ubiquitousness index of the gene.'
    ),
)

# ================================== DOWNLOAD ==================================

# The "get tsv file" link of the site ("now a stable link", the site's news).
# The server answers requests from command line clients on its pages with 403
# but serves this file.
download = Download(
    url='https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz',
    filename='PanglaoDB_markers_27_Mar_2020.tsv.gz',
    subfolder='panglaodb',
    large=True,
    ext='gz',
)

# =================================== SCHEMA ===================================

CELL_TYPE = 'panglaodb_cell_type'
ORGAN = 'panglaodb_organ'
GERM_LAYER = 'panglaodb_germ_layer'

# Meaning of the values of a statement. The paper describes the sensitivity as
# a measure of how frequently the marker identifies the cell type uniquely and
# the ubiquitousness index as how often the gene is expressed in cell clusters
# (0 to 1, toward 1 for housekeeping genes); it does not define the specificity.
# The table gives sensitivity and specificity once for each species. The values
# carry the scale of the table, 0 to 1.
SENSITIVITY = 'panglaodb:sensitivity'
SPECIFICITY = 'panglaodb:specificity'
UBIQUITOUSNESS_INDEX = 'panglaodb:ubiquitousness_index'
# 1 where the marker is canonical (from the literature), 0 where the table says
# it is not; the table leaves the others NA, which is no statement.
CANONICAL_MARKER = 'panglaodb:canonical_marker'
CANONICAL = {'1': 'true', '0': 'false'}

f = FieldConfig()


def _measure(field: str):
    def value(row) -> Measurement | None:
        return measurement(row.get(field), source_field=field)

    return value


def _canonical_marker(row) -> str | None:
    flag = row.get('canonical_marker')
    if flag is None:
        return None
    if flag not in CANONICAL:
        raise ValueError(f'Unknown PanglaoDB canonical marker flag: {flag!r}')
    return CANONICAL[flag]


# The gene is the subject: the table has gene symbols. The table writes them all
# in upper case; the parser gives the statement of a mouse gene the symbol in
# the mouse form (``Try4``), see ``parsers/panglaodb.py``.
marker_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=Gene,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f('gene')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=lambda row: f"NCBITaxon:{row['taxon']}"),
        ),
    ),
    # A marker of a cell type is a gene expressed in it: ``expressed_in`` is
    # defined for a gene or gene product and a cell, an anatomical entity (the
    # pinned Biolink model, RO:0002206). The table does not say how specific
    # the expression is; that is the sensitivity and specificity.
    predicate=slots.expressed_in,
    object=EntityBuilder(
        entity_type=OntologyClass,
        identifiers=IdentifiersBuilder(
            CV(term=CELL_TYPE, value=f('cell_type')),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=SENSITIVITY, value=_measure('sensitivity')),
        CV(term=SPECIFICITY, value=_measure('specificity')),
        CV(term=UBIQUITOUSNESS_INDEX, value=_measure('ubiquitousness_index')),
        CV(term=CANONICAL_MARKER, value=_canonical_marker),
    ),
)


def exclusion_reason(row: dict) -> str | None:
    """Why a parsed row gives no statement; None for a row that is mapped."""
    if row.get('taxon') is None:
        return 'species not stated'
    return None


def map_marker(row: dict) -> Relation | None:
    """A marker statement; None where the row states no species."""
    if exclusion_reason(row):
        return None
    return marker_schema(row)


def _term(namespace: str, relations=None) -> EntityBuilder:
    return EntityBuilder(
        entity_type=OntologyClass,
        identifiers=IdentifiersBuilder(
            CV(term=namespace, value=f('id')),
            CV(term=Namespace.NAME, value=f('id')),
        ),
        ontology_relations=relations,
    )


def _cell_type_relations(row) -> list[OntologyRelation]:
    relations = []
    # The table gives each cell type its organ: the cell type is found in, and
    # is a part of, that organ. It is not a kind of the organ (acinar cells are
    # not a pancreas), so the relation is ``part_of``, as BTO and Uberon write
    # it between their terms, and not ``subclass_of``.
    if row.get('organ'):
        relations.append(
            OntologyRelation(
                predicate=slots.part_of,
                object=EntityRef(OntologyClass, ORGAN, row['organ']),
            )
        )
    # The germ layer is a second attribute of the cell type, not a level of
    # the hierarchy: a germ layer holds several organs and an organ cells of
    # two germ layers.
    if row.get('germ_layer'):
        relations.append(
            OntologyRelation(
                predicate=slots.associated_with,
                object=EntityRef(OntologyClass, GERM_LAYER, row['germ_layer']),
            )
        )
    return relations


TERM_SCHEMAS = {
    'cell_type': _term(CELL_TYPE, _cell_type_relations),
    'organ': _term(ORGAN),
    'germ_layer': _term(GERM_LAYER),
}


def map_term(row: dict):
    """The entity of a cell type, an organ or a germ layer."""
    kind = row.get('kind')
    if kind not in TERM_SCHEMAS:
        raise ValueError(f'Unknown PanglaoDB term kind: {kind!r}')
    return TERM_SCHEMAS[kind](row)


# ================================= RESOURCE ===================================

resource = Resource(
    config,
    terms=Dataset(
        download=download,
        mapper=map_term,
        raw_parser=iter_terms,
    ),
    markers=Dataset(
        download=download,
        mapper=map_marker,
        raw_parser=iter_markers,
    ),
)

# ================================= REFERENCE ==================================

# species  official gene symbol  cell type  nicknames  ubiquitousness index
# product description  gene type  canonical marker  germ layer  organ
# sensitivity_human  sensitivity_mouse  specificity_human  specificity_mouse
# e.g. Mm Hs  CTRB1  Acinar cells  CTRB  0.017  chymotrypsinogen B1
#      protein-coding gene  1  Endoderm  Pancreas  1.0  0.957143  0.000628931  0.0159201
