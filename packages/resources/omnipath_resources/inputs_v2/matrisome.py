"""
Parse the matrisome master lists and emit Entity and Relation records.

The matrisome is the ensemble of the genes that encode the extracellular matrix
(ECM) proteins and the proteins that are associated with it. Naba et al. defined
it in silico from the domain organisation of the proteins and then by hand, for
human and mouse, and the Matrisome Portal publishes it as one master list for each
species. A gene is in a division (the core matrisome, or the matrisome-associated
proteins) and in a category of it:

    Core matrisome          Collagens, ECM Glycoproteins, Proteoglycans
    Matrisome-associated    ECM Regulators, ECM-affiliated Proteins, Secreted Factors

The divisions and categories are terms of the matrisome (the ``terms`` dataset), with
the category below its division. Each gene of the human and mouse lists (datasets
``human`` and ``mouse``) is associated with the category the list gives it; the
comment of the list on the gene (``Semaphorin``, ``Basement Membrane``) is an
attribute of that association. The lists are the Google Sheets that the portal links
as the download of each species.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.matrisome import iter_genes, iter_terms
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    EntityRef,
    Identifier,
    OntologyRelation,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.MATRISOME,
    name='Matrisome',
    full='The Matrisome Project',
    url='https://sites.google.com/uic.edu/matrisome/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='22159717',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary=(
        'Human and mouse genes of the matrisome (extracellular matrix proteins '
        'and associated proteins) with their division and category.'
    ),
    description=(
        'The matrisome is the ensemble of genes encoding the extracellular '
        'matrix (ECM) proteins and the proteins that are associated with the '
        'ECM or can modify it. Naba et al. defined it in silico, from the '
        'domain-based organisation of the ECM proteins, and classified it into '
        'the core matrisome (collagens, ECM glycoproteins and proteoglycans) '
        'and the matrisome-associated proteins (ECM regulators, ECM-affiliated '
        'proteins and secreted factors). The Matrisome Portal of the Naba lab '
        'publishes the master lists of the species.'
    ),
)

DIVISION_NAMESPACE = 'matrisome_division'
CATEGORY_NAMESPACE = 'matrisome_category'
NOTE = 'matrisome:note'
# The mouse list names its genes by MGI id (``MGI:1918645``).


def map_term(row: dict[str, Any]) -> Entity:
    """A division or category, a category below its division."""
    namespaces = {
        'division': DIVISION_NAMESPACE,
        'category': CATEGORY_NAMESPACE,
    }
    namespace = namespaces[row['kind']]
    relations = None
    if row['parent']:
        relations = [
            OntologyRelation(
                predicate=slots.subclass_of,
                object=EntityRef(model.OntologyClass, DIVISION_NAMESPACE, row['parent']),
                ontology_id='matrisome',
            )
        ]
    return Entity(
        type=model.OntologyClass,
        identifiers=[
            Identifier(type=namespace, value=row['term']),
            Identifier(type=Namespace.NAME, value=row['term']),
        ],
        ontology_relations=relations,
    )


def map_gene(row: dict[str, Any]) -> Relation:
    """A gene associated with its category; the list's comment is an attribute."""
    gene_namespace = Namespace.MGI if row['taxon'] == '10090' else Namespace.HGNC
    identifiers = [
        Identifier(type=gene_namespace, value=row['gene_id']),
        Identifier(type=Namespace.GENESYMBOL, value=row['symbol']),
    ]
    if row['name']:
        identifiers.append(Identifier(type=Namespace.NAME, value=row['name']))
    identifiers += [
        Identifier(type=Namespace.SYNONYM, value=synonym)
        for synonym in row['synonyms']
    ]
    annotations = []
    if row['note']:
        annotations.append(Annotation(term=NOTE, value=row['note']))
    return Relation(
        subject=Entity(
            type=model.Gene,
            identifiers=identifiers,
            annotations=[
                Annotation(term=slots.in_taxon, value=f'NCBITaxon:{row["taxon"]}')
            ],
        ),
        predicate=slots.associated_with,
        object=EntityRef(model.OntologyClass, CATEGORY_NAMESPACE, row['category']),
        annotations=annotations,
    )


# The sheets that the Matrisome Portal links under "Download the complete ...
# matrisome list (rev. 2014)" on its pages for Homo sapiens
# (sites.google.com/uic.edu/matrisome/matrisome-annotations/homo-sapiens) and Mus
# musculus. The legacy URLs (matrisomeproject.mit.edu) no longer answer; the rescued
# copy has the human list of 2014 only, with pseudogenes the portal has since removed.
SHEET_URL = 'https://docs.google.com/spreadsheets/d/%s/export?format=xlsx'
human_download = Download(
    url=SHEET_URL % '12OslOlOhQQtHMOUKAFFgAQXL2zYX-g9q',
    filename='Hs_Matrisome_Masterlist_Nabaetal_2012.xlsx',
    subfolder='matrisome',
    large=True,
    encoding=None,
    default_mode='rb',
)
mouse_download = Download(
    url=SHEET_URL % '15ijmTLYuCXCc19ZBfNoCIxbxcs9tCsgY',
    filename='Mm_Matrisome_Masterlist_Nabaetal_2012.xlsx',
    subfolder='matrisome',
    large=True,
    encoding=None,
    default_mode='rb',
)

resource = Resource(
    config,
    terms=Dataset(
        download=human_download,
        mapper=map_term,
        raw_parser=functools.partial(iter_terms, mouse_download=mouse_download),
    ),
    human=Dataset(
        download=human_download,
        mapper=map_gene,
        raw_parser=functools.partial(iter_genes, species='human'),
    ),
    mouse=Dataset(
        download=mouse_download,
        mapper=map_gene,
        raw_parser=functools.partial(iter_genes, species='mouse'),
    ),
)
