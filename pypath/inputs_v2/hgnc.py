"""
Parse the HGNC gene groups and emit gene families with their members.

The HUGO Gene Nomenclature Committee (HGNC) names human genes and groups them in
gene groups, also called gene families: sets of genes that share structure,
function or origin, such as the immunoglobulin like domain containing genes or the
collagens. A group has member genes and can have subgroups. Each group becomes a
``GeneFamily`` with the HGNC group id, its name, its description and references,
and has as members its genes (by HGNC id and approved symbol) and its subgroups,
nested as the database nests them. A gene is listed with the groups it belongs to
directly, as the database lists it; some genes of a subgroup are also listed in the
group above, and are members of both.

The HGNC complete set, which is a table of the gene identifiers, is not read; the
approved symbols come from a two-column download.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.hgnc import iter_groups
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Association,
    Entity,
    EntityRef,
    Identifier,
    Membership,
)

config = ResourceConfig(
    id=ResourceCv.HGNC,
    name='HGNC',
    full='HUGO Gene Nomenclature Committee',
    url='https://www.genenames.org/',
    # https://www.genenames.org/about/license/: "all data is released under the
    # Creative Commons Public Domain (CC0) License".
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='36243972',
    topics=(Topic.MOLECULAR_CLASSIFICATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    # The 2023 article (NAR): the Symbol Reports "contain data curated by HGNC
    # nomenclature advisors"; the gene groups carry the PubMed references of
    # their curation (read here as publications).
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='HGNC gene groups (gene families) with their member genes and subgroups.',
    description=(
        'The HUGO Gene Nomenclature Committee (HGNC) is the authority that '
        'approves the names and symbols of human genes. It also curates gene '
        'groups, or families, which collect genes by shared structure, function '
        'or evolutionary origin, with a description, references and a hierarchy '
        'of groups and subgroups.'
    ),
)

GROUP = 'hgnc_group'
TAXON = 'NCBITaxon:9606'


def _gene(member: dict[str, Any]) -> Entity:
    identifiers = [Identifier(type=Namespace.HGNC, value=member['hgnc_id'])]
    if member['symbol']:
        identifiers.append(
            Identifier(type=Namespace.GENESYMBOL, value=member['symbol'])
        )
    return Entity(
        type=model.Gene,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def _subgroup(group_id: str) -> Entity:
    return Entity(
        type=model.GeneFamily,
        identifiers=[Identifier(type=GROUP, value=group_id)],
    )


def map_group(row: dict[str, Any]) -> Entity:
    """A gene group with its genes and subgroups as members."""
    identifiers = [
        Identifier(type=GROUP, value=row['id']),
        Identifier(type=Namespace.NAME, value=row['name']),
    ]
    synonyms = [row['abbreviation'], *row['aliases']]
    identifiers += [
        Identifier(type=Namespace.SYNONYM, value=value)
        for value in dict.fromkeys(synonyms)
        if value
    ]
    annotations = []
    if row['description']:
        annotations.append(Annotation(term=slots.description, value=row['description']))
    annotations += [
        Annotation(term=slots.publications, value=pmid_curie(pmid))
        for pmid in row['pmids']
    ]
    associations = None
    if row['go']:
        associations = [
            Association(
                object=EntityRef(model.OntologyClass, Namespace.GO, row['go']),
                predicate=slots.associated_with,
            )
        ]
    return Entity(
        type=model.GeneFamily,
        identifiers=identifiers,
        annotations=annotations or None,
        associations=associations,
        membership=[
            Membership(member=_gene(member), predicate=slots.has_member)
            for member in row['members']
        ]
        + [
            Membership(member=_subgroup(child), predicate=slots.has_member)
            for child in row['children']
        ],
    )


BASE = 'https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/'


def _table(name: str) -> Download:
    return Download(
        url=f'{BASE}{name}.csv',
        filename=f'hgnc_{name}.csv',
        subfolder='hgnc',
        large=True,
    )


# The legacy module read the web download interface of HGNC, which names the groups
# of each gene but has neither group ids nor the hierarchy. The group tables are
# from the HGNC download page ("gene group data files"); the symbols are the
# HGNC ID and Approved symbol columns of the custom download.
symbols_download = Download(
    url=(
        'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&'
        'col=gd_app_sym&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&'
        'format=text&submit=submit'
    ),
    filename='hgnc_approved_symbols.txt',
    subfolder='hgnc',
    large=True,
)

resource = Resource(
    config,
    families=Dataset(
        download=_table('family'),
        mapper=map_group,
        raw_parser=functools.partial(
            iter_groups,
            alias_download=_table('family_alias'),
            hierarchy_download=_table('hierarchy'),
            members_download=_table('gene_has_family'),
            symbols_download=symbols_download,
        ),
    ),
)
