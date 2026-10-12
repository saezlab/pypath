"""
Parse the complexes of the AlzPathway map and emit Entity records.

AlzPathway is a CellDesigner map of the signaling pathways of Alzheimer's
disease: a reaction network in the process description notation, published
as Additional file 1 of Mizuno et al. 2012 (PMID 22647208). alzpathway.org
describes the map and does not serve it; the file is read from the PubMed
Central open access dataset.

The legacy input read ``alzpw-ppi.csv``, a table that the pypath authors
derived from the map: the proteins of each complex expanded to pairs, with
UniProt accessions translated from the names. It is not a file of the
source. This module reads the same information from the map itself: each
complex species with its components, as the map draws them, and nothing
expanded to pairs. The reactions of the map and their catalysts and
inhibitors are not ported.

Identity of the components. The authors name genes, RNAs and proteins by
HGNC symbols; a name written like a symbol is given as a gene symbol, any
other name (``AmyloidB``, ``mutant PSEN``) as a name. A UniProt accession or
PubChem identifier is given when the map annotates the molecule with exactly
one; many molecules have several (``WNT``) or none.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import (
    ChemicalEntity,
    Gene,
    MacromolecularComplex,
    Protein,
    RNAProduct,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.alzpathway import iter_complexes
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)

config = ResourceConfig(
    id=ResourceCv.ALZPATHWAY,
    name='AlzPathway',
    url='http://alzpathway.org/',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='22647208',
    topics=(Topic.SIGNALING, Topic.DISEASE_PHENOTYPE),
    content=(Content.PROTEIN_COMPLEXES,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Complexes of the AlzPathway map of Alzheimer\'s disease signaling, '
        'with their components and publications.'
    ),
    description=(
        'AlzPathway is a manually curated map of the signaling pathways of '
        'Alzheimer\'s disease, built with CellDesigner from more than 100 '
        'review articles. It describes about 1,350 molecules and 1,070 '
        'reactions in neurons, the blood-brain barrier, synapses, astrocytes '
        'and microglia, with the publications that support each reaction.'
    ),
)

download_map = Download(
    # Additional file 1 (alzpathway_sbml_map.xml) of PMID 22647208, in the
    # PubMed Central open access dataset.
    url='https://pmc-oa-opendata.s3.amazonaws.com/PMC3411424.1/1752-0509-6-52-S1.xml',
    filename='1752-0509-6-52-S1.xml',
    subfolder='alzpathway',
    encoding=None,
    default_mode='rb',
)

# The map describes human pathways: HGNC symbols and human UniProt entries.
TAXON = 'NCBITaxon:9606'

# CellDesigner class of a component -> entity type. A component of class
# UNKNOWN is a placeholder (``Unknown receptor``, ``unknown protein``): it
# names no molecule and gives no member.
MEMBER_TYPES = {
    'PROTEIN': Protein,
    'GENE': Gene,
    'RNA': RNAProduct,
    'SIMPLE_MOLECULE': ChemicalEntity,
    'ION': ChemicalEntity,
    'COMPLEX': MacromolecularComplex,
    'UNKNOWN': None,
}
_GENE_PRODUCTS = (Protein, Gene, RNAProduct)

# The form of an HGNC symbol.
_symbol = re.compile(r'[A-Z][A-Z0-9]*(?:orf[0-9]+)?(?:-[A-Z0-9]+)*')
_uniprot = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
)
_number = re.compile(r'[1-9]\d*')
# A complex without a name has its species id (``s78``) as its name.
_species_id = re.compile(r's\d+')


def _single(member: dict[str, Any], collection: str, pattern: re.Pattern) -> str | None:
    """The identifier of a collection, if the map gives the molecule exactly one.

    The annotations also hold search queries (``?query=ribosomal+protein``).
    """
    identifiers = [
        identifier
        for identifier in (member.get('xrefs') or {}).get(collection, ())
        if pattern.fullmatch(identifier)
    ]
    return identifiers[0] if len(identifiers) == 1 else None


def _name_identifier(name: str, symbol: bool) -> list[Identifier]:
    name = (name or '').strip()
    if not name or _species_id.fullmatch(name):
        return []
    if symbol and _symbol.fullmatch(name):
        return [Identifier(type=Namespace.GENESYMBOL, value=name)]
    return [Identifier(type=Namespace.NAME, value=name)]


def _member_entity(member: dict[str, Any]) -> Entity | None:
    cell_designer_class = member.get('class')
    if cell_designer_class not in MEMBER_TYPES:
        raise ValueError(
            f'Unknown CellDesigner class of a complex component: {cell_designer_class!r}'
        )
    entity_type = MEMBER_TYPES[cell_designer_class]
    if entity_type is None:
        return None
    if entity_type is MacromolecularComplex:
        return _complex_entity(member)
    identifiers = []
    if entity_type is Protein:
        accession = _single(member, 'uniprot', _uniprot)
        if accession:
            identifiers.append(Identifier(type=Namespace.UNIPROT, value=accession))
    elif entity_type is ChemicalEntity:
        for collection, namespace in (
            ('pubchem.compound', Namespace.PUBCHEM),
            ('pubchem.substance', Namespace.PUBCHEM_SUBSTANCE),
        ):
            number = _single(member, collection, _number)
            if number:
                identifiers.append(Identifier(type=namespace, value=number))
    gene_product = entity_type in _GENE_PRODUCTS
    identifiers.extend(_name_identifier(member.get('name'), symbol=gene_product))
    if not identifiers:
        return None
    return Entity(
        type=entity_type,
        identifiers=identifiers,
        annotations=(
            [Annotation(term=slots.in_taxon, value=TAXON)] if gene_product else None
        ),
    )


def _memberships(row: dict[str, Any]) -> list[Membership]:
    memberships = []
    for member in row.get('members') or ():
        entity = _member_entity(member)
        if entity is None:
            continue
        copies = member.get('homodimer')
        memberships.append(
            Membership(
                member=entity,
                # CellDesigner draws a multimer with its number of copies.
                annotations=(
                    [Annotation(term=slots.stoichiometry, value=str(copies))]
                    if copies
                    else None
                ),
            )
        )
    return memberships


def _complex_entity(
    row: dict[str, Any],
    annotations: list[Annotation] | None = None,
) -> Entity | None:
    memberships = _memberships(row)
    if not memberships:
        return None
    return Entity(
        type=MacromolecularComplex,
        identifiers=_name_identifier(row.get('name'), symbol=False),
        annotations=annotations or None,
        membership=memberships,
    )


def map_complex(row: dict[str, Any]) -> Entity | None:
    """A complex of the map, with the publications of the reactions that form it."""
    pmids: list[str] = []
    for reaction in row.get('producing_reactions') or ():
        for pmid in (reaction.get('xrefs') or {}).get('pubmed', ()):
            if is_pmid(pmid) and pmid not in pmids:
                pmids.append(pmid)
    return _complex_entity(
        row,
        [Annotation(term=slots.publications, value=f'PMID:{pmid}') for pmid in pmids],
    )


resource = Resource(
    config,
    complexes=Dataset(
        download=download_map,
        mapper=map_complex,
        raw_parser=iter_complexes,
    ),
)
