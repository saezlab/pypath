"""
Parse Membranome and emit the membranes that single-pass membrane proteins are in.

Membranome is a database of single-pass (bitopic) transmembrane proteins of 18
species, classified into superfamilies and families, with the membrane each
protein is in and the topology of the protein in it. This module imports the
membranes (the source's own vocabulary of 24) as the ``terms`` dataset and, for each protein, a
``located_in`` relation to every membrane the database lists for it, with the
publication of the evidence and the side of the membrane the protein's topology
is on, where the database gives it.

The legacy reader took one membrane for each protein, the one of the protein
table, and the side: the ``topology_in`` side of the membrane where the table's
``topology_show_in`` flag is true and the ``topology_out`` side where it is
false or empty. Here the side is only given where the flag is set.
"""

from __future__ import annotations

from functools import partial
from typing import Any

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.membranome import iter_localizations, iter_membranes
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Relation

config = ResourceConfig(
    id=ResourceCv.MEMBRANOME,
    name='Membranome',
    url='https://biomembhub.org/membranome/',
    # The pages of the site end with "Copyright 2024-2026 by the Im Lab | License", and
    # the link leads to https://creativecommons.org/licenses/by/3.0/ . (The
    # description of the data API names the Apache License 2.0.)
    license=LicenseCV.CC_BY_3_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='35481632',
    topics=(Topic.SUBCELLULAR_LOCALIZATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED, KnowledgeOrigin.COMPUTATIONALLY_PREDICTED),
    summary=(
        'Membranes of single-pass transmembrane proteins of 18 species, with '
        'the side of the membrane their topology is on.'
    ),
    description=(
        'Membranome is a database of single-pass (bitopic) transmembrane '
        'proteins of human, model organisms, pathogens and other species. It '
        'classifies the proteins into superfamilies and families and gives '
        'for each its membrane, the topology of the protein in it, the domain '
        'organization, the thickness and tilt of the transmembrane helix '
        'computed in the membrane, and links to structures, interactions, '
        'complexes and pathways.'
    ),
)

API = 'https://opm-back.cc.lehigh.edu/membranome-backend'
MEMBRANE_NAMESPACE = 'membranome_membrane'

# Species of the database -> NCBI taxon. The names are those of the species
# table of the API (18); a species that is not here raises.
TAXA = {
    'Arabidopsis thaliana': 'NCBITaxon:3702',
    'Bacillus subtilis': 'NCBITaxon:1423',
    'Caenorhabditis elegans': 'NCBITaxon:6239',
    'Dictyostelium discoideum': 'NCBITaxon:44689',
    'Drosophila melanogaster': 'NCBITaxon:7227',
    'Escherichia coli': 'NCBITaxon:562',
    'Homo sapiens': 'NCBITaxon:9606',
    'Klebsiella pneumoniae': 'NCBITaxon:573',
    'Methanocaldococcus jannaschii': 'NCBITaxon:2190',
    'Monosiga brevicollis': 'NCBITaxon:81824',
    'Mus musculus': 'NCBITaxon:10090',
    'Mycobacterium tuberculosis': 'NCBITaxon:1773',
    'Plasmodium falciparum': 'NCBITaxon:5833',
    'Pseudomonas aeruginosa': 'NCBITaxon:287',
    'Saccharomyces cerevisiae': 'NCBITaxon:4932',
    'Staphylococcus aureus': 'NCBITaxon:1280',
    'Synechocystis sp.': 'NCBITaxon:1143',
    'Trichoplax adhaerens': 'NCBITaxon:10228',
}

# Source attributes of the statement. The side is the name the membrane gives
# to the side the protein's topology is shown on ("cytoplasmic side").
_SIDE = 'membranome:topology_side'
# Whether the row is marked as the protein's main membrane ('true' or 'false').
# The parser reads the file's `true` / `false` into a boolean and raises on any
# other value; the crosswalk to the attribute value is exact.
_MAIN = 'membranome:main_membrane'
_FLAG = {True: 'true', False: 'false'}


def _protein(row: dict[str, Any]) -> Entity:
    identifiers = [Identifier(type=Namespace.UNIPROT, value=row['accession'])]
    if row['entry_name']:
        identifiers.append(Identifier(type=Namespace.UNIPROT_ENTRY, value=row['entry_name']))
    if row['gene']:
        identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=row['gene']))
    identifiers.append(Identifier(type=Namespace.NAME, value=row['name']))
    if row['species'] not in TAXA:
        raise ValueError(f'Unknown Membranome species: {row["species"]!r}')
    return Entity(
        type=Protein,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXA[row['species']])],
    )


def map_localization(row: dict[str, Any]) -> Relation:
    """A protein located in a membrane, on a side of it where the topology is given."""
    annotations = [Annotation(term=_MAIN, value=_FLAG[row['main']])]
    if row['side']:
        annotations.append(Annotation(term=_SIDE, value=row['side']))
    if row['pmid']:
        curie = pmid_curie(row['pmid'])
        if curie is None:
            raise ValueError(f'Not a PubMed identifier: {row["pmid"]!r}')
        annotations.append(Annotation(term=slots.publications, value=curie))
    return Relation(
        subject=_protein(row),
        predicate=slots.located_in,
        object=Entity(
            type=OntologyClass,
            identifiers=[Identifier(type=MEMBRANE_NAMESPACE, value=row['membrane'])],
        ),
        annotations=annotations,
    )


def map_membrane(row: dict[str, Any]) -> Entity:
    """A membrane of the database."""
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=MEMBRANE_NAMESPACE, value=row['name']),
            Identifier(type=Namespace.NAME, value=row['name']),
            Identifier(type=Namespace.SYNONYM, value=row['abbreviation']),
        ],
    )


membrane_download = Download(
    url=f'{API}/membranes',
    filename='membranome_membranes.json',
    subfolder='membranome',
    large=True,
)
protein_download = Download(
    url=f'{API}/proteins?fileFormat=csv',
    filename='membranome_proteins.csv',
    subfolder='membranome',
    large=True,
)
localization_download = Download(
    url=f'{API}/protein_membranes?fileFormat=csv',
    filename='membranome_protein_membranes.csv',
    subfolder='membranome',
    large=True,
)

resource = Resource(
    config,
    terms=Dataset(
        download=membrane_download,
        mapper=map_membrane,
        raw_parser=iter_membranes,
    ),
    localizations=Dataset(
        download=localization_download,
        mapper=map_localization,
        raw_parser=partial(
            iter_localizations,
            membrane_download=membrane_download,
            protein_download=protein_download,
        ),
    ),
)
