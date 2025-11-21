#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""
Parse Reactome BioPAX data and emit Entity records.

This module downloads Reactome BioPAX (OWL) data and converts it into Entity
records using the declarative schema pattern from tabular_builder.
"""

from __future__ import annotations

import pickle
from collections.abc import Generator
from pathlib import Path

from rdflib import Graph, Namespace
from rdflib.namespace import RDF

from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceAnnotationCv,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    Column,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    Map,
    Member,
    MembershipBuilder,
    MembersFromList,
)
from pypath.share.downloads import download_and_open
from pypath.share import settings


# BioPAX namespace
BP = Namespace("http://www.biopax.org/release/biopax-level3.owl#")

# Module-level cache for parsed data (keyed by data_type)
_DATA_CACHE: dict[str, list[dict]] = {}

# Mapping of BioPAX PhysicalEntity types to EntityTypeCv terms
PHYSICAL_ENTITY_TYPE_MAP = {
    'smallmolecule': EntityTypeCv.SMALL_MOLECULE,
    'protein': EntityTypeCv.PROTEIN,
    'gene': EntityTypeCv.GENE,
    'complex': EntityTypeCv.PROTEIN_COMPLEX,
    'complexassembly': EntityTypeCv.PROTEIN_COMPLEX,
    'dna': EntityTypeCv.DNA,
    'dnaregion': EntityTypeCv.DNA,
    'rna': EntityTypeCv.RNA,
    'rnaregion': EntityTypeCv.RNA,
    'physicalentity': EntityTypeCv.PHYSICAL_ENTITY,
}

# Mapping of BioPAX EntityReference subtypes to EntityTypeCv terms
ENTITY_REFERENCE_TYPE_MAP = {
    'proteinreference': EntityTypeCv.PROTEIN,
    'smallmoleculereference': EntityTypeCv.SMALL_MOLECULE,
    'dnareference': EntityTypeCv.DNA,
    'rnareference': EntityTypeCv.RNA,
}

# EntityReference types to process
ENTITY_REFERENCE_TYPES = [
    'ProteinReference',
    'SmallMoleculeReference',
    'DnaReference',
    'RnaReference',
]


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing Reactome metadata.
    """
    yield Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.REACTOME),
            Identifier(type=IdentifierNamespaceCv.NAME, value='Reactome'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='34788843'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://reactome.org/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'Reactome is a free, open-source, curated and peer-reviewed '
                'pathway database. It provides intuitive bioinformatics tools '
                'for the visualization, interpretation and analysis of pathway '
                'knowledge to support basic and clinical research, genome analysis, '
                'modeling, systems biology and education.'
            )),
        ],
    )


def _get_cache_path(data_type: str) -> Path:
    """Get the path for the cached pickle file."""
    cache_dir = Path(settings.get('cachedir'))
    return cache_dir / 'reactome' / f'{data_type}.pkl'


def _load_cached_data(data_type: str, force_refresh: bool = False) -> list[dict] | None:
    """
    Load cached parsed data from pickle file.

    Args:
        data_type: Type of data (e.g., 'entity_references', 'reactions')
        force_refresh: If True, ignore cache

    Returns:
        List of dictionaries, or None if not cached
    """
    if force_refresh:
        return None

    # Check memory cache first
    if data_type in _DATA_CACHE:
        return _DATA_CACHE[data_type]

    # Try loading from pickle
    pickle_path = _get_cache_path(data_type)
    if pickle_path.exists():
        try:
            with open(pickle_path, 'rb') as f:
                data = pickle.load(f)
            _DATA_CACHE[data_type] = data
            return data
        except Exception:
            pass

    return None


def _save_cached_data(data_type: str, data: list[dict]) -> None:
    """Save parsed data to pickle cache."""
    _DATA_CACHE[data_type] = data

    pickle_path = _get_cache_path(data_type)
    try:
        pickle_path.parent.mkdir(parents=True, exist_ok=True)
        with open(pickle_path, 'wb') as f:
            pickle.dump(data, f)
    except Exception:
        pass  # Caching is best-effort


def _ensure_all_caches_populated(species: str = 'Homo_sapiens', force_refresh: bool = False) -> bool:
    """
    Ensure all data caches are populated.

    Parses the RDF graph once and populates all caches (entity_references,
    reactions, pathways, controls) in a single pass.

    Args:
        species: Species name for the BioPAX file
        force_refresh: If True, force re-parsing

    Returns:
        True if caches are populated, False if loading failed
    """
    data_types = ['entity_references', 'reactions', 'pathways', 'controls']

    # Check if all caches are already populated
    if not force_refresh:
        all_cached = all(_load_cached_data(dt) is not None for dt in data_types)
        if all_cached:
            return True

    # Need to parse the RDF graph
    g = _load_biopax_graph(species, force_refresh)
    if g is None:
        return False

    # Parse all data types and cache them
    if _load_cached_data('entity_references', force_refresh) is None:
        data = list(_iterate_entity_references(g, max_records=None))
        _save_cached_data('entity_references', data)

    if _load_cached_data('reactions', force_refresh) is None:
        data = list(_iterate_reactions(g, max_records=None))
        _save_cached_data('reactions', data)

    if _load_cached_data('pathways', force_refresh) is None:
        data = list(_iterate_pathways(g, max_records=None))
        _save_cached_data('pathways', data)

    if _load_cached_data('controls', force_refresh) is None:
        data = list(_iterate_controls(g, max_records=None))
        _save_cached_data('controls', data)

    return True


def _load_biopax_graph(
    species: str = 'Homo_sapiens',
    force_refresh: bool = False,
) -> Graph | None:
    """
    Download and load a Reactome BioPAX file into an RDF graph.

    Args:
        species: Species name for the BioPAX file (e.g., 'Homo_sapiens')
        force_refresh: If True, force redownload of the data

    Returns:
        RDF Graph containing the BioPAX data, or None if loading fails
    """
    url = 'https://reactome.org/download/current/biopax.zip'
    opener = download_and_open(
        url,
        filename='reactome_biopax.zip',
        subfolder='reactome',
        large=True,
        ext='zip',
        default_mode='rb',
        force=force_refresh,
    )

    if not opener or not opener.result:
        return None

    # Find the OWL file for the specified species
    owl_file = None
    target_filename = f'{species}.owl'
    for filename, file_handle in opener.result.items():
        if target_filename in filename:
            owl_file = file_handle
            break

    if not owl_file:
        return None

    # Parse the OWL file into an RDF graph
    g = Graph()
    g.parse(owl_file, format='xml')
    return g


def _extract_xrefs(g: Graph, uri, bp_ns: Namespace) -> dict[str, list[str]]:
    """
    Extract cross-references from a BioPAX resource.

    Returns a dict mapping database names to lists of IDs.
    """
    xrefs: dict[str, list[str]] = {}

    for xref in g.objects(uri, bp_ns.xref):
        db = g.value(xref, bp_ns.db)
        id_value = g.value(xref, bp_ns.id)

        if db and id_value:
            db_str = str(db).lower()
            id_str = str(id_value)

            if 'reactome' in db_str and 'pubmed' not in db_str:
                key = 'reactome_stable_id' if ('stable' in db_str or 'R-' in id_str) else 'reactome_id'
                xrefs.setdefault(key, []).append(id_str)
            elif 'uniprot' in db_str:
                xrefs.setdefault('uniprot', []).append(id_str)
            elif 'chebi' in db_str:
                xrefs.setdefault('chebi', []).append(id_str)
            elif 'pubchem' in db_str or 'compound' in db_str:
                xrefs.setdefault('pubchem_compound', []).append(id_str)
            elif 'kegg' in db_str:
                xrefs.setdefault('kegg', []).append(id_str)
            elif 'pubmed' in db_str:
                xrefs.setdefault('pubmed', []).append(id_str)
            elif 'gene ontology' in db_str or id_str.startswith('GO:'):
                xrefs.setdefault('go', []).append(id_str)
            elif 'taxonomy' in db_str:
                xrefs.setdefault('ncbi_taxonomy', []).append(id_str)

    return xrefs


def _extract_names(g: Graph, uri, bp_ns: Namespace) -> dict[str, str | list[str]]:
    """
    Extract names from a BioPAX resource.
    """
    names: dict[str, str | list[str]] = {}

    display_name = g.value(uri, bp_ns.displayName)
    if display_name:
        names['display_name'] = str(display_name)

    standard_name = g.value(uri, bp_ns.standardName)
    if standard_name:
        names['standard_name'] = str(standard_name)

    synonyms = [str(name) for name in g.objects(uri, bp_ns.name)]
    if synonyms:
        names['synonyms'] = synonyms

    return names


def _extract_organism(g: Graph, uri, bp_ns: Namespace) -> dict[str, str]:
    """
    Extract organism information from a BioPAX resource.
    """
    organism_info: dict[str, str] = {}
    organism = g.value(uri, bp_ns.organism)

    if organism:
        organism_name = g.value(organism, bp_ns.name)
        if organism_name:
            organism_info['organism_name'] = str(organism_name)

        for xref in g.objects(organism, bp_ns.xref):
            db = g.value(xref, bp_ns.db)
            id_value = g.value(xref, bp_ns.id)
            if db and id_value and 'taxonomy' in str(db).lower():
                organism_info['ncbi_taxonomy_id'] = str(id_value)
                break

    return organism_info


def _extract_cellular_location(g: Graph, uri, bp_ns: Namespace) -> str | None:
    """
    Extract cellular location from a BioPAX resource.
    """
    location = g.value(uri, bp_ns.cellularLocation)
    if location:
        location_term = g.value(location, bp_ns.term)
        if location_term:
            return str(location_term)
    return None


def _iterate_entity_references(
    g: Graph,
    max_records: int | None = None,
) -> Generator[dict[str, str | list[str]], None, None]:
    """
    Iterate through BioPAX EntityReference elements and yield dictionaries.

    Args:
        g: RDF Graph containing BioPAX data
        max_records: Maximum number of records to parse

    Yields:
        Dictionary representations of entity references
    """
    count = 0

    for reference_type in ENTITY_REFERENCE_TYPES:
        if max_records is not None and count >= max_records:
            break

        query = f"""
        PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

        SELECT ?entity_ref
        WHERE {{
            ?entity_ref rdf:type bp:{reference_type} .
        }}
        """

        for row in g.query(query):
            if max_records is not None and count >= max_records:
                break

            uri = row.entity_ref
            ref_type_str = reference_type.lower()
            entity_type = ENTITY_REFERENCE_TYPE_MAP.get(ref_type_str, EntityTypeCv.PHYSICAL_ENTITY)

            # Extract all data
            names = _extract_names(g, uri, BP)
            xrefs = _extract_xrefs(g, uri, BP)
            organism = _extract_organism(g, uri, BP)

            # Build flat dictionary
            record = {
                'uri': str(uri),
                'entity_type': entity_type.value if hasattr(entity_type, 'value') else str(entity_type),
                'display_name': names.get('display_name', ''),
                'standard_name': names.get('standard_name', ''),
                'synonyms': ';'.join(names.get('synonyms', [])),
                'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
                'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
                'uniprot': ';'.join(xrefs.get('uniprot', [])),
                'chebi': ';'.join(xrefs.get('chebi', [])),
                'pubchem_compound': ';'.join(xrefs.get('pubchem_compound', [])),
                'kegg': ';'.join(xrefs.get('kegg', [])),
                'pubmed': ';'.join(xrefs.get('pubmed', [])),
                'go': ';'.join(xrefs.get('go', [])),
                'organism_name': organism.get('organism_name', ''),
                'ncbi_taxonomy_id': organism.get('ncbi_taxonomy_id', ''),
            }

            yield record
            count += 1


def _get_entity_type_cv(type_str: str) -> str:
    """Map entity type string to CV value."""
    mapping = {
        'protein': EntityTypeCv.PROTEIN.value,
        'small_molecule': EntityTypeCv.SMALL_MOLECULE.value,
        'dna': EntityTypeCv.DNA.value,
        'rna': EntityTypeCv.RNA.value,
    }
    return mapping.get(type_str, EntityTypeCv.PHYSICAL_ENTITY.value)


def reactome_entity_references(
    species: str = 'Homo_sapiens',
    max_records: int | None = None,
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
    """
    Download and parse Reactome EntityReference data as Entity records.

    EntityReferences represent abstract biological entities like proteins,
    small molecules, DNA, and RNA without specific cellular context.

    Args:
        species: Species name (e.g., 'Homo_sapiens')
        max_records: Maximum number of records to parse
        force_refresh: If True, force redownload of the data

    Yields:
        Entity records representing biological entities
    """
    # Ensure all caches are populated (single RDF parse for all data types)
    if not _ensure_all_caches_populated(species, force_refresh):
        return

    cached_data = _DATA_CACHE['entity_references']

    # Define the schema for entity references
    schema = EntityBuilder(
        entity_type=Map(col=Column('entity_type'), extract=[_get_entity_type_cv]),
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.NAME, value=Column('display_name')),
            CV(term=IdentifierNamespaceCv.SYNONYM, value=Column('synonyms', delimiter=';')),
            CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=Column('reactome_stable_id', delimiter=';')),
            CV(term=IdentifierNamespaceCv.REACTOME_ID, value=Column('reactome_id', delimiter=';')),
            CV(term=IdentifierNamespaceCv.UNIPROT, value=Column('uniprot', delimiter=';')),
            CV(term=IdentifierNamespaceCv.CHEBI, value=Column('chebi', delimiter=';')),
            CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=Column('pubchem_compound', delimiter=';')),
            CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=Column('kegg', delimiter=';')),
            CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=Column('go', delimiter=';')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=IdentifierNamespaceCv.PUBMED, value=Column('pubmed', delimiter=';')),
            CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=Column('ncbi_taxonomy_id')),
            CV(term='organism', value=Column('organism_name')),
        ),
    )

    # Apply max_records limit when iterating
    for i, record in enumerate(cached_data):
        if max_records is not None and i >= max_records:
            break
        yield schema(record)


def _extract_reaction_participants(
    g: Graph,
    reaction_uri,
    bp_ns: Namespace,
    role: str,
    predicate,
) -> list[dict]:
    """
    Extract reaction participants (reactants or products).
    """
    participants = []

    # Build stoichiometry map
    stoich_map = {}
    for stoich in g.objects(reaction_uri, bp_ns.participantStoichiometry):
        entity = g.value(stoich, bp_ns.physicalEntity)
        coeff = g.value(stoich, bp_ns.stoichiometricCoefficient)
        if entity and coeff:
            stoich_map[str(entity)] = str(coeff)

    for molecule in g.objects(reaction_uri, predicate):
        names = _extract_names(g, molecule, bp_ns)
        xrefs = _extract_xrefs(g, molecule, bp_ns)

        # Determine entity type
        entity_type_uri = g.value(molecule, RDF.type)
        entity_type_str = str(entity_type_uri).split('#')[-1].lower() if entity_type_uri else 'physicalentity'
        entity_type = PHYSICAL_ENTITY_TYPE_MAP.get(entity_type_str, EntityTypeCv.PHYSICAL_ENTITY)

        participant = {
            'role': role,
            'entity_type': entity_type.value if hasattr(entity_type, 'value') else str(entity_type),
            'display_name': names.get('display_name', ''),
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'uniprot': ';'.join(xrefs.get('uniprot', [])),
            'chebi': ';'.join(xrefs.get('chebi', [])),
            'stoichiometry': stoich_map.get(str(molecule), ''),
        }
        participants.append(participant)

    return participants


def _iterate_reactions(
    g: Graph,
    max_records: int | None = None,
) -> Generator[dict, None, None]:
    """
    Iterate through BioPAX reaction elements and yield dictionaries.

    Args:
        g: RDF Graph containing BioPAX data
        max_records: Maximum number of records to parse

    Yields:
        Dictionary representations of reactions
    """
    reaction_types = [
        ('BiochemicalReaction', EntityTypeCv.REACTION),
        ('Degradation', EntityTypeCv.DEGRADATION),
        ('TemplateReaction', EntityTypeCv.REACTION),
    ]

    count = 0

    for reaction_type, entity_type_cv in reaction_types:
        if max_records is not None and count >= max_records:
            break

        query = f"""
        PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

        SELECT ?reaction
        WHERE {{
            ?reaction rdf:type bp:{reaction_type} .
        }}
        """

        for row in g.query(query):
            if max_records is not None and count >= max_records:
                break

            uri = row.reaction
            names = _extract_names(g, uri, BP)
            xrefs = _extract_xrefs(g, uri, BP)

            # Get EC number for biochemical reactions
            ec_number = g.value(uri, BP.eCNumber)

            # Get conversion direction
            direction = g.value(uri, BP.conversionDirection)

            # Extract participants
            reactants = _extract_reaction_participants(g, uri, BP, 'reactant', BP.left)
            products = _extract_reaction_participants(g, uri, BP, 'product', BP.right)

            # For TemplateReaction, also get template
            templates = []
            if reaction_type == 'TemplateReaction':
                template = g.value(uri, BP.template)
                if template:
                    template_names = _extract_names(g, template, BP)
                    template_xrefs = _extract_xrefs(g, template, BP)
                    templates.append({
                        'role': 'template',
                        'display_name': template_names.get('display_name', ''),
                        'reactome_stable_id': ';'.join(template_xrefs.get('reactome_stable_id', [])),
                    })

            record = {
                'uri': str(uri),
                'reaction_type': reaction_type,
                'entity_type': entity_type_cv.value,
                'display_name': names.get('display_name', ''),
                'synonyms': ';'.join(names.get('synonyms', [])),
                'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
                'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
                'pubmed': ';'.join(xrefs.get('pubmed', [])),
                'ec_number': str(ec_number) if ec_number else '',
                'direction': str(direction) if direction else '',
                'reactants': reactants,
                'products': products,
                'templates': templates,
            }

            yield record
            count += 1


def _build_member_entities(participants: list[dict]) -> list[Entity]:
    """
    Build Entity objects for reaction participants.
    """
    entities = []
    for p in participants:
        entity_type_str = p.get('entity_type', 'physical_entity')
        # Map string back to CV
        type_mapping = {v.value: v for v in EntityTypeCv}
        entity_type = type_mapping.get(entity_type_str, EntityTypeCv.PHYSICAL_ENTITY)

        identifiers = []
        if p.get('display_name'):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=p['display_name']))
        if p.get('reactome_stable_id'):
            for rid in p['reactome_stable_id'].split(';'):
                if rid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))
        if p.get('uniprot'):
            for uid in p['uniprot'].split(';'):
                if uid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.UNIPROT, value=uid))
        if p.get('chebi'):
            for cid in p['chebi'].split(';'):
                if cid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.CHEBI, value=cid))

        annotations = [Annotation(term='role', value=p.get('role', ''))]
        if p.get('stoichiometry'):
            annotations.append(Annotation(term='stoichiometric_coefficient', value=p['stoichiometry']))

        entities.append((
            Entity(type=entity_type, identifiers=identifiers if identifiers else None),
            annotations,
        ))

    return entities


def reactome_reactions(
    species: str = 'Homo_sapiens',
    max_records: int | None = None,
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
    """
    Download and parse Reactome reaction data as Entity records.

    This includes BiochemicalReactions, Degradations, and TemplateReactions.

    Args:
        species: Species name (e.g., 'Homo_sapiens')
        max_records: Maximum number of records to parse
        force_refresh: If True, force redownload of the data

    Yields:
        Entity records representing reactions with their participants
    """
    from pypath.internals.silver_schema import Membership

    # Ensure all caches are populated (single RDF parse for all data types)
    if not _ensure_all_caches_populated(species, force_refresh):
        return

    cached_data = _DATA_CACHE['reactions']

    for i, record in enumerate(cached_data):
        if max_records is not None and i >= max_records:
            break
        # Build identifiers
        identifiers = []
        if record.get('display_name'):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=record['display_name']))
        if record.get('synonyms'):
            for syn in record['synonyms'].split(';'):
                if syn:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))
        if record.get('reactome_stable_id'):
            for rid in record['reactome_stable_id'].split(';'):
                if rid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))
        if record.get('reactome_id'):
            for rid in record['reactome_id'].split(';'):
                if rid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_ID, value=rid))

        # Build annotations
        annotations = []
        if record.get('pubmed'):
            for pmid in record['pubmed'].split(';'):
                if pmid:
                    annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED, value=pmid))
        if record.get('ec_number'):
            annotations.append(Annotation(term='ec_number', value=record['ec_number']))
        if record.get('direction'):
            annotations.append(Annotation(term='conversion_direction', value=record['direction']))

        # Build membership from participants
        membership = []
        all_participants = record.get('reactants', []) + record.get('products', []) + record.get('templates', [])

        for entity, member_annotations in _build_member_entities(all_participants):
            membership.append(Membership(
                member=entity,
                is_parent=False,
                annotations=member_annotations,
            ))

        # Get entity type
        type_mapping = {v.value: v for v in EntityTypeCv}
        entity_type = type_mapping.get(record.get('entity_type', ''), EntityTypeCv.REACTION)

        yield Entity(
            type=entity_type,
            identifiers=identifiers if identifiers else None,
            annotations=annotations if annotations else None,
            membership=membership if membership else None,
        )


def _iterate_pathways(
    g: Graph,
    max_records: int | None = None,
) -> Generator[dict, None, None]:
    """
    Iterate through BioPAX Pathway elements and yield dictionaries.

    Args:
        g: RDF Graph containing BioPAX data
        max_records: Maximum number of records to parse

    Yields:
        Dictionary representations of pathways
    """
    query = """
    PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

    SELECT ?pathway
    WHERE {
        ?pathway rdf:type bp:Pathway .
    }
    """

    count = 0
    for row in g.query(query):
        if max_records is not None and count >= max_records:
            break

        uri = row.pathway
        names = _extract_names(g, uri, BP)
        xrefs = _extract_xrefs(g, uri, BP)
        organism = _extract_organism(g, uri, BP)

        # Extract comments
        comments = []
        descriptions = []
        for comment in g.objects(uri, BP.comment):
            comment_str = str(comment)
            if "Reactome DB_ID:" not in comment_str:
                if any(comment_str.startswith(prefix) for prefix in ['Reviewed:', 'Authored:', 'Edited:']):
                    comments.append(comment_str)
                else:
                    descriptions.append(comment_str)

        # Build step order map
        step_order_map = {}
        step_index = 0
        for step in g.objects(uri, BP.pathwayOrder):
            step_process = g.value(step, BP.stepProcess)
            if step_process:
                step_order_map[str(step_process)] = step_index
                step_index += 1

        # Extract pathway components
        components = []
        for component in g.objects(uri, BP.pathwayComponent):
            component_type = g.value(component, RDF.type)
            if component_type:
                type_str = str(component_type).split('#')[-1]
                component_name = g.value(component, BP.displayName)
                component_xrefs = _extract_xrefs(g, component, BP)

                comp_dict = {
                    'type': type_str,
                    'display_name': str(component_name) if component_name else '',
                    'reactome_stable_id': ';'.join(component_xrefs.get('reactome_stable_id', [])),
                    'step_order': step_order_map.get(str(component), None),
                }
                components.append(comp_dict)

        record = {
            'uri': str(uri),
            'display_name': names.get('display_name', ''),
            'synonyms': ';'.join(names.get('synonyms', [])),
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
            'pubmed': ';'.join(xrefs.get('pubmed', [])),
            'go': ';'.join(xrefs.get('go', [])),
            'organism_name': organism.get('organism_name', ''),
            'ncbi_taxonomy_id': organism.get('ncbi_taxonomy_id', ''),
            'description': ' '.join(descriptions),
            'comments': ';'.join(comments),
            'components': components,
        }

        yield record
        count += 1


def reactome_pathways(
    species: str = 'Homo_sapiens',
    max_records: int | None = None,
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
    """
    Download and parse Reactome pathway data as Entity records.

    Args:
        species: Species name (e.g., 'Homo_sapiens')
        max_records: Maximum number of records to parse
        force_refresh: If True, force redownload of the data

    Yields:
        Entity records representing pathways with their components
    """
    from pypath.internals.silver_schema import Membership

    # Ensure all caches are populated (single RDF parse for all data types)
    if not _ensure_all_caches_populated(species, force_refresh):
        return

    cached_data = _DATA_CACHE['pathways']

    for i, record in enumerate(cached_data):
        if max_records is not None and i >= max_records:
            break
        # Build identifiers
        identifiers = []
        if record.get('display_name'):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=record['display_name']))
        if record.get('synonyms'):
            for syn in record['synonyms'].split(';'):
                if syn:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))
        if record.get('reactome_stable_id'):
            for rid in record['reactome_stable_id'].split(';'):
                if rid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))
        if record.get('reactome_id'):
            for rid in record['reactome_id'].split(';'):
                if rid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_ID, value=rid))
        if record.get('go'):
            for goid in record['go'].split(';'):
                if goid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

        # Build annotations
        annotations = []
        if record.get('pubmed'):
            for pmid in record['pubmed'].split(';'):
                if pmid:
                    annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED, value=pmid))
        if record.get('organism_name'):
            annotations.append(Annotation(term='organism', value=record['organism_name']))
        if record.get('ncbi_taxonomy_id'):
            annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=record['ncbi_taxonomy_id']))
        if record.get('description'):
            annotations.append(Annotation(term='description', value=record['description']))

        # Build membership from components
        membership = []
        for comp in record.get('components', []):
            # Determine component entity type
            type_str = comp.get('type', '')
            if 'Pathway' in type_str:
                comp_type = EntityTypeCv.PATHWAY
            elif 'BiochemicalReaction' in type_str:
                comp_type = EntityTypeCv.REACTION
            elif 'Catalysis' in type_str:
                comp_type = EntityTypeCv.CATALYSIS
            elif 'Control' in type_str:
                comp_type = EntityTypeCv.CONTROL
            else:
                comp_type = EntityTypeCv.INTERACTION

            comp_identifiers = []
            if comp.get('display_name'):
                comp_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=comp['display_name']))
            if comp.get('reactome_stable_id'):
                for rid in comp['reactome_stable_id'].split(';'):
                    if rid:
                        comp_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

            comp_annotations = [Annotation(term='role', value='pathway_component')]
            if comp.get('step_order') is not None:
                comp_annotations.append(Annotation(term='step_order', value=str(comp['step_order'])))

            membership.append(Membership(
                member=Entity(type=comp_type, identifiers=comp_identifiers if comp_identifiers else None),
                is_parent=False,
                annotations=comp_annotations,
            ))

        yield Entity(
            type=EntityTypeCv.PATHWAY,
            identifiers=identifiers if identifiers else None,
            annotations=annotations if annotations else None,
            membership=membership if membership else None,
        )


def _iterate_controls(
    g: Graph,
    max_records: int | None = None,
) -> Generator[dict, None, None]:
    """
    Iterate through BioPAX Control and Catalysis elements and yield dictionaries.

    Args:
        g: RDF Graph containing BioPAX data
        max_records: Maximum number of records to parse

    Yields:
        Dictionary representations of controls/catalysis
    """
    control_types = [
        ('Catalysis', EntityTypeCv.CATALYSIS),
        ('Control', EntityTypeCv.CONTROL),
    ]

    count = 0

    for control_type, entity_type_cv in control_types:
        if max_records is not None and count >= max_records:
            break

        query = f"""
        PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

        SELECT ?control
        WHERE {{
            ?control rdf:type bp:{control_type} .
        }}
        """

        for row in g.query(query):
            if max_records is not None and count >= max_records:
                break

            uri = row.control
            names = _extract_names(g, uri, BP)
            xrefs = _extract_xrefs(g, uri, BP)

            # Get control type annotation
            control_type_val = g.value(uri, BP.controlType)

            # Get controller
            controller = g.value(uri, BP.controller)
            controller_info = {}
            if controller:
                controller_names = _extract_names(g, controller, BP)
                controller_xrefs = _extract_xrefs(g, controller, BP)
                controller_type_uri = g.value(controller, RDF.type)
                controller_type_str = str(controller_type_uri).split('#')[-1].lower() if controller_type_uri else 'physicalentity'

                controller_info = {
                    'display_name': controller_names.get('display_name', ''),
                    'entity_type': PHYSICAL_ENTITY_TYPE_MAP.get(controller_type_str, EntityTypeCv.PHYSICAL_ENTITY).value,
                    'reactome_stable_id': ';'.join(controller_xrefs.get('reactome_stable_id', [])),
                    'uniprot': ';'.join(controller_xrefs.get('uniprot', [])),
                }

            # Get controlled process
            controlled = g.value(uri, BP.controlled)
            controlled_info = {}
            if controlled:
                controlled_names = _extract_names(g, controlled, BP)
                controlled_xrefs = _extract_xrefs(g, controlled, BP)
                controlled_type_uri = g.value(controlled, RDF.type)
                controlled_type_str = str(controlled_type_uri).split('#')[-1] if controlled_type_uri else ''

                controlled_info = {
                    'display_name': controlled_names.get('display_name', ''),
                    'type': controlled_type_str,
                    'reactome_stable_id': ';'.join(controlled_xrefs.get('reactome_stable_id', [])),
                }

            record = {
                'uri': str(uri),
                'control_class': control_type,
                'entity_type': entity_type_cv.value,
                'display_name': names.get('display_name', ''),
                'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
                'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
                'go': ';'.join(xrefs.get('go', [])),
                'control_type': str(control_type_val) if control_type_val else '',
                'controller': controller_info,
                'controlled': controlled_info,
            }

            yield record
            count += 1


def reactome_controls(
    species: str = 'Homo_sapiens',
    max_records: int | None = None,
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
    """
    Download and parse Reactome control/catalysis data as Entity records.

    This includes Catalysis and Control entities that represent regulatory
    relationships.

    Args:
        species: Species name (e.g., 'Homo_sapiens')
        max_records: Maximum number of records to parse
        force_refresh: If True, force redownload of the data

    Yields:
        Entity records representing controls/catalysis with their participants
    """
    from pypath.internals.silver_schema import Membership

    # Ensure all caches are populated (single RDF parse for all data types)
    if not _ensure_all_caches_populated(species, force_refresh):
        return

    cached_data = _DATA_CACHE['controls']

    for i, record in enumerate(cached_data):
        if max_records is not None and i >= max_records:
            break
        # Build identifiers
        identifiers = []
        if record.get('display_name'):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=record['display_name']))
        elif record.get('controller', {}).get('display_name'):
            # Use controller name if no display name
            identifiers.append(Identifier(
                type=IdentifierNamespaceCv.NAME,
                value=f"{record['control_class']} by {record['controller']['display_name']}",
            ))

        if record.get('reactome_stable_id'):
            for rid in record['reactome_stable_id'].split(';'):
                if rid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))
        if record.get('reactome_id'):
            for rid in record['reactome_id'].split(';'):
                if rid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_ID, value=rid))
        if record.get('go'):
            for goid in record['go'].split(';'):
                if goid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

        # Build annotations
        annotations = []
        if record.get('control_type'):
            annotations.append(Annotation(term='control_type', value=record['control_type']))

        # Build membership
        membership = []

        # Add controller
        controller = record.get('controller', {})
        if controller:
            controller_identifiers = []
            if controller.get('display_name'):
                controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=controller['display_name']))
            if controller.get('reactome_stable_id'):
                for rid in controller['reactome_stable_id'].split(';'):
                    if rid:
                        controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))
            if controller.get('uniprot'):
                for uid in controller['uniprot'].split(';'):
                    if uid:
                        controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.UNIPROT, value=uid))

            type_mapping = {v.value: v for v in EntityTypeCv}
            controller_type = type_mapping.get(controller.get('entity_type', ''), EntityTypeCv.PHYSICAL_ENTITY)

            if controller_identifiers:
                membership.append(Membership(
                    member=Entity(type=controller_type, identifiers=controller_identifiers),
                    is_parent=False,
                    annotations=[Annotation(term='role', value='controller')],
                ))

        # Add controlled
        controlled = record.get('controlled', {})
        if controlled:
            controlled_identifiers = []
            if controlled.get('display_name'):
                controlled_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=controlled['display_name']))
            if controlled.get('reactome_stable_id'):
                for rid in controlled['reactome_stable_id'].split(';'):
                    if rid:
                        controlled_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

            # Determine controlled type
            controlled_type_str = controlled.get('type', '')
            if 'BiochemicalReaction' in controlled_type_str:
                controlled_type = EntityTypeCv.REACTION
            elif 'Pathway' in controlled_type_str:
                controlled_type = EntityTypeCv.PATHWAY
            else:
                controlled_type = EntityTypeCv.INTERACTION

            if controlled_identifiers:
                membership.append(Membership(
                    member=Entity(type=controlled_type, identifiers=controlled_identifiers),
                    is_parent=False,
                    annotations=[Annotation(term='role', value='controlled')],
                ))

        # Get entity type
        type_mapping = {v.value: v for v in EntityTypeCv}
        entity_type = type_mapping.get(record.get('entity_type', ''), EntityTypeCv.CONTROL)

        yield Entity(
            type=entity_type,
            identifiers=identifiers if identifiers else None,
            annotations=annotations if annotations else None,
            membership=membership if membership else None,
        )
