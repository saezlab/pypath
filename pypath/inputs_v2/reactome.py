
"""
Parse Reactome BioPAX data and emit Entity records.

This module downloads Reactome BioPAX (OWL) data and converts it into Entity
records using the declarative schema pattern from tabular_builder.

OPTIMIZED VERSION:
- Replaced SPARQL with native generator iteration.
- Implemented global XRef caching to minimize graph traversal.
- Batched property fetching.
"""

from __future__ import annotations

import pickle
from collections import defaultdict
from collections.abc import Generator
from pathlib import Path

from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF

from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    MembershipRoleCv,
    MoleculeAnnotationsCv,
    ParticipantMetadataCv,
    ResourceAnnotationCv,
    ResourceCv,
    UpdateCategoryCV,
)

from pypath.share.downloads import download_and_open, DATA_DIR


# BioPAX namespace
BP = Namespace("http://www.biopax.org/release/biopax-level3.owl#")

# Module-level cache for parsed data (keyed by data_type)
_DATA_CACHE: dict[str, list[dict]] = {}

# Mapping of BioPAX PhysicalEntity types to EntityTypeCv terms
PHYSICAL_ENTITY_TYPE_MAP = {
    'smallmolecule': EntityTypeCv.SMALL_MOLECULE,
    'protein': EntityTypeCv.PROTEIN,
    'gene': EntityTypeCv.GENE,
    'complex': EntityTypeCv.COMPLEX,
    'complexassembly': EntityTypeCv.COMPLEX,
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
ENTITY_REFERENCE_TYPES = {
    BP.ProteinReference,
    BP.SmallMoleculeReference,
    BP.DnaReference,
    BP.RnaReference,
}


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.
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
                'pathway database.'
            )),
        ],
    )


# --------------------------------------------------------------------------- #
# Caching Infrastructure
# --------------------------------------------------------------------------- #

def _get_cache_path(data_type: str) -> Path:
    cache_dir = DATA_DIR / 'reactome'
    return cache_dir / f'{data_type}.pkl'


def _load_cached_data(data_type: str, force_refresh: bool = False) -> list[dict] | None:
    if force_refresh:
        return None
    if data_type in _DATA_CACHE:
        return _DATA_CACHE[data_type]
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
    _DATA_CACHE[data_type] = data
    pickle_path = _get_cache_path(data_type)
    try:
        pickle_path.parent.mkdir(parents=True, exist_ok=True)
        with open(pickle_path, 'wb') as f:
            pickle.dump(data, f)
    except Exception:
        pass


def _load_biopax_graph(species: str = 'Homo_sapiens', force_refresh: bool = False) -> Graph | None:
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

    owl_file = None
    target_filename = f'{species}.owl'
    for filename, file_handle in opener.result.items():
        if target_filename in filename:
            owl_file = file_handle
            break

    if not owl_file:
        return None

    g = Graph()
    # Using lxml/xml parser via rdflib
    g.parse(owl_file, format='xml')
    return g


# --------------------------------------------------------------------------- #
# Optimization Helpers (Graph Traversal)
# --------------------------------------------------------------------------- #

def _build_xref_cache(g: Graph, bp_ns: Namespace) -> dict[str, dict[str, str]]:
    """
    Pre-scan all Xref nodes (UnificationXref, RelationshipXref, etc.) in the graph.
    Returns a dict mapping the XRef URI string to {'db': ..., 'id': ...}.
    This turns O(N) graph queries into O(1) dict lookups.
    """
    cache = defaultdict(dict)
    
    # Scan all 'db' properties
    for s, o in g.subject_objects(bp_ns.db):
        cache[str(s)]['db'] = str(o)
        
    # Scan all 'id' properties
    for s, o in g.subject_objects(bp_ns.id):
        if str(s) in cache:
            cache[str(s)]['id'] = str(o)
            
    return cache


def _get_entity_props(g: Graph, uri: URIRef) -> dict[URIRef, list[URIRef | str]]:
    """
    Fetch all properties of a subject in one pass to avoid repetitive graph access.
    """
    props = defaultdict(list)
    for p, o in g.predicate_objects(uri):
        props[p].append(o)
    return props


def _extract_xrefs_from_props(
    props: dict, 
    xref_cache: dict[str, dict], 
    bp_ns: Namespace
) -> dict[str, list[str]]:
    """
    Extract cross-references using the property dict and global xref cache.
    """
    xrefs: dict[str, list[str]] = {}
    
    # Props contains the URIs of the Xref nodes
    xref_uris = props.get(bp_ns.xref, [])
    
    for xref_uri in xref_uris:
        info = xref_cache.get(str(xref_uri))
        if not info or 'db' not in info or 'id' not in info:
            continue
            
        db = info['db'].lower()
        id_str = info['id']

        if 'reactome' in db and 'pubmed' not in db:
            key = 'reactome_stable_id' if ('stable' in db or 'R-' in id_str) else 'reactome_id'
            xrefs.setdefault(key, []).append(id_str)
        elif 'uniprot' in db:
            xrefs.setdefault('uniprot', []).append(id_str)
        elif 'chebi' in db:
            xrefs.setdefault('chebi', []).append(id_str)
        elif 'pubchem' in db or 'compound' in db:
            xrefs.setdefault('pubchem_compound', []).append(id_str)
        elif 'kegg' in db:
            xrefs.setdefault('kegg', []).append(id_str)
        elif 'pubmed' in db:
            xrefs.setdefault('pubmed', []).append(id_str)
        elif 'gene ontology' in db or id_str.startswith('GO:'):
            xrefs.setdefault('go', []).append(id_str)
        elif 'taxonomy' in db:
            xrefs.setdefault('ncbi_taxonomy', []).append(id_str)

    return xrefs


def _extract_names_from_props(props: dict, bp_ns: Namespace) -> dict[str, str | list[str]]:
    names = {}
    
    display_name = props.get(bp_ns.displayName)
    if display_name:
        names['display_name'] = str(display_name[0])

    standard_name = props.get(bp_ns.standardName)
    if standard_name:
        names['standard_name'] = str(standard_name[0])

    synonyms = props.get(bp_ns.name, [])
    if synonyms:
        names['synonyms'] = [str(s) for s in synonyms]

    return names


def _extract_organism_from_props(props: dict, xref_cache: dict, bp_ns: Namespace) -> str:
    """Returns NCBI Tax ID or empty string."""
    organisms = props.get(bp_ns.organism, [])
    if not organisms:
        return ''
    
    # Usually one organism per entity
    # The organism node usually has an xref to taxonomy
    # We can't use the props dict of the entity here, we need to look up the organism node
    # However, to avoid passing 'g', we can check if we can resolve it via cache.
    # Actually, in BioPAX, organism is a BioSource. We might need to look it up.
    # For speed, let's assume organism node URI might be in our xref cache if it WAS an xref
    # but it's not. 
    # Optimized compromise: Since there are few unique organisms, we won't optimize 
    # this deeply, or we can rely on the caller to handle it.
    # Let's return the URI string and handle lookup in the main loop if needed, 
    # OR just assume we have to look it up. 
    
    # Given we removed 'g' from arguments, we will skip deep traversal here 
    # and rely on the calling function to have fetched the organism props if needed.
    # OR, we pass 'g' just for this rare lookup. 
    pass 
    return ''

def _get_organism_tax_id(g: Graph, organism_uri: URIRef, xref_cache: dict, bp_ns: Namespace) -> str:
    if not organism_uri:
        return ''
    
    # Fast lookup: get props of organism node
    # We do this "on demand" because there are few organism nodes compared to entities
    org_xrefs = list(g.objects(organism_uri, bp_ns.xref))
    
    for xref in org_xrefs:
        info = xref_cache.get(str(xref))
        if info and 'taxonomy' in info.get('db', '').lower():
            return info['id']
    return ''


# --------------------------------------------------------------------------- #
# Index Builders
# --------------------------------------------------------------------------- #

def _load_entity_reference_index(g: Graph, xref_cache: dict[str, dict]) -> dict[str, dict]:
    """
    Build an index of EntityReference URIs to their full Entity data.
    """
    reference_index = {}

    # Iterate all subjects that define a type
    # Check if that type is one of our target reference types
    for s, o in g.subject_objects(RDF.type):
        if o not in ENTITY_REFERENCE_TYPES:
            continue

        # Fetch all properties once
        props = _get_entity_props(g, s)

        # Determine Type
        type_str = str(o).split('#')[-1].lower()
        entity_type = ENTITY_REFERENCE_TYPE_MAP.get(type_str, EntityTypeCv.PHYSICAL_ENTITY)

        # Extract Data
        names = _extract_names_from_props(props, BP)
        xrefs = _extract_xrefs_from_props(props, xref_cache, BP)
        
        # Organism
        ncbi_tax_id = ''
        org_uris = props.get(BP.organism, [])
        if org_uris:
            ncbi_tax_id = _get_organism_tax_id(g, org_uris[0], xref_cache, BP)

        # Build identifiers
        identifiers = []
        if names.get('display_name'):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=names['display_name']))
        if names.get('standard_name'):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=names['standard_name']))
        for synonym in names.get('synonyms', []):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=synonym))

        for reactome_id in xrefs.get('reactome_stable_id', []):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=reactome_id))
        for reactome_id in xrefs.get('reactome_id', []):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_ID, value=reactome_id))
        for uniprot_id in xrefs.get('uniprot', []):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.UNIPROT, value=uniprot_id))
        for chebi_id in xrefs.get('chebi', []):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.CHEBI, value=chebi_id))
        for pubchem_id in xrefs.get('pubchem_compound', []):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=pubchem_id))
        for kegg_id in xrefs.get('kegg', []):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.KEGG_COMPOUND, value=kegg_id))
        for go_id in xrefs.get('go', []):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=go_id))

        # Build annotations
        annotations = []
        for pubmed_id in xrefs.get('pubmed', []):
            annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED, value=pubmed_id))
        if ncbi_tax_id:
            annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=ncbi_tax_id))

        # Create the full Entity
        entity = Entity(
            type=entity_type,
            identifiers=identifiers if identifiers else None,
            annotations=annotations if annotations else None,
        )

        # Store in index
        reference_index[str(s)] = {
            'type': entity_type,
            'primary_name': names.get('display_name', ''),
            'reactome_identifier': ';'.join(xrefs.get('reactome_stable_id', [])),
            'entity': entity,
        }

    return reference_index


# --------------------------------------------------------------------------- #
# Data Iterators
# --------------------------------------------------------------------------- #

def _extract_participant_data(
    g: Graph,
    molecule_uri: URIRef,
    role: str,
    entity_reference_index: dict,
    xref_cache: dict,
    stoich_map: dict[str, str]
) -> dict | list[dict]:
    """
    Extract data for a single participant molecule.

    Returns either a single dict for simple entities, or a list with:
    - First element: the parent set/family entity
    - Remaining elements: member entities
    """
    props = _get_entity_props(g, molecule_uri)
    names = _extract_names_from_props(props, BP)
    xrefs = _extract_xrefs_from_props(props, xref_cache, BP)

    # Entity Type
    type_uris = props.get(RDF.type, [])
    entity_type_str = str(type_uris[0]).split('#')[-1].lower() if type_uris else 'physicalentity'
    entity_type = PHYSICAL_ENTITY_TYPE_MAP.get(entity_type_str, EntityTypeCv.PHYSICAL_ENTITY)

    # Check if this is a set/family (has memberPhysicalEntity but no entityReference)
    refs = props.get(BP.entityReference, [])
    members = props.get(BP.memberPhysicalEntity, [])

    # If it has a direct entityReference, treat as simple entity
    if refs:
        participant = {
            'role': role,
            'entity_type': entity_type.value if hasattr(entity_type, 'value') else str(entity_type),
            'display_name': names.get('display_name', ''),
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'uniprot': ';'.join(xrefs.get('uniprot', [])),
            'chebi': ';'.join(xrefs.get('chebi', [])),
            'stoichiometry': stoich_map.get(str(molecule_uri), ''),
            'synonyms': '',
            'pubchem_compound': '',
            'kegg': '',
            'go': '',
            'ncbi_tax_id': '',
        }

        # Enrich from entityReference
        ref_entry = entity_reference_index.get(str(refs[0]))
        if ref_entry and ref_entry.get('entity'):
            ref_entity = ref_entry['entity']

            existing_uniprot = set(participant['uniprot'].split(';')) if participant['uniprot'] else set()
            existing_chebi = set(participant['chebi'].split(';')) if participant['chebi'] else set()
            all_synonyms = []
            all_pubchem = []
            all_kegg = []
            all_go = []

            for identifier in ref_entity.identifiers or []:
                if identifier.type == IdentifierNamespaceCv.SYNONYM:
                    all_synonyms.append(identifier.value)
                elif identifier.type == IdentifierNamespaceCv.UNIPROT and identifier.value not in existing_uniprot:
                    existing_uniprot.add(identifier.value)
                elif identifier.type == IdentifierNamespaceCv.CHEBI and identifier.value not in existing_chebi:
                    existing_chebi.add(identifier.value)
                elif identifier.type == IdentifierNamespaceCv.PUBCHEM_COMPOUND:
                    all_pubchem.append(identifier.value)
                elif identifier.type == IdentifierNamespaceCv.KEGG_COMPOUND:
                    all_kegg.append(identifier.value)
                elif identifier.type == IdentifierNamespaceCv.CV_TERM_ACCESSION:
                    all_go.append(identifier.value)

            for annotation in ref_entity.annotations or []:
                if annotation.term == IdentifierNamespaceCv.NCBI_TAX_ID and annotation.value:
                    participant['ncbi_tax_id'] = annotation.value

            participant['uniprot'] = ';'.join(sorted(existing_uniprot - {''}))
            participant['chebi'] = ';'.join(sorted(existing_chebi - {''}))
            participant['synonyms'] = ';'.join(all_synonyms)
            participant['pubchem_compound'] = ';'.join(all_pubchem)
            participant['kegg'] = ';'.join(all_kegg)
            participant['go'] = ';'.join(all_go)

        return participant

    # If it has members but no direct reference, it's a set/family
    elif members:
        # Determine family type based on base entity type
        if entity_type == EntityTypeCv.PROTEIN:
            family_type = EntityTypeCv.PROTEIN_FAMILY
        else:
            # For other types (complexes, small molecules), keep original type
            family_type = entity_type

        # Create parent family entity
        family = {
            'role': role,
            'entity_type': family_type.value if hasattr(family_type, 'value') else str(family_type),
            'display_name': names.get('display_name', ''),
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'uniprot': '',
            'chebi': '',
            'stoichiometry': stoich_map.get(str(molecule_uri), ''),
            'synonyms': '',
            'pubchem_compound': '',
            'kegg': '',
            'go': '',
            'ncbi_tax_id': '',
            'is_family': True,
            'members': [],  # Will be populated with member data
        }

        # Extract member entities
        member_list = []
        for member_uri in members:
            member_props = _get_entity_props(g, member_uri)
            member_names = _extract_names_from_props(member_props, BP)
            member_xrefs = _extract_xrefs_from_props(member_props, xref_cache, BP)

            member_data = {
                'role': 'member',  # Special role for family members
                'entity_type': entity_type.value if hasattr(entity_type, 'value') else str(entity_type),
                'display_name': member_names.get('display_name', ''),
                'reactome_stable_id': ';'.join(member_xrefs.get('reactome_stable_id', [])),
                'uniprot': ';'.join(member_xrefs.get('uniprot', [])),
                'chebi': ';'.join(member_xrefs.get('chebi', [])),
                'stoichiometry': '',
                'synonyms': '',
                'pubchem_compound': '',
                'kegg': '',
                'go': '',
                'ncbi_tax_id': '',
            }

            # Enrich member from its entityReference
            member_refs = member_props.get(BP.entityReference, [])
            if member_refs:
                ref_entry = entity_reference_index.get(str(member_refs[0]))
                if ref_entry and ref_entry.get('entity'):
                    ref_entity = ref_entry['entity']

                    existing_uniprot = set(member_data['uniprot'].split(';')) if member_data['uniprot'] else set()
                    existing_chebi = set(member_data['chebi'].split(';')) if member_data['chebi'] else set()
                    all_synonyms = []
                    all_pubchem = []
                    all_kegg = []
                    all_go = []

                    for identifier in ref_entity.identifiers or []:
                        if identifier.type == IdentifierNamespaceCv.SYNONYM:
                            all_synonyms.append(identifier.value)
                        elif identifier.type == IdentifierNamespaceCv.UNIPROT and identifier.value not in existing_uniprot:
                            existing_uniprot.add(identifier.value)
                        elif identifier.type == IdentifierNamespaceCv.CHEBI and identifier.value not in existing_chebi:
                            existing_chebi.add(identifier.value)
                        elif identifier.type == IdentifierNamespaceCv.PUBCHEM_COMPOUND:
                            all_pubchem.append(identifier.value)
                        elif identifier.type == IdentifierNamespaceCv.KEGG_COMPOUND:
                            all_kegg.append(identifier.value)
                        elif identifier.type == IdentifierNamespaceCv.CV_TERM_ACCESSION:
                            all_go.append(identifier.value)

                    for annotation in ref_entity.annotations or []:
                        if annotation.term == IdentifierNamespaceCv.NCBI_TAX_ID and annotation.value:
                            member_data['ncbi_tax_id'] = annotation.value

                    member_data['uniprot'] = ';'.join(sorted(existing_uniprot - {''}))
                    member_data['chebi'] = ';'.join(sorted(existing_chebi - {''}))
                    member_data['synonyms'] = ';'.join(all_synonyms)
                    member_data['pubchem_compound'] = ';'.join(all_pubchem)
                    member_data['kegg'] = ';'.join(all_kegg)
                    member_data['go'] = ';'.join(all_go)

            member_list.append(member_data)

        family['members'] = member_list
        return family

    # Fallback: no reference and no members (shouldn't happen often)
    else:
        return {
            'role': role,
            'entity_type': entity_type.value if hasattr(entity_type, 'value') else str(entity_type),
            'display_name': names.get('display_name', ''),
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'uniprot': ';'.join(xrefs.get('uniprot', [])),
            'chebi': ';'.join(xrefs.get('chebi', [])),
            'stoichiometry': stoich_map.get(str(molecule_uri), ''),
            'synonyms': '',
            'pubchem_compound': '',
            'kegg': '',
            'go': '',
            'ncbi_tax_id': '',
        }


def _iterate_reactions(
    g: Graph,
    xref_cache: dict,
    entity_reference_index: dict,
    max_records: int | None = None,
) -> Generator[dict, None, None]:
    
    reaction_targets = {
        BP.BiochemicalReaction: EntityTypeCv.REACTION,
        BP.Degradation: EntityTypeCv.DEGRADATION,
        BP.TemplateReaction: EntityTypeCv.REACTION,
    }

    count = 0

    # Iterate through all subjects that have a type in our target list
    for s, o in g.subject_objects(RDF.type):
        if o not in reaction_targets:
            continue
        
        if max_records is not None and count >= max_records:
            break

        reaction_uri = s
        entity_type_cv = reaction_targets[o]
        reaction_type_str = str(o).split('#')[-1]

        props = _get_entity_props(g, reaction_uri)
        names = _extract_names_from_props(props, BP)
        xrefs = _extract_xrefs_from_props(props, xref_cache, BP)

        ec_number = str(props[BP.eCNumber][0]) if props.get(BP.eCNumber) else ''
        direction = str(props[BP.conversionDirection][0]) if props.get(BP.conversionDirection) else ''

        # Build Stoichiometry Map (Requires traversing stoichiometric objects)
        stoich_map = {}
        for stoich_node in props.get(BP.participantStoichiometry, []):
            # Fetch props for the stoichiometry node (small overhead)
            s_props = _get_entity_props(g, stoich_node)
            pe = s_props.get(BP.physicalEntity)
            coeff = s_props.get(BP.stoichiometricCoefficient)
            if pe and coeff:
                stoich_map[str(pe[0])] = str(coeff[0])

        # Extract Reactants and Products
        reactants = []
        for mol in props.get(BP.left, []):
            reactants.append(_extract_participant_data(g, mol, 'reactant', entity_reference_index, xref_cache, stoich_map))

        products = []
        for mol in props.get(BP.right, []):
            products.append(_extract_participant_data(g, mol, 'product', entity_reference_index, xref_cache, stoich_map))

        # Templates
        templates = []
        if reaction_type_str == 'TemplateReaction':
            for templ in props.get(BP.template, []):
                t_props = _get_entity_props(g, templ)
                t_names = _extract_names_from_props(t_props, BP)
                t_xrefs = _extract_xrefs_from_props(t_props, xref_cache, BP)
                templates.append({
                    'role': 'template',
                    'display_name': t_names.get('display_name', ''),
                    'reactome_stable_id': ';'.join(t_xrefs.get('reactome_stable_id', [])),
                })

        yield {
            'uri': str(reaction_uri),
            'reaction_type': reaction_type_str,
            'entity_type': entity_type_cv.value,
            'display_name': names.get('display_name', ''),
            'synonyms': ';'.join(names.get('synonyms', [])),
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
            'pubmed': ';'.join(xrefs.get('pubmed', [])),
            'ec_number': ec_number,
            'direction': direction,
            'reactants': reactants,
            'products': products,
            'templates': templates,
        }
        count += 1


def _iterate_pathways(
    g: Graph,
    xref_cache: dict,
    max_records: int | None = None,
) -> Generator[dict, None, None]:

    count = 0
    for s in g.subjects(RDF.type, BP.Pathway):
        if max_records is not None and count >= max_records:
            break
        
        props = _get_entity_props(g, s)
        names = _extract_names_from_props(props, BP)
        xrefs = _extract_xrefs_from_props(props, xref_cache, BP)
        
        ncbi_tax_id = ''
        orgs = props.get(BP.organism, [])
        if orgs:
            ncbi_tax_id = _get_organism_tax_id(g, orgs[0], xref_cache, BP)

        comments = []
        descriptions = []
        for comment in props.get(BP.comment, []):
            c_str = str(comment)
            if "Reactome DB_ID:" not in c_str:
                if any(c_str.startswith(p) for p in ['Reviewed:', 'Authored:', 'Edited:']):
                    comments.append(c_str)
                else:
                    descriptions.append(c_str)

        # Step Order (Requires parsing PathwayStep)
        step_order_map = {}
        step_index = 0
        for step in props.get(BP.pathwayOrder, []):
            # steps are nodes
            step_props = _get_entity_props(g, step)
            for process in step_props.get(BP.stepProcess, []):
                step_order_map[str(process)] = step_index
            step_index += 1

        components = []
        for comp in props.get(BP.pathwayComponent, []):
            c_props = _get_entity_props(g, comp)
            c_types = c_props.get(RDF.type, [])
            type_str = str(c_types[0]).split('#')[-1] if c_types else ''
            
            c_names = _extract_names_from_props(c_props, BP)
            c_xrefs = _extract_xrefs_from_props(c_props, xref_cache, BP)
            
            components.append({
                'type': type_str,
                'display_name': c_names.get('display_name', ''),
                'reactome_stable_id': ';'.join(c_xrefs.get('reactome_stable_id', [])),
                'step_order': step_order_map.get(str(comp), None),
            })

        yield {
            'uri': str(s),
            'display_name': names.get('display_name', ''),
            'synonyms': ';'.join(names.get('synonyms', [])),
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
            'pubmed': ';'.join(xrefs.get('pubmed', [])),
            'go': ';'.join(xrefs.get('go', [])),
            'ncbi_tax_id': ncbi_tax_id,
            'description': ' '.join(descriptions),
            'comments': ';'.join(comments),
            'components': components,
        }
        count += 1


def _iterate_controls(
    g: Graph,
    xref_cache: dict,
    entity_reference_index: dict,
    max_records: int | None = None,
) -> Generator[dict, None, None]:
    
    control_targets = {
        BP.Catalysis: EntityTypeCv.CATALYSIS,
        BP.Control: EntityTypeCv.CONTROL,
    }
    
    count = 0
    for s, o in g.subject_objects(RDF.type):
        if o not in control_targets:
            continue
        
        if max_records is not None and count >= max_records:
            break

        entity_type_cv = control_targets[o]
        control_type_cls = str(o).split('#')[-1]
        
        props = _get_entity_props(g, s)
        names = _extract_names_from_props(props, BP)
        xrefs = _extract_xrefs_from_props(props, xref_cache, BP)
        
        control_type_val = str(props[BP.controlType][0]) if props.get(BP.controlType) else ''

        # Controller
        controller_info = {}
        controllers = props.get(BP.controller, [])
        if controllers:
            c_uri = controllers[0]
            c_props = _get_entity_props(g, c_uri)
            c_names = _extract_names_from_props(c_props, BP)
            c_xrefs = _extract_xrefs_from_props(c_props, xref_cache, BP)
            
            c_types = c_props.get(RDF.type, [])
            c_type_str = str(c_types[0]).split('#')[-1].lower() if c_types else 'physicalentity'
            
            controller_info = {
                'display_name': c_names.get('display_name', ''),
                'entity_type': PHYSICAL_ENTITY_TYPE_MAP.get(c_type_str, EntityTypeCv.PHYSICAL_ENTITY).value,
                'reactome_stable_id': ';'.join(c_xrefs.get('reactome_stable_id', [])),
                'uniprot': ';'.join(c_xrefs.get('uniprot', [])),
                'synonyms': '',
                'pubchem_compound': '',
                'kegg': '',
                'go': '',
                'ncbi_tax_id': '',
            }
            
            # Enrich - collect entity references to merge
            controller_refs_to_merge = []

            # Try direct entityReference
            c_refs = c_props.get(BP.entityReference, [])
            if c_refs:
                controller_refs_to_merge.append(str(c_refs[0]))

            # If no direct reference, check memberPhysicalEntity (for sets/complexes)
            if not controller_refs_to_merge:
                c_members = c_props.get(BP.memberPhysicalEntity, [])
                for member_uri in c_members:
                    member_props = _get_entity_props(g, member_uri)
                    member_refs = member_props.get(BP.entityReference, [])
                    if member_refs:
                        controller_refs_to_merge.append(str(member_refs[0]))

            # Merge identifiers from all found entity references
            if controller_refs_to_merge:
                existing_uniprot = set(controller_info['uniprot'].split(';')) if controller_info['uniprot'] else set()
                all_synonyms = []
                all_pubchem = []
                all_kegg = []
                all_go = []

                for ref_key in controller_refs_to_merge:
                    ref_entry = entity_reference_index.get(ref_key)
                    if ref_entry and ref_entry.get('entity'):
                        ref_entity = ref_entry['entity']

                        for identifier in ref_entity.identifiers or []:
                            if identifier.type == IdentifierNamespaceCv.SYNONYM:
                                all_synonyms.append(identifier.value)
                            elif identifier.type == IdentifierNamespaceCv.UNIPROT:
                                existing_uniprot.add(identifier.value)
                            elif identifier.type == IdentifierNamespaceCv.PUBCHEM_COMPOUND:
                                all_pubchem.append(identifier.value)
                            elif identifier.type == IdentifierNamespaceCv.KEGG_COMPOUND:
                                all_kegg.append(identifier.value)
                            elif identifier.type == IdentifierNamespaceCv.CV_TERM_ACCESSION:
                                all_go.append(identifier.value)

                        for annotation in ref_entity.annotations or []:
                            if annotation.term == IdentifierNamespaceCv.NCBI_TAX_ID and annotation.value:
                                controller_info['ncbi_tax_id'] = annotation.value

                controller_info['uniprot'] = ';'.join(sorted(existing_uniprot - {''}))
                controller_info['synonyms'] = ';'.join(all_synonyms)
                controller_info['pubchem_compound'] = ';'.join(all_pubchem)
                controller_info['kegg'] = ';'.join(all_kegg)
                controller_info['go'] = ';'.join(all_go)

        # Controlled
        controlled_info = {}
        controlled_list = props.get(BP.controlled, [])
        if controlled_list:
            cd_uri = controlled_list[0]
            cd_props = _get_entity_props(g, cd_uri)
            cd_names = _extract_names_from_props(cd_props, BP)
            cd_xrefs = _extract_xrefs_from_props(cd_props, xref_cache, BP)
            cd_types = cd_props.get(RDF.type, [])
            cd_type_str = str(cd_types[0]).split('#')[-1] if cd_types else ''
            
            controlled_info = {
                'display_name': cd_names.get('display_name', ''),
                'type': cd_type_str,
                'reactome_stable_id': ';'.join(cd_xrefs.get('reactome_stable_id', [])),
            }

        yield {
            'uri': str(s),
            'control_class': control_type_cls,
            'entity_type': entity_type_cv.value,
            'display_name': names.get('display_name', ''),
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
            'go': ';'.join(xrefs.get('go', [])),
            'control_type': control_type_val,
            'controller': controller_info,
            'controlled': controlled_info,
        }
        count += 1


# --------------------------------------------------------------------------- #
# Main Logic
# --------------------------------------------------------------------------- #

def _ensure_all_caches_populated(species: str = 'Homo_sapiens', force_refresh: bool = False) -> bool:
    data_types = ['reactions', 'pathways', 'controls']

    if not force_refresh:
        all_cached = all(_load_cached_data(dt) is not None for dt in data_types)
        if all_cached:
            return True

    g = _load_biopax_graph(species, force_refresh)
    if g is None:
        return False

    # 1. Build Global XRef Cache (The big performance win)
    xref_cache = _build_xref_cache(g, BP)

    # 2. Build EntityReference Index (Optimized)
    ref_index = _load_entity_reference_index(g, xref_cache)

    # 3. Parse Data
    if _load_cached_data('reactions', force_refresh) is None:
        data = list(_iterate_reactions(g, xref_cache, ref_index, max_records=None))
        _save_cached_data('reactions', data)

    if _load_cached_data('pathways', force_refresh) is None:
        data = list(_iterate_pathways(g, xref_cache, max_records=None))
        _save_cached_data('pathways', data)

    if _load_cached_data('controls', force_refresh) is None:
        data = list(_iterate_controls(g, xref_cache, ref_index, max_records=None))
        _save_cached_data('controls', data)

    return True


def _build_member_entities(participants: list[dict]) -> list[Entity]:
    """
    Build Entity objects for reaction participants.
    Handles both simple entities and families with members.
    """
    from pypath.internals.silver_schema import Membership

    role_mapping = {
        'reactant': BiologicalRoleCv.REACTANT,
        'product': BiologicalRoleCv.PRODUCT,
        'template': BiologicalRoleCv.TEMPLATE,
        'controller': BiologicalRoleCv.CONTROLLER,
        'controlled': BiologicalRoleCv.CONTROLLED,
        'pathway_component': BiologicalRoleCv.PATHWAY_COMPONENT,
        'member': MembershipRoleCv.MEMBER_OF,  # For family members
    }

    entities = []
    for p in participants:
        entity_type_str = p.get('entity_type', 'physical_entity')
        type_mapping = {v.value: v for v in EntityTypeCv}
        entity_type = type_mapping.get(entity_type_str, EntityTypeCv.PHYSICAL_ENTITY)

        identifiers = []
        annotations = []

        if p.get('display_name'):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=p['display_name']))

        if p.get('synonyms'):
            for syn in p['synonyms'].split(';'):
                if syn:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))

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

        if p.get('pubchem_compound'):
            for pcid in p['pubchem_compound'].split(';'):
                if pcid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=pcid))

        if p.get('kegg'):
            for kid in p['kegg'].split(';'):
                if kid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.KEGG_COMPOUND, value=kid))

        if p.get('go'):
            for goid in p['go'].split(';'):
                if goid:
                    identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

        role_str = p.get('role', '')
        role_cv = role_mapping.get(role_str)
        if role_cv:
            annotations.append(Annotation(term=role_cv))

        if p.get('stoichiometry'):
            annotations.append(Annotation(term=ParticipantMetadataCv.STOICHIOMETRY, value=p['stoichiometry']))

        if p.get('ncbi_tax_id'):
            annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=p['ncbi_tax_id']))

        # Check if this is a family with members
        if p.get('is_family') and p.get('members'):
            from pypath.internals.silver_schema import Membership

            # Build member entities first
            family_membership = []
            for member in p['members']:
                member_identifiers = []
                member_annotations = []

                if member.get('display_name'):
                    member_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=member['display_name']))

                if member.get('synonyms'):
                    for syn in member['synonyms'].split(';'):
                        if syn:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))

                if member.get('reactome_stable_id'):
                    for rid in member['reactome_stable_id'].split(';'):
                        if rid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

                if member.get('uniprot'):
                    for uid in member['uniprot'].split(';'):
                        if uid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.UNIPROT, value=uid))

                if member.get('chebi'):
                    for cid in member['chebi'].split(';'):
                        if cid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.CHEBI, value=cid))

                if member.get('pubchem_compound'):
                    for pcid in member['pubchem_compound'].split(';'):
                        if pcid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=pcid))

                if member.get('kegg'):
                    for kid in member['kegg'].split(';'):
                        if kid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.KEGG_COMPOUND, value=kid))

                if member.get('go'):
                    for goid in member['go'].split(';'):
                        if goid:
                            member_identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

                # Add member role
                member_annotations.append(Annotation(term=MembershipRoleCv.MEMBER_OF))

                if member.get('ncbi_tax_id'):
                    member_annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=member['ncbi_tax_id']))

                # Get member entity type
                member_type_str = member.get('entity_type', 'physical_entity')
                member_type = type_mapping.get(member_type_str, EntityTypeCv.PHYSICAL_ENTITY)

                # Create member entity
                member_entity = Entity(
                    type=member_type,
                    identifiers=member_identifiers if member_identifiers else None,
                    annotations=member_annotations if member_annotations else None
                )

                family_membership.append(Membership(
                    member=member_entity,
                    is_parent=False,
                    annotations=member_annotations if member_annotations else None
                ))

            # Create family entity with members
            family_entity = Entity(
                type=entity_type,
                identifiers=identifiers if identifiers else None,
                annotations=annotations if annotations else None,
                membership=family_membership if family_membership else None
            )

            entities.append((family_entity, annotations if annotations else None))
        else:
            # Simple entity without members
            entities.append((
                Entity(type=entity_type, identifiers=identifiers if identifiers else None, annotations=annotations if annotations else None),
                annotations if annotations else None,
            ))

    return entities


def reactome_reactions(
    species: str = 'Homo_sapiens',
    max_records: int | None = None,
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
    """
    Download and parse Reactome reaction data as Entity records.
    """
    from pypath.internals.silver_schema import Membership

    if not _ensure_all_caches_populated(species, force_refresh):
        return

    cached_data = _DATA_CACHE['reactions']

    for i, record in enumerate(cached_data):
        if max_records is not None and i >= max_records:
            break
        
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

        annotations = []
        if record.get('pubmed'):
            for pmid in record['pubmed'].split(';'):
                if pmid:
                    annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED, value=pmid))
        if record.get('ec_number'):
            annotations.append(Annotation(term=MoleculeAnnotationsCv.EC_NUMBER, value=record['ec_number']))
        if record.get('direction'):
            annotations.append(Annotation(term=InteractionMetadataCv.CONVERSION_DIRECTION, value=record['direction']))

        membership = []
        all_participants = record.get('reactants', []) + record.get('products', []) + record.get('templates', [])

        for entity, member_annotations in _build_member_entities(all_participants):
            membership.append(Membership(
                member=entity,
                is_parent=False,
                annotations=member_annotations,
            ))

        type_mapping = {v.value: v for v in EntityTypeCv}
        entity_type = type_mapping.get(record.get('entity_type', ''), EntityTypeCv.REACTION)

        yield Entity(
            type=entity_type,
            identifiers=identifiers if identifiers else None,
            annotations=annotations if annotations else None,
            membership=membership if membership else None,
        )


def reactome_pathways(
    species: str = 'Homo_sapiens',
    max_records: int | None = None,
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
    """
    Download and parse Reactome pathway data as Entity records.
    """
    from pypath.internals.silver_schema import Membership

    if not _ensure_all_caches_populated(species, force_refresh):
        return

    cached_data = _DATA_CACHE['pathways']

    for i, record in enumerate(cached_data):
        if max_records is not None and i >= max_records:
            break
        
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

        annotations = []
        if record.get('pubmed'):
            for pmid in record['pubmed'].split(';'):
                if pmid:
                    annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED, value=pmid))
        if record.get('ncbi_tax_id'):
            annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=record['ncbi_tax_id']))
        if record.get('description'):
            annotations.append(Annotation(term=MoleculeAnnotationsCv.DESCRIPTION, value=record['description']))

        membership = []
        for comp in record.get('components', []):
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

            comp_annotations = [Annotation(term=BiologicalRoleCv.PATHWAY_COMPONENT)]
            if comp.get('step_order') is not None:
                comp_annotations.append(Annotation(term=ParticipantMetadataCv.STEP_ORDER, value=str(comp['step_order'])))

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


def reactome_controls(
    species: str = 'Homo_sapiens',
    max_records: int | None = None,
    force_refresh: bool = False,
) -> Generator[Entity, None, None]:
    """
    Download and parse Reactome control/catalysis data as Entity records.
    """
    from pypath.internals.silver_schema import Membership

    if not _ensure_all_caches_populated(species, force_refresh):
        return

    cached_data = _DATA_CACHE['controls']

    for i, record in enumerate(cached_data):
        if max_records is not None and i >= max_records:
            break
        
        identifiers = []
        if record.get('display_name'):
            identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=record['display_name']))
        elif record.get('controller', {}).get('display_name'):
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

        annotations = []
        if record.get('control_type'):
            annotations.append(Annotation(term=InteractionMetadataCv.CONTROL_TYPE, value=record['control_type']))

        membership = []

        controller = record.get('controller', {})
        if controller:
            controller_identifiers = []
            controller_annotations = [Annotation(term=BiologicalRoleCv.CONTROLLER)]

            if controller.get('display_name'):
                controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=controller['display_name']))

            if controller.get('synonyms'):
                for syn in controller['synonyms'].split(';'):
                    if syn:
                        controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=syn))

            if controller.get('reactome_stable_id'):
                for rid in controller['reactome_stable_id'].split(';'):
                    if rid:
                        controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

            if controller.get('uniprot'):
                for uid in controller['uniprot'].split(';'):
                    if uid:
                        controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.UNIPROT, value=uid))

            if controller.get('pubchem_compound'):
                for pcid in controller['pubchem_compound'].split(';'):
                    if pcid:
                        controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=pcid))

            if controller.get('kegg'):
                for kid in controller['kegg'].split(';'):
                    if kid:
                        controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.KEGG_COMPOUND, value=kid))

            if controller.get('go'):
                for goid in controller['go'].split(';'):
                    if goid:
                        controller_identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=goid))

            if controller.get('ncbi_tax_id'):
                controller_annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=controller['ncbi_tax_id']))

            type_mapping = {v.value: v for v in EntityTypeCv}
            controller_type = type_mapping.get(controller.get('entity_type', ''), EntityTypeCv.PHYSICAL_ENTITY)

            if controller_identifiers:
                membership.append(Membership(
                    member=Entity(type=controller_type, identifiers=controller_identifiers, annotations=controller_annotations if controller_annotations else None),
                    is_parent=False,
                    annotations=[Annotation(term=BiologicalRoleCv.CONTROLLER)],
                ))

        controlled = record.get('controlled', {})
        if controlled:
            controlled_identifiers = []
            if controlled.get('display_name'):
                controlled_identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=controlled['display_name']))
            if controlled.get('reactome_stable_id'):
                for rid in controlled['reactome_stable_id'].split(';'):
                    if rid:
                        controlled_identifiers.append(Identifier(type=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=rid))

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
                    annotations=[Annotation(term=BiologicalRoleCv.CONTROLLED)],
                ))

        type_mapping = {v.value: v for v in EntityTypeCv}
        entity_type = type_mapping.get(record.get('entity_type', ''), EntityTypeCv.CONTROL)

        yield Entity(
            type=entity_type,
            identifiers=identifiers if identifiers else None,
            annotations=annotations if annotations else None,
            membership=membership if membership else None,
        )