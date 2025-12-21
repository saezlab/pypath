"""
Reactome BioPAX/RDF graph parser.

Parses Reactome BioPAX (OWL) data files using RDF graph traversal with optimizations:
- Global XRef caching to minimize graph queries
- Batched property fetching
- EntityReference indexing for fast lookup
- Pickle-based caching for parsed data
"""

from __future__ import annotations

import pickle
from collections import defaultdict
from collections.abc import Generator
from pathlib import Path

from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF

from pypath.internals.cv_terms import EntityTypeCv
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
    g.parse(owl_file, format='xml')
    return g


# --------------------------------------------------------------------------- #
# Optimization Helpers (Graph Traversal)
# --------------------------------------------------------------------------- #

def _build_xref_cache(g: Graph, bp_ns: Namespace) -> dict[str, dict[str, str]]:
    """
    Pre-scan all Xref nodes in the graph.
    Returns a dict mapping XRef URI to {'db': ..., 'id': ...}.
    """
    cache = defaultdict(dict)

    for s, o in g.subject_objects(bp_ns.db):
        cache[str(s)]['db'] = str(o)

    for s, o in g.subject_objects(bp_ns.id):
        if str(s) in cache:
            cache[str(s)]['id'] = str(o)

    return cache


def _get_entity_props(g: Graph, uri: URIRef) -> dict[URIRef, list[URIRef | str]]:
    """Fetch all properties of a subject in one pass."""
    props = defaultdict(list)
    for p, o in g.predicate_objects(uri):
        props[p].append(o)
    return props


def _extract_xrefs_from_props(
    props: dict,
    xref_cache: dict[str, dict],
    bp_ns: Namespace
) -> dict[str, list[str]]:
    """Extract cross-references using the property dict and global xref cache."""
    xrefs: dict[str, list[str]] = {}

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


def _get_organism_tax_id(g: Graph, organism_uri: URIRef, xref_cache: dict, bp_ns: Namespace) -> str:
    if not organism_uri:
        return ''

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
    """Build an index of EntityReference URIs to their full Entity data."""
    from pypath.internals.silver_schema import Entity, Identifier, Annotation
    from pypath.internals.cv_terms import IdentifierNamespaceCv

    reference_index = {}

    for s, o in g.subject_objects(RDF.type):
        if o not in ENTITY_REFERENCE_TYPES:
            continue

        props = _get_entity_props(g, s)

        type_str = str(o).split('#')[-1].lower()
        entity_type = ENTITY_REFERENCE_TYPE_MAP.get(type_str, EntityTypeCv.PHYSICAL_ENTITY)

        names = _extract_names_from_props(props, BP)
        xrefs = _extract_xrefs_from_props(props, xref_cache, BP)

        ncbi_tax_id = ''
        org_uris = props.get(BP.organism, [])
        if org_uris:
            ncbi_tax_id = _get_organism_tax_id(g, org_uris[0], xref_cache, BP)

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

        annotations = []
        for pubmed_id in xrefs.get('pubmed', []):
            annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED, value=pubmed_id))
        if ncbi_tax_id:
            annotations.append(Annotation(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=ncbi_tax_id))

        entity = Entity(
            type=entity_type,
            identifiers=identifiers if identifiers else None,
            annotations=annotations if annotations else None,
        )

        reference_index[str(s)] = {
            'type': entity_type,
            'primary_name': names.get('display_name', ''),
            'reactome_identifier': ';'.join(xrefs.get('reactome_stable_id', [])),
            'entity': entity,
        }

    return reference_index


# --------------------------------------------------------------------------- #
# Data Extraction
# --------------------------------------------------------------------------- #

def _extract_participant_data(
    g: Graph,
    molecule_uri: URIRef,
    role: str,
    entity_reference_index: dict,
    xref_cache: dict,
    stoich_map: dict[str, str]
) -> dict | list[dict]:
    """Extract data for a single participant molecule."""
    from pypath.internals.cv_terms import IdentifierNamespaceCv

    props = _get_entity_props(g, molecule_uri)
    names = _extract_names_from_props(props, BP)
    xrefs = _extract_xrefs_from_props(props, xref_cache, BP)

    type_uris = props.get(RDF.type, [])
    entity_type_str = str(type_uris[0]).split('#')[-1].lower() if type_uris else 'physicalentity'
    entity_type = PHYSICAL_ENTITY_TYPE_MAP.get(entity_type_str, EntityTypeCv.PHYSICAL_ENTITY)

    refs = props.get(BP.entityReference, [])
    members = props.get(BP.memberPhysicalEntity, [])

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

    elif members:
        if entity_type == EntityTypeCv.PROTEIN:
            family_type = EntityTypeCv.PROTEIN_FAMILY
        else:
            family_type = entity_type

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
            'members': [],
        }

        member_list = []
        for member_uri in members:
            member_props = _get_entity_props(g, member_uri)
            member_names = _extract_names_from_props(member_props, BP)
            member_xrefs = _extract_xrefs_from_props(member_props, xref_cache, BP)

            member_data = {
                'role': 'member',
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


# --------------------------------------------------------------------------- #
# Data Iterators
# --------------------------------------------------------------------------- #

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

        stoich_map = {}
        for stoich_node in props.get(BP.participantStoichiometry, []):
            s_props = _get_entity_props(g, stoich_node)
            pe = s_props.get(BP.physicalEntity)
            coeff = s_props.get(BP.stoichiometricCoefficient)
            if pe and coeff:
                stoich_map[str(pe[0])] = str(coeff[0])

        reactants = []
        for mol in props.get(BP.left, []):
            reactants.append(_extract_participant_data(g, mol, 'reactant', entity_reference_index, xref_cache, stoich_map))

        products = []
        for mol in props.get(BP.right, []):
            products.append(_extract_participant_data(g, mol, 'product', entity_reference_index, xref_cache, stoich_map))

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

        step_order_map = {}
        step_index = 0
        for step in props.get(BP.pathwayOrder, []):
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

    from pypath.internals.cv_terms import IdentifierNamespaceCv

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

            controller_refs_to_merge = []

            c_refs = c_props.get(BP.entityReference, [])
            if c_refs:
                controller_refs_to_merge.append(str(c_refs[0]))

            if not controller_refs_to_merge:
                c_members = c_props.get(BP.memberPhysicalEntity, [])
                for member_uri in c_members:
                    member_props = _get_entity_props(g, member_uri)
                    member_refs = member_props.get(BP.entityReference, [])
                    if member_refs:
                        controller_refs_to_merge.append(str(member_refs[0]))

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
# Main Parser Function
# --------------------------------------------------------------------------- #

def _ensure_all_caches_populated(species: str = 'Homo_sapiens', force_refresh: bool = False) -> bool:
    """Ensure all data types are cached."""
    data_types = ['reactions', 'pathways', 'controls']

    if not force_refresh:
        all_cached = all(_load_cached_data(dt) is not None for dt in data_types)
        if all_cached:
            return True

    g = _load_biopax_graph(species, force_refresh)
    if g is None:
        return False

    xref_cache = _build_xref_cache(g, BP)
    ref_index = _load_entity_reference_index(g, xref_cache)

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


def _raw(
    opener,
    data_type: str,
    species: str = 'Homo_sapiens',
    max_records: int | None = None,
    force_refresh: bool = False,
    **_kwargs: object,
):
    """
    Parse Reactome BioPAX data and yield records.

    Args:
        opener: File opener (not used, data is loaded via download_and_open internally)
        data_type: One of 'reactions', 'pathways', or 'controls'
        species: Species name (default: 'Homo_sapiens')
        max_records: Optional limit on records to yield
        force_refresh: Force re-parsing and cache refresh

    Yields:
        Dictionary for each record
    """
    if not _ensure_all_caches_populated(species, force_refresh):
        return
    cached_data = _DATA_CACHE[data_type]
    for i, record in enumerate(cached_data):
        if max_records is not None and i >= max_records:
            break
        yield record
