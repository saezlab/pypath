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
from pypath.share.downloads import DATA_DIR


# BioPAX namespace
BP = Namespace("http://www.biopax.org/release/biopax-level3.owl#")

# Module-level cache for parsed data (keyed by data_type)
_DATA_CACHE: dict[str, list[dict]] = {}

# Cache version to invalidate older pickled formats
_CACHE_VERSION = 5

# Delimiter used for list-of-participants and list-of-components fields
_LIST_DELIMITER = "||"
_MISSING_VALUE = "__MISSING__"

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
    return cache_dir / f'{data_type}_v{_CACHE_VERSION}.pkl'


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


def _load_biopax_graph(opener, species: str = 'Homo_sapiens') -> Graph | None:
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


_PARTICIPANT_FIELDS = [
    'role',
    'entity_type',
    'display_name',
    'synonyms',
    'reactome_stable_id',
    'uniprot',
    'chebi',
    'pubchem_compound',
    'kegg',
    'go',
    'ncbi_tax_id',
    'stoichiometry',
    'pathway_term_accession',
]

_CONTROLLER_MEMBER_FIELDS = [
    'entity_type',
    'display_name',
    'synonyms',
    'reactome_stable_id',
    'uniprot',
    'chebi',
    'pubchem_compound',
    'kegg',
    'go',
    'ncbi_tax_id',
    'pathway_term_accession',
]


def _join_list(values: list[str]) -> str:
    return _LIST_DELIMITER.join(values)


def _flatten_participants(participants: list[dict], prefix: str = 'participant') -> dict[str, str]:
    data: dict[str, str] = {}
    for field in _PARTICIPANT_FIELDS:
        items = []
        for participant in participants:
            value = participant.get(field, '')
            if value in (None, ''):
                value = _MISSING_VALUE
            items.append(str(value))
        data[f'{prefix}_{field}'] = _join_list(items)
    return data


def _flatten_child_pathways(child_pathways: list[dict], prefix: str = 'child_pathway') -> dict[str, str]:
    fields = ['display_name', 'reactome_stable_id', 'uri', 'step_order']
    data: dict[str, str] = {}
    for field in fields:
        items = []
        for child_pathway in child_pathways:
            value = child_pathway.get(field, '')
            if value in (None, ''):
                value = _MISSING_VALUE
            items.append(str(value))
        data[f'{prefix}_{field}'] = _join_list(items)
    return data


def _flatten_controller_members(members: list[dict], prefix: str = 'controller_member') -> dict[str, str]:
    data: dict[str, str] = {}
    for field in _CONTROLLER_MEMBER_FIELDS:
        items = []
        for member in members:
            value = member.get(field, '')
            if value in (None, ''):
                value = _MISSING_VALUE
            items.append(str(value))
        data[f'{prefix}_{field}'] = _join_list(items)
    return data


def _join_unique_values(values: list[str]) -> str:
    seen: list[str] = []
    for value in values:
        if value and value not in seen:
            seen.append(value)
    return ';'.join(seen)


def _split_field(value: str, delimiter: str = _LIST_DELIMITER) -> list[str]:
    if not value:
        return []
    return [item for item in value.split(delimiter) if item and item != _MISSING_VALUE]



def _build_pathway_terms(pathway_records: list[dict]) -> list[dict]:
    terms: dict[str, dict] = {}
    parent_map: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)

    for record in pathway_records:
        stable_id = record.get('reactome_stable_id', '').split(';')[0]
        if not stable_id:
            continue

        term = terms.setdefault(
            stable_id,
            {
                'id': stable_id,
                'name': record.get('display_name', ''),
                'definition': record.get('definition', ''),
                'synonyms': [],
                'comments': [],
                'xrefs': [],
                'parent_ids': [],
                'parent_names': [],
            },
        )
        term['synonyms'].extend([s for s in record.get('synonyms', '').split(';') if s])
        term['comments'].extend([c for c in record.get('comments', '').split(';') if c])

        reactome_id = record.get('reactome_id', '')
        if reactome_id:
            term['xrefs'].append(f'Reactome:{reactome_id}')
        term['xrefs'].extend([g for g in record.get('go', '').split(';') if g])
        taxon_id = record.get('ncbi_tax_id', '')
        if taxon_id:
            term['xrefs'].append(f'NCBITaxon:{taxon_id}')

        child_ids = _split_field(record.get('child_pathway_reactome_stable_id', ''))
        child_names = _split_field(record.get('child_pathway_display_name', ''))
        for idx, child_id_field in enumerate(child_ids):
            child_id = child_id_field.split(';')[0]
            if not child_id:
                continue
            child_name = child_names[idx] if idx < len(child_names) else ''
            parent_map[child_id].append((stable_id, record.get('display_name', '')))
            terms.setdefault(
                child_id,
                {
                    'id': child_id,
                    'name': child_name,
                    'definition': '',
                    'synonyms': [],
                    'comments': [],
                    'xrefs': [],
                    'parent_ids': [],
                    'parent_names': [],
                },
            )

    for term_id, parents in parent_map.items():
        seen: set[tuple[str, str]] = set()
        for parent_id, parent_name in parents:
            key = (parent_id, parent_name)
            if key in seen:
                continue
            seen.add(key)
            terms[term_id]['parent_ids'].append(parent_id)
            terms[term_id]['parent_names'].append(parent_name)

    for term in terms.values():
        term['synonyms'] = sorted(set(term['synonyms']))
        term['comments'] = list(dict.fromkeys(term['comments']))
        term['xrefs'] = list(dict.fromkeys(term['xrefs']))
        term['synonyms'] = ';'.join(term['synonyms'])
        term['comments'] = ';'.join(term['comments'])
        term['xrefs'] = ';'.join(term['xrefs'])
        term['parent_ids'] = ';'.join(term['parent_ids'])
        term['parent_names'] = ';'.join(term['parent_names'])

    return sorted(terms.values(), key=lambda t: t['id'])



def _build_pathway_membership_index(g: Graph, xref_cache: dict[str, dict]) -> dict[str, list[dict[str, str]]]:
    entity_to_pathways: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    pathway_info_by_uri: dict[str, dict[str, str]] = {}
    pathway_parent_map: defaultdict[str, set[str]] = defaultdict(set)

    for pathway_uri in g.subjects(RDF.type, BP.Pathway):
        props = _get_entity_props(g, pathway_uri)
        names = _extract_names_from_props(props, BP)
        xrefs = _extract_xrefs_from_props(props, xref_cache, BP)
        pathway_uri_str = str(pathway_uri)
        pathway_info = {
            'uri': pathway_uri_str,
            'display_name': names.get('display_name', ''),
            'reactome_stable_id': _join_unique_values(xrefs.get('reactome_stable_id', [])),
        }
        pathway_info_by_uri[pathway_uri_str] = pathway_info

        for component_uri in props.get(BP.pathwayComponent, []):
            entity_to_pathways[str(component_uri)].append(pathway_info)
            component_props = _get_entity_props(g, component_uri)
            component_types = {str(t).split('#')[-1].lower() for t in component_props.get(RDF.type, [])}
            if 'pathway' in component_types:
                pathway_parent_map[str(component_uri)].add(pathway_uri_str)

        for step_uri in props.get(BP.pathwayOrder, []):
            step_props = _get_entity_props(g, step_uri)
            for process_uri in step_props.get(BP.stepProcess, []):
                entity_to_pathways[str(process_uri)].append(pathway_info)
                process_props = _get_entity_props(g, process_uri)
                process_types = {str(t).split('#')[-1].lower() for t in process_props.get(RDF.type, [])}
                if 'pathway' in process_types:
                    pathway_parent_map[str(process_uri)].add(pathway_uri_str)

    def ancestors(pathway_uri_str: str, seen: set[str] | None = None) -> list[dict[str, str]]:
        seen = seen or set()
        result: list[dict[str, str]] = []
        for parent_uri in pathway_parent_map.get(pathway_uri_str, set()):
            if parent_uri in seen:
                continue
            seen.add(parent_uri)
            parent_info = pathway_info_by_uri.get(parent_uri)
            if parent_info:
                result.append(parent_info)
            result.extend(ancestors(parent_uri, seen))
        return result

    deduped: dict[str, list[dict[str, str]]] = {}
    for entity_uri, pathways in entity_to_pathways.items():
        expanded = list(pathways)
        for pathway in list(pathways):
            pathway_uri_str = pathway.get('uri', '')
            if pathway_uri_str:
                expanded.extend(ancestors(pathway_uri_str))

        seen: set[tuple[str, str]] = set()
        deduped[entity_uri] = []
        for pathway in expanded:
            key = (pathway.get('reactome_stable_id', ''), pathway.get('display_name', ''))
            if key in seen:
                continue
            seen.add(key)
            deduped[entity_uri].append(pathway)

    return deduped


def _pathway_term_accessions(pathway_index: dict[str, list[dict[str, str]]], entity_uri: URIRef | str) -> str:
    return _join_unique_values([
        pathway.get('reactome_stable_id', '')
        for pathway in pathway_index.get(str(entity_uri), [])
    ])


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
    pathway_index: dict[str, list[dict[str, str]]],
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
            participant = _extract_participant_data(
                g,
                mol,
                'reactant',
                entity_reference_index,
                xref_cache,
                stoich_map,
            )
            if isinstance(participant, list):
                reactants.extend(participant)
            else:
                participant.pop('members', None)
                participant.pop('is_family', None)
                reactants.append(participant)

        products = []
        for mol in props.get(BP.right, []):
            participant = _extract_participant_data(
                g,
                mol,
                'product',
                entity_reference_index,
                xref_cache,
                stoich_map,
            )
            if isinstance(participant, list):
                products.extend(participant)
            else:
                participant.pop('members', None)
                participant.pop('is_family', None)
                products.append(participant)

        templates = []
        if reaction_type_str == 'TemplateReaction':
            for templ in props.get(BP.template, []):
                t_props = _get_entity_props(g, templ)
                t_names = _extract_names_from_props(t_props, BP)
                t_xrefs = _extract_xrefs_from_props(t_props, xref_cache, BP)
                templates.append({
                    'role': 'template',
                    'entity_type': '',
                    'display_name': t_names.get('display_name', ''),
                    'synonyms': '',
                    'reactome_stable_id': ';'.join(t_xrefs.get('reactome_stable_id', [])),
                    'uniprot': '',
                    'chebi': '',
                    'pubchem_compound': '',
                    'kegg': '',
                    'go': '',
                    'ncbi_tax_id': '',
                    'stoichiometry': '',
                })

        pathway_term_accession = _pathway_term_accessions(pathway_index, reaction_uri)

        participants = reactants + products + templates
        for participant in participants:
            participant['pathway_term_accession'] = pathway_term_accession
        participant_data = _flatten_participants(participants, prefix='participant')

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
            'pathway_term_accession': pathway_term_accession,
            **participant_data,
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

        child_pathways = []
        seen_child_pathways: set[str] = set()
        for comp in props.get(BP.pathwayComponent, []):
            c_props = _get_entity_props(g, comp)
            c_types = c_props.get(RDF.type, [])
            type_str = str(c_types[0]).split('#')[-1] if c_types else ''
            if 'pathway' not in type_str.lower():
                continue
            if str(comp) in seen_child_pathways:
                continue
            seen_child_pathways.add(str(comp))

            c_names = _extract_names_from_props(c_props, BP)
            c_xrefs = _extract_xrefs_from_props(c_props, xref_cache, BP)

            child_pathways.append({
                'display_name': c_names.get('display_name', ''),
                'reactome_stable_id': ';'.join(c_xrefs.get('reactome_stable_id', [])),
                'uri': str(comp),
                'step_order': step_order_map.get(str(comp), None),
            })

        for step in props.get(BP.pathwayOrder, []):
            step_props = _get_entity_props(g, step)
            for process in step_props.get(BP.stepProcess, []):
                process_props = _get_entity_props(g, process)
                process_types = process_props.get(RDF.type, [])
                type_str = str(process_types[0]).split('#')[-1] if process_types else ''
                if 'pathway' not in type_str.lower():
                    continue
                if str(process) in seen_child_pathways:
                    continue
                seen_child_pathways.add(str(process))

                process_names = _extract_names_from_props(process_props, BP)
                process_xrefs = _extract_xrefs_from_props(process_props, xref_cache, BP)
                child_pathways.append({
                    'display_name': process_names.get('display_name', ''),
                    'reactome_stable_id': ';'.join(process_xrefs.get('reactome_stable_id', [])),
                    'uri': str(process),
                    'step_order': step_order_map.get(str(process), None),
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
            'definition': ' '.join(descriptions),
            'comments': ';'.join(comments),
            **_flatten_child_pathways(child_pathways, prefix='child_pathway'),
        }
        count += 1


def _classify_group_controller_entity_type(
    controller_type: EntityTypeCv,
    member_types: list[EntityTypeCv],
) -> EntityTypeCv:
    """Classify a Reactome controller assembled from memberPhysicalEntity.

    We only keep COMPLEX when the BioPAX controller itself is complex-like.
    All other grouped controllers are modeled conservatively as PHYSICAL_ENTITY,
    because Reactome memberPhysicalEntity sets often mean alternatives / sets /
    grouped active forms, not a curated protein family in the biological sense.
    """
    member_type_set = set(member_types)

    if controller_type == EntityTypeCv.COMPLEX or EntityTypeCv.COMPLEX in member_type_set:
        return EntityTypeCv.COMPLEX

    return EntityTypeCv.PHYSICAL_ENTITY


def _iterate_controls(
    g: Graph,
    xref_cache: dict,
    entity_reference_index: dict,
    pathway_index: dict[str, list[dict[str, str]]],
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
        pathway_term_accession = _pathway_term_accessions(pathway_index, s)

        controller_info = {}
        controller_members: list[dict[str, str]] = []
        controllers = props.get(BP.controller, [])
        if controllers:
            c_uri = controllers[0]
            c_props = _get_entity_props(g, c_uri)
            c_names = _extract_names_from_props(c_props, BP)
            c_xrefs = _extract_xrefs_from_props(c_props, xref_cache, BP)

            c_types = c_props.get(RDF.type, [])
            c_type_str = str(c_types[0]).split('#')[-1].lower() if c_types else 'physicalentity'
            controller_entity_type = PHYSICAL_ENTITY_TYPE_MAP.get(c_type_str, EntityTypeCv.PHYSICAL_ENTITY)

            controller_info = {
                'role': 'controller',
                'display_name': c_names.get('display_name', ''),
                'entity_type': controller_entity_type.value,
                'reactome_stable_id': ';'.join(c_xrefs.get('reactome_stable_id', [])),
                'uniprot': ';'.join(c_xrefs.get('uniprot', [])),
                'chebi': '',
                'synonyms': '',
                'pubchem_compound': '',
                'kegg': '',
                'go': '',
                'ncbi_tax_id': '',
                'stoichiometry': '',
                'pathway_term_accession': pathway_term_accession,
            }

            controller_refs_to_merge: list[str] = []
            controller_member_types: list[EntityTypeCv] = []
            has_member_physical_entities = False

            c_refs = c_props.get(BP.entityReference, [])
            if c_refs:
                controller_refs_to_merge.append(str(c_refs[0]))

            if not controller_refs_to_merge:
                c_members = c_props.get(BP.memberPhysicalEntity, [])
                has_member_physical_entities = bool(c_members)
                for member_uri in c_members:
                    member_props = _get_entity_props(g, member_uri)
                    member_names = _extract_names_from_props(member_props, BP)
                    member_xrefs = _extract_xrefs_from_props(member_props, xref_cache, BP)

                    member_types = member_props.get(RDF.type, [])
                    member_type_str = str(member_types[0]).split('#')[-1].lower() if member_types else 'physicalentity'
                    member_entity_type = PHYSICAL_ENTITY_TYPE_MAP.get(member_type_str, EntityTypeCv.PHYSICAL_ENTITY)
                    controller_member_types.append(member_entity_type)

                    member_data = {
                        'entity_type': member_entity_type.value,
                        'display_name': member_names.get('display_name', ''),
                        'reactome_stable_id': ';'.join(member_xrefs.get('reactome_stable_id', [])),
                        'uniprot': ';'.join(member_xrefs.get('uniprot', [])),
                        'chebi': ';'.join(member_xrefs.get('chebi', [])),
                        'synonyms': '',
                        'pubchem_compound': '',
                        'kegg': '',
                        'go': '',
                        'ncbi_tax_id': '',
                        'pathway_term_accession': pathway_term_accession,
                    }

                    member_refs = member_props.get(BP.entityReference, [])
                    if member_refs:
                        ref_key = str(member_refs[0])
                        controller_refs_to_merge.append(ref_key)
                        ref_entry = entity_reference_index.get(ref_key)
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
                                elif identifier.type == IdentifierNamespaceCv.UNIPROT:
                                    existing_uniprot.add(identifier.value)
                                elif identifier.type == IdentifierNamespaceCv.CHEBI:
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

                    controller_members.append(member_data)

            if controller_refs_to_merge:
                existing_uniprot = set(controller_info['uniprot'].split(';')) if controller_info['uniprot'] else set()
                existing_chebi = set(controller_info['chebi'].split(';')) if controller_info['chebi'] else set()
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
                            elif identifier.type == IdentifierNamespaceCv.CHEBI:
                                existing_chebi.add(identifier.value)
                            elif identifier.type == IdentifierNamespaceCv.PUBCHEM_COMPOUND:
                                all_pubchem.append(identifier.value)
                            elif identifier.type == IdentifierNamespaceCv.KEGG_COMPOUND:
                                all_kegg.append(identifier.value)
                            elif identifier.type == IdentifierNamespaceCv.CV_TERM_ACCESSION:
                                all_go.append(identifier.value)

                        for annotation in ref_entity.annotations or []:
                            if annotation.term == IdentifierNamespaceCv.NCBI_TAX_ID and annotation.value:
                                controller_info['ncbi_tax_id'] = annotation.value

                if has_member_physical_entities:
                    controller_entity_type = _classify_group_controller_entity_type(
                        controller_entity_type,
                        controller_member_types,
                    )
                    controller_info['entity_type'] = controller_entity_type.value

                    if controller_entity_type in {EntityTypeCv.PROTEIN_FAMILY, EntityTypeCv.PHYSICAL_ENTITY, EntityTypeCv.COMPLEX}:
                        controller_info['uniprot'] = ''
                        controller_info['chebi'] = ''
                    else:
                        controller_info['uniprot'] = ';'.join(sorted(existing_uniprot - {''}))
                        controller_info['chebi'] = ';'.join(sorted(existing_chebi - {''}))
                else:
                    controller_info['uniprot'] = ';'.join(sorted(existing_uniprot - {''}))
                    controller_info['chebi'] = ';'.join(sorted(existing_chebi - {''}))

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
            cd_type_str_lower = cd_type_str.lower()
            if 'biochemicalreaction' in cd_type_str_lower:
                controlled_entity_type = EntityTypeCv.REACTION
            elif 'pathway' in cd_type_str_lower:
                controlled_entity_type = EntityTypeCv.CV_TERM
            else:
                controlled_entity_type = EntityTypeCv.INTERACTION

            controlled_info = {
                'role': 'controlled',
                'display_name': cd_names.get('display_name', ''),
                'entity_type': controlled_entity_type.value,
                'reactome_stable_id': ';'.join(cd_xrefs.get('reactome_stable_id', [])),
                'uniprot': '',
                'chebi': '',
                'synonyms': '',
                'pubchem_compound': '',
                'kegg': '',
                'go': '',
                'ncbi_tax_id': '',
                'stoichiometry': '',
                'pathway_term_accession': pathway_term_accession,
            }

        control_display_name = names.get('display_name', '')
        if not control_display_name and controller_info.get('display_name'):
            control_display_name = f"{control_type_cls} by {controller_info['display_name']}"

        yield {
            'uri': str(s),
            'control_class': control_type_cls,
            'entity_type': entity_type_cv.value,
            'display_name': control_display_name,
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
            'go': ';'.join(xrefs.get('go', [])),
            'control_type': control_type_val,
            'pathway_term_accession': pathway_term_accession,
            'controller_entity_type': controller_info.get('entity_type', ''),
            'controller_display_name': controller_info.get('display_name', ''),
            'controller_synonyms': controller_info.get('synonyms', ''),
            'controller_reactome_stable_id': controller_info.get('reactome_stable_id', ''),
            'controller_uniprot': controller_info.get('uniprot', ''),
            'controller_chebi': controller_info.get('chebi', ''),
            'controller_pubchem_compound': controller_info.get('pubchem_compound', ''),
            'controller_kegg': controller_info.get('kegg', ''),
            'controller_go': controller_info.get('go', ''),
            'controller_ncbi_tax_id': controller_info.get('ncbi_tax_id', ''),
            'controller_pathway_term_accession': controller_info.get('pathway_term_accession', ''),
            'controlled_entity_type': controlled_info.get('entity_type', ''),
            'controlled_display_name': controlled_info.get('display_name', ''),
            'controlled_synonyms': controlled_info.get('synonyms', ''),
            'controlled_reactome_stable_id': controlled_info.get('reactome_stable_id', ''),
            'controlled_uniprot': controlled_info.get('uniprot', ''),
            'controlled_chebi': controlled_info.get('chebi', ''),
            'controlled_pubchem_compound': controlled_info.get('pubchem_compound', ''),
            'controlled_kegg': controlled_info.get('kegg', ''),
            'controlled_go': controlled_info.get('go', ''),
            'controlled_ncbi_tax_id': controlled_info.get('ncbi_tax_id', ''),
            'controlled_pathway_term_accession': controlled_info.get('pathway_term_accession', ''),
            **_flatten_controller_members(controller_members, prefix='controller_member'),
        }
        count += 1


def _iterate_control_groups(
    g: Graph,
    xref_cache: dict,
    entity_reference_index: dict,
    pathway_index: dict[str, list[dict[str, str]]],
    max_records: int | None = None,
) -> Generator[dict, None, None]:
    seen: set[tuple[str, str, str, str]] = set()
    count = 0

    for record in _iterate_controls(g, xref_cache, entity_reference_index, pathway_index, max_records=None):
        member_types = record.get('controller_member_entity_type', '')
        member_names = record.get('controller_member_display_name', '')
        member_stable_ids = record.get('controller_member_reactome_stable_id', '')
        member_uniprots = record.get('controller_member_uniprot', '')
        member_chebis = record.get('controller_member_chebi', '')
        member_pubchem = record.get('controller_member_pubchem_compound', '')
        member_kegg = record.get('controller_member_kegg', '')
        member_go = record.get('controller_member_go', '')
        member_tax = record.get('controller_member_ncbi_tax_id', '')

        has_members = any(
            value and value != _MISSING_VALUE
            for value in [member_types, member_names, member_stable_ids, member_uniprots, member_chebis]
        )
        if not has_members:
            continue

        key = (
            str(record.get('controller_entity_type', '')),
            str(record.get('controller_display_name', '')),
            str(record.get('controller_reactome_stable_id', '')),
            str(record.get('controller_uniprot', '')),
        )
        if key in seen:
            continue
        seen.add(key)

        if max_records is not None and count >= max_records:
            break

        yield {
            'controller_entity_type': record.get('controller_entity_type', ''),
            'controller_display_name': record.get('controller_display_name', ''),
            'controller_synonyms': record.get('controller_synonyms', ''),
            'controller_reactome_stable_id': record.get('controller_reactome_stable_id', ''),
            'controller_uniprot': record.get('controller_uniprot', ''),
            'controller_chebi': record.get('controller_chebi', ''),
            'controller_pubchem_compound': record.get('controller_pubchem_compound', ''),
            'controller_kegg': record.get('controller_kegg', ''),
            'controller_go': record.get('controller_go', ''),
            'controller_ncbi_tax_id': record.get('controller_ncbi_tax_id', ''),
            'pathway_term_accession': record.get('pathway_term_accession', ''),
            'controller_pathway_term_accession': record.get('controller_pathway_term_accession', ''),
            'controller_member_entity_type': member_types,
            'controller_member_display_name': member_names,
            'controller_member_synonyms': record.get('controller_member_synonyms', ''),
            'controller_member_reactome_stable_id': member_stable_ids,
            'controller_member_uniprot': member_uniprots,
            'controller_member_chebi': member_chebis,
            'controller_member_pubchem_compound': member_pubchem,
            'controller_member_kegg': member_kegg,
            'controller_member_go': member_go,
            'controller_member_ncbi_tax_id': member_tax,
            'controller_member_pathway_term_accession': record.get('controller_member_pathway_term_accession', ''),
        }
        count += 1


# --------------------------------------------------------------------------- #
# Main Parser Function
# --------------------------------------------------------------------------- #

def _ensure_all_caches_populated(
    opener,
    species: str = 'Homo_sapiens',
    force_refresh: bool = False,
) -> bool:
    """Ensure all data types are cached."""
    data_types = ['reactions', 'pathways', 'pathway_terms', 'controls', 'control_groups']

    if not force_refresh:
        all_cached = all(_load_cached_data(dt) is not None for dt in data_types)
        if all_cached:
            return True

    g = _load_biopax_graph(opener, species)
    if g is None:
        return False

    xref_cache = _build_xref_cache(g, BP)
    ref_index = _load_entity_reference_index(g, xref_cache)
    pathway_index = _build_pathway_membership_index(g, xref_cache)

    if _load_cached_data('reactions', force_refresh) is None:
        data = list(_iterate_reactions(g, xref_cache, ref_index, pathway_index, max_records=None))
        _save_cached_data('reactions', data)

    if _load_cached_data('pathways', force_refresh) is None:
        data = list(_iterate_pathways(g, xref_cache, max_records=None))
        _save_cached_data('pathways', data)

    if _load_cached_data('pathway_terms', force_refresh) is None:
        pathway_records = _load_cached_data('pathways', force_refresh) or []
        data = _build_pathway_terms(pathway_records)
        _save_cached_data('pathway_terms', data)

    if _load_cached_data('controls', force_refresh) is None:
        data = list(_iterate_controls(g, xref_cache, ref_index, pathway_index, max_records=None))
        _save_cached_data('controls', data)

    if _load_cached_data('control_groups', force_refresh) is None:
        data = list(_iterate_control_groups(g, xref_cache, ref_index, pathway_index, max_records=None))
        _save_cached_data('control_groups', data)

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
    if not _ensure_all_caches_populated(opener, species, force_refresh):
        return
    cached_data = _DATA_CACHE[data_type]
    for i, record in enumerate(cached_data):
        if max_records is not None and i >= max_records:
            break
        yield record
