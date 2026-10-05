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
import json
import re
from collections import defaultdict
from collections.abc import Generator
from pathlib import Path

from rdflib import Graph, Namespace, URIRef
from rdflib.namespace import RDF

from pypath.share.downloads import DATA_DIR
from omnipath_core.molecular_forms import molecular_form_from_identifiers, normalize_molecular_form
from pypath.inputs_v2._molecular_forms import combine_forms, sequence_form


# BioPAX namespace
BP = Namespace("http://www.biopax.org/release/biopax-level3.owl#")

# Module-level cache for parsed data (keyed by data_type)
_DATA_CACHE: dict[str, list[dict]] = {}

# Cache version to invalidate older pickled formats
_CACHE_VERSION = 13  # Includes occurrence molecular forms before feature flattening.

# Delimiter used for list-of-participants and list-of-components fields
_LIST_DELIMITER = "||"
_MISSING_VALUE = "__MISSING__"

# Mapping of BioPAX PhysicalEntity types to str terms
PHYSICAL_ENTITY_TYPE_MAP = {
    'smallmolecule': 'chemical',
    'protein': 'protein',
    'gene': 'gene',
    'complex': 'complex',
    'complexassembly': 'complex',
    'dna': 'dna',
    'dnaregion': 'dna',
    'rna': 'rna',
    'rnaregion': 'rna',
    'physicalentity': 'physical_entity',
}

# Mapping of BioPAX EntityReference subtypes to str terms
ENTITY_REFERENCE_TYPE_MAP = {
    'proteinreference': 'protein',
    'smallmoleculereference': 'chemical',
    'dnareference': 'dna',
    'rnareference': 'rna',
    'rnaregionreference': 'rna',
    'dnaregionreference': 'dna',
}

# EntityReference types to process
ENTITY_REFERENCE_TYPES = {
    BP.ProteinReference,
    BP.SmallMoleculeReference,
    BP.DnaReference,
    BP.RnaReference,
    BP.RnaRegionReference,
    BP.DnaRegionReference,
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


def _prepared_cache_available(
    *,
    data_type: str,
    force_refresh: bool = False,
    **_kwargs: object,
) -> bool:
    return not force_refresh and _load_cached_data(data_type) is not None


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
        elif db in {'refseq', 'refseq protein', 'refseq rna'}:
            xrefs.setdefault('refseq', []).append(id_str)
        elif db in {'ensembl', 'ensembl protein', 'ensembl transcript'}:
            xrefs.setdefault('ensembl', []).append(id_str)
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
    'source_physical_entity', 'compartment', 'modification',
    'molecular_form', 'feature_context', 'refseq', 'ensembl',
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
    'molecular_form', 'feature_context', 'refseq', 'ensembl', 'source_physical_entity',
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
            if field in {'molecular_form', 'feature_context'} and value:
                # Preserve form structure through the source's tabular transport.
                # Escape the participant delimiter even in source descriptions.
                value = json.dumps(value, sort_keys=True).replace('|', '\\u007c')
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


def _flatten_parent_pathways(parent_pathways: list[dict], prefix: str = 'parent_pathway') -> dict[str, str]:
    fields = ['display_name', 'reactome_stable_id', 'uri']
    data: dict[str, str] = {}
    for field in fields:
        items = []
        for parent_pathway in parent_pathways:
            value = parent_pathway.get(field, '')
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
            if field in {'molecular_form', 'feature_context'} and value:
                value = json.dumps(value, sort_keys=True).replace('|', '\\u007c')
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

def _build_degradation_index(
    g: Graph,
    xref_cache: dict,
    entity_reference_index: dict,
) -> dict[str, dict]:
    """Build an index mapping Degradation reaction URIs to their reactant."""
    index: dict[str, dict] = {}

    for s, o in g.subject_objects(RDF.type):
        if o != BP.Degradation:
            continue

        props = _get_entity_props(g, s)

        stoich_map = {}
        for stoich_node in props.get(BP.participantStoichiometry, []):
            s_props = _get_entity_props(g, stoich_node)
            pe = s_props.get(BP.physicalEntity)
            coeff = s_props.get(BP.stoichiometricCoefficient)
            if pe and coeff:
                stoich_map[str(pe[0])] = str(coeff[0])

        for mol in [*props.get(BP.left, []), *props.get(BP.template, [])]:
            participant = _extract_participant_data(
                g,
                mol,
                'template' if mol in props.get(BP.template, []) else 'reactant',
                entity_reference_index,
                xref_cache,
                stoich_map,
            )
            if isinstance(participant, list):
                participant = participant[0] if participant else {}
            else:
                participant.pop('members', None)
                participant.pop('is_family', None)

            if participant:
                index[str(s)] = participant
                break  # Only take the first reactant

    return index


def _load_entity_reference_index(g: Graph, xref_cache: dict[str, dict]) -> dict[str, dict]:
    """Build an index of EntityReference URIs to their full Entity data."""
    from pypath.internals.silver_schema import Entity, Identifier, Annotation

    reference_index = {}

    for s, o in g.subject_objects(RDF.type):
        if o not in ENTITY_REFERENCE_TYPES:
            continue

        props = _get_entity_props(g, s)

        type_str = str(o).split('#')[-1].lower()
        entity_type = ENTITY_REFERENCE_TYPE_MAP.get(type_str, 'physical_entity')

        names = _extract_names_from_props(props, BP)
        xrefs = _extract_xrefs_from_props(props, xref_cache, BP)

        ncbi_tax_id = ''
        org_uris = props.get(BP.organism, [])
        if org_uris:
            ncbi_tax_id = _get_organism_tax_id(g, org_uris[0], xref_cache, BP)

        identifiers = []
        if names.get('display_name'):
            identifiers.append(Identifier(type='name', value=names['display_name']))
        if names.get('standard_name'):
            identifiers.append(Identifier(type='synonym', value=names['standard_name']))
        for synonym in names.get('synonyms', []):
            identifiers.append(Identifier(type='synonym', value=synonym))

        for reactome_id in xrefs.get('reactome_stable_id', []):
            identifiers.append(Identifier(type='reactome', value=reactome_id))
        for reactome_id in xrefs.get('reactome_id', []):
            identifiers.append(Identifier(type='reactome_id', value=reactome_id))
        for uniprot_id in xrefs.get('uniprot', []):
            identifiers.append(Identifier(type='uniprot', value=uniprot_id))
        for chebi_id in xrefs.get('chebi', []):
            identifiers.append(Identifier(type='chebi', value=chebi_id))
        for pubchem_id in xrefs.get('pubchem_compound', []):
            identifiers.append(Identifier(type='pubchem', value=pubchem_id))
        for kegg_id in xrefs.get('kegg', []):
            identifiers.append(Identifier(type='kegg', value=kegg_id))
        for go_id in xrefs.get('go', []):
            identifiers.append(Identifier(type='go', value=go_id))

        for namespace in ('refseq', 'ensembl'):
            for accession in xrefs.get(namespace, []):
                identifiers.append(Identifier(type=namespace, value=accession))

        annotations = []
        for pubmed_id in xrefs.get('pubmed', []):
            annotations.append(Annotation(term='publications', value=pubmed_id))
        if ncbi_tax_id:
            annotations.append(Annotation(term='in_taxon', value=ncbi_tax_id))

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
            'sequence': str(next(g.objects(s, BP.sequence), '')),
        }

    return reference_index


# --------------------------------------------------------------------------- #
# Data Extraction
# --------------------------------------------------------------------------- #

def _participant_molecular_form(g, molecule_uri, participant):
    """Read explicit BioPAX features without assigning unknown sequence positions."""
    source_type = participant.get('entity_type')
    if source_type not in {'protein', 'rna', 'dna'}:
        return None
    form = molecular_form_from_identifiers([
        {'ns': 'uniprot', 'id': accession}
        for accession in str(participant.get('uniprot') or '').split(';') if accession
    ] + [
        {'ns': namespace, 'id': accession}
        for namespace in ('refseq', 'ensembl')
        for accession in str(participant.get(namespace) or '').split(';') if accession
    ]) or {}
    form = combine_forms(form, sequence_form(participant.get('sequence'),
        system='protein' if source_type == 'protein' else 'transcript' if source_type == 'rna' else 'genomic')) or {}
    sequences = form.get('sequence_identifiers') or []
    explicit_sequence = sequence_form(participant.get('sequence'), system='protein' if source_type == 'protein' else 'transcript' if source_type == 'rna' else 'genomic')
    coordinate = {
        'identifier': explicit_sequence['sequence_identifiers'][0] if explicit_sequence else sequences[0] if len(sequences) == 1 else None,
        'coordinate_system': 'protein' if source_type == 'protein' else 'transcript' if source_type == 'rna' else 'genomic',
        'position_base': 1,
    }
    modifications = []

    def exact_position(location):
        if location is None:
            return None
        statuses = {str(status) for status in g.objects(location, BP.positionStatus)}
        # A missing status is not an assertion of equality.
        if statuses != {'EQUAL'}:
            return None
        positions = list(g.objects(location, BP.sequencePosition))
        if len(positions) != 1 or not str(positions[0]).isdigit():
            return None
        return int(positions[0]) or None

    for feature in g.objects(molecule_uri, BP.feature):
        vocabularies = list(g.objects(feature, BP.modificationType))
        if not vocabularies:
            continue
        labels = [str(term) for vocabulary in vocabularies for term in g.objects(vocabulary, BP['term'])]
        mod_accessions = []
        for vocabulary in vocabularies:
            for xref in g.objects(vocabulary, BP.xref):
                databases = {str(db).lower() for db in g.objects(xref, BP.db)}
                if not databases.intersection({'mod', 'psi-mod'}):
                    continue
                for identifier in g.objects(xref, BP.id):
                    match = re.fullmatch(r'(?:MOD:)?(\d+)', str(identifier))
                    if match and f'MOD:{match[1]}' not in mod_accessions:
                        mod_accessions.append(f'MOD:{match[1]}')
        locations = list(g.objects(feature, BP.featureLocation)) or [None]
        for location in locations:
            starts = list(g.objects(location, BP.sequenceIntervalBegin)) if location else []
            ends = list(g.objects(location, BP.sequenceIntervalEnd)) if location else []
            start = exact_position(starts[0]) if len(starts) == 1 else exact_position(location)
            end = exact_position(ends[0]) if len(ends) == 1 else start
            source_locations = starts + ends if starts or ends else [location]
            source_ranges = [
                {'position': [str(p) for p in g.objects(source_location, BP.sequencePosition)],
                 'status': [str(s) for s in g.objects(source_location, BP.positionStatus)]}
                for source_location in source_locations if source_location is not None
            ]
            modifications.append({
                'term': mod_accessions[0] if len(mod_accessions) == 1 else '; '.join(labels) or None,
                'position': start, 'end_position': end,
                'coordinate_reference': coordinate,
                'description': json.dumps({'source_feature': str(feature),
                                           'source_terms': labels,
                                           'mod_accessions': mod_accessions,
                                           'source_ranges': source_ranges}, sort_keys=True),
            })
    form['modifications'] = modifications
    return normalize_molecular_form(form, allow_resolved=False)


def _attach_molecular_context(g, uri, participant, reference_index, xref_cache):
    """Keep one physical entity's form and uninterpreted source feature context."""
    props = _get_entity_props(g, uri)
    xrefs = _extract_xrefs_from_props(props, xref_cache, BP)
    references = [
        reference_index.get(str(ref), {})
        for ref in props.get(BP.entityReference, [])
    ]
    for namespace in ('refseq', 'ensembl'):
        values = list(xrefs.get(namespace, []))
        for reference in references:
            for identifier in (
                getattr(reference.get('entity'), 'identifiers', None) or []
            ):
                if (
                    identifier.type == namespace
                    and identifier.value not in values
                ):
                    values.append(identifier.value)
        participant[namespace] = ';'.join(values)
    sequences = list(
        dict.fromkeys(
            reference.get('sequence')
            for reference in references
            if reference.get('sequence')
        )
    )
    participant['sequence'] = sequences[0] if len(sequences) == 1 else ''
    participant['source_physical_entity'] = str(uri)
    context = []
    for predicate, present in ((BP.feature, True), (BP.notFeature, False)):
        for feature in props.get(predicate, []):
            feature_props = _get_entity_props(g, feature)
            context.append(
                {
                    'source_feature': str(feature),
                    'present': present,
                    'properties': {
                        str(key): [str(value) for value in values]
                        for key, values in feature_props.items()
                    },
                    'locations': [
                        {
                            'properties': {
                                str(key): [str(value) for value in values]
                                for key, values in _get_entity_props(
                                    g, location
                                ).items()
                            },
                            'boundaries': [
                                {
                                    str(key): [str(value) for value in values]
                                    for key, values in _get_entity_props(
                                        g, boundary
                                    ).items()
                                }
                                for predicate in (
                                    BP.sequenceIntervalBegin,
                                    BP.sequenceIntervalEnd,
                                )
                                for boundary in g.objects(location, predicate)
                            ],
                        }
                        for location in feature_props.get(
                            BP.featureLocation, []
                        )
                    ],
                    'modification_types': [
                        {
                            str(key): [str(value) for value in values]
                            for key, values in _get_entity_props(
                                g, vocabulary
                            ).items()
                        }
                        for vocabulary in feature_props.get(
                            BP.modificationType, []
                        )
                    ],
                }
            )
    participant['feature_context'] = context or None
    participant['molecular_form'] = _participant_molecular_form(
        g, uri, participant
    )
    return participant


def _extract_participant_data(g, molecule_uri, role, entity_reference_index, xref_cache, stoich_map):
    """Preserve physical-state context independently of the reference identity."""
    result = _extract_participant_core(g, molecule_uri, role, entity_reference_index, xref_cache, stoich_map)
    props = _get_entity_props(g, molecule_uri)
    compartments = []
    for location in props.get(BP.cellularLocation, []):
        terms = list(g.objects(location, BP['term']))
        compartments.extend(str(term) for term in terms)
    features = []
    for feature in props.get(BP.feature, []):
        labels = [str(term) for vocabulary in g.objects(feature, BP.modificationType)
                  for term in g.objects(vocabulary, BP['term'])]
        positions = [str(position) for location in g.objects(feature, BP.featureLocation)
                     for position in g.objects(location, BP.sequencePosition)]
        features.append('; '.join([*labels, *[f'position {p}' for p in positions]]) or str(feature))
    for participant in result if isinstance(result, list) else [result]:
        participant['source_physical_entity'] = str(molecule_uri)
        participant['compartment'] = '; '.join(compartments)
        participant['modification'] = '; '.join(features)
        _attach_molecular_context(g, molecule_uri, participant, entity_reference_index, xref_cache)
    return result


def _extract_participant_core(
    g: Graph,
    molecule_uri: URIRef,
    role: str,
    entity_reference_index: dict,
    xref_cache: dict,
    stoich_map: dict[str, str]
) -> dict | list[dict]:
    """Extract data for a single participant molecule."""

    props = _get_entity_props(g, molecule_uri)
    names = _extract_names_from_props(props, BP)
    xrefs = _extract_xrefs_from_props(props, xref_cache, BP)

    type_uris = props.get(RDF.type, [])
    entity_type_str = str(type_uris[0]).split('#')[-1].lower() if type_uris else 'physicalentity'
    entity_type = PHYSICAL_ENTITY_TYPE_MAP.get(entity_type_str, 'physical_entity')

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
                if identifier.type == 'synonym':
                    all_synonyms.append(identifier.value)
                elif identifier.type == 'uniprot' and identifier.value not in existing_uniprot:
                    existing_uniprot.add(identifier.value)
                elif identifier.type == 'chebi' and identifier.value not in existing_chebi:
                    existing_chebi.add(identifier.value)
                elif identifier.type == 'pubchem':
                    all_pubchem.append(identifier.value)
                elif identifier.type == 'kegg':
                    all_kegg.append(identifier.value)
                elif identifier.type == 'go':
                    all_go.append(identifier.value)

            for annotation in ref_entity.annotations or []:
                if annotation.term == 'in_taxon' and annotation.value:
                    participant['ncbi_tax_id'] = annotation.value

            participant['uniprot'] = ';'.join(sorted(existing_uniprot - {''}))
            participant['chebi'] = ';'.join(sorted(existing_chebi - {''}))
            participant['synonyms'] = ';'.join(all_synonyms)
            participant['pubchem_compound'] = ';'.join(all_pubchem)
            participant['kegg'] = ';'.join(all_kegg)
            participant['go'] = ';'.join(all_go)

        return participant

    elif members:
        # memberPhysicalEntity denotes alternatives, not a physical assembly
        # or an evolutionary protein family. The children keep their own types.
        family_type = 'complex' if entity_type == 'complex' else 'physical_entity'

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
                        if identifier.type == 'synonym':
                            all_synonyms.append(identifier.value)
                        elif identifier.type == 'uniprot' and identifier.value not in existing_uniprot:
                            existing_uniprot.add(identifier.value)
                        elif identifier.type == 'chebi' and identifier.value not in existing_chebi:
                            existing_chebi.add(identifier.value)
                        elif identifier.type == 'pubchem':
                            all_pubchem.append(identifier.value)
                        elif identifier.type == 'kegg':
                            all_kegg.append(identifier.value)
                        elif identifier.type == 'go':
                            all_go.append(identifier.value)

                    for annotation in ref_entity.annotations or []:
                        if annotation.term == 'in_taxon' and annotation.value:
                            member_data['ncbi_tax_id'] = annotation.value

                    member_data['uniprot'] = ';'.join(sorted(existing_uniprot - {''}))
                    member_data['chebi'] = ';'.join(sorted(existing_chebi - {''}))
                    member_data['synonyms'] = ';'.join(all_synonyms)
                    member_data['pubchem_compound'] = ';'.join(all_pubchem)
                    member_data['kegg'] = ';'.join(all_kegg)
                    member_data['go'] = ';'.join(all_go)

            _attach_molecular_context(g, member_uri, member_data, entity_reference_index, xref_cache)
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

_NUCLEIC_ACID_TYPES = {str('dna'), str('rna')}
_PROTEIN_OR_RNA_TYPES = {str('protein'), str('rna')}


def _is_transcription_or_translation(
    reactants: list[dict],
    products: list[dict],
) -> bool:
    """Skip reactions that represent DNA/RNA → RNA/Protein information flow."""
    if not reactants or not products:
        return False

    reactant_types = {r.get('entity_type', '') for r in reactants if r.get('entity_type')}
    product_types = {p.get('entity_type', '') for p in products if p.get('entity_type')}

    if not reactant_types or not product_types:
        return False

    # All reactants must be DNA or RNA
    if not reactant_types.issubset(_NUCLEIC_ACID_TYPES):
        return False

    # All products must be RNA or Protein
    if not product_types.issubset(_PROTEIN_OR_RNA_TYPES):
        return False

    return True


def _iterate_reactions(
    g: Graph,
    xref_cache: dict,
    entity_reference_index: dict,
    pathway_index: dict[str, list[dict[str, str]]],
    max_records: int | None = None,
) -> Generator[dict, None, None]:

    reaction_targets = {
        BP.BiochemicalReaction: 'reaction',
        BP.Degradation: 'degradation',
        BP.Transport: 'reaction',
        BP.TransportWithBiochemicalReaction: 'reaction',
        BP.ComplexAssembly: 'reaction',
        BP.TemplateReaction: 'reaction',
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
        for mol in [*props.get(BP.left, []), *props.get(BP.template, [])]:
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
        for mol in [*props.get(BP.right, []), *props.get(BP.product, [])]:
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

        # Keep explicitly reported transcription/translation participants too.

        pathway_term_accession = _pathway_term_accessions(pathway_index, reaction_uri)

        participants = reactants + products
        for participant in participants:
            participant['pathway_term_accession'] = pathway_term_accession
        participant_data = _flatten_participants(participants, prefix='participant')

        yield {
            'uri': str(reaction_uri),
            'reaction_type': reaction_type_str,
            'entity_type': entity_type_cv,
            'display_name': names.get('display_name', ''),
            'synonyms': ';'.join(names.get('synonyms', [])),
            'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
            'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
            'pubmed': ';'.join(xrefs.get('pubmed', [])),
            'ec_number': ec_number,
            'direction': direction,
            'template_direction': str(next(iter(props.get(BP.templateDirection, [])), '')),
            'pathway_term_accession': pathway_term_accession,
            **participant_data,
        }
        count += 1


def _pathway_info(
    g: Graph,
    pathway_uri: URIRef,
    xref_cache: dict,
) -> dict[str, str]:
    props = _get_entity_props(g, pathway_uri)
    names = _extract_names_from_props(props, BP)
    xrefs = _extract_xrefs_from_props(props, xref_cache, BP)
    return {
        'display_name': names.get('display_name', ''),
        'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
        'uri': str(pathway_uri),
    }


def _is_pathway_uri(g: Graph, uri: URIRef) -> bool:
    props = _get_entity_props(g, uri)
    types = {str(type_).split('#')[-1].lower() for type_ in props.get(RDF.type, [])}
    return any('pathway' in type_ for type_ in types)


def _build_pathway_parent_index(
    g: Graph,
    xref_cache: dict,
) -> dict[str, list[dict[str, str]]]:
    parent_index: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    seen: defaultdict[str, set[str]] = defaultdict(set)

    for parent_uri in g.subjects(RDF.type, BP.Pathway):
        parent_info = _pathway_info(g, parent_uri, xref_cache)
        parent_key = str(parent_uri)
        props = _get_entity_props(g, parent_uri)

        for child_uri in props.get(BP.pathwayComponent, []):
            if not _is_pathway_uri(g, child_uri):
                continue
            child_key = str(child_uri)
            if parent_key in seen[child_key]:
                continue
            seen[child_key].add(parent_key)
            parent_index[child_key].append(parent_info)

        for step_uri in props.get(BP.pathwayOrder, []):
            step_props = _get_entity_props(g, step_uri)
            for process_uri in step_props.get(BP.stepProcess, []):
                if not _is_pathway_uri(g, process_uri):
                    continue
                child_key = str(process_uri)
                if parent_key in seen[child_key]:
                    continue
                seen[child_key].add(parent_key)
                parent_index[child_key].append(parent_info)

    return parent_index


def _iterate_pathways(
    g: Graph,
    xref_cache: dict,
    max_records: int | None = None,
) -> Generator[dict, None, None]:

    parent_index = _build_pathway_parent_index(g, xref_cache)
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
            **_flatten_parent_pathways(parent_index.get(str(s), []), prefix='parent_pathway'),
        }
        count += 1


def _classify_group_controller_entity_type(
    controller_type: str,
    member_types: list[str],
) -> str:
    """Classify a Reactome controller assembled from memberPhysicalEntity.

    We only keep COMPLEX when the BioPAX controller itself is complex-like.
    All other grouped controllers are modeled conservatively as PHYSICAL_ENTITY,
    because Reactome memberPhysicalEntity sets often mean alternatives / sets /
    grouped active forms, not a curated protein family in the biological sense.
    """
    if controller_type == 'complex':
        return 'complex'

    return 'physical_entity'


def _causal_statement_for_control(control_type_val: str, is_degradation: bool = False) -> str | None:
    """Preserve source control type; the resource module chooses biological semantics."""
    return control_type_val or None


def _controller_set_info(control, controllers):
    return {
        'entity_type': 'logical_control_set',
        'display_name': 'Controllers of ' + str(control),
        'source_physical_entity': str(control) + '#controllers',
        'control_set': {
            'source_control': str(control),
            'logic': 'AND',
            'controllers': [str(uri) for uri in controllers],
        },
    }


def _controller_set_members(g, controllers, reference_index, xref_cache):
    members = []
    for controller in controllers:
        values = _extract_participant_data(
            g, controller, 'member', reference_index, xref_cache, {}
        )
        members.extend(values if isinstance(values, list) else [values])
    return members


def _iterate_controls(
    g: Graph,
    xref_cache: dict,
    entity_reference_index: dict,
    pathway_index: dict[str, list[dict[str, str]]],
    max_records: int | None = None,
    degradation_index: dict[str, dict] | None = None,
) -> Generator[dict, None, None]:


    control_targets = {
        BP.Catalysis: 'catalysis',
        BP.Control: 'control',
    }

    degradation_index = degradation_index or {}

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
        if len(controllers) == 1:
            c_uri = controllers[0]
            c_props = _get_entity_props(g, c_uri)
            c_names = _extract_names_from_props(c_props, BP)
            c_xrefs = _extract_xrefs_from_props(c_props, xref_cache, BP)

            c_types = c_props.get(RDF.type, [])
            c_type_str = str(c_types[0]).split('#')[-1].lower() if c_types else 'physicalentity'
            controller_entity_type = PHYSICAL_ENTITY_TYPE_MAP.get(c_type_str, 'physical_entity')

            controller_info = {
                'role': 'controller',
                'display_name': c_names.get('display_name', ''),
                'entity_type': controller_entity_type,
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
            controller_member_types: list[str] = []
            has_member_physical_entities = False

            c_refs = c_props.get(BP.entityReference, [])
            if c_refs:
                controller_refs_to_merge.append(str(c_refs[0]))

            if not controller_refs_to_merge:
                c_members = list(dict.fromkeys([*c_props.get(BP.memberPhysicalEntity, []), *c_props.get(BP.component, [])]))
                has_member_physical_entities = bool(c_members)
                for member_uri in c_members:
                    member_props = _get_entity_props(g, member_uri)
                    member_names = _extract_names_from_props(member_props, BP)
                    member_xrefs = _extract_xrefs_from_props(member_props, xref_cache, BP)

                    member_types = member_props.get(RDF.type, [])
                    member_type_str = str(member_types[0]).split('#')[-1].lower() if member_types else 'physicalentity'
                    member_entity_type = PHYSICAL_ENTITY_TYPE_MAP.get(member_type_str, 'physical_entity')
                    controller_member_types.append(member_entity_type)

                    member_data = {
                        'entity_type': member_entity_type,
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
                                if identifier.type == 'synonym':
                                    all_synonyms.append(identifier.value)
                                elif identifier.type == 'uniprot':
                                    existing_uniprot.add(identifier.value)
                                elif identifier.type == 'chebi':
                                    existing_chebi.add(identifier.value)
                                elif identifier.type == 'pubchem':
                                    all_pubchem.append(identifier.value)
                                elif identifier.type == 'kegg':
                                    all_kegg.append(identifier.value)
                                elif identifier.type == 'go':
                                    all_go.append(identifier.value)

                            for annotation in ref_entity.annotations or []:
                                if annotation.term == 'in_taxon' and annotation.value:
                                    member_data['ncbi_tax_id'] = annotation.value

                            member_data['uniprot'] = ';'.join(sorted(existing_uniprot - {''}))
                            member_data['chebi'] = ';'.join(sorted(existing_chebi - {''}))
                            member_data['synonyms'] = ';'.join(all_synonyms)
                            member_data['pubchem_compound'] = ';'.join(all_pubchem)
                            member_data['kegg'] = ';'.join(all_kegg)
                            member_data['go'] = ';'.join(all_go)

                    _attach_molecular_context(g, member_uri, member_data, entity_reference_index, xref_cache)
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
                            if identifier.type == 'synonym':
                                all_synonyms.append(identifier.value)
                            elif identifier.type == 'uniprot':
                                existing_uniprot.add(identifier.value)
                            elif identifier.type == 'chebi':
                                existing_chebi.add(identifier.value)
                            elif identifier.type == 'pubchem':
                                all_pubchem.append(identifier.value)
                            elif identifier.type == 'kegg':
                                all_kegg.append(identifier.value)
                            elif identifier.type == 'go':
                                all_go.append(identifier.value)

                        for annotation in ref_entity.annotations or []:
                            if annotation.term == 'in_taxon' and annotation.value:
                                controller_info['ncbi_tax_id'] = annotation.value

                if has_member_physical_entities:
                    controller_entity_type = _classify_group_controller_entity_type(
                        controller_entity_type,
                        controller_member_types,
                    )
                    controller_info['entity_type'] = controller_entity_type

                    if controller_entity_type in {'protein_family', 'physical_entity', 'complex'}:
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

        if len(controllers) == 1:
            _attach_molecular_context(g, c_uri, controller_info, entity_reference_index, xref_cache)
        elif len(controllers) > 1:
            controller_info = _controller_set_info(s, controllers)
            controller_members = _controller_set_members(g, controllers, entity_reference_index, xref_cache)

        for cd_uri in props.get(BP.controlled, []) or [None]:
            controlled_info = {}
            is_degradation_controlled = False
            if cd_uri is not None:
                cd_uri_str = str(cd_uri)
                cd_props = _get_entity_props(g, cd_uri)
                cd_names = _extract_names_from_props(cd_props, BP)
                cd_xrefs = _extract_xrefs_from_props(cd_props, xref_cache, BP)
                cd_types = cd_props.get(RDF.type, [])
                cd_type_str = str(cd_types[0]).split('#')[-1] if cd_types else ''
                cd_type_str_lower = cd_type_str.lower()

                if cd_type_str_lower == 'degradation':
                    # Preserve the controlled process, not an inferred effect on its substrate.
                    is_degradation_controlled = True
                    controlled_info = {
                        'role': 'controlled',
                        'display_name': cd_names.get('display_name', ''),
                        'entity_type': 'degradation',
                        'reactome_stable_id': ';'.join(cd_xrefs.get('reactome_stable_id', [])),
                        'pathway_term_accession': pathway_term_accession,
                    }
                elif cd_type_str_lower in {'biochemicalreaction', 'transport', 'transportwithbiochemicalreaction', 'complexassembly', 'templatereaction'}:
                    controlled_entity_type = 'reaction'
                    controlled_info = {
                        'role': 'controlled',
                        'display_name': cd_names.get('display_name', ''),
                        'entity_type': controlled_entity_type,
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
                elif cd_type_str_lower == 'pathway':
                    controlled_entity_type = 'pathway'
                    controlled_info = {
                        'role': 'controlled',
                        'display_name': cd_names.get('display_name', ''),
                        'entity_type': controlled_entity_type,
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
                else:
                    controlled_entity_type = 'interaction'
                    controlled_info = {
                        'role': 'controlled',
                        'display_name': cd_names.get('display_name', ''),
                        'entity_type': controlled_entity_type,
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
                'controlled_source_uri': str(cd_uri) if cd_uri is not None else '',
                'uri': str(s),
                'control_class': control_type_cls,
                'entity_type': entity_type_cv,
                'display_name': control_display_name,
                'reactome_stable_id': ';'.join(xrefs.get('reactome_stable_id', [])),
                'reactome_id': ';'.join(xrefs.get('reactome_id', [])),
                'go': ';'.join(xrefs.get('go', [])),
                'control_type': control_type_val,
                'causal_statement': _causal_statement_for_control(
                    control_type_val,
                    is_degradation=is_degradation_controlled,
                ) or '',
                'pathway_term_accession': pathway_term_accession,
                **{f'controller_{field}': (json.dumps(controller_info.get(field), sort_keys=True).replace('|', '\\u007c')
                                          if field in {'molecular_form', 'feature_context', 'control_set'} and controller_info.get(field)
                                          else controller_info.get(field, ''))
                   for field in ('molecular_form', 'feature_context', 'refseq', 'ensembl', 'source_physical_entity', 'control_set')},
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
            str(record.get('controller_source_physical_entity', '')),
            str(record.get('controller_member_molecular_form', '')),
        )
        if key in seen:
            continue
        seen.add(key)

        if max_records is not None and count >= max_records:
            break

        yield {
            **{key: value for key, value in record.items()
               if key.startswith('controller_')},
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

def _iterate_physical_groups(g, xref_cache, reference_index, max_records=None):
    """Preserve member/component forms even for groups outside control records."""
    for control in dict.fromkeys(g.subjects(BP.controller, None)):
        controllers = list(g.objects(control, BP.controller))
        if len(controllers) > 1:
            info = _controller_set_info(control, controllers)
            yield {
                **{
                    'controller_' + key: value
                    for key, value in info.items()
                    if key != 'control_set'
                },
                'controller_control_set': json.dumps(
                    info['control_set'], sort_keys=True
                ),
                **_flatten_controller_members(
                    _controller_set_members(
                        g, controllers, reference_index, xref_cache
                    )
                ),
            }
    count = 0
    uris = dict.fromkeys(
        [
            *g.subjects(BP.memberPhysicalEntity, None),
            *g.subjects(BP.component, None),
        ]
    )
    for uri in uris:
        if max_records is not None and count >= max_records:
            return
        props = _get_entity_props(g, uri)
        xrefs = _extract_xrefs_from_props(props, xref_cache, BP)
        names = _extract_names_from_props(props, BP)
        members = []
        for member in dict.fromkeys(
            [
                *props.get(BP.memberPhysicalEntity, []),
                *props.get(BP.component, []),
            ]
        ):
            values = _extract_participant_data(
                g, member, 'member', reference_index, xref_cache, {}
            )
            members.extend(values if isinstance(values, list) else [values])
        root = _extract_participant_data(
            g, uri, 'member', reference_index, xref_cache, {}
        )
        yield {
            'controller_feature_context': json.dumps(
                root.get('feature_context'), sort_keys=True
            )
            if root.get('feature_context')
            else '',
            'controller_entity_type': 'complex'
            if BP.Complex in props.get(RDF.type, [])
            else 'physical_entity',
            'controller_reactome_stable_id': ';'.join(
                xrefs.get('reactome_stable_id', [])
            ),
            'controller_display_name': names.get('display_name', ''),
            'controller_source_physical_entity': str(uri),
            **_flatten_controller_members(members),
        }
        count += 1


def _ensure_all_caches_populated(
    opener,
    species: str = 'Homo_sapiens',
    force_refresh: bool = False,
) -> bool:
    """Ensure all data types are cached."""
    data_types = ['reactions', 'pathways', 'controls', 'control_groups', 'physical_groups']

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
    degradation_index = _build_degradation_index(g, xref_cache, ref_index)

    if _load_cached_data('reactions', force_refresh) is None:
        data = list(_iterate_reactions(g, xref_cache, ref_index, pathway_index, max_records=None))
        _save_cached_data('reactions', data)

    if _load_cached_data('pathways', force_refresh) is None:
        data = list(_iterate_pathways(g, xref_cache, max_records=None))
        _save_cached_data('pathways', data)

    if _load_cached_data('controls', force_refresh) is None:
        data = list(_iterate_controls(
            g, xref_cache, ref_index, pathway_index,
            max_records=None,
            degradation_index=degradation_index,
        ))
        _save_cached_data('controls', data)

    if _load_cached_data('control_groups', force_refresh) is None:
        data = list(_iterate_control_groups(g, xref_cache, ref_index, pathway_index, max_records=None))
        _save_cached_data('control_groups', data)

    if _load_cached_data('physical_groups', force_refresh) is None:
        data = list(_iterate_physical_groups(g, xref_cache, ref_index))
        _save_cached_data('physical_groups', data)

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


_raw.prepared_cache_available = _prepared_cache_available
