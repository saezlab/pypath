"""ChEBI OBO parser for molecule entities and ontology terms."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Generator, Iterable
import re
from typing import Any


_DATA_CACHE: dict[str, list[dict[str, Any]]] = {}

_CHEMICAL_ENTITY_ROOT = 'CHEBI:24431'
_EXCLUDED_ROOTS = {
    'CHEBI:50906',  # role
    'CHEBI:33232',  # application
}

_TRUE_VALUES = {'1', 'true', 'yes'}
_STRUCTURAL_PROPERTY_KEYS = {
    'chemrof:smiles_string': 'smiles',
    'chemrof:inchi_string': 'inchi',
    'chemrof:inchi_key_string': 'inchikey',
    'chemrof:generalized_empirical_formula': 'formula',
    'chemrof:mass': 'mass',
    'chemrof:monoisotopic_mass': 'monoisotopic_mass',
    'chemrof:charge': 'charge',
}
_STRUCTURE_SIGNAL_KEYS = {'smiles', 'inchi', 'inchikey', 'formula'}

_PROPERTY_VALUE_RE = re.compile(r'^(\S+)\s+(".*?"|\S+)(?:\s+\S+)?$')
_KEGG_RE = re.compile(r'(?:kegg(?:\s+compound)?|kegg\.compound)[:\s]+(C\d+)', re.IGNORECASE)
_PUBCHEM_RE = re.compile(r'(?:pubchem(?:\s+compound)?|cid)[:\s]+(\d+)', re.IGNORECASE)
_HMDB_RE = re.compile(r'(?:^|\b)(HMDB\d+)(?:\b|$)', re.IGNORECASE)
_LIPIDMAPS_RE = re.compile(r'(?:lipid\s*maps|lipidmaps)[:\s]+([A-Z0-9]+)', re.IGNORECASE)
_CAS_RE = re.compile(r'\b(\d{2,7}-\d{2}-\d)\b')


def _raw(
    opener,
    data_type: str,
    force_refresh: bool = False,
    max_records: int | None = None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield normalized ChEBI records for the requested dataset."""

    records = _load_records(opener, force_refresh=force_refresh)
    rows = records.get(data_type, [])

    if max_records is not None:
        rows = rows[:max_records]

    yield from rows


def _load_records(opener, force_refresh: bool = False) -> dict[str, list[dict[str, Any]]]:

    if _DATA_CACHE and not force_refresh:
        return _DATA_CACHE

    terms = _parse_terms(_first_handle(opener))
    records = _build_records(terms)

    _DATA_CACHE.clear()
    _DATA_CACHE.update(records)

    return records


def _first_handle(opener):
    if not opener or not opener.result:
        return None
    if isinstance(opener.result, dict):
        return next(iter(opener.result.values()), None)
    return opener.result


def _parse_terms(handle) -> dict[str, dict[str, Any]]:
    if not handle:
        return {}

    if hasattr(handle, 'seek'):
        handle.seek(0)

    terms: dict[str, dict[str, Any]] = {}
    current: dict[str, Any] | None = None

    for raw_line in handle:
        if isinstance(raw_line, bytes):
            raw_line = raw_line.decode('utf-8', 'ignore')

        stripped = raw_line.rstrip('\n').strip()

        if not stripped or stripped.startswith('!'):
            continue

        if stripped == '[Term]':
            _store_term(current, terms)
            current = _new_term()
            continue

        if stripped.startswith('['):
            _store_term(current, terms)
            current = None
            continue

        if current is None or ':' not in stripped:
            continue

        tag, value = stripped.split(':', 1)
        _apply_tag(current, tag.strip(), value.strip())

    _store_term(current, terms)
    return terms


def _new_term() -> dict[str, Any]:
    return {
        'id': '',
        'name': '',
        'definition': '',
        'synonyms': [],
        'comments': [],
        'xrefs': [],
        'is_a': [],
        'relationships': [],
        'properties': {},
        'is_obsolete': False,
    }


def _store_term(term: dict[str, Any] | None, terms: dict[str, dict[str, Any]]) -> None:
    if not term:
        return

    term_id = term.get('id')
    if not term_id:
        return

    properties = {
        key: value
        for key, value in term.get('properties', {}).items()
        if value not in (None, '')
    }

    terms[term_id] = {
        **term,
        'synonyms': _dedupe(term['synonyms']),
        'comments': _dedupe(term['comments']),
        'xrefs': _dedupe(term['xrefs']),
        'is_a': _dedupe(term['is_a']),
        'relationships': _dedupe_relationships(term['relationships']),
        'properties': properties,
    }


def _apply_tag(term: dict[str, Any], tag: str, value: str) -> None:
    if tag == 'id':
        term['id'] = value
    elif tag == 'name':
        term['name'] = value
    elif tag == 'def':
        term['definition'] = _extract_quoted(value)
    elif tag == 'synonym':
        synonym = _extract_quoted(value)
        if synonym:
            term['synonyms'].append(synonym)
    elif tag == 'comment':
        if value:
            term['comments'].append(value)
    elif tag == 'xref':
        xref = _strip_inline_comment(value)
        if xref:
            term['xrefs'].append(xref)
    elif tag == 'is_a':
        parent = _strip_inline_comment(value)
        if parent:
            term['is_a'].append(parent)
    elif tag == 'relationship':
        relationship = _parse_relationship(value)
        if relationship:
            term['relationships'].append(relationship)
    elif tag == 'property_value':
        key, property_value = _parse_property_value(value)
        if key and property_value not in (None, ''):
            term['properties'][key] = property_value
    elif tag == 'is_obsolete':
        term['is_obsolete'] = value.strip().lower() in _TRUE_VALUES


def _strip_inline_comment(value: str) -> str:
    return value.split(' ! ', 1)[0].strip()


def _extract_quoted(value: str) -> str:
    match = re.search(r'"((?:\\.|[^"\\])*)"', value)
    if not match:
        return value.strip()
    return bytes(match.group(1), 'utf-8').decode('unicode_escape').strip()


def _parse_relationship(value: str) -> dict[str, str] | None:
    base, _, target_name = value.partition(' ! ')
    parts = base.strip().split(None, 2)
    if len(parts) < 2:
        return None
    return {
        'type': parts[0],
        'target': parts[1],
        'target_name': target_name.strip() or '',
    }


def _parse_property_value(value: str) -> tuple[str | None, str | None]:
    match = _PROPERTY_VALUE_RE.match(value)
    if not match:
        return None, None

    raw_key, raw_value = match.groups()
    key = _STRUCTURAL_PROPERTY_KEYS.get(raw_key)
    if not key:
        return None, None

    if raw_value.startswith('"'):
        raw_value = _extract_quoted(raw_value)

    return key, raw_value.strip()


def _dedupe(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def _dedupe_relationships(values: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    seen: set[tuple[str, str, str]] = set()
    result: list[dict[str, str]] = []
    for value in values:
        key = (value.get('type', ''), value.get('target', ''), value.get('target_name', ''))
        if not key[0] or not key[1] or key in seen:
            continue
        seen.add(key)
        result.append(value)
    return result


def _build_records(terms: dict[str, dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    parent_map = {
        term_id: [parent for parent in term.get('is_a', []) if parent in terms]
        for term_id, term in terms.items()
    }
    ancestor_map = {
        term_id: _ancestor_closure(term_id, parent_map)
        for term_id in terms
    }

    included_terms = {
        term_id
        for term_id, term in terms.items()
        if _is_included_chemical_term(term_id, term, ancestor_map)
    }

    child_map: dict[str, set[str]] = defaultdict(set)
    for child_id, parents in parent_map.items():
        if child_id not in included_terms:
            continue
        for parent_id in parents:
            if parent_id in included_terms:
                child_map[parent_id].add(child_id)

    structured_terms = {
        term_id
        for term_id in included_terms
        if _has_structure(terms[term_id])
    }

    structured_leaf_terms = {
        term_id
        for term_id in structured_terms
        if not child_map.get(term_id)
    }

    ontology_term_ids = included_terms - structured_leaf_terms

    ontology_terms = [
        {
            'id': term_id,
            'name': terms[term_id].get('name', ''),
            'definition': terms[term_id].get('definition', ''),
            'synonyms': terms[term_id].get('synonyms', []),
            'comments': terms[term_id].get('comments', []),
            'xrefs': terms[term_id].get('xrefs', []),
            'is_a': [parent for parent in terms[term_id].get('is_a', []) if parent in ontology_term_ids],
            'relationships': [
                relationship
                for relationship in terms[term_id].get('relationships', [])
                if relationship.get('target') in ontology_term_ids
            ],
            'is_obsolete': terms[term_id].get('is_obsolete', False),
        }
        for term_id in sorted(ontology_term_ids)
    ]

    molecules: list[dict[str, Any]] = []
    for term_id in sorted(structured_terms):
        term = terms[term_id]
        xrefs = _extract_identifier_xrefs(term.get('xrefs', []))
        properties = term.get('properties', {})
        molecules.append(
            {
                'chebi_id': term_id,
                'name': term.get('name', ''),
                'synonyms': term.get('synonyms', []),
                'definition': term.get('definition', ''),
                'ancestor_terms': sorted(ancestor_map.get(term_id, set()) & ontology_term_ids),
                'kegg_compound': sorted(xrefs['kegg_compound']),
                'pubchem_compound': sorted(xrefs['pubchem_compound']),
                'hmdb': sorted(xrefs['hmdb']),
                'lipidmaps': sorted(xrefs['lipidmaps']),
                'cas': sorted(xrefs['cas']),
                'smiles': properties.get('smiles', ''),
                'inchi': properties.get('inchi', ''),
                'inchikey': properties.get('inchikey', ''),
                'formula': properties.get('formula', ''),
                'mass': properties.get('mass', ''),
                'monoisotopic_mass': properties.get('monoisotopic_mass', ''),
                'charge': properties.get('charge', ''),
            }
        )

    return {
        'ontology_terms': ontology_terms,
        'molecules': molecules,
    }


def _ancestor_closure(term_id: str, parent_map: dict[str, list[str]]) -> set[str]:
    ancestors: set[str] = set()
    stack = list(parent_map.get(term_id, []))

    while stack:
        parent = stack.pop()
        if parent in ancestors:
            continue
        ancestors.add(parent)
        stack.extend(parent_map.get(parent, []))

    return ancestors


def _is_included_chemical_term(term_id: str, term: dict[str, Any], ancestor_map: dict[str, set[str]]) -> bool:
    if term.get('is_obsolete'):
        return False
    if not term.get('name'):
        return False

    ancestors = ancestor_map.get(term_id, set())

    if term_id != _CHEMICAL_ENTITY_ROOT and _CHEMICAL_ENTITY_ROOT not in ancestors:
        return False

    if any(root == term_id or root in ancestors for root in _EXCLUDED_ROOTS):
        return False

    return True


def _has_structure(term: dict[str, Any]) -> bool:
    properties = term.get('properties', {})
    return any(properties.get(key) for key in _STRUCTURE_SIGNAL_KEYS)


def _extract_identifier_xrefs(xrefs: Iterable[str]) -> dict[str, set[str]]:
    result: dict[str, set[str]] = defaultdict(set)

    for xref in xrefs:
        text = xref.strip()
        lower = text.lower()

        match = _KEGG_RE.search(text)
        if match:
            result['kegg_compound'].add(match.group(1).upper())

        match = _PUBCHEM_RE.search(text)
        if match:
            result['pubchem_compound'].add(match.group(1))

        match = _HMDB_RE.search(text)
        if match and 'hmdb' in lower:
            result['hmdb'].add(match.group(1).upper())

        match = _LIPIDMAPS_RE.search(text)
        if match:
            result['lipidmaps'].add(match.group(1).upper())

        match = _CAS_RE.search(text)
        if match and ('cas' in lower or 'registry' in lower):
            result['cas'].add(match.group(1))

    return result
