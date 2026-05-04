"""Generic OBO 1.2 parser for ontology datasets."""

from __future__ import annotations

from collections.abc import Generator, Iterable
import re
from typing import Any

from pypath.internals.ontology_schema import OntologyRelationship, OntologyTerm

_TRUE_VALUES = {'1', 'true', 'yes'}


def iter_obo(opener, *, preprocess_text=None, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield normalized ``[Term]`` records from an OBO opener.

    The parser keeps commonly used OBO tags and ignores non-term stanzas.
    ``preprocess_text`` can be used for source-specific text cleanup before
    parsing (for example PSI-MI date-line cleanup).
    """
    text = _read_text(opener)
    if preprocess_text is not None:
        text = preprocess_text(text)
    yield from parse_obo_text(text)


def parse_obo_text(text: str) -> Generator[dict[str, Any], None, None]:
    ontology_id: str | None = None
    current: dict[str, Any] | None = None

    for raw_line in text.splitlines():
        stripped = raw_line.strip()

        if not stripped or stripped.startswith('!'):
            if not stripped and current is not None:
                yield _finalize_term(current)
                current = None
            continue

        if current is None and ontology_id is None and stripped.startswith('ontology:'):
            ontology_id = stripped.partition(':')[2].strip() or None
            continue

        if stripped == '[Term]':
            if current is not None:
                yield _finalize_term(current)
            current = _new_term(ontology_id)
            continue

        if stripped.startswith('['):
            if current is not None:
                yield _finalize_term(current)
                current = None
            continue

        if current is None:
            continue

        tag, sep, value = stripped.partition(':')
        if not sep:
            continue
        _apply_tag(current, tag.strip(), value.strip())

    if current is not None:
        yield _finalize_term(current)


def obo_record_to_term(row: dict[str, Any]) -> OntologyTerm | None:
    """Map a parsed OBO record to the canonical ``OntologyTerm`` schema."""
    term_id = row.get('id')
    name = row.get('name') or term_id
    if not term_id or not name:
        return None

    relationships = [
        OntologyRelationship(
            type=relationship['type'],
            target=relationship['target'],
            target_name=relationship.get('target_name') or None,
        )
        for relationship in row.get('relationships', [])
        if relationship.get('type') and relationship.get('target')
    ]

    return OntologyTerm(
        id=term_id,
        name=name,
        definition=row.get('definition') or None,
        synonyms=row.get('synonyms') or None,
        comments=row.get('comments') or None,
        xrefs=row.get('xrefs') or None,
        is_a=row.get('is_a') or None,
        relationships=relationships or None,
        is_obsolete=bool(row.get('is_obsolete')),
        alt_ids=row.get('alt_ids') or None,
        namespace=row.get('namespace') or None,
    )


def _read_text(opener) -> str:
    handle = None
    if opener and getattr(opener, 'result', None):
        handle = next(iter(opener.result.values()), None) if isinstance(opener.result, dict) else opener.result
    if handle is None:
        return ''
    if hasattr(handle, 'seek'):
        handle.seek(0)
    if hasattr(handle, 'read'):
        content = handle.read()
        return content.decode('utf-8', 'ignore') if isinstance(content, bytes) else str(content)
    chunks = []
    for chunk in handle:
        chunks.append(chunk.decode('utf-8', 'ignore') if isinstance(chunk, bytes) else str(chunk))
    return ''.join(chunks)


def _new_term(ontology_id: str | None) -> dict[str, Any]:
    return {
        'id': '',
        'name': '',
        'definition': '',
        'synonyms': [],
        'alt_ids': [],
        'namespace': None,
        'comments': [],
        'xrefs': [],
        'is_a': [],
        'relationships': [],
        'is_obsolete': False,
        'ontology': ontology_id,
    }


def _finalize_term(term: dict[str, Any]) -> dict[str, Any]:
    return {
        **term,
        'synonyms': _dedupe(term.get('synonyms', [])),
        'alt_ids': _dedupe(term.get('alt_ids', [])),
        'comments': _dedupe(term.get('comments', [])),
        'xrefs': _dedupe(term.get('xrefs', [])),
        'is_a': _dedupe(term.get('is_a', [])),
        'relationships': _dedupe_relationships(term.get('relationships', [])),
    }


def _apply_tag(term: dict[str, Any], tag: str, value: str) -> None:
    if tag == 'id' and not term['id']:
        term['id'] = value
    elif tag == 'name' and not term['name']:
        term['name'] = value
    elif tag == 'def' and not term['definition']:
        term['definition'] = _extract_quoted(value)
    elif tag == 'synonym':
        synonym = _extract_quoted(value)
        if synonym:
            term['synonyms'].append(synonym)
    elif tag == 'alt_id':
        term['alt_ids'].append(_strip_inline_comment(value))
    elif tag == 'namespace':
        term['namespace'] = value or None
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
    elif tag == 'is_obsolete':
        term['is_obsolete'] = value.strip().lower() in _TRUE_VALUES


def _extract_quoted(value: str) -> str:
    match = re.search(r'"((?:\\.|[^"\\])*)"', value)
    if match:
        return bytes(match.group(1), 'utf-8').decode('unicode_escape').strip()
    return _strip_inline_comment(value)


def _strip_inline_comment(value: str) -> str:
    return value.split(' ! ', 1)[0].strip()


def _parse_relationship(value: str) -> dict[str, str] | None:
    base, _, target_name = value.partition(' ! ')
    parts = base.strip().split(None, 2)
    if len(parts) < 2:
        return None
    return {'type': parts[0], 'target': parts[1], 'target_name': target_name.strip()}


def _dedupe(values: Iterable[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value and value not in seen:
            seen.add(value)
            out.append(value)
    return out


def _dedupe_relationships(values: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    seen = set()
    out = []
    for relationship in values:
        key = (relationship.get('type'), relationship.get('target'), relationship.get('target_name'))
        if relationship.get('type') and relationship.get('target') and key not in seen:
            seen.add(key)
            out.append(relationship)
    return out
