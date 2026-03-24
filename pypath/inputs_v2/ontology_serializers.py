"""Serializers for ontology datasets."""
from __future__ import annotations

from datetime import UTC, datetime

from pypath.internals.ontology_schema import OntologyDocument, OntologyTerm

__all__ = ['escape_obo_string', 'format_obo']


def escape_obo_string(value: str) -> str:
    r"""Escape OBO text while preserving already-valid OBO escape sequences.

    Some upstream parsers pass through source values that are already OBO-escaped
    (for example ``\{`` / ``\}`` in ChEBI names). The old implementation
    doubled every backslash, turning valid escapes into invalid ones like
    ``\\{``. Preserve known escape sequences first, then escape any remaining
    reserved characters.
    """
    placeholders = {
        '\\{': '__OBO_ESCAPED_OPEN__',
        '\\}': '__OBO_ESCAPED_CLOSE__',
        '\\\\': '__OBO_ESCAPED_BACKSLASH__',
        '\\"': '__OBO_ESCAPED_QUOTE__',
        '\\n': '__OBO_ESCAPED_NEWLINE__',
    }

    for old, marker in placeholders.items():
        value = value.replace(old, marker)

    value = (
        value
        .replace('\\', '\\\\')
        .replace('{', '\\{')
        .replace('}', '\\}')
        .replace('"', '\\"')
        .replace('\n', ' ')
    )

    for old, marker in reversed(list(placeholders.items())):
        value = value.replace(marker, old)

    return value


def format_obo(document: OntologyDocument, terms: list[OntologyTerm]) -> str:
    date_str = datetime.now(UTC).strftime('%d:%m:%Y %H:%M')
    lines = [
        'format-version: 1.2',
        f'ontology: {document.ontology}',
        f'date: {date_str}',
    ]
    if document.default_namespace:
        lines.append(f'default-namespace: {document.default_namespace}')
    if document.remark:
        lines.append(f'remark: {document.remark}')
    lines.append('')

    for typedef in document.typedefs or []:
        lines.extend([
            '[Typedef]',
            f'id: {typedef.id}',
            f'name: {typedef.name}',
            '',
        ])

    for term in sorted(terms, key=lambda item: item.id):
        lines.append('[Term]')
        lines.append(f'id: {term.id}')
        lines.append(f'name: {escape_obo_string(term.name or term.id)}')

        if term.definition:
            lines.append(f'def: "{escape_obo_string(term.definition)}" []')

        for parent in term.is_a or []:
            lines.append(f'is_a: {parent}')

        for synonym in term.synonyms or []:
            lines.append(f'synonym: "{escape_obo_string(synonym)}" EXACT []')

        for comment in term.comments or []:
            lines.append(f'comment: {escape_obo_string(comment)}')

        for xref in term.xrefs or []:
            lines.append(f'xref: {xref}')

        for relationship in term.relationships or []:
            suffix = (
                f' ! {escape_obo_string(relationship.target_name)}'
                if relationship.target_name else ''
            )
            lines.append(
                f'relationship: {relationship.type} {relationship.target}{suffix}'
            )

        if term.is_obsolete:
            lines.append('is_obsolete: true')

        lines.append('')

    return '\n'.join(lines)
