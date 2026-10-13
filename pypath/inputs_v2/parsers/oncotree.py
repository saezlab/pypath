"""Parser of the OncoTree tumor type tree (``/api/tumorTypes/tree``)."""

from __future__ import annotations

from collections.abc import Generator
import json
import re
from typing import Any

from pypath.inputs_v2.parsers.base import read_text
from pypath.internals.ontology_schema import OntologyTerm

# The databases of ``externalReferences`` and how the ontology modules write them.
XREF_PREFIXES = {'UMLS': 'UMLS', 'NCI': 'NCIT'}
# A previous code of a tumor type is an identifier of it when it is a plain code.
_PLAIN_CODE = re.compile(r'[A-Za-z0-9][A-Za-z0-9_/-]*')


def _flatten(code: str, node: dict[str, Any], parent: str | None, out: list[dict[str, Any]]) -> None:
    if node.get('code') != code:
        raise ValueError(f'OncoTree node {code!r} has the code {node.get("code")!r}')
    if node.get('parent') != parent:
        raise ValueError(
            f'OncoTree node {code!r} is under {parent!r} but names the parent {node.get("parent")!r}'
        )
    name = node.get('name')
    if not isinstance(name, str) or not name.strip():
        raise ValueError(f'OncoTree node {code!r} has no name')
    references = node.get('externalReferences') or {}
    unknown = set(references) - set(XREF_PREFIXES)
    if unknown:
        raise ValueError(f'OncoTree node {code!r} has references to {sorted(unknown)}')
    out.append({
        'code': code,
        'name': name.strip(),
        'parent': parent,
        'level': node.get('level'),
        'main_type': node.get('mainType') or None,
        'tissue': node.get('tissue') or None,
        'xrefs': {database: list(references.get(database) or []) for database in XREF_PREFIXES},
        'history': list(node.get('history') or []),
    })
    for child_code, child in (node.get('children') or {}).items():
        _flatten(child_code, child, code, out)


def iter_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the tumor types of the tree, each after its parent.

    The tree (one JSON object: the root's code to the root node, every node
    holding its ``children`` by code) is the one document that has the root,
    ``TISSUE``, which the flat list of ``/api/tumorTypes`` lacks. A row has
    ``code``, ``name``, ``parent`` (None for the root), ``level``, ``main_type``
    and ``tissue`` (None where the source has none), ``xrefs`` (the codes of
    ``UMLS`` and ``NCI``, as lists) and ``history`` (the previous codes).

    A code that is repeated, a node without a name or a previous code that is
    the code of a tumor type raises.
    """
    rows: list[dict[str, Any]] = []
    text = read_text(opener)
    if not text.strip():
        raise ValueError('OncoTree download is empty')
    for code, node in json.loads(text).items():
        _flatten(code, node, None, rows)
    codes = [row['code'] for row in rows]
    if len(set(codes)) != len(codes):
        repeated = sorted({code for code in codes if codes.count(code) > 1})
        raise ValueError(f'OncoTree codes are repeated: {repeated}')
    live = set(codes)
    for row in rows:
        reused = live.intersection(row['history'])
        if reused:
            raise ValueError(f'OncoTree {row["code"]} lists the live codes {sorted(reused)} as history')
    yield from rows


def oncotree_record_to_term(row: dict[str, Any]) -> OntologyTerm | None:
    """The canonical ``OntologyTerm`` of a tumor type: its xrefs, parent and previous codes."""
    if not row.get('code') or not row.get('name'):
        return None
    return OntologyTerm(
        id=row['code'],
        name=row['name'],
        xrefs=[
            f'{XREF_PREFIXES[database]}:{code}'
            for database, codes in (row.get('xrefs') or {}).items()
            for code in codes
        ]
        or None,
        is_a=[row['parent']] if row.get('parent') else None,
        alt_ids=[code for code in row.get('history') or [] if _PLAIN_CODE.fullmatch(code)] or None,
    )
