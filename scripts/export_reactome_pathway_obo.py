#!/usr/bin/env python
"""Export Reactome pathways as an OBO ontology.

Builds a Reactome pathway ontology from the BioPAX pathway records parsed by
``pypath.inputs_v2.parsers.reactome`` and writes an OBO file suitable for the
ontology service.

Usage:
    uv run python pypath/scripts/export_reactome_pathway_obo.py [output_path]
"""

from __future__ import annotations

import sys
from collections import defaultdict
from datetime import datetime, UTC
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from pypath.inputs_v2.reactome import download
from pypath.inputs_v2.parsers.reactome import _raw


def escape_obo_string(value: str) -> str:
    return value.replace('\\', '\\\\').replace('"', '\\"').replace('\n', ' ')


def _split_field(value: str, delimiter: str = '||') -> list[str]:
    if not value:
        return []
    return [item for item in value.split(delimiter) if item and item != '__MISSING__']


def load_pathway_records(force_refresh: bool = False) -> list[dict]:
    opener = download.open(force_refresh=force_refresh)
    return list(_raw(opener, data_type='pathways', force_refresh=force_refresh))


def build_terms(records: list[dict]) -> list[dict]:
    terms: dict[str, dict] = {}
    parent_map: defaultdict[str, list[tuple[str, str]]] = defaultdict(list)

    for record in records:
        stable_id = record.get('reactome_stable_id', '').split(';')[0]
        if not stable_id:
            continue

        terms.setdefault(
            stable_id,
            {
                'id': stable_id,
                'name': record.get('display_name', ''),
                'definition': record.get('definition', ''),
                'synonyms': [s for s in record.get('synonyms', '').split(';') if s],
                'comments': [c for c in record.get('comments', '').split(';') if c],
                'xrefs': [],
            },
        )

        reactome_id = record.get('reactome_id', '')
        if reactome_id:
            terms[stable_id]['xrefs'].append(f'Reactome:{reactome_id}')

        for go_id in [g for g in record.get('go', '').split(';') if g]:
            terms[stable_id]['xrefs'].append(go_id)

        taxon_id = record.get('ncbi_tax_id', '')
        if taxon_id:
            terms[stable_id]['xrefs'].append(f'NCBITaxon:{taxon_id}')

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
                },
            )

    for term_id, parents in parent_map.items():
        seen = set()
        deduped = []
        for parent_id, parent_name in parents:
            key = (parent_id, parent_name)
            if key in seen:
                continue
            seen.add(key)
            deduped.append((parent_id, parent_name))
        terms[term_id]['parents'] = deduped

    for term in terms.values():
        term.setdefault('parents', [])
        term['synonyms'] = sorted(set(term['synonyms']))
        term['comments'] = list(dict.fromkeys(term['comments']))
        term['xrefs'] = list(dict.fromkeys(term['xrefs']))

    return sorted(terms.values(), key=lambda t: t['id'])


def format_obo(terms: list[dict]) -> str:
    date_str = datetime.now(UTC).strftime('%d:%m:%Y %H:%M')
    lines = [
        'format-version: 1.2',
        'ontology: reactome_pathways',
        f'date: {date_str}',
        'default-namespace: reactome_pathways',
        'remark: Reactome pathway ontology exported from Reactome BioPAX via pypath.',
        '',
        '[Typedef]',
        'id: part_of',
        'name: part_of',
        '',
    ]

    for term in terms:
        lines.append('[Term]')
        lines.append(f"id: {term['id']}")
        lines.append(f"name: {escape_obo_string(term['name'] or term['id'])}")

        if term['definition']:
            lines.append(f'def: "{escape_obo_string(term["definition"])}" []')

        for synonym in term['synonyms']:
            lines.append(f'synonym: "{escape_obo_string(synonym)}" EXACT []')

        for comment in term['comments']:
            lines.append(f'comment: {escape_obo_string(comment)}')

        for xref in term['xrefs']:
            lines.append(f'xref: {xref}')

        for parent_id, parent_name in term['parents']:
            suffix = f' ! {escape_obo_string(parent_name)}' if parent_name else ''
            lines.append(f'relationship: part_of {parent_id}{suffix}')

        lines.append('')

    return '\n'.join(lines)


def main(output_path: str = 'reactome_pathways.obo', force_refresh: bool = False) -> None:
    records = load_pathway_records(force_refresh=force_refresh)
    terms = build_terms(records)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(format_obo(terms), encoding='utf-8')
    print(f'Wrote {len(terms)} Reactome pathway terms to {output.resolve()}')


if __name__ == '__main__':
    output = sys.argv[1] if len(sys.argv) > 1 else 'reactome_pathways.obo'
    refresh = '--force-refresh' in sys.argv[2:] or '--force-refresh' in sys.argv[1:]
    main(output, force_refresh=refresh)
